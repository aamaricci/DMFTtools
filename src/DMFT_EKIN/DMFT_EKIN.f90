MODULE DMFT_EKIN
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS,   only:free_unit,str,reg
  USE SF_ARRAYS,    only:arange
  USE SF_TIMER
  USE SF_LINALG,    only:inv,eye,eigh,diag
  USE SF_MISC,      only:assert_shape
  USE SF_SPECIAL,   only:fermi
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  interface dmft_kinetic_energy
     !NsoNso-Shape ([Nlat])[Nso][Nso][L] for Sigma & Self
     module procedure :: dmft_kinetic_energy_normal_Nso_main
     module procedure :: dmft_kinetic_energy_normal_Nso_dos
     module procedure :: dmft_kinetic_energy_normal_Nso_lattice
     module procedure :: dmft_kinetic_energy_superc_Nso_main
     module procedure :: dmft_kinetic_energy_superc_Nso_dos
     module procedure :: dmft_kinetic_energy_superc_Nso_lattice
     !NNNN-Shape ([Nlat])[Nspin][Nspin][Norb][Norb][L] for Sigma & Self
     module procedure :: dmft_kinetic_energy_normal_NN
     module procedure :: dmft_kinetic_energy_normal_NN_dos
     module procedure :: dmft_kinetic_energy_normal_NN_lattice
     module procedure :: dmft_kinetic_energy_superc_NN
     module procedure :: dmft_kinetic_energy_superc_NN_dos
     module procedure :: dmft_kinetic_energy_superc_NN_lattice
  end interface dmft_kinetic_energy


  interface print_Hloc
     module procedure print_Hloc_2
     module procedure print_Hloc_4
  end interface print_Hloc



  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block


  real(8),dimension(:),allocatable :: wm
  integer                          :: mpi_ierr
  integer                          :: mpi_rank
  integer                          :: mpi_size
  logical                          :: mpi_master

  !PUBLIC in DMFT
  public :: dmft_kinetic_energy

contains 






  !----------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy of the lattice model. 
  ! input: 
  ! 1. the Hamiltonian matrix H(k)
  ! 2. DMFT normal self-energy Sigma.
  !
  ! The main routines accept self-energy as:
  ! - Sigma: [Nspin*Norb][Nspin*Norb][L]
  ! - Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  !----------------------------------------------------------------------------------------
  subroutine dmft_kinetic_energy_normal_Nso_main(Hk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)                 :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:)                 :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    real(8),dimension(size(Hk,1)),optional      :: Ekin,Eloc
    !
    integer                                     :: Lk,Nso,Liw
    integer                                     :: i,ik,iso
    !
    integer                                     :: Norb,Nporb
    integer                                     :: Nspin  
    real(8)                                     :: beta
    real(8)                                     :: xmu
    !
    real(8),dimension(size(Hk,1),size(Hk,1))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Ak,Bk,Ck,Dk,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Gk,Tk
    !
    real(8),dimension(size(Hk,1))               :: Tail0,Tail1
    real(8),dimension(size(Hk,1))               :: Lail0,Lail1
    real(8)                                     :: spin_degeneracy
    !
    real(8),dimension(size(Hk,1))               :: H0,Hl
    real(8),dimension(size(Hk,1))               :: H0tmp,Hltmp
    !
    real(8),dimension(size(Hk,1))               :: Ekin_,Eloc_
    !
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nso = size(Hk,1)
    Lk  = size(Hk,3)
    Liw = size(Sigma,3)
    !Testing:
    !if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_main_mpi: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_normal_Nso_main","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_Nso_main","Sigma")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
    if(size(Hloc,1)<16)then
       if(mpi_master)call print_hloc(Hloc)
    else
       if(mpi_master)call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer()
    H0=0d0
    Hl=0d0
    H0tmp= 0d0
    Hltmp= 0d0
    do ik=1,Lk
       Ak = Hk(:,:,ik) - Hloc(:,:)
       Bk =-Hk(:,:,ik) - Sigma_HF(:,:)
       do i=1+mpi_rank,Liw,mpi_size
          Gk = (xi*wm(i)+xmu)*eye(Nso) - Sigma(:,:,i) - Hk(:,:,ik) 
          select case(Nso)
          case default
             call inv(Gk)
          case(1)
             Gk = 1d0/Gk
          end select
          Tk = eye(Nso)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
          Ck = matmul(Ak  ,Gk - Tk)
          Dk = matmul(Hloc,Gk - Tk)
          do iso=1,Nso
             H0tmp(iso) = H0tmp(iso) + Ck(iso,iso)/dble(Lk)
             Hltmp(iso) = Hltmp(iso) + Dk(iso,iso)/dble(Lk)
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer()
    spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2*spin_degeneracy
    Hl=Hl/beta*2*spin_degeneracy
    !
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       Ck= matmul(Ak,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Dk= matmul(Hloc,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       do iso=1,Nso
          Tail0(iso) = Tail0(iso) + 0.5d0*Ak(iso,iso)/dble(Lk)
          Tail1(iso) = Tail1(iso) + 0.25d0*Ck(iso,iso)/dble(Lk)
          Lail0(iso) = Lail0(iso) + 0.5d0*Hloc(iso,iso)/dble(Lk)
          Lail1(iso) = Lail1(iso) + 0.25d0*Dk(iso,iso)/dble(Lk)
       enddo
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin_=H0+Tail0+Tail1
    Eloc_=Hl+Lail0+Lail1
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_normal_Nso_main


  subroutine dmft_kinetic_energy_normal_Nso_dos(Ebands,Dbands,Hloc,Sigma,Ekin,Eloc)
    real(8),dimension(:,:),intent(in)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:)                                 :: Sigma   ![Nspin*Norb][Nspin*Norb][L]
    real(8),dimension(size(Ebands,1)),optional                  :: Ekin,Eloc
    !
    integer                                                     :: Lk,Nso,Liw
    integer                                                     :: i,ik,iso
    !
    integer                                                     :: Norb,Nporb
    integer                                                     :: Nspin  
    real(8)                                                     :: beta
    real(8)                                                     :: xmu
    !
    real(8),dimension(size(Ebands,1),size(Ebands,1))            :: Sigma_HF
    !
    complex(8)                                                  :: Ak,Bk,Ck,Dk
    complex(8)                                                  :: Gk,Tk
    !
    real(8),dimension(size(Ebands,1))                           :: Tail0,Tail1
    real(8),dimension(size(Ebands,1))                           :: Lail0,Lail1
    real(8)                                                     :: spin_degeneracy
    !
    real(8),dimension(size(Ebands,1))                           :: H0,Hl
    real(8),dimension(size(Ebands,1))                           :: H0tmp,Hltmp
    !
    real(8),dimension(size(Ebands,1))                           :: Ekin_,Eloc_
    !
    real(8),dimension(:),allocatable                            :: wm
    !
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nso = size(Ebands,1)
    Lk  = size(Ebands,2)
    Liw = size(Sigma,3)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_Nso_dos: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_Nso_dos","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    !
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer()
    H0=0d0
    Hl=0d0
    H0tmp= 0d0
    Hltmp= 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       do iso=1,Nso
          Ak = Ebands(iso,ik)
          Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso) 
          do i=1+mpi_rank,Liw,mpi_size
             Gk = (xi*wm(i)+xmu) - Sigma(iso,iso,i) - Ebands(iso,ik) - Hloc(iso)
             Gk = 1d0/Gk
             Tk = 1d0/(xi*wm(i)) - Bk/(xi*wm(i))**2
             Ck = Ak*(Gk - Tk)
             Dk = Hloc(iso)*(Gk - Tk)
             H0tmp(iso) = H0tmp(iso) + Dbands(iso,ik)*Ck
             Hltmp(iso) = Hltmp(iso) + Dbands(iso,ik)*Dk
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer()
    spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2*spin_degeneracy
    Hl=Hl/beta*2*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       do iso=1,Nso
          Ak = Ebands(iso,ik)
          Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
          Ck= Ak*Bk
          Dk= Hloc(iso)*Bk
          Tail0(iso) = Tail0(iso) + 0.5d0*Dbands(iso,ik)*Ak
          Tail1(iso) = Tail1(iso) + 0.25d0*Dbands(iso,ik)*Ck
          Lail0(iso) = Lail0(iso) + 0.5d0*Dbands(iso,ik)*Hloc(iso)
          Lail1(iso) = Lail1(iso) + 0.25d0*Dbands(iso,ik)*Dk
       enddo
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin_=H0+Tail0+Tail1
    Eloc_=Hl+Lail0+Lail1
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_normal_Nso_dos


  subroutine dmft_kinetic_energy_normal_Nso_lattice(Hk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma     ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    real(8),dimension(size(Hk,1)),optional                          :: Ekin,Eloc
    !aux
    integer                                                         :: Lk,Nlso,Liw,Nlat,Nso
    integer                                                         :: ik,iso
    integer                                                         :: i,iorb,ilat,ispin,io,is
    integer                                                         :: j,jorb,jlat,jspin,jo,js
    !
    integer                                                         :: Norb,Nporb
    integer                                                         :: Nspin  
    real(8)                                                         :: beta
    real(8)                                                         :: xmu
    !
    complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Ak,Bk,Ck,Dk,Hloc,Hloc_tmp
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Gk
    !
    real(8),dimension(size(Hk,1))                                   :: Tail0,Tail1
    real(8),dimension(size(Hk,1))                                   :: Lail0,Lail1
    real(8)                                                         :: spin_degeneracy
    !
    real(8),dimension(size(Hk,1))                                   :: H0,Hl
    real(8),dimension(size(Hk,1))                                   :: H0tmp,Hltmp
    !
    real(8),dimension(size(Hk,1))                                   :: Ekin_,Eloc_
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Nlat = size(Sigma,1)
    Nso  = size(Sigma,2)
    Liw  = size(Sigma,4)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_Nso_lattice: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_normal_Nso_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_normal_Nso_lattice","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF(:,:,:) = dreal(Sigma(:,:,:,Liw))
    !
    ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
    Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
    Hloc=0d0
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hloc(is,js)=Hloc_tmp(is,js) 
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
    if(size(Hloc,1)<16)then
       if(mpi_master)call print_hloc(Hloc)
    else
       if(mpi_master)call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    !Start the timer:
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer
    H0 = 0d0
    Hl = 0d0
    H0tmp= 0d0
    Hltmp= 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    =  Hk(:,:,ik) - Hloc(:,:)
       Bk    = -Hk(:,:,ik) - blocks_to_matrix(Sigma_HF(:,:,:),Nlat,Nso) !Sigma_HF [Nlat,Nso,Nso]--> [Nlso,Nslo]
       do i=1+mpi_rank,Liw,mpi_size
          Gk = (xi*wm(i)+xmu)*eye(Nlso) - blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso) - Hk(:,:,ik) !Sigma [Nlat,Nso,Nso,*]--> [Nlso,Nslo]
          call inv(Gk(:,:))
          Tk = eye(Nlso)/(xi*wm(i)) - Bk/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          Dk = matmul(Hloc,Gk(:,:) - Tk)
          do is=1,Nlso
             H0tmp(is) = H0tmp(is) + Ck(is,is)/dble(Lk)
             Hltmp(is) = Hltmp(is) + Dk(is,is)/dble(Lk)
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer
    spin_degeneracy=3d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    Hl = Hl/beta*2d0*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    do ik=1,Lk
       Ak    =  Hk(:,:,ik) - Hloc(:,:)
       Bk    = -Hk(:,:,ik) - blocks_to_matrix(Sigma_HF(:,:,:),Nlat,Nso) !Sigma_HF [Nlat,Nso,Nso]--> [Nlso,Nslo]
       Ck= matmul(Ak,Bk)
       Dk= matmul(Hloc,Bk)
       do is=1,Nlso
          Tail0(is) = Tail0(is) + 0.5d0*Ak(is,is)/dble(Lk)
          Tail1(is) = Tail1(is) + 0.25d0*Ck(is,is)/dble(Lk)
          Lail0(is) = Lail0(is) + 0.5d0*Hloc(is,is)/dble(Lk)
          Lail1(is) = Lail1(is) + 0.25d0*Dk(is,is)/dble(Lk)
       enddo
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin_=H0+Tail0+Tail1
    Eloc_=Hl+Lail0+Lail1
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_,Nlat,Nso)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_normal_Nso_lattice









  subroutine dmft_kinetic_energy_superc_Nso_main(Hk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)                     :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:)                     :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:)                     :: Self ![Nspin*Norb][Nspin*Norb][L]
    real(8),dimension(size(Hk,1)),optional          :: Ekin,Eloc
    !
    integer                                         :: Lk,Nso,Liw
    integer                                         :: i,is,ik
    !
    integer                                         :: Norb,Nporb
    integer                                         :: Nspin  
    real(8)                                         :: beta
    real(8)                                         :: xmu
    !
    real(8),dimension(size(Hk,1),size(Hk,1))        :: Sigma_HF,Self_HF
    complex(8),dimension(size(Hk,1),size(Hk,1))     :: Ak,Bk,Ck,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,1))     :: Gk,Tk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,1)) :: GkNambu,HkNambu,HlocNambu,AkNambu
    real(8),dimension(2*size(Hk,1))                 :: NkNambu
    complex(8),dimension(2*size(Hk,1),2*size(Hk,1)) :: Evec
    real(8),dimension(2*size(Hk,1))                 :: Eval,Coef
    real(8)                                         :: spin_degeneracy
    !
    real(8),dimension(size(Hk,1))                   :: H0,Hl
    real(8),dimension(size(Hk,1))                   :: H0tmp,Hltmp
    real(8),dimension(size(Hk,1))                   :: H0free,Hlfree
    !
    real(8),dimension(size(Hk,1))                   :: Ekin_,Eloc_

    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nso = size(Hk,1)
    Lk  = size(Hk,3)
    Liw = size(Sigma,3)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_Nso_main: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_superc_Nso_main","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_Nso_main","Sigma")
    call assert_shape(Self,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_Nso_main","Self")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,Liw))
    Self_HF  = dreal(Self(:,:,Liw))
    !
    !Get the local Hamiltonian
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
    if(size(Hloc,1)<16)then
       if(mpi_master)call print_hloc(Hloc)
    else
       if(mpi_master)call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer()
    H0=0d0
    Hl=0d0
    H0tmp= 0d0
    Hltmp= 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak= Hk(:,:,ik)-Hloc
       do i=1+mpi_rank,Liw,mpi_size
          Gknambu=zero
          Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)&
               - Sigma(:,:,i)        - Hk(:,:,ik)
          Gknambu(1:Nso,Nso+1:2*Nso)       = - Self(:,:,i)
          Gknambu(Nso+1:2*Nso,1:Nso)       = - Self(:,:,i)
          Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)&
               + conjg(Sigma(:,:,i)) + Hk(:,:,ik)
          call inv(Gknambu)
          Gk  = Gknambu(1:Nso,1:Nso) !Gk(iw)
          !
          Gknambu=zero
          Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)&
               - Sigma_hf  - Hk(:,:,ik)
          Gknambu(1:Nso,Nso+1:2*Nso)       = - Self_hf
          Gknambu(Nso+1:2*Nso,1:Nso)       = - Self_hf
          Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)&
               + Sigma_hf  + Hk(:,:,ik)
          call inv(Gknambu)
          Tk = Gknambu(1:Nso,1:Nso) !G0k(iw)
          !
          Bk = matmul(Ak, Gk - Tk) !
          Ck = matmul(Hloc,Gk - Tk)
          do is=1,Nso
             H0tmp(is) = H0tmp(is) + Bk(is,is)/dble(Lk)
             Hltmp(is) = Hltmp(is) + Ck(is,is)/dble(Lk)
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer
    spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    Hl=Hl/beta*2.d0*spin_degeneracy
    !
    !
    H0free=0d0
    Hlfree=0d0
    HkNambu=zero
    HlocNambu=zero
    do ik=1,Lk
       HkNambu(1:Nso,1:Nso)               =  Hk(:,:,ik)-Hloc
       HkNambu(Nso+1:2*Nso,Nso+1:2*Nso)   = -Hk(:,:,ik)+Hloc
       HlocNambu(1:Nso,1:Nso)             =  Hloc
       HlocNambu(Nso+1:2*Nso,Nso+1:2*Nso) = -Hloc
       Evec(1:Nso,1:Nso)                 =  Hk(:,:,ik) + Sigma_hf
       Evec(1:Nso,Nso+1:2*Nso)           =             + Self_hf
       Evec(Nso+1:2*Nso,1:Nso)           =             + Self_hf
       Evec(Nso+1:2*Nso,Nso+1:2*Nso)     = -Hk(:,:,ik) - Sigma_hf
       call eigh(Evec,Eval)
       do is=1,2*Nso
          Coef = Evec(:,is)*conjg(Evec(:,is))
          NkNambu(is) = dot_product(Coef,fermi(Eval,beta))
       enddo
       AkNambu = matmul( diag(NkNambu) , HkNambu )
       do is=1,Nso
          H0free(is) = H0free(is) + AkNambu(is,is)/dble(Lk) !Take only the 11 part
       enddo
       AkNambu = matmul( diag(NkNambu) , HlocNambu )
       do is=1,Nso
          Hlfree(is) = Hlfree(is) + AkNambu(is,is)/dble(Lk)
       enddo
    enddo
    H0free=spin_degeneracy*H0free
    Hlfree=spin_degeneracy*Hlfree
    !
    Ekin_=H0+H0free
    Eloc_=Hl+Hlfree
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_superc_Nso_main


  !
  !DOS
  subroutine dmft_kinetic_energy_superc_Nso_dos(Ebands,Dbands,Hloc,Sigma,Self,Ekin,Eloc)
    real(8),dimension(:,:),intent(in)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:)                                 :: Sigma   ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:)                                 :: Self    ![Nspin*Norb][Nspin*Norb][L]
    real(8),dimension(size(Ebands,1)),optional                  :: Ekin,Eloc
    !
    integer                                                     :: Lk,Nso,Liw
    integer                                                     :: i,is,ik
    !
    integer                                                     :: Norb,Nporb
    integer                                                     :: Nspin  
    real(8)                                                     :: beta
    real(8)                                                     :: xmu
    !
    real(8),dimension(size(Ebands,1),size(Ebands,1))            :: Sigma_HF,Self_HF
    complex(8)                                                  :: Ak,Bk,Ck
    complex(8)                                                  :: Gk,Tk
    !
    complex(8),dimension(2,2)                                   :: GkNambu
    complex(8),dimension(2,2)                                   :: HkNambu,HlocNambu,AkNambu
    real(8),dimension(2)                                        :: NkNambu
    complex(8),dimension(2,2)                                   :: Evec
    real(8),dimension(2)                                        :: Eval,Coef
    real(8)                                                     :: spin_degeneracy
    !
    real(8),dimension(size(Ebands,1))                           :: H0,Hl
    real(8),dimension(size(Ebands,1))                           :: H0tmp,Hltmp
    real(8),dimension(size(Ebands,1))                           :: H0free,Hlfree
    real(8),dimension(size(Ebands,1))                           :: Ekin_,Eloc_
    !
    real(8),dimension(:),allocatable                            :: wm
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nso = size(Ebands,1)
    Lk  = size(Ebands,2)
    Liw = size(Sigma,3)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_dos: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_dos","Sigma")
    call assert_shape(Self,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_dos","Self")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,Liw))
    Self_HF  = dreal(Self(:,:,Liw))
    !
    !
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer()
    H0=0d0
    Hl=0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1+mpi_rank,Lk,mpi_size
       do is=1,Nso
          Ak = Ebands(is,ik)
          do i=1,Liw
             Gknambu=zero
             Gknambu(1,1) = (xi*wm(i) + xmu)  - Sigma(is,is,i)        - Ebands(is,ik) - Hloc(is)
             Gknambu(1,2) =                   - Self(is,is,i)
             Gknambu(2,1) =                   - Self(is,is,i)
             Gknambu(2,2) = (xi*wm(i) - xmu)  + conjg(Sigma(is,is,i)) + Ebands(is,ik) + Hloc(is)
             call inv(Gknambu)
             Gk  = Gknambu(1,1) !Gk(iw)
             !
             Gknambu=zero
             Gknambu(1,1) = (xi*wm(i) + xmu)  - Sigma_hf(is,is)  - Ebands(is,ik) - Hloc(is)
             Gknambu(1,2) =                   - Self_hf(is,is)
             Gknambu(2,1) =                   - Self_hf(is,is)
             Gknambu(2,2) = (xi*wm(i) - xmu)  + Sigma_hf(is,is)  + Ebands(is,ik) + Hloc(is)
             call inv(Gknambu)
             Tk = Gknambu(1,1) !G0k(iw)
             !
             Bk = Ak*(Gk - Tk)
             Ck = Hloc(is)*(Gk - Tk)
             !
             H0tmp(is) = H0tmp(is) + Dbands(is,ik)*Bk
             Hltmp(is) = Hltmp(is) + Dbands(is,ik)*Ck
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer
    spin_degeneracy=3d0-dble(Nspin) !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2.d0*spin_degeneracy
    Hl=Hl/beta*2.d0*spin_degeneracy
    !
    !
    H0free=0d0
    Hlfree=0d0
    do ik=1,Lk
       do is=1,Nso
          HkNambu=zero
          HlocNambu=zero
          HkNambu(1,1)   =  Ebands(is,ik)
          HkNambu(2,2)   = -Ebands(is,ik)
          HlocNambu(1,1) =  Hloc(is)
          HlocNambu(2,2) = -Hloc(is)
          Evec(1,1)      =  Ebands(is,ik) + Sigma_hf(is,is)
          Evec(1,2)      =                 + Self_hf(is,is)
          Evec(2,1)      =                 + Self_hf(is,is)
          Evec(2,2)      = -EBands(is,ik) - Sigma_hf(is,is)
          !
          call eigh(Evec,Eval)
          Coef = Evec(:,is)*conjg(Evec(:,is))
          NkNambu(is) = dot_product(Coef,fermi(Eval,beta))
          !
          AkNambu = matmul( diag(NkNambu) , HkNambu )
          H0free(is) = H0free(is) + Dbands(is,ik)*AkNambu(is,is) !Take only the 11 part
          !
          AkNambu = matmul( diag(NkNambu) , HlocNambu )
          Hlfree(is) = Hlfree(is) + Dbands(is,ik)*AkNambu(is,is)
          !
       enddo
    enddo
    H0free=spin_degeneracy*H0free
    Hlfree=spin_degeneracy*Hlfree
    !
    Ekin_=H0+H0free
    Eloc_=Hl+Hlfree
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_superc_Nso_dos
  !
  !
  !
  subroutine dmft_kinetic_energy_superc_Nso_lattice(Hk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)                                     :: Hk ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:,:)                                   :: Self ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    real(8),dimension(size(Hk,1)),optional                          :: Ekin,Eloc
    !
    integer                                                         :: Lk,Nlso,Liw,Nso,Nlat
    integer                                                         :: ik
    integer                                                         :: i,iorb,ilat,ispin,io,is
    integer                                                         :: j,jorb,jlat,jspin,jo,js
    !
    integer                                                         :: mpi_ierr
    integer                                                         :: mpi_rank
    integer                                                         :: mpi_size
    logical                                                         :: mpi_master
    !
    integer                                                         :: Norb,Nporb
    integer                                                         :: Nspin  
    real(8)                                                         :: beta
    real(8)                                                         :: xmu
    !
    complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Sigma_HF ![Nlat][Nso][Nso]
    complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Self_HF  ![Nlat][Nso][Nso]
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Ak,Bk,Ck,Hloc,Hloc_tmp
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Gk,Tk
    complex(8),dimension(2*size(Hk,1),2*size(Hk,2))                 :: Gknambu,HkNambu,HlocNambu,AkNambu
    real(8),dimension(2*size(Hk,1))                                 :: NkNambu
    complex(8),dimension(2*size(Hk,1),2*size(Hk,2))                 :: Evec
    real(8),dimension(2*size(Hk,1))                                 :: Eval,Coef
    real(8)                                                         :: spin_degeneracy
    !
    real(8),dimension(size(Hk,1))                                   :: H0,Hl
    real(8),dimension(size(Hk,1))                                   :: H0free,Hlfree
    real(8),dimension(size(Hk,1))                                   :: H0tmp,Hltmp
    real(8),dimension(size(Hk,1))                                   :: Ekin_,Eloc_

    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(Norb,"NORB")
    call get_ctrl_var(Nspin,"NSPIN")
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !Get generalized Lattice-Spin-Orbital index
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Nlat = size(Sigma,1)
    Nso  = size(Sigma,2)
    Liw  = size(Sigma,4)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_lattice: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_superc_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_lattice","Sigma")
    call assert_shape(Self,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_lattice","Self")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,:,Liw))![Nlat,Nso,Nso]
    Self_HF  = dreal(Self(:,:,:,Liw)) ![Nlat,Nso,Nso]
    !
    ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
    Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
    Hloc=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hloc(is,js)=Hloc_tmp(is,js) 
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
    if(size(Hloc,1)<16)then
       if(mpi_master)call print_hloc(Hloc)
    else
       if(mpi_master)call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    !Start the timer:
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer
    H0    = 0d0
    Hl    = 0d0
    H0tmp = 0d0
    Hltmp = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       do i=1+mpi_rank,Liw,mpi_size
          Gknambu=zero
          Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - &
               blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso)  - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        = - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        = - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + &
               conjg(blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso)) + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          Gk = Gknambu(1:Nlso,1:Nlso)
          !
          Gknambu=zero
          Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - &
               blocks_to_matrix(Sigma_HF,Nlat,Nso)  - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        =  - blocks_to_matrix(Self_HF,Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        =  - blocks_to_matrix(Self_HF,Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + &
               blocks_to_matrix(Sigma_HF,Nlat,Nso)  + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          Tk = Gknambu(1:Nlso,1:Nlso)
          !
          Bk = matmul(Ak, Gk - Tk)
          Ck = matmul(Hloc, Gk - Tk)
          do is=1,Nlso
             H0tmp(is) = H0tmp(is) + Bk(is,is)/dble(Lk)
             Hltmp(is) = Hltmp(is) + Ck(is,is)/dble(Lk)
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer
    spin_degeneracy=3d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    Hl = Hl/beta*2d0*spin_degeneracy
    !
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    H0free=0d0
    Hlfree=0d0
    HkNambu=zero
    HlocNambu=zero
    do ik=1,Lk
       HkNambu(1:Nlso,1:Nlso)                 =  Hk(:,:,ik)-Hloc
       HkNambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso)   = -Hk(:,:,ik)+Hloc
       HlocNambu(1:Nlso,1:Nlso)               =  Hloc
       HlocNambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = -Hloc
       Evec(1:Nlso,1:Nlso)                    =  Hk(:,:,ik) +  blocks_to_matrix(Sigma_HF,Nlat,Nso)
       Evec(1:Nlso,Nlso+1:2*Nlso)             =             +  blocks_to_matrix(Self_HF,Nlat,Nso)
       Evec(Nlso+1:2*Nlso,1:Nlso)             =             +  blocks_to_matrix(Self_HF,Nlat,Nso)
       Evec(Nlso+1:2*Nlso,Nlso+1:2*Nlso)      = -Hk(:,:,ik) -  blocks_to_matrix(Sigma_HF,Nlat,Nso)
       call eigh(Evec,Eval)
       NkNambu = fermi(Eval,beta)
       GkNambu = matmul(Evec,matmul(diag(NkNambu),conjg(transpose(Evec))))
       AkNambu = matmul(HkNambu  , GkNambu)
       do is=1,Nlso
          H0free(is) = H0free(is) + AkNambu(is,is)/dble(Lk) !Take only the 11 part
       enddo
       AkNambu = matmul(HlocNambu, GkNambu)
       do is=1,Nlso
          Hlfree(is) = Hlfree(is) + AkNambu(is,is)/dble(Lk)
       enddo
    enddo
    H0free=spin_degeneracy*H0free
    Hlfree=spin_degeneracy*Hlfree
    !
    Ekin_=H0+H0free
    Eloc_=Hl+Hlfree
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_,Nlat,Nso)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_superc_Nso_lattice








  !####################################################################
  !                    ADDITIONAL INTERFACES 
  !####################################################################

  !********************************************************************
  !                             NN FORM
  ! - Sigma & Self shape:       [Nspin][Nspin][Norb][Norb][L]
  ! - Sigma & Self shape: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !********************************************************************
  subroutine dmft_kinetic_energy_normal_NN(Hk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Hk,1)),optional  :: Ekin,Eloc
    real(8),dimension(size(Hk,1))           :: Ekin_,Eloc_
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN","Sigma")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_normal_NN","Hk")
    allocate(Sigma_(Nso,Nso,Lmats))
    Sigma_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_Nso_main(Hk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN

  subroutine dmft_kinetic_energy_normal_NN_dos(Ebands,Dbands,Hloc,Sigma,Ekin,Eloc)
    real(8),dimension(:,:)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1))                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:)                  :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                          :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Ebands,1)),optional       :: Ekin,Eloc
    real(8),dimension(size(Ebands,1))                :: Ekin_,Eloc_
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN","Sigma")
    allocate(Sigma_(Nso,Nso,Lmats))
    Sigma_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_Nso_dos(Ebands,Dbands,Hloc,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_dos

  subroutine dmft_kinetic_energy_normal_NN_lattice(Hk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1)),optional    :: Ekin,Eloc
    real(8),dimension(size(Hk,1))             :: Ekin_,Eloc_
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"dmft_kinetic_energy_normal_NN_lattice","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"dmft_kinetic_energy_normal_NN_lattice","Sigma")
    allocate(Sigma_(Nlat,Nso,Nso,Liw))
    Sigma_ = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Sigma_(:,is,js,:) = Sigma(:,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_Nso_lattice(Hk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_lattice

  subroutine dmft_kinetic_energy_superc_NN(Hk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)         :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Hk,1)),optional  :: Ekin,Eloc
    real(8),dimension(size(Hk,1))           :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Self")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_superc_NN","Hk")
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) =  Sigma(ispin,jspin,iorb,jorb,:)
                Self_(is,js,:)= Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_Nso_main(Hk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN

  subroutine dmft_kinetic_energy_superc_NN_dos(Ebands,Dbands,Hloc,Sigma,Self,Ekin,Eloc)
    real(8),dimension(:,:)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1))                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:)                  :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)                  :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                          :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Ebands,1)),optional       :: Ekin,Eloc
    real(8),dimension(size(Ebands,1))                :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Self")
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) =  Sigma(ispin,jspin,iorb,jorb,:)
                Self_(is,js,:)= Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_Nso_dos(Ebands,Dbands,Hloc,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_dos

  subroutine dmft_kinetic_energy_superc_NN_lattice(Hk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)         :: Self   ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Self_  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                   :: i,iorb,ilat,ispin,is
    integer                                   :: j,jorb,jlat,jspin,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1)),optional    :: Ekin,Eloc
    real(8),dimension(size(Hk,1))             :: Ekin_,Eloc_
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    !Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"kinetic_energy_lattice_superc_2","Hk") !implictly check Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Sigma")
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Self")
    allocate(Sigma_(Nlat,Nso,Nso,Liw))
    allocate(Self_(Nlat,Nso,Nso,Liw))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(:,is,js,:) =  Sigma(:,ispin,jspin,iorb,jorb,:)
                Self_(:,is,js,:)= Self(:,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_Nso_lattice(Hk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_lattice









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8)                    :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !+-----------------------------------------------------------------------------+!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix,Nlat,Nso) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !+-----------------------------------------------------------------------------+!
  function select_block_Nlso(ip,Matrix,Nlat,Nso) result(Vblock)
    integer                                 :: ip
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nso,Nso)           :: Vblock
    integer                                 :: i,j
    Vblock=zero
    i = 1+(ip-1)*Nso
    j =       ip*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ip,Matrix,Nlat,Nspin,Norb) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block_nnn



  !+------------------------------------------------------------------+
  !PURPOSE: Print the local part of the Non-interacting Hamiltonian
  !+------------------------------------------------------------------+
  subroutine print_Hloc_2(hloc,file)
    complex(8),dimension(:,:) :: hloc
    character(len=*),optional :: file
    integer                   :: Ni,Nj,iorb,jorb,unit
    unit=6;
    if(present(file))then
       write(unit,"(A)")"print_Hloc on file :"//reg(file)
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Ni = size(hloc,1)
    Nj = size(hloc,2)
    if(present(file))then
       do iorb=1,Ni
          write(unit,"(10000F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       do iorb=1,Ni
          write(unit,"(10000F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       close(unit)
    else
       do iorb=1,Ni
          write(unit,"(10000(A1,F7.3,A1,F7.3,A1,2x))")&
               ('(',dreal(Hloc(iorb,jorb)),',',dimag(Hloc(iorb,jorb)),')',jorb =1,Nj)
       enddo
       write(unit,*)""
    endif
  end subroutine print_Hloc_2
  !
  subroutine print_Hloc_4(hloc,file)
    complex(8),dimension(:,:,:,:) :: hloc ![Nspin][Nspin][Norb][Norb]
    character(len=*),optional     :: file
    integer                       :: Nspin,Norb
    integer                       :: iorb,jorb,ispin,jspin,unit
    unit=6;
    if(present(file))then
       write(unit,"(A)")"print_Hloc on file :"//reg(file)
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Nspin = size(hloc,1)
    Norb  = size(hloc,3)
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"print_Hloc","Hloc")
    if(present(file))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(10000F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(10000F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    else
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(10000(A1,F7.3,A1,F7.3,A1,2x))")&
                  (&
                  (&
                  '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
                  jorb =1,Norb),&
                  jspin=1,Nspin)
          enddo
       enddo
    endif
  end subroutine print_Hloc_4





  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_info()
    !to be removed
  end subroutine write_kinetic_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_value(Ekin,Eloc,Nlat,Nso)
    real(8),dimension(:)               :: Ekin
    real(8),dimension(size(Ekin))      :: Eloc
    real(8),dimension(:,:),allocatable :: Ekin_,Eloc_
    integer,optional                   :: Nlat,Nso
    integer                            :: Nlso
    integer                            :: i,iso,ilat,unit
    if(.not.present(Nlat))then
       !
       Nlso = size(Ekin)
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nlso),&
            (str(2+Nlso+i)//"<Eloc"//str(i)//">",i=1,Nlso)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.dat")
       write(unit,"(90F15.9)")sum(Ekin),sum(Eloc),(Ekin(i),i=1,Nlso),(Eloc(i),i=1,Nlso)
       close(unit)
       !
    else
       !
       if(.not.present(Nso))stop "ERROR write_kinetic_value: Nlat present but Nso not present."
       !
       Nlso = size(Ekin)
       if(Nlso /= Nlat*Nso)stop "Error write_kinetic_value: Nlso != Nlat*Nso" 
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nso),&
            (str(2+Nso+i)//"<Eloc"//str(i)//">",i=1,Nso)
       close(unit)
       !
       !
       allocate(Ekin_(Nlat,Nso))
       allocate(Eloc_(Nlat,Nso))
       do ilat=1,Nlat
          do iso=1,Nso
             i = iso + (ilat-1)*Nso
             Ekin_(ilat,iso) = Ekin(i)
             Eloc_(ilat,iso) = Eloc(i)
          enddo
       enddo
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.dat")       
       write(unit,"(90F15.9)")sum(Ekin_)/Nlat,sum(Eloc_)/Nlat
       do ilat=1,Nlat
          write(unit,"(100000F15.9)")sum(Ekin_(ilat,:)),sum(Eloc_(ilat,:)),&
               (Ekin_(ilat,i),i=1,Nso),&
               (Eloc_(ilat,i),i=1,Nso)
       enddo
       close(unit)
    endif
  end subroutine write_kinetic_value









#ifdef _MPI
  function MPI_Get_size(comm) result(size)
    integer :: comm
    integer :: size,ierr
    call MPI_Comm_size(comm,size,ierr)
  end function MPI_Get_size



  function MPI_Get_rank(comm) result(rank)
    integer :: comm
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
  end function MPI_Get_rank



  function MPI_Get_master(comm) result(master)
    integer :: comm
    logical :: master
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    master=.false.
    if(rank==0)master=.true.
  end function MPI_Get_master
#endif



end MODULE DMFT_EKIN


