MODULE DMFT_EKIN
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS,   only:free_unit,reg,txtfy
  USE SF_ARRAYS,    only:arange
  USE SF_TIMER
  USE SF_LINALG,    only:inv,eye,eigh,diag
  USE SF_MISC,      only:assert_shape
  USE SF_SPECIAL,   only:fermi
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private


  interface dmft_kinetic_energy
     module procedure :: dmft_kinetic_energy_normal_main
     module procedure :: dmft_kinetic_energy_normal_main_lattice
     module procedure :: dmft_kinetic_energy_superc_main
     module procedure :: dmft_kinetic_energy_superc_main_lattice
#ifdef _MPI
     module procedure :: dmft_kinetic_energy_normal_main_mpi
     module procedure :: dmft_kinetic_energy_normal_main_lattice_mpi
     module procedure :: dmft_kinetic_energy_superc_main_mpi
     module procedure :: dmft_kinetic_energy_superc_main_lattice_mpi
#endif
     !ADDITIONAL INTERFACES:
     !1-Band ([Nlat])[L] for Sigma & Self
     module procedure :: dmft_kinetic_energy_normal_1band
     module procedure :: dmft_kinetic_energy_normal_1band_lattice
     module procedure :: dmft_kinetic_energy_superc_1band
     module procedure :: dmft_kinetic_energy_superc_1band_lattice
#ifdef _MPI
     module procedure :: dmft_kinetic_energy_normal_1band_mpi
     module procedure :: dmft_kinetic_energy_normal_1band_lattice_mpi
     module procedure :: dmft_kinetic_energy_superc_1band_mpi
     module procedure :: dmft_kinetic_energy_superc_1band_lattice_mpi
#endif
     !NN-Shape ([Nlat])[Nspin][Nspin][Norb][Norb][L] for Sigma & Self
     module procedure :: dmft_kinetic_energy_normal_NN
     module procedure :: dmft_kinetic_energy_normal_NN_lattice
     module procedure :: dmft_kinetic_energy_superc_NN
     module procedure :: dmft_kinetic_energy_superc_NN_lattice
#ifdef _MPI
     module procedure :: dmft_kinetic_energy_normal_NN_mpi
     module procedure :: dmft_kinetic_energy_normal_NN_lattice_mpi
     module procedure :: dmft_kinetic_energy_superc_NN_mpi
     module procedure :: dmft_kinetic_energy_superc_NN_lattice_mpi
#endif
  end interface dmft_kinetic_energy


  interface print_Hloc
     module procedure print_Hloc_2
     module procedure print_Hloc_4
  end interface print_Hloc



  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block


  real(8),dimension(:),allocatable        :: wm



  !PUBLIC in DMFT
  public :: dmft_kinetic_energy

contains 










  !----------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy of the lattice model. (NORMAL)
  ! input: 
  ! 1. the Hamiltonian matrix H(k)
  ! 2. DMFT normal self-energy Sigma.
  !
  ! The main routines accept self-energy as:
  ! - Sigma: [Nspin*Norb][Nspin*Norb][L]
  ! - Sigma: [Nlat][Nspin*Norb][Nspin*Norb][L]
  !
  ! Additional interfaces with different shapes of the Sigma are
  ! appended to this module.
  !----------------------------------------------------------------------------------------
  function dmft_kinetic_energy_normal_main(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:,:,:)                 :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))               :: Wtk   ![Lk]
    complex(8),dimension(:,:,:)                 :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    !
    integer                                     :: Lk,Nso,Liw
    integer                                     :: i,ik
    !
    integer                                     :: Norb,Nporb
    integer                                     :: Nspin  
    real(8)                                     :: beta
    real(8)                                     :: xmu
    !
    real(8),dimension(size(Hk,1),size(Hk,1))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Ak,Bk,Ck,Dk,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Gk,Tk
    real(8)                                     :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                     :: H0,Hl
    real(8)                                     :: Ekin,Eloc
    real(8)                                     :: Eout(2)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_main: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_normal_main","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_main","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    !Get the local Hamiltonian
    Hloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
    if(size(Hloc,1)<16)then
       call print_hloc(Hloc)
    else
       call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    write(*,"(A)") "Kinetic energy computation"
    call start_timer()
    H0=0d0
    Hl=0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak = Hk(:,:,ik) - Hloc(:,:)
       Bk =-Hk(:,:,ik) - Sigma_HF(:,:)
       do i=1,Liw
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
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,Nso)
          Hl = Hl + Wtk(ik)*trace_matrix(Dk,Nso)
       enddo
       call eta(ik,Lk)
    enddo
    call stop_timer()
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
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       Ck= matmul(Ak,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Dk= matmul(Hloc,(-Hk(:,:,ik)-Sigma_HF(:,:)))
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nso)
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin=H0+Tail0+Tail1
    Eloc=Hl+Lail0+Lail1
    Eout = [Ekin,Eloc]
    !
    call write_kinetic_info()
    call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_normal_main


  function dmft_kinetic_energy_normal_main_lattice(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))                                   :: Wtk       ! [Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma     ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    !aux
    integer                                                         :: Lk,Nlso,Liw,Nlat,Nso
    integer                                                         :: ik
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
    real(8)                                                         :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                                         :: H0,Hl
    real(8)                                                         :: H0k,Hlk
    real(8)                                                         :: H0ktmp,Hlktmp 
    real(8)                                                         :: Ekin_lattice,Eloc_lattice
    real(8)                                                         :: Eout(2)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_main_lattice: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_normal_main_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_normal_main_lattice","Sigma")
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
       call print_hloc(Hloc)
    else
       call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    !Start the timer:
    write(*,"(A)") "Kinetic energy computation"
    call start_timer
    H0           = 0d0
    Hl           = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    =  Hk(:,:,ik) - Hloc(:,:)
       Bk    = -Hk(:,:,ik) - blocks_to_matrix(Sigma_HF(:,:,:),Nlat,Nso) !Sigma_HF [Nlat,Nso,Nso]--> [Nlso,Nslo]
       do i=1,Liw
          Gk = (xi*wm(i)+xmu)*eye(Nlso) - blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso) - Hk(:,:,ik) !Sigma [Nlat,Nso,Nso,*]--> [Nlso,Nslo]
          call inv(Gk(:,:))
          Tk = eye(Nlso)/(xi*wm(i)) - Bk/(xi*wm(i))**2
          Ck = matmul(Ak,Gk(:,:) - Tk)
          Dk = matmul(Hloc,Gk(:,:) - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Ck,Nlso)
          Hl = Hl + Wtk(ik)*trace_matrix(Dk,Nlso)
       enddo
       call eta(ik,Lk)
    enddo
    call stop_timer
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
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nlso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nlso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nlso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nlso)
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin_lattice=H0+Tail0+Tail1
    Eloc_lattice=Hl+Lail0+Lail1
    Ekin_lattice=Ekin_lattice/dble(Nlat)
    Eloc_lattice=Eloc_lattice/dble(Nlat)
    Eout = [Ekin_lattice,Eloc_lattice]
    !
    call write_kinetic_info()
    call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_normal_main_lattice



#ifdef _MPI
  function dmft_kinetic_energy_normal_main_mpi(MpiComm,Hk,Wtk,Sigma) result(Eout)
    integer                                     :: MpiComm
    complex(8),dimension(:,:,:)                 :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))               :: Wtk   ![Lk]
    complex(8),dimension(:,:,:)                 :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    !
    integer                                     :: Lk,Nso,Liw
    integer                                     :: i,ik
    !
    integer                                     :: mpi_ierr
    integer                                     :: mpi_rank
    integer                                     :: mpi_size
    logical                                     :: mpi_master
    !
    integer                                     :: Norb,Nporb
    integer                                     :: Nspin  
    real(8)                                     :: beta
    real(8)                                     :: xmu
    !
    real(8),dimension(size(Hk,1),size(Hk,1))    :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Ak,Bk,Ck,Dk,Hloc
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Gk,Tk
    real(8)                                     :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                     :: H0,Hl
    real(8)                                     :: H0tmp,Hltmp
    real(8)                                     :: Ekin,Eloc
    real(8)                                     :: Eout(2)
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_main_mpi: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_normal_main_mpi","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_main_mpi","Sigma")
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
          H0tmp = H0tmp + Wtk(ik)*trace_matrix(Ck,Nso)
          Hltmp = Hltmp + Wtk(ik)*trace_matrix(Dk,Nso)
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
    call MPI_ALLREDUCE(H0tmp, H0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
    call MPI_ALLREDUCE(Hltmp, Hl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
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
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nso)
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin=H0+Tail0+Tail1
    Eloc=Hl+Lail0+Lail1
    Eout = [Ekin,Eloc]
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_normal_main_mpi


  function dmft_kinetic_energy_normal_main_lattice_mpi(MpiComm,Hk,Wtk,Sigma) result(Eout)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))                                   :: Wtk       ! [Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma     ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    !aux
    integer                                                         :: Lk,Nlso,Liw,Nlat,Nso
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
    complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Sigma_HF
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Ak,Bk,Ck,Dk,Hloc,Hloc_tmp
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Tk
    complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Gk
    real(8)                                                         :: Tail0,Tail1,Lail0,Lail1,spin_degeneracy
    !
    real(8)                                                         :: H0,Hl
    real(8)                                                         :: H0tmp,Hltmp
    real(8)                                                         :: Ekin_lattice,Eloc_lattice
    real(8)                                                         :: Eout(2)
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_lattice_normal_main: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_lattice_normal_main","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_lattice_normal_main","Sigma")
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
          H0tmp = H0tmp + Wtk(ik)*trace_matrix(Ck,Nlso)
          Hltmp = Hltmp + Wtk(ik)*trace_matrix(Dk,Nlso)
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
    call MPI_ALLREDUCE(H0tmp, H0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
    call MPI_ALLREDUCE(Hltmp, Hl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
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
       Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,Nlso)
       Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,Nlso)
       Lail0 = Lail0 + 0.5d0*Wtk(ik)*trace_matrix(Hloc(:,:),Nlso)
       Lail1 = Lail1 + 0.25d0*Wtk(ik)*trace_matrix(Dk,Nlso)
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin_lattice=H0+Tail0+Tail1
    Eloc_lattice=Hl+Lail0+Lail1
    Ekin_lattice=Ekin_lattice/dble(Nlat)
    Eloc_lattice=Eloc_lattice/dble(Nlat)
    Eout = [Ekin_lattice,Eloc_lattice]
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_normal_main_lattice_mpi
#endif





  !----------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy of the lattice model. (SUPERC)
  ! input: 
  ! 1. the Hamiltonian matrix H(k)
  ! 2. DMFT normal&anomalous self-energy Sigma & Self.
  !
  ! The main routines accept self-energy as:
  ! - Sigma & Self: [Nspin*Norb][Nspin*Norb][L]
  ! - Sigma & Self: [Nlat][Nspin*Norb][Nspin*Norb][L]
  !
  ! Additional interfaces with different shapes of the Sigma are
  ! appended to this module.
  !----------------------------------------------------------------------------------------
  function dmft_kinetic_energy_superc_main(Hk,Wtk,Sigma,Self) result(Eout)
    complex(8),dimension(:,:,:)                     :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))                   :: Wtk   ![Lk]
    complex(8),dimension(:,:,:)                     :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:)                     :: Self ![Nspin*Norb][Nspin*Norb][L]
    !
    integer                                         :: Lk,Nso,Liw
    integer                                         :: i,is,ik
    !
    integer                                     :: Norb,Nporb
    integer                                     :: Nspin  
    real(8)                                     :: beta
    real(8)                                     :: xmu
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
    real(8)                                         :: H0,Hl,H0free,Hlfree
    real(8)                                         :: Ekin,Eloc
    real(8)                                         :: Eout(2)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_main: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_superc_main","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_main","Sigma")
    call assert_shape(Self,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_main","Self")
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
       call print_hloc(Hloc)
    else
       call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    write(*,"(A)") "Kinetic energy computation"
    call start_timer()
    H0=0d0
    Hl=0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak= Hk(:,:,ik)-Hloc
       do i=1,Liw
          Gknambu=zero
          Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)  - Sigma(:,:,i)        - Hk(:,:,ik)
          Gknambu(1:Nso,Nso+1:2*Nso)       =                            - Self(:,:,i)
          Gknambu(Nso+1:2*Nso,1:Nso)       =                            - Self(:,:,i)
          Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)  + conjg(Sigma(:,:,i)) + Hk(:,:,ik)
          call inv(Gknambu)
          Gk  = Gknambu(1:Nso,1:Nso) !Gk(iw)
          !
          Gknambu=zero
          Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)  - Sigma_hf  - Hk(:,:,ik)
          Gknambu(1:Nso,Nso+1:2*Nso)       =                            - Self_hf
          Gknambu(Nso+1:2*Nso,1:Nso)       =                            - Self_hf
          Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)  + Sigma_hf  + Hk(:,:,ik)
          call inv(Gknambu)
          Tk = Gknambu(1:Nso,1:Nso) !G0k(iw)
          !
          Bk = matmul(Ak, Gk - Tk) !
          H0 = H0 + Wtk(ik)*trace_matrix(Bk,Nso)
          Ck = matmul(Hloc,Gk - Tk)
          Hl = Hl + Wtk(ik)*trace_matrix(Ck,Nso)
       enddo
       call eta(ik,Lk)
    enddo
    call stop_timer
    spin_degeneracy=3d0-dble(Nspin) !2 if Nspin=1, 1 if Nspin=2
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
       H0free = H0free + Wtk(ik)*trace_matrix(AkNambu(:Nso,:Nso),Nso) !Take only the 11 part
       AkNambu = matmul( diag(NkNambu) , HlocNambu )
       Hlfree = Hlfree + Wtk(ik)*trace_matrix(AkNambu(:Nso,:Nso),Nso)
    enddo
    H0free=spin_degeneracy*H0free
    Hlfree=spin_degeneracy*Hlfree
    !
    Ekin=H0+H0free
    Eloc=Hl+Hlfree
    Eout = [Ekin,Eloc]
    !
    call write_kinetic_info()
    call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_superc_main


  function dmft_kinetic_energy_superc_main_lattice(Hk,Wtk,Sigma,Self) result(Eout)
    complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))                                   :: Wtk       ! [Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:,:)                                   :: Self   ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    !
    integer                                                         :: Lk,Nlso,Liw,Nso,Nlat
    integer                                                         :: ik
    integer                                                         :: i,iorb,ilat,ispin,io,is
    integer                                                         :: j,jorb,jlat,jspin,jo,js
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
    real(8)                                                         :: H0,Hl
    real(8)                                                         :: H0free,Hlfree
    real(8)                                                         :: Ekin_lattice,Eloc_lattice
    real(8)                                                         :: Eout(2)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_main_lattice: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_superc_main_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_main_lattice","Sigma")
    call assert_shape(Self,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_main_lattice","Self")
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
       call print_hloc(Hloc)
    else
       call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    endif
    !
    !Start the timer:
    write(*,"(A)") "Kinetic energy computation"
    call start_timer
    H0              = 0d0
    Hl              = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       do i=1,Liw
          Gknambu=zero
          Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - &
               blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso) - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        =                            - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        =                            - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + &
               conjg(blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso)) + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          Gk = Gknambu(1:Nlso,1:Nlso)
          !
          Gknambu=zero
          Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - blocks_to_matrix(Sigma_HF,Nlat,Nso)  - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        =                            - blocks_to_matrix(Self_HF,Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        =                            - blocks_to_matrix(Self_HF,Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + blocks_to_matrix(Sigma_HF,Nlat,Nso)  + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          Tk = Gknambu(1:Nlso,1:Nlso)
          !
          Bk = matmul(Ak, Gk - Tk)
          H0 = H0 + Wtk(ik)*trace_matrix(Bk,Nlso)
          Ck = matmul(Hloc, Gk - Tk)
          Hl = Hl + Wtk(ik)*trace_matrix(Ck,Nlso)
       enddo
       call eta(ik,Lk)
    enddo
    call stop_timer
    spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy!;print*,"Ekin_=",H0/Nlat
    Hl = Hl/beta*2d0*spin_degeneracy!;print*,"Eloc_=",Hl/Nlat
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
       H0free  = H0free + Wtk(ik)*trace_matrix( AkNambu(:Nlso,:Nlso) , Nlso) !take only the 11 part
       AkNambu = matmul(HlocNambu, GkNambu)
       Hlfree  = Hlfree + Wtk(ik)*trace_matrix( AkNambu(:Nlso,:Nlso) , Nlso)
    enddo
    H0free=spin_degeneracy*H0free!;print*,"Efree=",H0free/Nlat
    Hlfree=spin_degeneracy*Hlfree!;print*,"Efree_loc=",Hlfree/Nlat
    !
    Ekin_lattice=H0+H0free
    Eloc_lattice=Hl+Hlfree
    Ekin_lattice=Ekin_lattice/dble(Nlat)
    Eloc_lattice=Eloc_lattice/dble(Nlat)
    Eout = [Ekin_lattice,Eloc_lattice]
    !
    call write_kinetic_info()
    call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_superc_main_lattice











#ifdef _MPI
  function dmft_kinetic_energy_superc_main_mpi(MpiComm,Hk,Wtk,Sigma,Self) result(Eout)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:)                     :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))                   :: Wtk   ![Lk]
    complex(8),dimension(:,:,:)                     :: Sigma ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:)                     :: Self ![Nspin*Norb][Nspin*Norb][L]
    !
    integer                                         :: Lk,Nso,Liw
    integer                                         :: i,is,ik
    !
    integer                                         :: mpi_ierr
    integer                                         :: mpi_rank
    integer                                         :: mpi_size
    logical                                         :: mpi_master
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
    real(8)                                         :: H0,Hl,H0free,Hlfree
    real(8)                                         :: H0tmp,Hltmp
    real(8)                                         :: Ekin,Eloc
    real(8)                                         :: Eout(2)
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_main_mpi: Nso != Norb*Nspin [from Hk]"
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_superc_main_mpi","Hk")
    call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_main_mpi","Sigma")
    call assert_shape(Self,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_main_mpi","Self")
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
          Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)  - Sigma(:,:,i)        - Hk(:,:,ik)
          Gknambu(1:Nso,Nso+1:2*Nso)       =                            - Self(:,:,i)
          Gknambu(Nso+1:2*Nso,1:Nso)       =                            - Self(:,:,i)
          Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)  + conjg(Sigma(:,:,i)) + Hk(:,:,ik)
          call inv(Gknambu)
          Gk  = Gknambu(1:Nso,1:Nso) !Gk(iw)
          !
          Gknambu=zero
          Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)  - Sigma_hf  - Hk(:,:,ik)
          Gknambu(1:Nso,Nso+1:2*Nso)       =                            - Self_hf
          Gknambu(Nso+1:2*Nso,1:Nso)       =                            - Self_hf
          Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)  + Sigma_hf  + Hk(:,:,ik)
          call inv(Gknambu)
          Tk = Gknambu(1:Nso,1:Nso) !G0k(iw)
          !
          Bk = matmul(Ak, Gk - Tk) !
          H0tmp = H0tmp + Wtk(ik)*trace_matrix(Bk,Nso)
          Ck = matmul(Hloc,Gk - Tk)
          Hltmp = Hltmp + Wtk(ik)*trace_matrix(Ck,Nso)
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
    call MPI_ALLREDUCE(H0tmp, H0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
    call MPI_ALLREDUCE(Hltmp, Hl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
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
       H0free = H0free + Wtk(ik)*trace_matrix(AkNambu(:Nso,:Nso),Nso) !Take only the 11 part
       AkNambu = matmul( diag(NkNambu) , HlocNambu )
       Hlfree = Hlfree + Wtk(ik)*trace_matrix(AkNambu(:Nso,:Nso),Nso)
    enddo
    H0free=spin_degeneracy*H0free
    Hlfree=spin_degeneracy*Hlfree
    !
    Ekin=H0+H0free
    Eloc=Hl+Hlfree
    Eout = [Ekin,Eloc]
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_superc_main_mpi


  function dmft_kinetic_energy_superc_main_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Self) result(Eout)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3))                                   :: Wtk       ! [Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:,:)                                   :: Self   ! [Nlat][Nspin*Norb][Nspin*Norb][L]
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
    real(8)                                                         :: H0,Hl
    real(8)                                                         :: H0free,Hlfree
    real(8)                                                         :: H0tmp,Hltmp
    real(8)                                                         :: Ekin_lattice,Eloc_lattice
    real(8)                                                         :: Eout(2)
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
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
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_main_lattice: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_superc_main_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_main_lattice","Sigma")
    call assert_shape(Self,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_main_lattice","Self")
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
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        =                            - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        =                            - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + &
               conjg(blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso)) + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          Gk = Gknambu(1:Nlso,1:Nlso)
          !
          Gknambu=zero
          Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - blocks_to_matrix(Sigma_HF,Nlat,Nso)  - Hk(:,:,ik)
          Gknambu(1:Nlso,Nlso+1:2*Nlso)        =                            - blocks_to_matrix(Self_HF,Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,1:Nlso)        =                            - blocks_to_matrix(Self_HF,Nlat,Nso)
          Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + blocks_to_matrix(Sigma_HF,Nlat,Nso)  + Hk(:,:,ik)
          call inv(Gknambu(:,:))
          Tk = Gknambu(1:Nlso,1:Nlso)
          !
          Bk = matmul(Ak, Gk - Tk)
          H0tmp = H0tmp + Wtk(ik)*trace_matrix(Bk,Nlso)
          Ck = matmul(Hloc, Gk - Tk)
          Hltmp = Hltmp + Wtk(ik)*trace_matrix(Ck,Nlso)
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
    call MPI_ALLREDUCE(H0tmp, H0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
    call MPI_ALLREDUCE(Hltmp, Hl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MpiComm, mpi_ierr)
    if(mpi_master)call stop_timer
    spin_degeneracy=3d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy!;print*,"Ekin_=",H0/Nlat
    Hl = Hl/beta*2d0*spin_degeneracy!;print*,"Eloc_=",Hl/Nlat
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
       H0free  = H0free + Wtk(ik)*trace_matrix( AkNambu(:Nlso,:Nlso) , Nlso) !take only the 11 part
       AkNambu = matmul(HlocNambu, GkNambu)
       Hlfree  = Hlfree + Wtk(ik)*trace_matrix( AkNambu(:Nlso,:Nlso) , Nlso)
    enddo
    H0free=spin_degeneracy*H0free;print*,"Efree=",H0free/Nlat
    Hlfree=spin_degeneracy*Hlfree;print*,"Efree_loc=",Hlfree/Nlat
    !
    Ekin_lattice=H0+H0free
    Eloc_lattice=Hl+Hlfree
    Ekin_lattice=Ekin_lattice/Nlat
    Eloc_lattice=Eloc_lattice/Nlat
    Eout = [Ekin_lattice,Eloc_lattice]
    !
    if(mpi_master)call write_kinetic_info()
    if(mpi_master)call write_kinetic_value(Eout)
    !
    deallocate(wm)
  end function dmft_kinetic_energy_superc_main_lattice_mpi
#endif














  !####################################################################
  !                    ADDITIONAL INTERFACES 
  !####################################################################
  !********************************************************************
  !                           1 BAND
  ! - Sigma & Self shape:       [L] (no Nspin[=1], no Norb[=1])
  ! - Sigma & Self shape: [Nlat][L] (no Nspin[=1], no Norb[=1])
  !********************************************************************
  function dmft_kinetic_energy_normal_1band(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:)                   :: Hk
    real(8),dimension(size(Hk))               :: Wtk
    complex(8),dimension(:)                   :: Sigma
    complex(8),dimension(1,1,size(Sigma))     :: Sigma_
    complex(8),dimension(1,1,size(Hk))        :: Hk_
    real(8),dimension(2)                      :: Eout
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    Eout = dmft_kinetic_energy_normal_main(Hk_,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_1band
#ifdef _MPI
  function dmft_kinetic_energy_normal_1band_mpi(MpiComm,Hk,Wtk,Sigma) result(Eout)
    integer                                   :: MpiComm
    complex(8),dimension(:)                   :: Hk
    real(8),dimension(size(Hk))               :: Wtk
    complex(8),dimension(:)                   :: Sigma
    complex(8),dimension(1,1,size(Sigma))     :: Sigma_
    complex(8),dimension(1,1,size(Hk))        :: Hk_
    real(8),dimension(2)                      :: Eout
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    Eout = dmft_kinetic_energy_normal_main_mpi(MpiComm,Hk_,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_1band_mpi
#endif



  function dmft_kinetic_energy_normal_1band_lattice(Hk,Wtk,Sigma) result(Ekin_lattice)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*1][Nlat*1][Nk] !Nspin*Norb=1
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:)                 :: Sigma  ! [Nlat][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][1][1][L]
    integer                                   :: ilat
    integer                                   :: Nlat,Nlso,Nso,Lk
    real(8)                                   :: Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_normal_1band_lattice","Hk") !implictly test Nlat*1*=Nlso
    allocate(Sigma_(Nlat,1,1,size(Sigma,2)))
    Sigma_(:,1,1,:) = Sigma(:,:)
    Ekin_lattice = dmft_kinetic_energy_normal_main_lattice(Hk,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_1band_lattice
#ifdef _MPI
  function dmft_kinetic_energy_normal_1band_lattice_mpi(MpiComm,Hk,Wtk,Sigma) result(Ekin_lattice)
    integer                                   :: MpiComm
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*1][Nlat*1][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:)                 :: Sigma  ! [Nlat][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][1][1][L]
    integer                                   :: ilat
    integer                                   :: Nlat,Nlso,Nso,Lk
    real(8)                                   :: Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_normal_1band_lattice_mpi","Hk") !implictly test Nlat*1*=Nlso
    allocate(Sigma_(Nlat,1,1,size(Sigma,2)))
    Sigma_(:,1,1,:) = Sigma(:,:)
    Ekin_lattice = dmft_kinetic_energy_normal_main_lattice_mpi(MpiComm,Hk,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_1band_lattice_mpi
#endif




  function dmft_kinetic_energy_superc_1band(Hk,Wtk,Sigma,Self) result(Eout)
    real(8),dimension(:)                              :: Hk    ![Lk]
    real(8),dimension(size(Hk))                       :: Wtk   ![Lk]
    complex(8),dimension(:)                           :: Sigma ![L]
    complex(8),dimension(size(Sigma))                 :: Self  ![L]
    complex(8),dimension(1,1,size(Sigma))             :: Sigma_
    complex(8),dimension(1,1,size(Sigma))             :: Self_
    complex(8),dimension(1,1,size(Hk))                :: Hk_
    real(8)                                           :: Eout(2)
    Sigma_(1,1,:)  = Sigma(:)
    Self_(1,1,:)   = Self(:)
    Hk_(1,1,:)     = Hk
    Eout = dmft_kinetic_energy_superc_main(Hk_,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_1band
#ifdef _MPI
  function dmft_kinetic_energy_superc_1band_mpi(MpiComm,Hk,Wtk,Sigma,Self) result(Eout)
    integer                                           :: MpiComm
    real(8),dimension(:)                              :: Hk
    real(8),dimension(size(Hk))                       :: Wtk
    complex(8),dimension(:)                           :: Sigma
    complex(8),dimension(size(Sigma))                 :: Self
    complex(8),dimension(1,1,size(Sigma))             :: Sigma_
    complex(8),dimension(1,1,size(Sigma))             :: Self_
    complex(8),dimension(1,1,size(Hk))                :: Hk_
    real(8)                                           :: Eout(2)
    Sigma_(1,1,:)  = Sigma(:)
    Self_(1,1,:)   = Self(:)
    Hk_(1,1,:)     = Hk
    Eout = dmft_kinetic_energy_superc_main_mpi(MpiComm,Hk_,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_1band_mpi
#endif



  function dmft_kinetic_energy_superc_1band_lattice(Hk,Wtk,Sigma,Self) result(Ekin_lattice)
    complex(8),dimension(:,:,:)                       :: Hk     ! [Nlat*1][Nlat*1][Nk]
    real(8),dimension(size(Hk,3))                     :: wtk    ! [Nk]
    complex(8),dimension(:,:)                         :: Sigma  ! [Nlat][L]
    complex(8),dimension(size(Sigma,1),size(Sigma,2)) :: Self   ! [Nlat][L]  
    complex(8),dimension(:,:,:,:),allocatable         :: Sigma_ 
    complex(8),dimension(:,:,:,:),allocatable         :: Self_
    integer                                           :: i,iorb,ilat,ispin,io,is
    integer                                           :: j,jorb,jlat,jspin,jo,js
    integer                                           :: Nlat,Nlso,Nso,Lk,Liw
    real(8)                                           :: Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Liw  = size(Sigma,2)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_superc_1band_lattice","Hk") !implictly test Nlat*1*1=Nlso
    allocate(Sigma_(Nlat,1,1,Liw))
    allocate(Self_(Nlat,1,1,Liw))
    Sigma_(:,1,1,:)  = Sigma
    Self_(:,1,1,:) = Self
    Ekin_lattice = dmft_kinetic_energy_superc_main_lattice(Hk,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_1band_lattice
#ifdef _MPI
  function dmft_kinetic_energy_superc_1band_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Self) result(Ekin_lattice)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:)                       :: Hk     ! [Nlat*1][Nlat*1][Nk]
    real(8),dimension(size(Hk,3))                     :: wtk    ! [Nk]
    complex(8),dimension(:,:)                         :: Sigma  ! [Nlat][L]
    complex(8),dimension(size(Sigma,1),size(Sigma,2)) :: Self ! [Nlat][L]  
    complex(8),dimension(:,:,:,:),allocatable         :: Sigma_ 
    complex(8),dimension(:,:,:,:),allocatable         :: Self_
    integer                                           :: i,iorb,ilat,ispin,io,is
    integer                                           :: j,jorb,jlat,jspin,jo,js
    integer                                           :: Nlat,Nlso,Nso,Lk,Liw
    real(8)                                           :: Ekin_lattice(2)
    Nlat = size(Sigma,1)
    Liw  = size(Sigma,2)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_superc_1band_lattice_mpi","Hk") !implictly test Nlat*1*1=Nlso
    allocate(Sigma_(Nlat,1,1,Liw))
    allocate(Self_(Nlat,1,1,Liw))
    Sigma_(:,1,1,:)  = Sigma
    Self_(:,1,1,:) = Self
    Ekin_lattice = dmft_kinetic_energy_superc_main_lattice_mpi(MpiComm,Hk,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_1band_lattice_mpi
#endif




  
  !********************************************************************
  !                             NN FORM
  ! - Sigma & Self shape:       [Nspin][Nspin][Norb][Norb][L]
  ! - Sigma & Self shape: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !********************************************************************
  function dmft_kinetic_energy_normal_NN(Hk,Wtk,Sigma) result(Eout)
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(2)                    :: Eout
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
    Eout = dmft_kinetic_energy_normal_main(Hk,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_NN
#ifdef _MPI
  function dmft_kinetic_energy_normal_NN_mpi(MpiComm,Hk,Wtk,Sigma) result(Eout)
    integer                                 :: MpiComm
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(2)                    :: Eout
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN_mpi","Sigma")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_normal_NN_mpi","Hk")
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
    Eout = dmft_kinetic_energy_normal_main_mpi(MpiComm,Hk,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_NN_mpi
#endif







  function dmft_kinetic_energy_normal_NN_lattice(Hk,Wtk,Sigma) result(Ekin_lattice)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8)                                   :: Ekin_lattice(2)
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
    Ekin_lattice = dmft_kinetic_energy_normal_main_lattice(Hk,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_NN_lattice
#ifdef _MPI
  function dmft_kinetic_energy_normal_NN_lattice_mpi(MpiComm,Hk,Wtk,Sigma) result(Ekin_lattice)
    integer                                   :: MpiComm
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8)                                   :: Ekin_lattice(2)
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"dmft_kinetic_energy_normal_NN_lattice_mpi","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"dmft_kinetic_energy_normal_NN_lattice_mpi","Sigma")
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
    Ekin_lattice = dmft_kinetic_energy_normal_main_lattice_mpi(MpiComm,Hk,Wtk,Sigma_)
  end function dmft_kinetic_energy_normal_NN_lattice_mpi
#endif





  function dmft_kinetic_energy_superc_NN(Hk,Wtk,Sigma,Self) result(Eout)
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)         :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8)                                 :: Eout(2)
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
    Eout = dmft_kinetic_energy_superc_main(Hk,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_NN
#ifdef _MPI
  function dmft_kinetic_energy_superc_NN_mpi(MpiComm,Hk,Wtk,Sigma,Self) result(Eout)
    integer                                 :: MpiComm
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)         :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8)                                 :: Eout(2)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN_mpi","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN_mpi","Self")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_superc_NN_mpi","Hk")
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
    Eout = dmft_kinetic_energy_superc_main_mpi(MpiComm,Hk,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_NN_mpi
#endif



  function dmft_kinetic_energy_superc_NN_lattice(Hk,Wtk,Sigma,Self) result(Ekin_lattice)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)         :: Self   ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Self_  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                   :: i,iorb,ilat,ispin,is
    integer                                   :: j,jorb,jlat,jspin,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8)                                   :: Ekin_lattice(2)
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
    Ekin_lattice = dmft_kinetic_energy_superc_main_lattice(Hk,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_NN_lattice
#ifdef _MPI
  function dmft_kinetic_energy_superc_NN_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Self) result(Ekin_lattice)
    integer                                   :: MpiComm
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)         :: Self   ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Self_  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                   :: i,iorb,ilat,ispin,is
    integer                                   :: j,jorb,jlat,jspin,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8)                                   :: Ekin_lattice(2)
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
    Ekin_lattice = dmft_kinetic_energy_superc_main_lattice_mpi(MpiComm,Hk,Wtk,Sigma_,Self_)
  end function dmft_kinetic_energy_superc_NN_lattice_mpi
#endif








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
    integer :: unit
    unit = free_unit()
    open(unit,file="dmft_kinetic_energy.info")
    write(unit,"(A1,90(A14,1X))")"#",reg(txtfy(1))//"<K>",reg(txtfy(2))//"<Eloc>"
    close(unit)
  end subroutine write_kinetic_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_value(Eout)
    real(8) :: Eout(2),Ekin,Eloc
    integer :: unit
    unit = free_unit()
    Ekin=Eout(1)
    Eloc=Eout(2)
    open(unit,file="dmft_kinetic_energy.dat")
    write(unit,"(90F15.9)")Ekin,Eloc
    close(unit)
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


