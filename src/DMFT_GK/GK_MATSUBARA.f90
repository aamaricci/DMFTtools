module GK_MATSUBARA
  USE DMFT_CTRL_VARS  
  USE GK_COMMON
  implicit none
  private

  interface dmft_gk_matsubara
     module procedure :: dmft_get_gk_matsubara_normal_main
     module procedure :: dmft_get_gk_matsubara_normal_dos
     module procedure :: dmft_get_gk_matsubara_normal_ineq
     module procedure :: dmft_get_gk_matsubara_superc_main
     module procedure :: dmft_get_gk_matsubara_superc_dos
     module procedure :: dmft_get_gk_matsubara_superc_ineq
  end interface dmft_gk_matsubara

  public :: dmft_gk_matsubara

contains



  subroutine dmft_get_gk_matsubara_normal_main(Hk,Gkmats,Smats)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
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
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Nso   = Nspin*Norb    
    !Testing part:
    if(mpi_master)then
       call assert_shape(Hk,[Nso,Nso],"dmft_get_gk_matsubara_normal_main_mpi","Hk")
       call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_main_mpi","Smats")
       call assert_shape(Gkmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_main_mpi","Gkmats")
    endif
    !
    !Allocate and setup the Matsubara freq.
    allocate(zeta_mats(Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    call invert_gk_normal(zeta_mats,Hk,.false.,Gkmats)
  end subroutine dmft_get_gk_matsubara_normal_main


  subroutine dmft_get_gk_matsubara_normal_dos(Ebands,Dbands,Hloc,Gkmats,Smats)
    real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    !
    complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Smats
    !
    real(8)                                       :: beta
    real(8)                                       :: xmu,eps
    real(8)                                       :: wini,wfin
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
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Nso   = Nspin*Norb    
    !Testing part:
    if(mpi_master)then
       if(size(Ebands)/=Nso)stop "dmft_get_gk_matsubara_normal_dos_main_mpi ERROR: size(Ebands)!=Nso"
       call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_dos_main_mpi","Smats")
       call assert_shape(Gkmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_dos_main_mpi","Gkmats")
    endif
    !
    !Allocate and setup the Matsubara freq.
    allocate(zeta_mats(Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats))
    Gtmp=zero
    do i = 1+mpi_rank, Lmats, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             Gtmp(ispin,ispin,iorb,iorb,i) = Dbands(io)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io) )
          enddo
       enddo
    end do
    !
#ifdef _MPI    
    if(check_MPI())then
       Gkmats=zero
       call AllReduce_MPI(MPI_COMM_WORLD,Gtmp,Gkmats)
    else
       Gkmats=Gtmp
    endif
#else
    Gkmats=Gtmp
#endif
  end subroutine dmft_get_gk_matsubara_normal_dos


  subroutine dmft_get_gk_matsubara_normal_ineq(Hk,Gkmats,Smats,tridiag)
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
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
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    if(mpi_master)then
       call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Hk")
       call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Smats")
       call assert_shape(Gkmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Gkmats")
    endif
    !
    !
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(.not.tridiag_)then
       call invert_gk_normal_ineq(zeta_mats,Hk,.false.,Gkmats)
    else
       call invert_gk_normal_tridiag(zeta_mats,Hk,.false.,Gkmats)
    endif
  end subroutine dmft_get_gk_matsubara_normal_ineq












  subroutine dmft_get_gk_matsubara_superc_main(Hk,Gkmats,Smats)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
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
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Nso   = Nspin*Norb
    !Testing part:
    if(mpi_master)then
       call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_matsubara_superc_main_mpi',"Hk")
       call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main_mpi',"Smats")
       call assert_shape(Gkmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main_mpi',"Gkmats")
    endif
    !
    allocate(zeta_mats(2,2,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkmats=zero
    call invert_gk_superc(zeta_mats,Hk,.false.,Gkmats)
  end subroutine dmft_get_gk_matsubara_superc_main


  subroutine dmft_get_gk_matsubara_superc_dos(Ebands,Dbands,Hloc,Gkmats,Smats)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(size(Ebands,2)),intent(in)    :: Dbands    ![Nspin*Norb]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    !allocatable arrays
    complex(8)                                      :: gktmp(2),cdet
    complex(8)                                      :: zeta_11,zeta_12,zeta_22 
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
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
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Nso   = Nspin*Norb
    !Testing part:
    if(mpi_master)then
       call assert_shape(Ebands,[2,Nso],'dmft_get_gk_matsubara_superc_dos',"Ebands")
       call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Smats")
       call assert_shape(Gkmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Gkmats")
    endif
    !
    allocate(zeta_mats(2,2,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkmats=zero
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
    do i = 1+mpi_rank, Lmats, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_11 = zeta_mats(1,1,io,io,i)
             zeta_12 = zeta_mats(1,2,io,io,i)
             zeta_12 = zeta_mats(2,2,io,io,i)
             !
             cdet = (zeta_11-Hloc(1,io)-Ebands(1,io))*(zeta_22-Hloc(2,io)-Ebands(2,io)) - zeta_12**2
             gktmp(1)=-(zeta_22-Hloc(2,io)-Ebands(2,io))/cdet
             gktmp(2)=  zeta_12/cdet
             Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + gktmp(1)*Dbands(io)
             Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + gktmp(2)*Dbands(io)
          enddo
       enddo
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Gkmats=zero
       call AllReduce_MPI(MPI_COMM_WORLD,Gtmp,Gkmats)
    else
       Gkmats=Gtmp
    endif
#else
    Gkmats=Gtmp
#endif
  end subroutine dmft_get_gk_matsubara_superc_dos


  subroutine dmft_get_gk_matsubara_superc_ineq(Hk,Gkmats,Smats)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
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
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    if(mpi_master)then
       call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_matsubara_superc_ineq_main_mpi',"Hk")
       call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_ineq_main_mpi',"Smats")
       call assert_shape(Gkmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_ineq_main_mpi',"Gkmats")
    endif
    !
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))    
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
       !G22(iw) = -[G11[iw]]*
       !G21(iw) =   G12[w]
       do i=1,Lmats
          zeta_mats(1,1,ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,2,ilat,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb) )
       enddo
    enddo
    !
    Gkmats=zero
    call invert_gk_superc_ineq(zeta_mats,Hk,.false.,Gkmats)
  end subroutine dmft_get_gk_matsubara_superc_ineq





end module GK_MATSUBARA
