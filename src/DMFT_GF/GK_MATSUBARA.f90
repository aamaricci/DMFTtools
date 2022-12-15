module GK_MATSUBARA
  USE DMFT_CTRL_VARS  
  USE GF_COMMON_OLD
  implicit none
  private



  interface get_gk_matsubara
     module procedure :: dmft_get_gk_matsubara_normal_main
     module procedure :: dmft_get_gk_matsubara_normal_cluster
     module procedure :: dmft_get_gk_matsubara_normal_dos
     module procedure :: dmft_get_gk_matsubara_normal_ineq
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_matsubara_superc_main
     module procedure :: dmft_get_gk_matsubara_superc_dos
     module procedure :: dmft_get_gk_matsubara_superc_ineq
  end interface get_gk_matsubara




  interface dmft_gk_matsubara
     module procedure :: dmft_get_gk_matsubara_normal_main
     module procedure :: dmft_get_gk_matsubara_normal_cluster
     module procedure :: dmft_get_gk_matsubara_normal_dos
     module procedure :: dmft_get_gk_matsubara_normal_ineq
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_matsubara_superc_main
     module procedure :: dmft_get_gk_matsubara_superc_dos
     module procedure :: dmft_get_gk_matsubara_superc_ineq
  end interface dmft_gk_matsubara



  !PUBLIC IN DMFT:
  public :: get_gk_matsubara
  public :: dmft_gk_matsubara


contains


  subroutine dmft_get_gk_matsubara_normal_main(Hk,Gkmats,Smats,hk_symm)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                              :: hk_symm   !
    logical                                       :: hk_symm_  !
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
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso],"dmft_get_gk_matsubara_normal_main_mpi","Hk")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_main_mpi","Smats")
    call assert_shape(Gkmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_main_mpi","Gkmats")
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
    Gkmats=zero
    call invert_gk_normal_mpi(zeta_mats,Hk,hk_symm_,Gkmats)    
  end subroutine dmft_get_gk_matsubara_normal_main


  subroutine dmft_get_gk_matsubara_normal_cluster(Hk,Gkmats,Smats,hk_symm)
    complex(8),dimension(:,:),intent(in)              :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                                  :: hk_symm   !
    logical                                           :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable           :: zeta_mats ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
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
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Nlso   = Nlat*Nspin*Norb
    !
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],"dmft_get_gk_matsubara_normal_cluster_mpi","Hk")
    call assert_shape(Smats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_cluster_mpi","Smats")
    call assert_shape(Gkmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_cluster_mpi","Gkmats")
    !
    !Allocate and setup the Matsubara freq.
    allocate(zeta_mats(Nlso,Nlso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Smats(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkmats=zero
    call invert_gk_normal_cluster_mpi(zeta_mats,Hk,hk_symm_,Gkmats)      
  end subroutine dmft_get_gk_matsubara_normal_cluster


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
    if(size(Ebands)/=Nso)stop "dmft_get_gk_matsubara_normal_dos_main_mpi ERROR: size(Ebands)!=Nso"
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_dos_main_mpi","Smats")
    call assert_shape(Gkmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gk_matsubara_normal_dos_main_mpi","Gkmats")
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
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
    !
    do i = 1+mpi_rank, Lmats, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             Gtmp(ispin,ispin,iorb,iorb,i) = Dbands(io)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io) )
          enddo
       enddo
    end do
#ifdef _MPI    
    if(check_MPI())then
       Gkmats=zero
       call Mpi_AllReduce(Gtmp,Gkmats, size(Gkmats), MPI_Double_Complex, MPI_Sum, Mpi_Comm_World, MPI_ierr)
    else
       Gkmats = Gtmp
    endif
#else
    Gkmats = Gtmp
#endif
  end subroutine dmft_get_gk_matsubara_normal_dos


  subroutine dmft_get_gk_matsubara_normal_ineq(Hk,Gkmats,Smats,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,optional                                :: hk_symm   !
    logical                                         :: hk_symm_  !
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Smats")
    call assert_shape(Gkmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Gkmats")
    !
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
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
    Gkmats=zero
    if(.not.tridiag_)then
       call invert_gk_normal_ineq_mpi(zeta_mats,Hk,hk_symm_,Gkmats)
    else
       call invert_gk_normal_tridiag_mpi(zeta_mats,Hk,hk_symm_,Gkmats)
    endif
  end subroutine dmft_get_gk_matsubara_normal_ineq


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gk_matsubara_normal_cluster_ineq(Hk,Gkmats,Smats,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)                :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Smats     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,optional                                    :: hk_symm   !
    logical                                             :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:),allocatable           :: zeta_mats ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
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
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nineq  = size(Smats,1)
    Nlat  = size(Smats,2)
    Nspin = size(Smats,4)
    Norb  = size(Smats,6)
    Lmats = size(Smats,8)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    Nilso = Nineq*Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nilso,Nilso],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Hk")
    call assert_shape(Smats,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Smats")
    call assert_shape(Gkmats,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main_mpi',"Gkmats")
    !
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(zeta_mats(Nineq,Nlso,Nlso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do iineq=1,Nineq
       do i=1,Lmats
          zeta_mats(iineq,:,:,i) = (xi*wm(i)+xmu)*eye(Nlso)     - nnn2lso_cluster_reshape(Smats(iineq,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    Gkmats=zero
    if(.not.tridiag_)then
       call invert_gk_normal_cluster_ineq_mpi(zeta_mats,Hk,hk_symm_,Gkmats)
    else
       call invert_gk_normal_cluster_ineq_tridiag_mpi(zeta_mats,Hk,hk_symm_,Gkmats)
    endif
  end subroutine dmft_get_gk_matsubara_normal_cluster_ineq
#endif












  subroutine dmft_get_gk_matsubara_superc_main(Hk,Gkmats,Smats,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                                :: hk_symm   !
    logical                                         :: hk_symm_  !
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
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_matsubara_superc_main_mpi',"Hk")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main_mpi',"Smats")
    call assert_shape(Gkmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main_mpi',"Gkmats")
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
    call invert_gk_superc_mpi(zeta_mats,Hk,hk_symm_,Gkmats)
  end subroutine dmft_get_gk_matsubara_superc_main


  subroutine dmft_get_gk_matsubara_superc_dos(Ebands,Dbands,Hloc,Gkmats,Smats)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(size(Ebands,2)),intent(in)    :: Dbands    ![Nspin*Norb]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)                      :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout)                   :: Gkmats     !as Smats
    ! arrays
    complex(8)                                                        :: gktmp(2),cdet
    complex(8)                                                        :: zeta_11,zeta_12,zeta_22 
    complex(8),dimension(2,2,size(Ebands),size(Ebands),size(Smats,6)) :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(2*size(Ebands),2*size(Ebands))               :: Gmatrix
    complex(8),dimension(:,:,:,:,:,:),allocatable                     :: Gtmp    !as Smats
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
    call assert_shape(Ebands,[2,Nso],'dmft_get_gk_matsubara_superc_dos',"Ebands")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Smats")
    call assert_shape(Gkmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Gkmats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    zeta_mats=zero
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) - diag(Hloc(1,:)) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                                      -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                                      -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) - diag(Hloc(2,:)) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
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
       call Mpi_AllReduce(Gtmp,Gkmats, size(Gkmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gkmats=Gtmp
    endif
#else
    Gkmats=Gtmp
#endif
  end subroutine dmft_get_gk_matsubara_superc_dos


  subroutine dmft_get_gk_matsubara_superc_ineq(Hk,Gkmats,Smats,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
    logical,optional                                  :: hk_symm   !
    logical                                           :: hk_symm_  !
    !allocatable arrays
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
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_matsubara_superc_ineq_main_mpi',"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_ineq_main_mpi',"Smats")
    call assert_shape(Gkmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_ineq_main_mpi',"Gkmats")
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
    call invert_gk_superc_ineq_mpi(zeta_mats,Hk,hk_symm_,Gkmats)
  end subroutine dmft_get_gk_matsubara_superc_ineq







end module GK_MATSUBARA
