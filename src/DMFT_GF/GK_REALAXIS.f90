module GK_REALAXIS
  USE DMFT_CTRL_VARS  
  USE GF_COMMON_OLD
  implicit none
  private


  interface get_gk_realaxis
     module procedure :: dmft_get_gk_realaxis_normal_main
     module procedure :: dmft_get_gk_realaxis_normal_cluster
     module procedure :: dmft_get_gk_realaxis_normal_dos
     module procedure :: dmft_get_gk_realaxis_normal_ineq
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_realaxis_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_realaxis_superc_main
     module procedure :: dmft_get_gk_realaxis_superc_dos
     module procedure :: dmft_get_gk_realaxis_superc_ineq
  end interface get_gk_realaxis



  interface dmft_gk_realaxis
     module procedure :: dmft_get_gk_realaxis_normal_main
     module procedure :: dmft_get_gk_realaxis_normal_cluster
     module procedure :: dmft_get_gk_realaxis_normal_dos
     module procedure :: dmft_get_gk_realaxis_normal_ineq
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_realaxis_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_realaxis_superc_main
     module procedure :: dmft_get_gk_realaxis_superc_dos
     module procedure :: dmft_get_gk_realaxis_superc_ineq
  end interface dmft_gk_realaxis


  !PUBLIC IN DMFT:
  public :: get_gk_realaxis
  public :: dmft_gk_realaxis


contains



  subroutine dmft_get_gk_realaxis_normal_main(hk,Gkreal,Sreal,hk_symm)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    logical,optional                              :: hk_symm   !
    logical                                       :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal] 
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
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso],'dmft_get_gk_realaxis_normal_main_mpi',"Hk")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_main_mpi',"Gkreal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkreal=zero
    call invert_gk_normal_mpi(zeta_real,Hk,hk_symm_,Gkreal)      
  end subroutine dmft_get_gk_realaxis_normal_main


  subroutine dmft_get_gk_realaxis_normal_cluster(hk,Gkreal,Sreal,hk_symm)
    complex(8),dimension(:,:),intent(in)              :: Hk    ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    logical,optional                                  :: hk_symm   !
    logical                                           :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable           :: zeta_real ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
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
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Nlso   = Nlat*Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_realaxis_normal_cluster_mpi',"Hk")
    call assert_shape(Sreal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_cluster_mpi',"Sreal")
    call assert_shape(Gkreal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_cluster_mpi',"Gkreal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(zeta_real(Nlso,Nlso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Sreal(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkreal=zero
    call invert_gk_normal_cluster_mpi(zeta_real,Hk,hk_symm_,Gkreal)      
  end subroutine dmft_get_gk_realaxis_normal_cluster


  subroutine dmft_get_gk_realaxis_normal_dos(Ebands,Dbands,Hloc,Gkreal,Sreal)
    real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
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
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb    
    !Testing part:
    !Testing part:
    if(size(Ebands)/=Nso)stop "dmft_get_gk_realaxis_normal_dos_main_mpi ERROR: size(Ebands)!=Nso"
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_dos_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_dos_main_mpi',"Gkreal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkreal=zero
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
    do i = 1+mpi_rank, Lreal, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             Gtmp(ispin,ispin,iorb,iorb,i) = Dbands(io)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io) )
          enddo
       enddo
    end do
#ifdef _MPI
    if(check_MPI())then
       call Mpi_AllReduce(Gtmp,Gkreal, size(Gkreal), MPI_Double_Complex, MPI_Sum, Mpi_Comm_World, MPI_ierr)
    else
       Gkreal=Gtmp
    endif
#else
    Gkreal=Gtmp
#endif    
  end subroutine dmft_get_gk_realaxis_normal_dos


  subroutine dmft_get_gk_realaxis_normal_ineq(hk,Gkreal,Sreal,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,optional                                :: hk_symm   !
    logical                                         :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
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
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Gkreal")
    !
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    Gkreal=zero
    if(.not.tridiag_)then
       call invert_gk_normal_ineq_mpi(zeta_real,Hk,hk_symm_,Gkreal)
    else
       call invert_gk_normal_tridiag_mpi(zeta_real,Hk,hk_symm_,Gkreal)
    endif
  end subroutine dmft_get_gk_realaxis_normal_ineq



#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gk_realaxis_normal_cluster_ineq(Hk,Gkreal,Sreal,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)                :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sreal     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,optional                                    :: hk_symm   !
    logical                                             :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable           :: zeta_real ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
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
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nineq  = size(Sreal,1)
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,4)
    Norb  = size(Sreal,6)
    Lreal = size(Sreal,8)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    Nilso = Nineq*Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nilso,Nilso],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Gkreal")
    !
    allocate(zeta_real(Nineq,Nlso,Nlso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do iineq=1,Nineq
       do i=1,Lreal
          zeta_real(iineq,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nlso)     - nnn2lso_cluster_reshape(Sreal(iineq,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    Gkreal=zero
    if(.not.tridiag_)then
       call invert_gk_normal_cluster_ineq_mpi(zeta_real,Hk,hk_symm_,Gkreal)
    else
       call invert_gk_normal_cluster_ineq_tridiag_mpi(zeta_real,Hk,hk_symm_,Gkreal)
    endif
  end subroutine dmft_get_gk_realaxis_normal_cluster_ineq
#endif








  subroutine dmft_get_gk_realaxis_superc_main(hk,Gkreal,Sreal,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !allocatable arrays
    logical,optional                                :: hk_symm   !
    logical                                         :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
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
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_realaxis_superc_main_mpi',"Hk")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main_mpi',"Sreal")
    call assert_shape(Gkreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main_mpi',"Gkreal")
    !
    allocate(zeta_real(2,2,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
            nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
            conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkreal=zero
    call invert_gk_superc_mpi(zeta_real,Hk,hk_symm_,Gkreal)
  end subroutine dmft_get_gk_realaxis_superc_main


  subroutine dmft_get_gk_realaxis_superc_dos(Ebands,Dbands,Hloc,Gkreal,Sreal)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(size(Ebands,2)),intent(in)    :: Dbands    ![Nspin*Norb]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    ! arrays
    complex(8)                                      :: gktmp(2),cdet
    complex(8)                                      :: zeta_11,zeta_12,zeta_22 
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp     ![2][Nspin][Nspin][Norb][Norb][Lreal]
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
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Ebands,[2,Nso],'dmft_get_gk_realaxis_superc_dos',"Ebands")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main',"Sreal")
    call assert_shape(Gkreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main',"Gkreal")
    !
    allocate(zeta_real(2,2,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
            nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
            conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    Gkreal=zero
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
    !
    do i = 1+mpi_rank, Lreal, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_11 = zeta_real(1,1,io,io,i)
             zeta_12 = zeta_real(1,2,io,io,i)
             zeta_12 = zeta_real(2,2,io,io,i)
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
       call Mpi_AllReduce(Gtmp,Gkreal, size(Gkreal), MPI_Double_Complex, MPI_Sum, Mpi_Comm_World, MPI_ierr)
    else
       Gkreal=Gtmp
    endif
#else
    Gkreal=Gtmp
#endif
  end subroutine dmft_get_gk_realaxis_superc_dos


  subroutine dmft_get_gk_realaxis_superc_ineq(hk,Gkreal,Sreal,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    logical,optional                                  :: hk_symm   !
    logical                                           :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
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
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_realaxis_superc_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_ineq_main_mpi',"Sreal")
    call assert_shape(Gkreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_ineq_main_mpi',"Gkreal")
    !
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       !SYMMETRIES in real-frequencies   [assuming a real order parameter]
       !G22(w)  = -[G11[-w]]*
       !G21(w)  =   G12[w]   
       do i=1,Lreal
          zeta_real(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
               nn2so_reshape(Sreal(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(1,2,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,1,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
               conjg( nn2so_reshape(Sreal(1,ilat,:,:,:,:,Lreal+1-i),Nspin,Norb) )
       enddo
    enddo
    !
    Gkreal=zero
    call invert_gk_superc_ineq_mpi(zeta_real,Hk,hk_symm_,Gkreal)
  end subroutine dmft_get_gk_realaxis_superc_ineq















end module GK_REALAXIS
