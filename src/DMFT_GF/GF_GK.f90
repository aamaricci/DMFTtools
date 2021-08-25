module GF_GK
  USE GF_COMMON
  implicit none
  private


  !##################################################################
  interface get_gk
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gk_normal_main
     module procedure :: dmft_get_gk_normal_dos
     module procedure :: dmft_get_gk_normal_ineq
     module procedure :: dmft_get_gk_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_superc_main
     module procedure :: dmft_get_gk_superc_dos
     module procedure :: dmft_get_gk_superc_ineq
  end interface get_gk

  interface get_gk_matsubara
     module procedure :: dmft_get_gk_matsubara_normal_main
     module procedure :: dmft_get_gk_matsubara_normal_dos
     module procedure :: dmft_get_gk_matsubara_normal_ineq
     module procedure :: dmft_get_gk_matsubara_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_matsubara_superc_main
     module procedure :: dmft_get_gk_matsubara_superc_dos
     module procedure :: dmft_get_gk_matsubara_superc_ineq
  end interface get_gk_matsubara

  interface get_gk_realaxis
     module procedure :: dmft_get_gk_realaxis_normal_main
     module procedure :: dmft_get_gk_realaxis_normal_dos
     module procedure :: dmft_get_gk_realaxis_normal_ineq
     module procedure :: dmft_get_gk_realaxis_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_realaxis_normal_cluster_ineq
#endif
     module procedure :: dmft_get_gk_realaxis_superc_main
     module procedure :: dmft_get_gk_realaxis_superc_dos
     module procedure :: dmft_get_gk_realaxis_superc_ineq
  end interface get_gk_realaxis


  !##################################################################



  !##################################################################
  interface dmft_gk
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gk_normal_main
     module procedure :: dmft_get_gk_normal_dos
     module procedure :: dmft_get_gk_normal_ineq
     module procedure :: dmft_get_gk_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_superc_main
     module procedure :: dmft_get_gk_superc_dos
     module procedure :: dmft_get_gk_superc_ineq
  end interface dmft_gk

  interface dmft_gk_matsubara
     module procedure :: dmft_get_gk_matsubara_normal_main
     module procedure :: dmft_get_gk_matsubara_normal_dos
     module procedure :: dmft_get_gk_matsubara_normal_ineq
     module procedure :: dmft_get_gk_matsubara_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_matsubara_superc_main
     module procedure :: dmft_get_gk_matsubara_superc_dos
     module procedure :: dmft_get_gk_matsubara_superc_ineq
  end interface dmft_gk_matsubara

  interface dmft_gk_realaxis
     module procedure :: dmft_get_gk_realaxis_normal_main
     module procedure :: dmft_get_gk_realaxis_normal_dos
     module procedure :: dmft_get_gk_realaxis_normal_ineq
     module procedure :: dmft_get_gk_realaxis_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_realaxis_normal_cluster_ineq
#endif
     module procedure :: dmft_get_gk_realaxis_superc_main
     module procedure :: dmft_get_gk_realaxis_superc_dos
     module procedure :: dmft_get_gk_realaxis_superc_ineq
  end interface dmft_gk_realaxis




  !##################################################################
  interface dmft_get_gk
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gk_normal_main
     module procedure :: dmft_get_gk_normal_dos
     module procedure :: dmft_get_gk_normal_ineq
     module procedure :: dmft_get_gk_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_superc_main
     module procedure :: dmft_get_gk_superc_dos
     module procedure :: dmft_get_gk_superc_ineq
  end interface dmft_get_gk

  interface dmft_get_gk_matsubara
     module procedure :: dmft_get_gk_matsubara_normal_main
     module procedure :: dmft_get_gk_matsubara_normal_dos
     module procedure :: dmft_get_gk_matsubara_normal_ineq
     module procedure :: dmft_get_gk_matsubara_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gk_matsubara_superc_main
     module procedure :: dmft_get_gk_matsubara_superc_dos
     module procedure :: dmft_get_gk_matsubara_superc_ineq
  end interface dmft_get_gk_matsubara

  interface dmft_get_gk_realaxis
     module procedure :: dmft_get_gk_realaxis_normal_main
     module procedure :: dmft_get_gk_realaxis_normal_dos
     module procedure :: dmft_get_gk_realaxis_normal_ineq
     module procedure :: dmft_get_gk_realaxis_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gk_realaxis_normal_cluster_ineq
#endif
     module procedure :: dmft_get_gk_realaxis_superc_main
     module procedure :: dmft_get_gk_realaxis_superc_dos
     module procedure :: dmft_get_gk_realaxis_superc_ineq
  end interface dmft_get_gk_realaxis
  !##################################################################




  !PUBLIC IN DMFT:
  public :: get_gk
  public :: get_gk_matsubara
  public :: get_gk_realaxis


  public :: dmft_gk
  public :: dmft_gk_matsubara
  public :: dmft_gk_realaxis

  public :: dmft_get_gk
  public :: dmft_get_gk_matsubara
  public :: dmft_get_gk_realaxis

contains




  subroutine dmft_get_gk_normal_main(Hk,Gk,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:),optional              :: zeta
    logical,optional                              :: hk_symm   !
    logical                                       :: hk_symm_  !
    complex(8),dimension(:,:,:),allocatable       :: csi ![Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso],"dmft_get_gk_normal_main_mpi","Hk")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gk_normal_main_mpi","Sigma")
    call assert_shape(Gk,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gk_normal_main_mpi","Gk")
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Nso,Nso,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Nso) - nn2so_reshape(Sigma(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Csi-Hk) for each k-point
    Gk=zero
    call invert_gk_normal_mpi(csi,Hk,hk_symm_,Gk)    
  end subroutine dmft_get_gk_normal_main




  subroutine dmft_get_gk_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,axis,zeta)
    real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb]
    real(8),dimension(size(Ebands)),intent(in)    :: Dbands    ![Nspin*Norb]
    real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:),optional              :: zeta
    !
    complex(8),dimension(:,:,:),allocatable       :: csi ![Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sigma
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)
    Nso   = Nspin*Norb    
    !Testing part:
    if(size(Ebands)/=Nso)stop "dmft_get_gk_normal_dos_main_mpi ERROR: size(Ebands)!=Nso"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gk_normal_dos_main_mpi","Sigma")
    call assert_shape(Gk,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gk_normal_dos_main_mpi","Gk")
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Nso,Nso,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Nso) - nn2so_reshape(Sigma(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
    !
    do i = 1+mpi_rank, Lfreq, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             Gtmp(ispin,ispin,iorb,iorb,i) = Dbands(io)/( csi(io,io,i)-Hloc(io)-Ebands(io) )
          enddo
       enddo
    end do
#ifdef _MPI    
    if(check_MPI())then
       Gk=zero
       call Mpi_AllReduce(Gtmp,Gk, size(Gk), MPI_Double_Complex, MPI_Sum, Mpi_Comm_World, MPI_ierr)
    else
       Gk = Gtmp
    endif
#else
    Gk = Gtmp
#endif
  end subroutine dmft_get_gk_normal_dos










  subroutine dmft_get_gk_normal_ineq(Hk,Gk,Sigma,axis,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                            :: axis
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,optional                                :: hk_symm   !
    logical                                         :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
    complex(8),dimension(:,:,:,:),allocatable       :: csi ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_normal_ineq_main_mpi',"Sigma")
    call assert_shape(Gk,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_normal_ineq_main_mpi',"Gk")
    !
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(csi(Nlat,Nso,Nso,Lfreq))
    !
    do ilat=1,Nlat
       do i=1,Lfreq
          csi(ilat,:,:,i) = (wfreq(i)+xmu)*eye(Nso)     - nn2so_reshape(Sigma(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Csi-Hk) for each k-point 
    Gk=zero
    if(.not.tridiag_)then
       call invert_gk_normal_ineq_mpi(csi,Hk,hk_symm_,Gk)
    else
       call invert_gk_normal_tridiag_mpi(csi,Hk,hk_symm_,Gk)
    endif
  end subroutine dmft_get_gk_normal_ineq










  subroutine dmft_get_gk_normal_cluster(Hk,Gk,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:),intent(in)              :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm   !
    logical                                           :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable           :: csi ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso   = Nlat*Nspin*Norb
    !
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],"dmft_get_gk_normal_cluster_mpi","Hk")
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gk_normal_cluster_mpi","Sigma")
    call assert_shape(Gk,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gk_normal_cluster_mpi","Gk")
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Nlso,Nlso,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Sigma(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    !invert (Csi-Hk) for each k-point
    Gk=zero
    call invert_gk_normal_cluster_mpi(csi,Hk,hk_symm_,Gk)      
  end subroutine dmft_get_gk_normal_cluster




#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gk_normal_cluster_ineq(Hk,Gk,Sigma,axis,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)                :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:),optional                    :: zeta
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,optional                                    :: hk_symm   !
    logical                                             :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
    complex(8),dimension(:,:,:,:),allocatable           :: csi ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nineq  = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    Nilso = Nineq*Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nilso,Nilso],'dmft_get_gk_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_normal_ineq_main_mpi',"Sigma")
    call assert_shape(Gk,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_normal_ineq_main_mpi',"Gk")
    !
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(csi(Nineq,Nlso,Nlso,Lfreq))
    do iineq=1,Nineq
       do i=1,Lfreq
          csi(iineq,:,:,i) = (wfreq(i)+xmu)*eye(Nlso)     - nnn2lso_cluster_reshape(Sigma(iineq,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Csi-Hk) for each k-point 
    Gk=zero
    if(.not.tridiag_)then
       call invert_gk_normal_cluster_ineq_mpi(csi,Hk,hk_symm_,Gk)
    else
       call invert_gk_normal_cluster_ineq_tridiag_mpi(csi,Hk,hk_symm_,Gk)
    endif
  end subroutine dmft_get_gk_normal_cluster_ineq
#endif












  subroutine dmft_get_gk_superc_main(Hk,Gk,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                            :: axis
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: hk_symm   !
    logical                                         :: hk_symm_  !
    complex(8),dimension(:,:,:,:,:),allocatable     :: csi ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
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
    !
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_superc_main_mpi',"Hk")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_superc_main_mpi',"Sigma")
    call assert_shape(Gk,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_superc_main_mpi',"Gk")
    !
    allocate(csi(2,2,Nso,Nso,Lfreq))
    select case(axis)
    case default; stop "dmft_get_gk_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso) -        nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) =                         -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) =                         -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb) )
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i) + xmu)*eye(Nso)                 - nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) =                                           - nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) =                                           - nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = -conjg(wfreq(Lfreq+1-i) + xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,Lfreq+1-i),Nspin,Norb) )
       enddo
    end select
    !invert (Csi-Hk) for each k-point
    Gk=zero
    call invert_gk_superc_mpi(csi,Hk,hk_symm_,Gk)
  end subroutine dmft_get_gk_superc_main


  subroutine dmft_get_gk_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,axis,zeta)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(:),intent(in)                 :: Dbands    ![Nspin*Norb]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:),optional                :: zeta
    ! arrays
    complex(8)                                      :: gktmp(2),cdet
    complex(8)                                      :: zeta_11,zeta_12,zeta_21,zeta_22 
    complex(8),dimension(:,:,:,:,:),allocatable     :: csi ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Ebands,[2,Nso],'dmft_get_gk_superc_dos',"Ebands")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_superc_main',"Sigma")
    call assert_shape(Gk,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_superc_main',"Gk")
    !
    allocate(csi(2,2,Nso,Nso,Lfreq))
    csi=zero
    do i=1,Lfreq
       csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso) - diag(Hloc(1,:)) -   nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
       csi(1,2,:,:,i) =                                      -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
       csi(2,1,:,:,i) =                                      -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
       csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Nso) - diag(Hloc(2,:)) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    select case(axis)
    case default; stop "dmft_get_gk_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso) - diag(Hloc(1,:)) -   nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) =                                      -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) =                                      -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Nso) - diag(Hloc(2,:)) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb) )
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso)                 - diag(Hloc(1,:)) - nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) = -nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) = -nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = -conjg( wfreq(Lfreq+1-i)+xmu )*eye(Nso) + diag(Hloc(2,:)) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,Lfreq+1-i),Nspin,Norb) )
       enddo
    end select
    !
    !invert (Csi-Hk) for each k-point
    Gk=zero
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
    do i = 1+mpi_rank, Lfreq, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_11 = csi(1,1,io,io,i)
             zeta_12 = csi(1,2,io,io,i)
             zeta_21 = csi(2,1,io,io,i)
             zeta_22 = csi(2,2,io,io,i)
             !
             cdet = (zeta_11-Hloc(1,io)-Ebands(1,io))*(zeta_22-Hloc(2,io)-Ebands(2,io)) - zeta_12*zeta_21
             gktmp(1)=-(zeta_22-Hloc(2,io)-Ebands(2,io))/cdet
             gktmp(2)=  zeta_12/cdet
             Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + gktmp(1)*Dbands(io)
             Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + gktmp(2)*Dbands(io)
          enddo
       enddo
    enddo
#ifdef _MPI    
    if(check_MPI())then
       call Mpi_AllReduce(Gtmp,Gk, size(Gk), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gk=Gtmp
    endif
#else
    Gk=Gtmp
#endif
  end subroutine dmft_get_gk_superc_dos


  subroutine dmft_get_gk_superc_ineq(Hk,Gk,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm   !
    logical                                           :: hk_symm_  !
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: csi ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_superc_ineq_main_mpi',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_superc_ineq_main_mpi',"Sigma")
    call assert_shape(Gk,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gk_superc_ineq_main_mpi',"Gk")
    !
    allocate(csi(2,2,Nlat,Nso,Nso,Lfreq));csi=zero
    select case(axis)
    case default; stop "dmft_get_gk_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do ilat=1,Nlat
          !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
          !G22(iw) = -[G11[iw]]*
          !G21(iw) =   G12[w]
          do i=1,Lfreq
             csi(1,1,ilat,:,:,i) = (wfreq(i)+xmu)*eye(Nso) -        nn2so_reshape(Sigma(1,ilat,:,:,:,:,i),Nspin,Norb)
             csi(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,2,ilat,:,:,i) = (wfreq(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,ilat,:,:,:,:,i),Nspin,Norb) )
          enddo
       enddo
    case("realaxis","real","Realaxis","Real")
       do ilat=1,Nlat
          !SYMMETRIES in real-frequencies   [assuming a real order parameter]
          !G22(w)  = -[G11[-w]]*
          !G21(w)  =   G12[w]   
          do i=1,Lfreq
             csi(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
                  nn2so_reshape(Sigma(1,ilat,:,:,:,:,i),Nspin,Norb)
             csi(1,2,ilat,:,:,i) =                                                  - &
                  nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,1,ilat,:,:,i) =                                                  - &
                  nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lfreq+1-i),eps)+xmu )*eye(Nso) + &
                  conjg( nn2so_reshape(Sigma(1,ilat,:,:,:,:,Lfreq+1-i),Nspin,Norb) )
          enddo
       enddo
    end select
    !
    Gk=zero
    call invert_gk_superc_ineq_mpi(csi,Hk,hk_symm_,Gk)
  end subroutine dmft_get_gk_superc_ineq













  !##################################################################
  !##################################################################
  !                       INTERFACES
  !##################################################################
  !##################################################################
  subroutine dmft_get_gk_matsubara_normal_main(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk         !as Sigma
    complex(8),dimension(:),optional              :: zeta
    logical,optional                              :: hk_symm   ![Lk]
    logical                                       :: hk_symm_  ![Lk]
    !
    hk_symm_  =.false. ; if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_main(Hk,Gk,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_main(Hk,Gk,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_matsubara_normal_main

  subroutine dmft_get_gk_realaxis_normal_main(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk         !as Sigma
    complex(8),dimension(:),optional              :: zeta
    logical,optional                              :: hk_symm   ![Lk]
    logical                                       :: hk_symm_  ![Lk]
    !
    hk_symm_  =.false. ; if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_main(Hk,Gk,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_main(Hk,Gk,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_realaxis_normal_main




  !##################################################################
  !##################################################################





  subroutine dmft_get_gk_matsubara_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,zeta)
    real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb]
    real(8),dimension(:),intent(in)               :: Dbands    ![Nspin*Norb /[1]
    real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk        !as Sigma
    complex(8),dimension(:),optional              :: zeta
    if(present(zeta))then
       call dmft_get_gk_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,"matsubara")
    else
       call dmft_get_gk_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,"matsubara",zeta)
    endif
  end subroutine dmft_get_gk_matsubara_normal_dos

  subroutine dmft_get_gk_realaxis_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,zeta)
    real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(:),intent(in)               :: Dbands    ![Nspin*Norb][Lk] /[1][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk         !as Sigma
    complex(8),dimension(:),optional              :: zeta
    if(present(zeta))then
       call dmft_get_gk_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,"realaxis")
    else
       call dmft_get_gk_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,"realaxis",zeta)
    endif
  end subroutine dmft_get_gk_realaxis_normal_dos



  !##################################################################
  !##################################################################


  subroutine dmft_get_gk_matsubara_normal_ineq(Hk,Gk,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)            :: Hk   ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,optional                                :: hk_symm
    logical                                         :: hk_symm_
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_ineq(Hk,Gk,Sigma,"matsubara",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_ineq(Hk,Gk,Sigma,"matsubara",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_matsubara_normal_ineq

  subroutine dmft_get_gk_realaxis_normal_ineq(Hk,Gk,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)            :: Hk   ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,optional                                :: hk_symm
    logical                                         :: hk_symm_
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_ineq(Hk,Gk,Sigma,"realaxis",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_ineq(Hk,Gk,Sigma,"realaxis",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_realaxis_normal_ineq






  !##################################################################
  !##################################################################



  subroutine dmft_get_gk_matsubara_normal_cluster(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:),intent(in)              :: Hk ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm   ![Lk]
    logical                                           :: hk_symm_  ![Lk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_cluster(Hk,Gk,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_cluster(Hk,Gk,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_matsubara_normal_cluster

  subroutine dmft_get_gk_realaxis_normal_cluster(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:),intent(in)              :: Hk ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm   ![Lk]
    logical                                           :: hk_symm_  ![Lk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_cluster(Hk,Gk,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_cluster(Hk,Gk,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_realaxis_normal_cluster





  !##################################################################
  !##################################################################

#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gk_matsubara_normal_cluster_ineq(Hk,Gk,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)                :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                    :: zeta
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,optional                                    :: hk_symm
    logical                                             :: hk_symm_
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_cluster_ineq(Hk,Gk,Sigma,"matsubara",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_cluster_ineq(Hk,Gk,Sigma,"matsubara",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_matsubara_normal_cluster_ineq

  subroutine dmft_get_gk_realaxis_normal_cluster_ineq(Hk,Gk,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:),intent(in)                :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                    :: zeta
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,optional                                    :: hk_symm
    logical                                             :: hk_symm_
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_normal_cluster_ineq(Hk,Gk,Sigma,"realaxis",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gk_normal_cluster_ineq(Hk,Gk,Sigma,"realaxis",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_realaxis_normal_cluster_ineq
#endif





  !##################################################################
  !##################################################################


  subroutine dmft_get_gk_matsubara_superc_main(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk     ![2][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma  ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: hk_symm
    logical                                         :: hk_symm_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_superc_main(Hk,Gk,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gk_superc_main(Hk,Gk,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_matsubara_superc_main

  subroutine dmft_get_gk_realaxis_superc_main(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk     ![2][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma  ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: hk_symm
    logical                                         :: hk_symm_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_superc_main(Hk,Gk,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gk_superc_main(Hk,Gk,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_realaxis_superc_main






  !##################################################################
  !##################################################################

  subroutine dmft_get_gk_matsubara_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,zeta)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(:),intent(in)                 :: Dbands    ![Nspin*Norb]/[1]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    !
    if(present(zeta))then
       call dmft_get_gk_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,"matsubara")
    else
       call dmft_get_gk_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,"matsubara",zeta)
    endif
  end subroutine dmft_get_gk_matsubara_superc_dos


  subroutine dmft_get_gk_realaxis_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,zeta)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(:),intent(in)                 :: Dbands    ![Nspin*Norb]/[1]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                :: zeta  !
    if(present(zeta))then
       call dmft_get_gk_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,"realaxis")
    else
       call dmft_get_gk_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,"realaxis",zeta)
    endif
  end subroutine dmft_get_gk_realaxis_superc_dos





  !##################################################################
  !##################################################################

  subroutine dmft_get_gk_matsubara_superc_ineq(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk   ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma![2][Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm
    logical                                           :: hk_symm_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_superc_ineq(Hk,Gk,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gk_superc_ineq(Hk,Gk,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_matsubara_superc_ineq



  subroutine dmft_get_gk_realaxis_superc_ineq(Hk,Gk,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk   ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma![2][Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk     !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm
    logical                                           :: hk_symm_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gk_superc_ineq(Hk,Gk,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gk_superc_ineq(Hk,Gk,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gk_realaxis_superc_ineq








end module GF_GK
