module LEGACY_GF_GK
  USE LEGACY_GF_COMMON
  implicit none
  private


  public :: get_gk_normal_main
  public :: get_gk_normal_tridiag
  public :: get_gk_normal_dos
  !
  public :: get_gk_superc_main
  public :: get_gk_superc_dos
  !
  public :: get_gk_normal_hk_rank4
  public :: get_gk_normal_dos_rank4
  public :: get_gk_normal_hk_rank5
  public :: get_gk_normal_tridiag_rank5
  public :: get_gk_normal_hk_rank5_6
  public :: get_gk_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gk_normal_hk_rank7
  public :: get_gk_normal_tridiag_rank7
#endif
  !
  public :: get_gk_superc_hk_rank4
  public :: get_gk_superc_dos_rank4
  public :: get_gk_superc_hk_rank5
  public :: get_gk_superc_hk_rank5_6



contains



  !##################################################################
  !##################################################################
  !                       NORMAL PHASE
  !##################################################################
  !##################################################################
  subroutine get_gk_normal_main(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)      :: Hk        ![N,N]
    complex(8),dimension(:,:,:),intent(in)    :: Sigma     ![N,N][Lfreq]
    complex(8),dimension(:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                          :: axis
    complex(8),dimension(:,:,:),allocatable   :: csi   ![N,N][Lfreq]
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
    Ntot  = size(Hk,1)
    Lfreq = size(Sigma,3)
    !Testing part:
    call assert_shape(Hk,[Ntot,Ntot],"get_gk_normal_main","Hk")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"get_gk_normal_main","Sigma")
    call assert_shape(Gk, [Ntot,Ntot,Lfreq],"get_gk_normal_main","Gk")
    !
    call build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Ntot,Ntot,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
    enddo
    !
    !invert (Csi-Hk) for each k-point
    Gk = zero
    call invert_gk_normal_mpi(csi,Hk,Gk)      
  end subroutine get_gk_normal_main



  !TRIDIAG Hk
  !Hk has a blocks tridiagonal form with blocks of size [Ncell,Ncell]. Ntot=Nsites*Ncell,
  !with: Nsites the number of sites in the super-cell
  !      Ncell the dimension of the unit cell 
  subroutine get_gk_normal_tridiag(Hk,Gk,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:),intent(in)      :: Hk      ![N,N]
    complex(8),dimension(:,:,:),intent(in)    :: Sigma   ![N,N,Lfreq]
    complex(8),dimension(:,:,:),intent(inout) :: Gk    !as Sigma
    character(len=*)                          :: axis
    integer,intent(in)                        :: Nsites,Ncell
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable   :: csi
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
    Ntot  = size(Hk,1)
    Lfreq = size(Sigma,3)
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop "get_gk_normal_tridiag ERROR: passed Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Ntot],'get_gk_normal_tridiag',"Hk")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],'get_gk_normal_tridiag',"Sigma")
    call assert_shape(Gk,[Ntot,Ntot,Lfreq],'get_gk_normal_tridiag',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(csi(Ntot,Ntot,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i) = (wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
    enddo
    !
    Gk=zero
    call invert_gk_normal_tridiag_mpi(csi,Hk,Gk,Nsites,Ncell)
  end subroutine get_gk_normal_tridiag


  !DOS case: Hk--> Es,DOSs
  subroutine get_gk_normal_dos(Ebands,Dbands,Hloc,Gk,Sigma,axis)
    real(8),dimension(:),intent(in)            :: Ebands    ![N]
    real(8),dimension(:),intent(in)            :: Dbands    ![N]
    real(8),dimension(size(Ebands)),intent(in) :: Hloc      ![N]
    complex(8),dimension(:,:,:),intent(in)     :: Sigma     ![N,N][Lmats]
    complex(8),dimension(:,:,:),intent(inout)  :: Gk      !as Sigma
    character(len=*)                           :: axis
    !
    complex(8),dimension(:,:,:),allocatable    :: csi ![N,N][Lmats]
    complex(8),dimension(:,:,:),allocatable    :: Gktmp      !as Sigma
    complex(8),dimension(:,:),allocatable      :: Gdos_tmp ![N][N]
    logical                                    :: dos_diag !1. T / 2. F
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
    Ntot  = size(Ebands)
    Lfreq = size(Sigma,3)
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    dos_diag = .not.(size(Dbands) < size(Ebands))
    !
    !Testing part:
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"dmft_get_gk_normal_dos","Sigma")
    call assert_shape(Gk,[Ntot,Ntot,Lfreq],"dmft_get_gk_normal_dos","Gk")
    !
    call build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Ntot,Ntot,Lfreq))
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
    enddo
    !
    !invert (Z-Hk) for each k-point
    allocate(Gktmp(Ntot,Ntot,Lfreq));Gktmp=zero
    !
    ! diagonal case
    if(dos_diag)then
       do i=1+mpi_rank, Lfreq, mpi_size
          do io=1,Ntot
             Gktmp(io,io,i) = Dbands(io)/( csi(io,io,i)-Hloc(io)-Ebands(io) )
          enddo
       enddo
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       do i = 1+mpi_rank, Lfreq, mpi_size
          Gdos_tmp = csi(:,:,i)-diag(Hloc(:))-diag(Ebands(:)) !G(e,w) = csi - Hloc - H(e)
          call inv(Gdos_tmp)
          Gktmp(:,:,i) = Dbands(1)*Gdos_tmp
       enddo
    end if
    !
#ifdef _MPI    
    if(check_MPI())then
       Gk=zero
       call Mpi_AllReduce(Gktmp,Gk, size(Gk), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gk=Gktmp
    endif
#else
    Gk=Gktmp
#endif
  end subroutine get_gk_normal_dos





















  subroutine get_gk_superc_main(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)                      :: Hk       ![2][N,N]
    complex(8),dimension(:,:,:,:),intent(in)                    :: Sigma    ![2][N,N,Lfreq]
    complex(8),dimension(:,:,:,:),intent(inout)                 :: Gk     !as Sigma
    character(len=*)                                            :: axis
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable                 :: csi     ![2][2][N][N][Lfreq]
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
    Ntot  = size(Hk,2)
    Lfreq = size(Sigma,4)
    !
    call assert_shape(Hk,[2,Ntot,Ntot],'get_gk_superc',"Hk")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gk_superc',"Sigma")
    call assert_shape(Gk,[2,Ntot,Ntot,Lfreq],'get_gk_superc',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(csi(2,2,Ntot,Ntot,Lfreq))
    !
    select case(axis)
    case default; stop "get_gk_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Ntot) -        Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                          -        Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                          -        Sigma(2,:,:,i)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Ntot) + conjg( Sigma(1,:,:,i))
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i) + xmu)*eye(Ntot)                - Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                                           - Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                                           - Sigma(2,:,:,i)
          csi(2,2,:,:,i) = -conjg(wfreq(Lfreq+1-i) + xmu)*eye(Ntot)  + conjg( Sigma(1,:,:,Lfreq+1-i) )
       enddo
    end select
    !
    Gk=zero
    call invert_gk_superc_mpi(csi,Hk,Gk)
  end subroutine get_gk_superc_main



  subroutine get_gk_superc_dos(Ebands,Dbands,Hloc,Gk,Sigma,axis)
    real(8),dimension(:,:),intent(in)              :: Ebands    ![2,N]
    real(8),dimension(:),intent(in)                :: Dbands    ![N]
    real(8),dimension(2,size(Ebands,2)),intent(in) :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:),intent(in)       :: Sigma     ![2,N,N,Lfreq]
    complex(8),dimension(:,:,:,:),intent(inout)    :: Gk      !as Sigma
    character(len=*)                               :: axis
    ! arrays
    complex(8),dimension(:,:,:,:,:),allocatable    :: csi ![2][2][N,N,Lfreq]
    complex(8),dimension(:,:),allocatable          :: Gmatrix 
    complex(8),dimension(:,:,:,:),allocatable      :: Gtmp    !as Sigma
    complex(8),dimension(:,:),allocatable          :: Gdos_tmp ![N][N]
    logical                                        :: dos_diag !1. T / 2. F
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
    Ntot  = size(Ebands,2)
    Lfreq = size(Sigma,4)
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[2,Ntot],'get_gk_superc_dos',"Ebands")
    call assert_shape(Dbands,[Ntot],'get_gk_superc_dos',"Dbands")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gk_superc_dos',"Sigma")
    call assert_shape(Gk,[2,Ntot,Ntot,Lfreq],'get_gk_superc_dos',"Gk")
    !
    call build_frequency_array(axis)
    !
    !
    allocate(csi(2,2,Ntot,Ntot,Lfreq));csi = zero
    select case(axis)
    case default; stop "get_gk_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Ntot) - diag(Hloc(1,:)) -        Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                                            -        Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                                            -        Sigma(2,:,:,i)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Ntot) + diag(Hloc(2,:)) + conjg( Sigma(1,:,:,i) )
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Ntot)                 - diag(Hloc(1,:)) -        Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                                                            -        Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                                                            -        Sigma(2,:,:,i)
          csi(2,2,:,:,i) = -conjg( wfreq(Lfreq+1-i)+xmu )*eye(Ntot) + diag(Hloc(2,:)) + conjg( Sigma(1,:,:,Lfreq+1-i) )
       enddo
    end select
    !
    !invert (Z-Hk) for each k-point
    allocate(Gtmp(2,Ntot,Ntot,Lfreq))
    allocate(Gmatrix(2*Ntot,2*Ntot))
    Gtmp=zero
    !
    do i=1+mpi_rank,Lfreq, mpi_size
       Gmatrix  = zero
       Gmatrix(1:Ntot,1:Ntot)               = csi(1,1,:,:,i) - diag(Ebands(1,:))
       Gmatrix(1:Ntot,Ntot+1:2*Ntot)        = csi(1,2,:,:,i)
       Gmatrix(Ntot+1:2*Ntot,1:Ntot)        = csi(2,1,:,:,i)
       Gmatrix(Ntot+1:2*Ntot,Ntot+1:2*Ntot) = csi(2,2,:,:,i) - diag(Ebands(2,:))
       call inv(Gmatrix)
       do io=1,Ntot
          Gtmp(1,io,io,i) = Gtmp(1,io,io,i) + Gmatrix(io,io)*Dbands(io)
          Gtmp(2,io,io,i) = Gtmp(2,io,io,i) + Gmatrix(io,Ntot+io)*Dbands(io)
       enddo
    enddo
    !
#ifdef _MPI    
    if(check_MPI())then
       Gk=zero
       call Mpi_AllReduce(Gtmp,Gk, size(Gk), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gk=Gtmp
    endif
#else
    Gk=Gtmp
#endif
  end subroutine get_gk_superc_dos














  !##################################################################
  !##################################################################
  !                       INTERFACES
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !                         NORMAL
  !##################################################################
  !##################################################################
  subroutine get_gk_normal_hk_rank4(Hk,Gk,Sigma,axis) !N=Nspin*Norb
    complex(8),dimension(:,:),intent(in)          :: Hk        ![N,N]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk        !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:,:,:),allocatable       :: SF,GF
    Ntot  = size(Hk,1)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Hk,[Ntot,Nso],"get_gk_normal_rank4","Hk")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank4","Sigma")
    call assert_shape(Gk, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank4","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    call get_gk_normal_main(Hk,GF,SF,axis)
    Gk = reshape_matrix_to_rank4(GF,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_hk_rank4

  subroutine get_gk_normal_dos_rank4(Ebands,Dbands,Hloc,Gk,Sigma,axis)!N=Nspin*Norb
    real(8),dimension(:),intent(in)             :: Ebands    ![N]
    real(8),dimension(:),intent(in)             :: Dbands    ![N]
    real(8),dimension(size(Ebands)),intent(in)  :: Hloc      ![N]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:,:,:),allocatable       :: SF,GF
    !
    Ntot  = size(Ebands,1)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank4","Sigma")
    call assert_shape(Gk, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank4","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    call get_gk_normal_dos(Ebands,Dbands,Hloc,GF,SF,axis)
    Gk = reshape_matrix_to_rank4(GF,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_dos_rank4

  subroutine get_gk_normal_hk_rank5(Hk,Gk,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:),intent(in)            :: Hk        ![N,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:),allocatable         :: SF,GF
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5","Sigma")
    call assert_shape(Gk, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gk_normal_main(Hk,GF,SF,axis)
    Gk = reshape_matrix_to_rank5(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_hk_rank5

  subroutine get_gk_normal_tridiag_rank5(Hk,Gk,Sigma,axis,Nsites,Ncell) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:),intent(in)            :: Hk        ![N,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                                :: axis
    integer,intent(in)                              :: Nsites,Ncell
    complex(8),dimension(:,:,:),allocatable         :: SF,GF
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5","Sigma")
    call assert_shape(Gk, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gk_normal_tridiag(Hk,GF,SF,axis,Nsites,Ncell)
    Gk = reshape_matrix_to_rank5(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_tridiag_rank5

  subroutine get_gk_normal_hk_rank5_6(Hk,Gk,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:),intent(in)              :: Hk        ![N,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     ![Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk        ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:),allocatable           :: SF,GF
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_rank5_6","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5_6","Sigma")
    call assert_shape(Gk, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5_6","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gk_normal_main(Hk,GF,SF,axis)
    Gk = reshape_matrix_to_rank6(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_hk_rank5_6


  subroutine get_gk_normal_hk_rank6(Hk,Gk,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:),intent(in)              :: Hk     ![N,N]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk   !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:),allocatable           :: SF,GF
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin]"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_rank6","Hk")
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank6","Sigma")
    call assert_shape(Gk, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank6","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank6_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gk_normal_main(Hk,GF,SF,axis)
    Gk = reshape_matrix_to_rank6(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_hk_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine get_gk_normal_hk_rank7(Hk,Gk,Sigma,axis) !N=Nineq*Nlat*Nspin*Norb
    complex(8),dimension(:,:),intent(in)                :: Hk     ![N,N]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk   !as Sigma
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:),allocatable             :: SF,GF
    Ntot  = size(Hk,1)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nlso  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin,Nineq]"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank7","Sigma")
    call assert_shape(Gk, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank7","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank7_to_matrix(Sigma,Nineq,Nlat,Nspin,Norb,Lfreq)
    call get_gk_normal_main(Hk,GF,SF,axis)
    Gk = reshape_matrix_to_rank7(GF,Nineq,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_hk_rank7

  subroutine get_gk_normal_tridiag_rank7(Hk,Gk,Sigma,axis,Nsites,Ncell) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:),intent(in)                :: Hk     ![N,N]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk   !as Sigma
    character(len=*)                                    :: axis
    integer,intent(in)                                  :: Nsites,Ncell
    complex(8),dimension(:,:,:),allocatable             :: SF,GF
    Ntot  = size(Hk,1)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nlso  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin,Nineq]"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank7","Sigma")
    call assert_shape(Gk, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank7","Gk")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank7_to_matrix(Sigma,Nineq,Nlat,Nspin,Norb,Lfreq)
    call get_gk_normal_tridiag(Hk,GF,SF,axis,Nsites,Ncell)
    Gk = reshape_matrix_to_rank7(GF,Nineq,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_normal_tridiag_rank7
#endif







  !##################################################################
  !##################################################################
  !                         SUPERC
  !##################################################################
  !##################################################################

  subroutine get_gk_superc_hk_rank4(Hk,Gk,Sigma,axis) !N=Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][N,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: SF,GF ![2][N,N,Lfreq]
    Ntot  = size(Hk,2)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Hk,[2,Ntot,Nso],"get_gk_superc_rank4","Hk")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank4","Sigma")
    call assert_shape(Gk, [2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank4","Gk")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank4_to_matrix(Sigma(1,:,:,:,:,:),Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank4_to_matrix(Sigma(2,:,:,:,:,:),Nspin,Norb,Lfreq)
    call get_gk_superc_main(Hk,GF,SF,axis)
    Gk(1,:,:,:,:,:) = reshape_matrix_to_rank4(GF(1,:,:,:),Nspin,Norb,Lfreq)
    Gk(2,:,:,:,:,:) = reshape_matrix_to_rank4(GF(2,:,:,:),Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_superc_hk_rank4

  subroutine get_gk_superc_dos_rank4(Ebands,Dbands,Hloc,Gk,Sigma,axis)!N=Nspin*Norb
    real(8),dimension(:,:),intent(in)               :: Ebands   ![2][N]
    real(8),dimension(:),intent(in)                 :: Dbands    ![N]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: SF,GF ![2][N,N,Lfreq]
    !
    Ntot  = size(Ebands,2)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank4","Sigma")
    call assert_shape(Gk, [2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank4","Gk")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank4_to_matrix(Sigma(1,:,:,:,:,:),Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank4_to_matrix(Sigma(2,:,:,:,:,:),Nspin,Norb,Lfreq)
    call get_gk_superc_dos(Ebands,Dbands,Hloc,GF,SF,axis)
    Gk(1,:,:,:,:,:) = reshape_matrix_to_rank4(GF(1,:,:,:),Nspin,Norb,Lfreq)
    Gk(2,:,:,:,:,:) = reshape_matrix_to_rank4(GF(2,:,:,:),Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_superc_dos_rank4


  subroutine get_gk_superc_hk_rank5(Hk,Gk,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2][N,N]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF ![2][N,N,Lfreq]
    Ntot  = size(Hk,2)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Hk,[2,Ntot,Nso],"get_gk_superc_rank5","Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank5","Sigma")
    call assert_shape(Gk, [2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank5","Gk")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank5_to_matrix(Sigma(1,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank5_to_matrix(Sigma(2,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    call get_gk_superc_main(Hk,GF,SF,axis)
    Gk(1,:,:,:,:,:,:) = reshape_matrix_to_rank5(GF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Gk(2,:,:,:,:,:,:) = reshape_matrix_to_rank5(GF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_superc_hk_rank5


  subroutine get_gk_superc_hk_rank5_6(Hk,Gk,Fk,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2,N,N]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk        ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Fk        ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][Lfreq]
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF
    Ntot  = size(Hk,2)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get k-dependent Green's function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[2,Ntot,Nlso],"get_gk_superc_rank5_6","Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank5_6","Sigma")
    call assert_shape(Gk, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank5_6","Gk")
    call assert_shape(Fk, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_rank5_6","Floc")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank5_to_matrix(Sigma(1,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank5_to_matrix(Sigma(2,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    call get_gk_superc_main(Hk,GF,SF,axis)
    Gk = reshape_matrix_to_rank6(GF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Fk = reshape_matrix_to_rank6(GF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gk_superc_hk_rank5_6






end module LEGACY_GF_GK
