module GF_COMMON
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,txtfy,splot,file_gzip,str
  USE SF_LINALG,    only: eye,inv,inv_sym,inv_tridiag,get_tridiag,diag
  USE SF_ARRAYS,    only: linspace,arange
  USE SF_MISC,      only: assert_shape
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  !
  USE DMFT_CTRL_VARS
  implicit none




  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block

  interface reshape_matrix_to_rank4
     module procedure :: d_matrix_TO_rank4
     module procedure :: c_matrix_TO_rank4
     module procedure :: d_matrixL_TO_rank4L
     module procedure :: c_matrixL_TO_rank4L
  end interface reshape_matrix_to_rank4

  interface reshape_rank4_to_matrix
     module procedure :: d_rank4_TO_matrix
     module procedure :: c_rank4_TO_matrix
     module procedure :: d_rank4L_TO_matrixL
     module procedure :: c_rank4L_TO_matrixL
  end interface reshape_rank4_to_matrix

  interface reshape_matrix_to_rank5
     module procedure :: d_matrix_TO_rank5
     module procedure :: c_matrix_TO_rank5
     module procedure :: d_matrixL_TO_rank5L
     module procedure :: c_matrixL_TO_rank5L
  end interface reshape_matrix_to_rank5

  interface reshape_rank5_to_matrix
     module procedure :: d_rank5_TO_matrix
     module procedure :: c_rank5_TO_matrix
     module procedure :: d_rank5L_TO_matrixL
     module procedure :: c_rank5L_TO_matrixL
  end interface reshape_rank5_to_matrix

  interface reshape_matrix_to_rank6
     module procedure :: d_matrix_TO_rank6
     module procedure :: c_matrix_TO_rank6
     module procedure :: d_matrixL_TO_rank6L
     module procedure :: c_matrixL_TO_rank6L
  end interface reshape_matrix_to_rank6

  interface reshape_rank6_to_matrix
     module procedure :: d_rank6_TO_matrix
     module procedure :: c_rank6_TO_matrix
     module procedure :: d_rank6L_TO_matrixL
     module procedure :: c_rank6L_TO_matrixL
  end interface reshape_rank6_to_matrix



  interface reshape_matrix_to_rank7
     module procedure :: d_matrix_TO_rank7
     module procedure :: c_matrix_TO_rank7
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: d_matrixL_TO_rank7L
     module procedure :: c_matrixL_TO_rank7L
#endif
  end interface reshape_matrix_to_rank7

  interface reshape_rank7_to_matrix
     module procedure :: d_rank7_TO_matrix
     module procedure :: c_rank7_TO_matrix
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: d_rank7L_TO_matrixL
     module procedure :: c_rank7L_TO_matrixL
#endif
  end interface reshape_rank7_to_matrix

  integer                                   :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats,Nineq,Nilso,Ntot
  integer                                   :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js,iineq
  !
  integer                                   :: mpi_ierr
  integer                                   :: mpi_rank
  integer                                   :: mpi_size
  logical                                   :: mpi_master
  !
  real(8)                                   :: beta
  real(8)                                   :: xmu,eps
  real(8)                                   :: wini,wfin 
  !
  real(8),dimension(:),allocatable    :: wm !Matsubara frequencies
  real(8),dimension(:),allocatable    :: wr !Real frequencies
  integer                             :: Lfreq
  complex(8),dimension(:),allocatable :: wfreq
  logical                             :: pushed=.false.


contains

  !TO BE MOVED DOWN:
  subroutine gf_push_zeta(zeta)
    complex(8),dimension(:) :: zeta
    Lfreq=size(zeta)
    if(allocated(wfreq))deallocate(wfreq)
    allocate(wfreq(Lfreq))
    wfreq=zeta
    pushed=.true.
  end subroutine gf_push_zeta

  subroutine build_frequency_array(axis)
    character(len=*) :: axis
    if(pushed)then
       if(size(wfreq)/=Lfreq)stop "build_frequency_array ERROR: pushed wfreq has wrong size"
    else
       if(allocated(wfreq))deallocate(wfreq)
       allocate(wfreq(Lfreq))
       select case(axis)
       case default;
          stop "build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
       case("matsubara","mats","Mats","Matsubara","M","m")
          call get_ctrl_var(beta,"BETA")
          wfreq = dcmplx(0d0,pi/beta*(2*arange(1,Lfreq)-1))
       case("realaxis","real","Realaxis","Real","R","r")
          call get_ctrl_var(wini,"WINI")
          call get_ctrl_var(wfin,"WFIN")
          call get_ctrl_var(eps,"EPS")
          wfreq = dcmplx(linspace(wini,wfin,Lfreq),eps)
       end select
    endif
    return
  end subroutine build_frequency_array







  ! INVERT_GK_NORMAL(_MPI)
  !SERIAL (OR PARALLEL ON K):
  subroutine invert_gk_normal(zeta,Hk,Gkout)
    complex(8),dimension(:,:,:),intent(in)    :: zeta    ![N,N,Lfreq]
    complex(8),dimension(:,:),intent(in)      :: Hk      ![N,N]
    complex(8),dimension(:,:,:),intent(inout) :: Gkout   ![N,N,Lfreq]
    complex(8)                                :: Gmatrix(size(Hk,1),size(Hk,2))
    integer                                   :: Ntot,Lfreq
    integer                                   :: i
    !
    Ntot  = size(Hk,1)
    Lfreq = size(zeta,3)
    do i=1,Lfreq
       Gmatrix  = zeta(:,:,i) - Hk
       call inv(Gmatrix)
       Gkout(:,:,i) = Gmatrix
    enddo
  end subroutine invert_gk_normal

  !PARALLEL ON FREQ:
  subroutine invert_gk_normal_mpi(zeta,Hk,Gkout)
    complex(8),dimension(:,:,:),intent(in)    :: zeta    ![N,N,Lfreq]
    complex(8),dimension(:,:),intent(in)      :: Hk      ![N,N]
    complex(8),dimension(:,:,:),intent(inout) :: Gkout   ![N,N,Lfreq]
    complex(8)                                :: Gktmp(size(Hk,1),size(Hk,2),size(zeta,3))   !as Gkout
    complex(8)                                :: Gmatrix(size(Hk,1),size(Hk,2))
    integer                                   :: Ntot,Lfreq
    integer                                   :: i
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
    Ntot  = size(Hk,1)
    Lfreq = size(zeta,3)
    !
    Gktmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = zeta(:,:,i) - Hk
       call inv(Gmatrix)
       Gktmp(:,:,i) = Gmatrix
    enddo
    !
#ifdef _MPI    
    if(check_MPI())then
       Gkout=zero
       call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Gkout=Gktmp
    endif
#else
    Gkout=Gktmp
#endif
  end subroutine invert_gk_normal_mpi




  !
  ! INVERT_GK_NORMAL_TRIDIAG(_MPI)
  !SERIAL (OR PARALLEL ON K)
  subroutine invert_gk_normal_tridiag(zeta,Hk,Gkout,Nsites,Ncell)
    complex(8),dimension(:,:,:),intent(in)     :: zeta    ![N,N,Lfreq]
    complex(8),dimension(:,:),intent(in)       :: Hk      ![N,N]
    complex(8),dimension(:,:,:),intent(inout)  :: Gkout   ![N,N,Lfreq]
    integer                                    :: Nsites,Ncell
    complex(8),dimension(Nsites,  Ncell,Ncell) :: D
    complex(8),dimension(Nsites-1,Ncell,Ncell) :: Sub
    complex(8),dimension(Nsites-1,Ncell,Ncell) :: Over
    complex(8),dimension(Nsites,Ncell,Ncell)   :: Gmatrix(Nsites,Ncell,Ncell)
    integer                                    :: Ntot,Lfreq
    integer                                    :: i,isite,ic,jc,io,jo
    !
    Ntot  = size(Hk,1)
    Lfreq = size(zeta,3)
    !
    Gkout=zero
    do i=1,Lfreq      
       call get_tridiag(Nsites,Ncell,(zeta(:,:,i)-Hk),Sub,D,Over)
       call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
       do isite=1,Nsites
          do ic=1,Ncell
             do jc=1,Ncell
                io = ic + (isite-1)*Ncell
                jo = jc + (isite-1)*Ncell
                Gkout(io,jo,i) = Gmatrix(isite,ic,jc)
             enddo
          enddo
       enddo
    enddo
  end subroutine invert_gk_normal_tridiag

  !PARALLEL ON FREQ:
  subroutine invert_gk_normal_tridiag_mpi(zeta,Hk,Gkout,Nsites,Ncell)
    complex(8),dimension(:,:,:),intent(in)     :: zeta    ![N,N,Lfreq]
    complex(8),dimension(:,:),intent(in)       :: Hk      ![N,N]
    complex(8),dimension(:,:,:),intent(inout)  :: Gkout   ![N,N,Lfreq]
    integer                                    :: Nsites,Ncell
    complex(8),dimension(Nsites,  Ncell,Ncell) :: D
    complex(8),dimension(Nsites-1,Ncell,Ncell) :: Sub
    complex(8),dimension(Nsites-1,Ncell,Ncell) :: Over
    complex(8)                                 :: Gktmp(size(Hk,1),size(Hk,2),size(zeta,3))   !as Gkout
    complex(8),dimension(Nsites,Ncell,Ncell)   :: Gmatrix(Nsites,Ncell,Ncell)
    integer                                    :: Ntot,Lfreq
    integer                                    :: i,isite,ic,jc,io,jo
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
    Lfreq = size(zeta,3)
    !
    Gktmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       call get_tridiag(Nsites,Ncell,(zeta(:,:,i)-Hk),Sub,D,Over)
       call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
       do isite=1,Nsites
          do ic=1,Ncell
             do jc=1,Ncell
                io = ic + (isite-1)*Ncell
                jo = jc + (isite-1)*Ncell
                Gktmp(io,jo,i) = Gmatrix(isite,ic,jc)
             enddo
          enddo
       enddo
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Gkout=zero
       call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Gkout=Gktmp
    endif
#else
    Gkout=Gktmp
#endif
  end subroutine invert_gk_normal_tridiag_mpi










  !
  ! INVERT_GK_SUPERC(_MPI)
  !
  subroutine invert_gk_superc(zeta,Hk,Gkout)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: zeta    ![2][2][N][N][Lfreq]
    complex(8),dimension(:,:,:),intent(in)          :: Hk      ![2][N][N]
    complex(8),dimension(:,:,:,:),intent(inout)     :: Gkout   ![2][N][N][Lfreq]
    complex(8),dimension(2*size(Hk,2),2*size(Hk,2)) :: Gmatrix ![2*N][2*N]
    integer                                         :: Ntot,Lfreq
    integer                                         :: i
    !
    Ntot = size(Hk,2)
    Lfreq = size(zeta,5)
    !
    Gkout = zero
    do i=1,Lfreq
       Gmatrix  = zero
       Gmatrix(1:Ntot,1:Ntot)               = zeta(1,1,:,:,i) - Hk(1,:,:)
       Gmatrix(1:Ntot,Ntot+1:2*Ntot)        = zeta(1,2,:,:,i)
       Gmatrix(Ntot+1:2*Ntot,1:Ntot)        = zeta(2,1,:,:,i)
       Gmatrix(Ntot+1:2*Ntot,Ntot+1:2*Ntot) = zeta(2,2,:,:,i) - Hk(2,:,:)
       call inv(Gmatrix)
       Gkout(1,:,:,i) = Gmatrix(1:Ntot,1:Ntot)
       Gkout(2,:,:,i) = Gmatrix(1:Ntot,Ntot+1:2*Ntot)
    enddo
  end subroutine invert_gk_superc

  subroutine invert_gk_superc_mpi(zeta,Hk,Gkout)
    complex(8),dimension(:,:,:,:,:),intent(in)                 :: zeta    ![2][2][N][N][Lfreq]
    complex(8),dimension(:,:,:),intent(in)                     :: Hk      ![2][N][N]
    complex(8),dimension(:,:,:,:),intent(inout)                :: Gkout   ![2][N][N][Lfreq]
    complex(8),dimension(2,size(Hk,2),size(Hk,2),size(zeta,5)) :: Gktmp   ![2][N][N][Lfreq]
    complex(8),dimension(2*size(Hk,2),2*size(Hk,2))            :: Gmatrix ![2*N][2*N]
    integer                                                    :: Nspin,Norb,Nso,Lfreq
    integer                                                    :: i,iorb,jorb,ispin,jspin,io,jo
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
    Ntot = size(Hk,2)
    Lfreq = size(zeta,5)
    !
    Gkout = zero
    Gktmp = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = zero
       Gmatrix(1:Ntot,1:Ntot)               = zeta(1,1,:,:,i) - Hk(1,:,:)
       Gmatrix(1:Ntot,Ntot+1:2*Ntot)        = zeta(1,2,:,:,i)
       Gmatrix(Ntot+1:2*Ntot,1:Ntot)        = zeta(2,1,:,:,i)
       Gmatrix(Ntot+1:2*Ntot,Ntot+1:2*Ntot) = zeta(2,2,:,:,i) - Hk(2,:,:)
       call inv(Gmatrix)
       Gktmp(1,:,:,i) = Gmatrix(1:Ntot,1:Ntot)
       Gktmp(2,:,:,i) = Gmatrix(1:Ntot,Ntot+1:2*Ntot)
    enddo
#ifdef _MPI
    if(check_MPI())then
       Gkout=zero
       call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Gkout=Gktmp
    endif
#else
    Gkout=Gktmp
#endif
  end subroutine invert_gk_superc_mpi











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





  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !MATRIX --> RANK4_7: DBLE
  function d_matrix_TO_rank4(Hmat,Nspin,Norb) result(Hrank)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hmat
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hrank
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hrank(ispin,jspin,iorb,jorb) = Hmat(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_matrix_TO_rank4


  function d_matrix_TO_rank5(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hrank
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hrank(ilat,ispin,jspin,iorb,jorb) = Hmat(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_matrix_TO_rank5


  function d_matrix_TO_rank6(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hrank
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hrank(ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_matrix_TO_rank6


  function d_matrix_TO_rank7(Hmat,Nineq,Nlat,Nspin,Norb) result(Hrank)
    integer                                                        :: Nineq,Nlat,Nspin,Norb
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb) :: Hmat
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: Hrank
    integer                                                        :: iineq
    integer                                                        :: ilat,jlat
    integer                                                        :: iorb,jorb
    integer                                                        :: ispin,jspin
    integer                                                        :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_matrix_TO_rank7



  !MATRIX --> RANK4_7: CMPLX
  function c_matrix_TO_rank4(Hmat,Nspin,Norb) result(Hrank)
    integer                                  :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hmat
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hrank
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hrank(ispin,jspin,iorb,jorb) = Hmat(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_matrix_TO_rank4

  function c_matrix_TO_rank5(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hrank
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hrank(ilat,ispin,jspin,iorb,jorb) = Hmat(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_matrix_TO_rank5

  function c_matrix_TO_rank6(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hrank
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hrank(ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_matrix_TO_rank6



  function c_matrix_TO_rank7(Hmat,Nineq,Nlat,Nspin,Norb) result(Hrank)
    integer                                                        :: Nineq,Nlat,Nspin,Norb
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb) :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: Hrank
    integer                                                        :: iineq
    integer                                                        :: ilat,jlat
    integer                                                        :: iorb,jorb
    integer                                                        :: ispin,jspin
    integer                                                        :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_matrix_TO_rank7





  !##################################################################
  !##################################################################
  !##################################################################




  function d_matrixL_TO_rank4L(Hmat,Nspin,Norb,L) result(Hrank)
    integer                                    :: Nspin,Norb,L
    real(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hmat
    real(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hrank
    integer                                    :: iorb,ispin,is
    integer                                    :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hrank(ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
             enddo
          enddo
       enddo
    enddo
  end function d_matrixL_TO_rank4L

  function d_matrixL_TO_rank5L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                              :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hrank
    integer                                              :: iorb,ispin,ilat,is
    integer                                              :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   Hrank(ilat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_matrixL_TO_rank5L

  function d_matrixL_TO_rank6L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                              :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L) :: Hrank
    integer                                              :: ilat,jlat
    integer                                              :: iorb,jorb
    integer                                              :: ispin,jspin
    integer                                              :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_matrixL_TO_rank6L

#if __GFORTRAN__ &&  __GNUC__ > 8
  function d_matrixL_TO_rank7L(Hmat,Nineq,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_matrixL_TO_rank7L
#endif

  function c_matrixL_TO_rank4L(Hmat,Nspin,Norb,L) result(Hrank)
    integer                                    :: Nspin,Norb,L
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hmat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hrank
    integer                                    :: iorb,ispin,is
    integer                                    :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hrank(ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
             enddo
          enddo
       enddo
    enddo
  end function c_matrixL_TO_rank4L

  function c_matrixL_TO_rank5L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                              :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hrank
    integer                                              :: iorb,ispin,ilat,is
    integer                                              :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   Hrank(ilat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_matrixL_TO_rank5L

  function c_matrixL_TO_rank6L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                 :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L) :: Hrank
    integer                                                 :: ilat,jlat
    integer                                                 :: iorb,jorb
    integer                                                 :: ispin,jspin
    integer                                                 :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_matrixL_TO_rank6L

#if __GFORTRAN__ &&  __GNUC__ > 8
  function c_matrixL_TO_rank7L(Hmat,Nineq,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_matrixL_TO_rank7L
#endif






  !##################################################################
  !##################################################################
  !##################################################################

  !RANK4_7 --> MATRIX

  !##################################################################
  !##################################################################
  !##################################################################

  function d_rank4_TO_matrix(Hrank,Nspin,Norb) result(Hmat)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hrank
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hmat
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hmat(is,js) = Hrank(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_rank4_TO_matrix


  function d_rank5_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hmat(is,js) = Hrank(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_rank5_TO_matrix

  function d_rank6_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hmat(is,js) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_rank6_TO_matrix



  function d_rank7_TO_matrix(Hrank,Nineq,Nlat,Nspin,Norb) result(Hmat)
    integer                                                        :: Nineq,Nlat,Nspin,Norb
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: Hrank
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb) :: Hmat
    integer                                                        :: iineq
    integer                                                        :: ilat,jlat
    integer                                                        :: iorb,jorb
    integer                                                        :: ispin,jspin
    integer                                                        :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hmat(is,js) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_rank7_TO_matrix



  !COMPLEX
  function c_rank4_TO_matrix(Hrank,Nspin,Norb) result(Hmat)
    integer                                  :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hmat
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hmat(is,js) = Hrank(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_rank4_TO_matrix

  function c_rank5_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hmat(is,js) = Hrank(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_rank5_TO_matrix

  function c_rank6_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hmat
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hmat(is,js) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_rank6_TO_matrix

  function c_rank7_TO_matrix(Hrank,Nineq,Nlat,Nspin,Norb) result(Hmat)
    integer                                                        :: Nineq,Nlat,Nspin,Norb
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: Hrank
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb) :: Hmat
    integer                                                        :: iineq
    integer                                                        :: ilat,jlat
    integer                                                        :: iorb,jorb
    integer                                                        :: ispin,jspin
    integer                                                        :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hmat(is,js) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_rank7_TO_matrix






  !##################################################################
  !##################################################################
  !##################################################################




  function d_rank4L_TO_matrixL(Hrank,Nspin,Norb,L) result(Hmat)
    integer                                    :: Nspin,Norb,L
    real(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hrank
    real(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hmat
    integer                                    :: iorb,ispin,is
    integer                                    :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hmat(is,js,:) = Hrank(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end function d_rank4L_TO_matrixL

  function d_rank5L_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                              :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    integer                                              :: iorb,ispin,ilat,is
    integer                                              :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   Hmat(is,js,:) = Hrank(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_rank5L_TO_matrixL

  function d_rank6L_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                              :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L) :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    integer                                              :: ilat,jlat
    integer                                              :: iorb,jorb
    integer                                              :: ispin,jspin
    integer                                              :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hmat(is,js,:) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_rank6L_TO_matrixL

#if __GFORTRAN__ &&  __GNUC__ > 8
  function d_rank7L_TO_matrixL(Hrank,Nineq,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hmat(is,js,:) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) 
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_rank7L_TO_matrixL
#endif

  function c_rank4L_TO_matrixL(Hrank,Nspin,Norb,L) result(Hmat)
    integer                                    :: Nspin,Norb,L
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hmat
    integer                                    :: iorb,ispin,is
    integer                                    :: jorb,jspin,js
    Hrank=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hmat(is,js,:) = Hrank(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
  end function c_rank4L_TO_matrixL

  function c_rank5L_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                              :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    integer                                              :: iorb,ispin,ilat,is
    integer                                              :: jorb,jspin,jlat,js
    Hrank=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
                   Hmat(is,js,:) = Hrank(ilat,ispin,jspin,iorb,jorb,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_rank5L_TO_matrixL

  function c_rank6L_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                              :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L) :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hmat
    integer                                              :: ilat,jlat
    integer                                              :: iorb,jorb
    integer                                              :: ispin,jspin
    integer                                              :: is,js
    Hrank=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      Hmat(is,js,:) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_rank6L_TO_matrixL

#if __GFORTRAN__ &&  __GNUC__ > 8
  function c_rank7L_TO_matrixL(Hrank,Nineq,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do iineq=1,Nineq
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
                         Hmat(is,js,:) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) 
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_rank7L_TO_matrixL
#endif








  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
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






  !--------------------------------------------------------------------!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !--------------------------------------------------------------------!
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






end module GF_COMMON













! !TO BE REMOVED:



! !+-----------------------------------------------------------------------------+!
! !PURPOSE: 
! ! reshape a matrix from the [Nlso][Nlso] shape
! ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
! ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
! ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
! !+-----------------------------------------------------------------------------+!
! function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
!   real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   integer                                            :: Nlat,Nspin,Norb
!   real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
!   integer                                            :: iorb,ispin,ilat,is
!   integer                                            :: jorb,jspin,jlat,js
!   Hnnn=zero
!   do ilat=1,Nlat
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function d_nlso2nnn
! function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
!   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   integer                                               :: Nlat,Nspin,Norb
!   complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
!   integer                                               :: iorb,ispin,ilat,is
!   integer                                               :: jorb,jspin,jlat,js
!   Hnnn=zero
!   do ilat=1,Nlat
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function c_nlso2nnn

! function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
!   real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
!   integer                                  :: Nspin,Norb
!   real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
!   integer                                  :: iorb,ispin,is
!   integer                                  :: jorb,jspin,js
!   Hnn=zero
!   do ispin=1,Nspin
!      do jspin=1,Nspin
!         do iorb=1,Norb
!            do jorb=1,Norb
!               is = iorb + (ispin-1)*Norb  !spin-orbit stride
!               js = jorb + (jspin-1)*Norb  !spin-orbit stride
!               Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
!            enddo
!         enddo
!      enddo
!   enddo
! end function d_nso2nn
! function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
!   complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
!   integer                                     :: Nspin,Norb
!   complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
!   integer                                     :: iorb,ispin,is
!   integer                                     :: jorb,jspin,js
!   Hnn=zero
!   do ispin=1,Nspin
!      do jspin=1,Nspin
!         do iorb=1,Norb
!            do jorb=1,Norb
!               is = iorb + (ispin-1)*Norb  !spin-orbit stride
!               js = jorb + (jspin-1)*Norb  !spin-orbit stride
!               Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
!            enddo
!         enddo
!      enddo
!   enddo
! end function c_nso2nn




! !+-----------------------------------------------------------------------------+!
! !PURPOSE: 
! ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
! ! from/to the [Nlso][Nlso] shape.
! ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
! ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
! !+-----------------------------------------------------------------------------+!
! function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
!   real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
!   integer                                            :: Nlat,Nspin,Norb
!   real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   integer                                            :: iorb,ispin,ilat,is
!   integer                                            :: jorb,jspin,jlat,js
!   Hlso=zero
!   do ilat=1,Nlat
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function d_nnn2nlso

! function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
!   complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
!   integer                                               :: Nlat,Nspin,Norb
!   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   integer                                               :: iorb,ispin,ilat,is
!   integer                                               :: jorb,jspin,jlat,js
!   Hlso=zero
!   do ilat=1,Nlat
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
!                  Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function c_nnn2nlso

! function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
!   real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
!   integer                                  :: Nspin,Norb
!   real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
!   integer                                  :: iorb,ispin,is
!   integer                                  :: jorb,jspin,js
!   Hso=zero
!   do ispin=1,Nspin
!      do jspin=1,Nspin
!         do iorb=1,Norb
!            do jorb=1,Norb
!               is = iorb + (ispin-1)*Norb  !spin-orbit stride
!               js = jorb + (jspin-1)*Norb  !spin-orbit stride
!               Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
!            enddo
!         enddo
!      enddo
!   enddo
! end function d_nn2nso

! function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
!   complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
!   integer                                     :: Nspin,Norb
!   complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
!   integer                                     :: iorb,ispin,is
!   integer                                     :: jorb,jspin,js
!   Hso=zero
!   do ispin=1,Nspin
!      do jspin=1,Nspin
!         do iorb=1,Norb
!            do jorb=1,Norb
!               is = iorb + (ispin-1)*Norb  !spin-orbit stride
!               js = jorb + (jspin-1)*Norb  !spin-orbit stride
!               Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
!            enddo
!         enddo
!      enddo
!   enddo
! end function c_nn2nso



! !+-----------------------------------------------------------------------------+!
! !PURPOSE: 
! ! reshape a matrix from the [Nlso][Nlso] shape
! ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nlat][Nspin][Nspin][Norb][Norb]  !
! !+-----------------------------------------------------------------------------+!

! function d_nlso2nnn_cluster(Hlso,Nlat,Nspin,Norb) result(Hnnn)
!   integer                                            :: Nlat,Nspin,Norb
!   real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
!   integer                                            :: ilat,jlat
!   integer                                            :: iorb,jorb
!   integer                                            :: ispin,jspin
!   integer                                            :: is,js
!   Hnnn=zero
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
!                     js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
!                     Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function d_nlso2nnn_cluster
! !
! function c_nlso2nnn_cluster(Hlso,Nlat,Nspin,Norb) result(Hnnn)
!   integer                                               :: Nlat,Nspin,Norb
!   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
!   integer                                               :: ilat,jlat
!   integer                                               :: iorb,jorb
!   integer                                               :: ispin,jspin
!   integer                                               :: is,js
!   Hnnn=zero
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
!                     js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
!                     Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function c_nlso2nnn_cluster


! !+-----------------------------------------------------------------------------+!
! !PURPOSE: 
! ! reshape a matrix from the [Nlat][Nlat][Nspin][Nspin][Norb][Norb] shape
! ! _nnn2nlso : from [Nlat][Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
! !+-----------------------------------------------------------------------------+!

! function d_nnn2nlso_cluster(Hnnn,Nlat,Nspin,Norb) result(Hlso)
!   integer                                            :: Nlat,Nspin,Norb
!   real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
!   real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   integer                                            :: ilat,jlat
!   integer                                            :: iorb,jorb
!   integer                                            :: ispin,jspin
!   integer                                            :: is,js
!   Hlso=zero
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
!                     js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
!                     Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function d_nnn2nlso_cluster
! !
! function c_nnn2nlso_cluster(Hnnn,Nlat,Nspin,Norb) result(Hlso)
!   integer                                               :: Nlat,Nspin,Norb
!   complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
!   complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
!   integer                                               :: ilat,jlat
!   integer                                               :: iorb,jorb
!   integer                                               :: ispin,jspin
!   integer                                               :: is,js
!   Hlso=zero
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
!                     js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
!                     Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function c_nnn2nlso_cluster












!   ! INVERT_GK_NORMAL_CLUSTER(_MPI)
!   !
!   !SERIAL (OR PARALLEL ON K):
!   subroutine invert_gk_normal_cluster(zeta,Hk,hk_symm,Gkout)
!     complex(8),dimension(:,:,:),intent(in)            :: zeta    ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     logical,intent(in)                                :: hk_symm
!     complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:),allocatable             :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     integer                                           :: Nspin,Norb,Nlso,Lfreq
!     integer                                           :: i,iorb,jorb,ispin,jspin,io,jo
!     !
!     Nlat  = size(Gkout,1)
!     Nspin = size(Gkout,3)
!     Norb  = size(Gkout,5)
!     Lfreq = size(zeta,3)
!     Nlso  = Nlat*Nspin*Norb
!     !testing
!     call assert_shape(zeta,[Nlso,Nlso,Lfreq],"invert_gk_normal_cluster","zeta")
!     call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_cluster","Hk")
!     call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_cluster","Gkout")
!     !
!     allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
!     allocate(Gmatrix(Nlso,Nlso))
!     Gktmp=zero
!     do i=1,Lfreq
!        Gmatrix  = zeta(:,:,i) - Hk
!        if(hk_symm) then
!           call inv_sym(Gmatrix)
!        else
!           call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
!        end if
!        !store the diagonal blocks directly into the tmp output 
!        Gktmp(:,:,:,:,:,:,i)=lso2nnn_cluster_reshape(Gmatrix,Nlat,Nspin,Norb)
!     enddo
!     Gkout = Gktmp
!   end subroutine invert_gk_normal_cluster

!   !PARALLEL ON FREQ:
!   subroutine invert_gk_normal_cluster_mpi(zeta,Hk,hk_symm,Gkout)
!     complex(8),dimension(:,:,:),intent(in)            :: zeta    ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     logical,intent(in)                                :: hk_symm
!     complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:),allocatable             :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     integer                                           :: Nspin,Norb,Nlso,Lfreq
!     integer                                           :: i,iorb,jorb,ispin,jspin,io,jo
!     !
!     !
!     !MPI setup:
! #ifdef _MPI    
!     if(check_MPI())then
!        mpi_size  = get_size_MPI()
!        mpi_rank =  get_rank_MPI()
!        mpi_master= get_master_MPI()
!     else
!        mpi_size=1
!        mpi_rank=0
!        mpi_master=.true.
!     endif
! #else
!     mpi_size=1
!     mpi_rank=0
!     mpi_master=.true.
! #endif
!     !
!     Nlat  = size(Gkout,1)
!     Nspin = size(Gkout,3)
!     Norb  = size(Gkout,5)
!     Lfreq = size(zeta,3)
!     Nlso   = Nlat*Nspin*Norb
!     !testing
!     call assert_shape(zeta,[Nlso,Nlso,Lfreq],"invert_gk_normal_mpi","zeta")
!     call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_mpi","Hk")
!     call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_mpi","Gkout")
!     !
!     allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
!     allocate(Gmatrix(Nlso,Nlso))
!     Gktmp=zero
!     do i=1+mpi_rank,Lfreq,mpi_size
!        Gmatrix  = zeta(:,:,i) - Hk
!        if(hk_symm) then
!           call inv_sym(Gmatrix)
!        else
!           call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
!        end if
!        !store the diagonal blocks directly into the tmp output 
!        Gktmp(:,:,:,:,:,:,i)=lso2nnn_cluster_reshape(Gmatrix,Nlat,Nspin,Norb)
!     enddo
! #ifdef _MPI    
!     if(check_MPI())then
!        Gkout=zero
!        call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
!     else
!        Gkout=Gktmp
!     endif
! #else
!     Gkout=Gktmp
! #endif
!   end subroutine invert_gk_normal_cluster_mpi





!   !
!   ! INVERT_GK_NORMAL_INEQ(_MPI)
!   !
!   !SERIAL (OR PARALLEL ON K)
!   subroutine invert_gk_normal_ineq(zeta,Hk,hk_symm,Gkout)
!     complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     logical,intent(in)                              :: hk_symm
!     complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     !
!     complex(8),dimension(:,:,:),allocatable         :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:),allocatable           :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
!     integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
!     !
!     Nlat  = size(zeta,1)
!     Nspin = size(Gkout,2)
!     Norb  = size(Gkout,4)
!     Lfreq = size(zeta,4)
!     Nso   = Nspin*Norb
!     Nlso  = Nlat*Nspin*Norb
!     !Testing:
!     call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_ineq","zeta")
!     call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_ineq","Hk")
!     call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_ineq","Gkout")
!     !
!     if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
!          stop "invert_gk_normal_ineq: lattice_Gamma_mats & lattice_Gamma_real both allocated"
!     if(allocated(lattice_gamma_mats))then
!        call assert_shape(lattice_gamma_mats,[Nlso,Nlso,Lfreq],"invert_gk_normal_ineq","lattice_gamma_mats")
!        allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_mats
!     endif
!     if(allocated(lattice_gamma_real))then
!        call assert_shape(lattice_gamma_real,[Nlso,Nlso,Lfreq],"invert_gk_normal_ineq","lattice_gamma_real")
!        allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_real
!     endif
!     !
!     allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
!     allocate(Gmatrix(Nlso,Nlso))
!     Gktmp=zero
!     do i=1,Lfreq
!        Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nlat,Nso) - Hk
!        if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
!        if(hk_symm) then
!           call inv_sym(Gmatrix)
!        else
!           call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
!        end if
!        !store the diagonal blocks directly into the tmp output 
!        do ilat=1,Nlat
!           do ispin=1,Nspin
!              do jspin=1,Nspin
!                 do iorb=1,Norb
!                    do jorb=1,Norb
!                       io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
!                       jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
!                       Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
!     Gkout = Gktmp
!   end subroutine invert_gk_normal_ineq

!   !PARALLEL ON FREQ:
!   subroutine invert_gk_normal_ineq_mpi(zeta,Hk,hk_symm,Gkout)
!     complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     logical,intent(in)                              :: hk_symm
!     complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     !
!     complex(8),dimension(:,:,:),allocatable         :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:),allocatable           :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
!     integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
!     !
!     !MPI setup:
! #ifdef _MPI    
!     if(check_MPI())then
!        mpi_size  = get_size_MPI()
!        mpi_rank =  get_rank_MPI()
!        mpi_master= get_master_MPI()
!     else
!        mpi_size=1
!        mpi_rank=0
!        mpi_master=.true.
!     endif
! #else
!     mpi_size=1
!     mpi_rank=0
!     mpi_master=.true.
! #endif
!     !
!     Nlat  = size(zeta,1)
!     Nspin = size(Gkout,2)
!     Norb  = size(Gkout,4)
!     Lfreq = size(zeta,4)
!     Nso   = Nspin*Norb
!     Nlso  = Nlat*Nspin*Norb
!     call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_ineq_mpi","zeta")
!     call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_ineq_mpi","Hk")
!     call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_ineq_mpi","Gkout")
!     !
!     if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
!          stop "invert_gk_normal_ineq_mpi: lattice_Gamma_mats & lattice_Gamma_real both allocated"
!     if(allocated(lattice_gamma_mats))then
!        call assert_shape(lattice_gamma_mats,[Nlso,Nlso,Lfreq],"invert_gk_normal_ineq_mpi","lattice_gamma_mats")
!        allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_mats
!     endif
!     if(allocated(lattice_gamma_real))then
!        call assert_shape(lattice_gamma_real,[Nlso,Nlso,Lfreq],"invert_gk_normal_ineq_mpi","lattice_gamma_real")
!        allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_real
!     endif
!     !
!     allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
!     allocate(Gmatrix(Nlso,Nlso))
!     Gktmp=zero
!     do i=1+mpi_rank,Lfreq,mpi_size
!        Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nlat,Nso) - Hk
!        if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
!        if(hk_symm) then
!           call inv_sym(Gmatrix)
!        else
!           call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
!        end if
!        !store the diagonal blocks directly into the tmp output 
!        do ilat=1,Nlat
!           do ispin=1,Nspin
!              do jspin=1,Nspin
!                 do iorb=1,Norb
!                    do jorb=1,Norb
!                       io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
!                       jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
!                       Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! #ifdef _MPI    
!     if(check_MPI())then
!        Gkout=zero
!        call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
!     else
!        Gkout=Gktmp
!     endif
! #else
!     Gkout=Gktmp
! #endif
!   end subroutine invert_gk_normal_ineq_mpi


!   !INVERT GK CLUSTER + REAL SPACE 
!   !SERIAL (OR PARALLEL ON K)
! #if __GFORTRAN__ &&  __GNUC__ > 8
!   subroutine invert_gk_normal_cluster_ineq(zeta,Hk,hk_symm,Gkout)
!     complex(8),dimension(:,:,:,:),intent(in)            :: zeta    ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:),intent(in)                :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     logical,intent(in)                                  :: hk_symm
!     complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     !
!     complex(8),dimension(:,:,:),allocatable             :: Gembed  ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:),allocatable               :: Gmatrix ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb]
!     integer                                             :: Nineq,Nlat,Nspin,Norb,Nso,Nlso,Nilso,Lfreq
!     integer                                             :: i,is,iineq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
!     !
!     Nineq = size(zeta,1)
!     Nlat  = size(Gkout,2)
!     Nspin = size(Gkout,4)
!     Norb  = size(Gkout,6)
!     Lfreq = size(zeta,4)
!     Nso   = Nspin*Norb
!     Nlso  = Nlat*Nspin*Norb
!     Nilso = Nineq*Nlat*Nspin*Norb
!     !Testing:
!     call assert_shape(zeta,[Nineq,Nlso,Nlso,Lfreq],"invert_gk_normal_ineq","zeta")
!     call assert_shape(Hk,[Nilso,Nilso],"invert_gk_normal_ineq","Hk")
!     call assert_shape(Gkout,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_ineq","Gkout")
!     !
!     if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
!          stop "invert_gk_normal_ineq: lattice_Gamma_mats & lattice_Gamma_real both allocated"
!     if(allocated(lattice_gamma_mats))then
!        call assert_shape(lattice_gamma_mats,[Nilso,Nilso,Lfreq],"invert_gk_normal_ineq","lattice_gamma_mats")
!        allocate(Gembed(Nilso,Nilso,Lfreq));Gembed=lattice_gamma_mats
!     endif
!     if(allocated(lattice_gamma_real))then
!        call assert_shape(lattice_gamma_real,[Nilso,Nilso,Lfreq],"invert_gk_normal_ineq","lattice_gamma_real")
!        allocate(Gembed(Nilso,Nilso,Lfreq));Gembed=lattice_gamma_real
!     endif
!     !
!     allocate(Gktmp(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
!     allocate(Gmatrix(Nilso,Nilso))
!     Gktmp=zero
!     do i=1,Lfreq
!        Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nineq,Nlso) - Hk
!        if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
!        if(hk_symm) then
!           call inv_sym(Gmatrix)
!        else
!           call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
!        end if
!        !store the diagonal blocks directly into the tmp output 
!        do iineq=1,Nineq
!           do ilat=1,Nlat
!              do jlat=1,Nlat
!                 do ispin=1,Nspin
!                    do jspin=1,Nspin
!                       do iorb=1,Norb
!                          do jorb=1,Norb
!                             io = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
!                             jo = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
!                             Gktmp(iineq,ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
!                          enddo
!                       enddo
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
!     Gkout = Gktmp
!   end subroutine invert_gk_normal_cluster_ineq


!   subroutine invert_gk_normal_cluster_ineq_mpi(zeta,Hk,hk_symm,Gkout)
!     complex(8),dimension(:,:,:,:),intent(in)            :: zeta    ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:),intent(in)                :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
!     logical,intent(in)                                  :: hk_symm
!     complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     !
!     complex(8),dimension(:,:,:),allocatable             :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
!     complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
!     complex(8),dimension(:,:),allocatable               :: Gmatrix ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb]
!     integer                                             :: Nineq,Nlat,Nspin,Norb,Nso,Nlso,Nilso,Lfreq
!     integer                                             :: i,is,iineq,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
!     !
!     !MPI setup:
! #ifdef _MPI    
!     if(check_MPI())then
!        mpi_size  = get_size_MPI()
!        mpi_rank =  get_rank_MPI()
!        mpi_master= get_master_MPI()
!     else
!        mpi_size=1
!        mpi_rank=0
!        mpi_master=.true.
!     endif
! #else
!     mpi_size=1
!     mpi_rank=0
!     mpi_master=.true.
! #endif
!     !
!     Nineq = size(zeta,1)
!     Nlat  = size(Gkout,2)
!     Nspin = size(Gkout,4)
!     Norb  = size(Gkout,6)
!     Lfreq = size(zeta,4)
!     Nso   = Nspin*Norb
!     Nlso  = Nlat*Nspin*Norb
!     Nilso = Nineq*Nlat*Nspin*Norb
!     !Testing:
!     call assert_shape(zeta,[Nineq,Nlso,Nlso,Lfreq],"invert_gk_normal_ineq","zeta")
!     call assert_shape(Hk,[Nilso,Nilso],"invert_gk_normal_ineq","Hk")
!     call assert_shape(Gkout,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_ineq","Gkout")
!     !
!     if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
!          stop "invert_gk_normal_ineq_mpi: lattice_Gamma_mats & lattice_Gamma_real both allocated"
!     if(allocated(lattice_gamma_mats))then
!        call assert_shape(lattice_gamma_mats,[Nilso,Nilso,Lfreq],"invert_gk_normal_ineq_mpi","lattice_gamma_mats")
!        allocate(Gembed(Nilso,Nilso,Lfreq));Gembed=lattice_gamma_mats
!     endif
!     if(allocated(lattice_gamma_real))then
!        call assert_shape(lattice_gamma_real,[Nilso,Nilso,Lfreq],"invert_gk_normal_ineq_mpi","lattice_gamma_real")
!        allocate(Gembed(Nilso,Nilso,Lfreq));Gembed=lattice_gamma_real
!     endif
!     !
!     allocate(Gktmp(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
!     allocate(Gmatrix(Nilso,Nilso))
!     Gktmp=zero
!     do i=1+mpi_rank,Lfreq,mpi_size
!        Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nineq,Nlso) - Hk
!        if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
!        if(hk_symm) then
!           call inv_sym(Gmatrix)
!        else
!           call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
!        end if
!        !store the diagonal blocks directly into the tmp output 
!        do iineq=1,Nineq
!           do ilat=1,Nlat
!              do jlat=1,Nlat
!                 do ispin=1,Nspin
!                    do jspin=1,Nspin
!                       do iorb=1,Norb
!                          do jorb=1,Norb
!                             io = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
!                             jo = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
!                             Gktmp(iineq,ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
!                          enddo
!                       enddo
!                    enddo
!                 enddo
!              enddo
!           enddo
!        enddo
!     enddo
! #ifdef _MPI    
!     if(check_MPI())then
!        Gkout=zero
!        call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
!     else
!        Gkout=Gktmp
!     endif
! #else
!     Gkout=Gktmp
! #endif
!   end subroutine invert_gk_normal_cluster_ineq_mpi

! #endif
