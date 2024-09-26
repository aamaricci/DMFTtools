module SC_COMMON
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,str,txtfy,to_lower
  USE SF_ARRAYS,    only:linspace,arange
  USE SF_LINALG,    only:eye,inv
  USE SF_MISC,      only:assert_shape
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
  interface get_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface get_block

  !-----------------------------------------!
  ! FROM Matrix=Rank2  TO rank3-7           !
  !-----------------------------------------!
  ![N,N] => [Nlat,Nso,Nso]
  interface to_rank3
     module procedure :: d_matrix_TO_rank3
     module procedure :: c_matrix_TO_rank3
     module procedure :: d_matrix_TO_rank3L
     module procedure :: c_matrix_TO_rank3L
  end interface to_rank3
  !
  ![N,N] =>  [Nspin,Nspin,Norb,Norb]
  interface to_rank4
     module procedure :: d_matrix_TO_rank4
     module procedure :: c_matrix_TO_rank4
     module procedure :: d_matrix_TO_rank4L
     module procedure :: c_matrix_TO_rank4L
  end interface to_rank4
  !
  ![N,N] => [Nlat,Nspin,Nspin,Norb,Norb]
  interface to_rank5
     module procedure :: d_matrix_TO_rank5
     module procedure :: c_matrix_TO_rank5
     module procedure :: d_matrix_TO_rank5L
     module procedure :: c_matrix_TO_rank5L
  end interface to_rank5
  !
  ![N,N] => [Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  interface to_rank6
     module procedure :: d_matrix_TO_rank6
     module procedure :: c_matrix_TO_rank6
     module procedure :: d_matrix_TO_rank6L
     module procedure :: c_matrix_TO_rank6L
  end interface to_rank6
  !
  ![N,N] => [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  interface to_rank7
     module procedure :: d_matrix_TO_rank7
     module procedure :: c_matrix_TO_rank7
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: d_matrix_TO_rank7L
     module procedure :: c_matrix_TO_rank7L
#endif
  end interface to_rank7
  !-----------------------------------------!


  !-----------------------------------------!
  ! FROM rank3-7 TO Matrix=Rank2            !
  !-----------------------------------------!
  ![Nlat,Nso,Nso] => [N,N]
  interface from_rank3
     module procedure :: d_rank3_TO_matrix
     module procedure :: c_rank3_TO_matrix
     module procedure :: d_rank3_TO_matrixL
     module procedure :: c_rank3_TO_matrixL
  end interface from_rank3
  !
  ![Nspin,Nspin,Norb,Norb] => [N,N]
  interface from_rank4
     module procedure :: d_rank4_TO_matrix
     module procedure :: c_rank4_TO_matrix
     module procedure :: d_rank4_TO_matrixL
     module procedure :: c_rank4_TO_matrixL
  end interface from_rank4
  !
  ![Nlat,Nspin,Nspin,Norb,Norb] => [N,N]
  interface from_rank5
     module procedure :: d_rank5_TO_matrix
     module procedure :: c_rank5_TO_matrix
     module procedure :: d_rank5_TO_matrixL
     module procedure :: c_rank5_TO_matrixL
  end interface from_rank5
  !
  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb] => [N,N]
  interface from_rank6
     module procedure :: d_rank6_TO_matrix
     module procedure :: c_rank6_TO_matrix
     module procedure :: d_rank6_TO_matrixL
     module procedure :: c_rank6_TO_matrixL
  end interface from_rank6
  !
  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb] => [N,N]
  interface from_rank7
     module procedure :: d_rank7_TO_matrix
     module procedure :: c_rank7_TO_matrix
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: d_rank7_TO_matrixL
     module procedure :: c_rank7_TO_matrixL
#endif
  end interface from_rank7
  !-----------------------------------------!



  interface gf_reduce
     module procedure  :: gf_reduce_rank2
     module procedure  :: gf_reduce_rank3
     module procedure  :: gf_reduce_rank4
     module procedure  :: gf_reduce_rank5
     module procedure  :: gf_reduce_rank6
     module procedure  :: gf_reduce_rank7
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure  :: gf_reduce_rank8
#endif
  end interface gf_reduce



  integer                             :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats,Nineq,Nilso,Ntot
  integer                             :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js,iineq
  !
  integer                             :: mpi_ierr
  integer                             :: mpi_rank
  integer                             :: mpi_size
  logical                             :: mpi_master
  !
  real(8)                             :: beta
  real(8)                             :: xmu,eps
  real(8)                             :: wini,wfin 
  !
  real(8),dimension(:),allocatable    :: wm !Matsubara frequencies
  real(8),dimension(:),allocatable    :: wr !Real frequencies
  integer                             :: Lfreq
  complex(8),dimension(:),allocatable :: wfreq
  logical                             :: pushed=.false.

  character(len=128)                  :: suffix
  character(len=128)                  :: weiss_suffix=".dat"


contains







  !####################################################################
  !####################################################################
  !
  !            GF REDUCTION MPI Rank2-->8
  !
  !####################################################################
  !####################################################################
  subroutine gf_reduce_rank2(G)
    complex(8),dimension(:,:)                       :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank2

  subroutine gf_reduce_rank3(G)
    complex(8),dimension(:,:,:) :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank3

  subroutine gf_reduce_rank4(G)
    complex(8),dimension(:,:,:,:) :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank4

  subroutine gf_reduce_rank5(G)
    complex(8),dimension(:,:,:,:,:) :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank5

  subroutine gf_reduce_rank6(G)
    complex(8),dimension(:,:,:,:,:,:) :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank6

  subroutine gf_reduce_rank7(G)
    complex(8),dimension(:,:,:,:,:,:,:) :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank7

#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine gf_reduce_rank8(G)
    complex(8),dimension(:,:,:,:,:,:,:,:) :: G
#ifdef _MPI    
    if(check_MPI())then
       call MPI_AllReduce(MPI_IN_PLACE, G, size(G), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    endif
#endif
  end subroutine gf_reduce_rank8
#endif




  !####################################################################
  !####################################################################
  !
  !        GENERATE FREQUENCY ARRAYS
  !
  !####################################################################
  !####################################################################
  subroutine gf_push_zeta(zeta,axis)
    complex(8),dimension(:) :: zeta
    character(len=*)        :: axis
    Lfreq=size(zeta)
    if(allocated(wfreq))deallocate(wfreq)
    allocate(wfreq(Lfreq))
    select case(to_lower(axis(1:1)))
    case default;
       stop "gf_push_zeta error: axis undefined. axis=[matsubara,realaxis]"
    case("m")
       wfreq = zeta
    case("r")
       wfreq = zeta
    end select
    pushed=.true.
  end subroutine gf_push_zeta


  subroutine build_frequency_array(axis)
    character(len=*) :: axis
    call get_ctrl_var(xmu,"XMU")
    if(pushed)then
       if(size(wfreq)/=Lfreq)stop "build_frequency_array ERROR: pushed wfreq has wrong size"
    else
       if(allocated(wfreq))deallocate(wfreq)
       allocate(wfreq(Lfreq))
       select case(to_lower(str(axis(1:1))))
       case default;
          stop "build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
       case("m")
          call get_ctrl_var(beta,"BETA")
          wfreq = dcmplx(0d0,pi/beta*(2*arange(1,Lfreq)-1))
       case("r")
          call get_ctrl_var(wini,"WINI")
          call get_ctrl_var(wfin,"WFIN")
          call get_ctrl_var(eps,"EPS")
          wfreq = dcmplx(linspace(wini,wfin,Lfreq),eps)
       end select
    endif
    return
  end subroutine build_frequency_array





  !####################################################################
  !####################################################################
  !
  !       BUILD DELTA FUNCTIONS NORMAL/SUPERC per Frequency
  !
  !####################################################################
  !####################################################################
  !                                           !
  ! [Delta] = [w-Hloc] - [Sigma] -  [Gloc]^-1 !
  !                                           !
  subroutine dmft_delta_normal(z,Gi,Si,Di)
    complex(8),dimension(:,:)                     :: z
    complex(8),dimension(size(z,1),size(z,2))     :: Gi
    complex(8),dimension(size(z,1),size(z,2))     :: Si
    complex(8),dimension(size(z,1),size(z,2))     :: Di
    complex(8),dimension(size(z,1),size(z,2))     :: invG
    !
    invG  = Gi
    !
    call inv(invG)
    !
    Di = z - Si - invG
    !
  end subroutine dmft_delta_normal

  subroutine dmft_delta_superc(z,Gi,Fi,nSi,aSi,Wi,Ti)
    complex(8),dimension(:,:)                     :: z
    complex(8),dimension(size(z,1),size(z,1))     :: Gi
    complex(8),dimension(size(z,1),size(z,1))     :: Fi
    complex(8),dimension(size(z,1),size(z,1))     :: nSi
    complex(8),dimension(size(z,1),size(z,1))     :: aSi
    complex(8),dimension(size(z,1),size(z,1))     :: Wi
    complex(8),dimension(size(z,1),size(z,1))     :: Ti
    integer                                       :: i,j,N
    complex(8),dimension(2*size(z,1),2*size(z,1)) :: invG
    !
    N=size(z,1)
    !
    do concurrent(i=1:N,j=1:N)
       invG(i,j)      =        Gi(i,j)
       invG(i,j+N)    =        Fi(i,j)
       invG(i+N,j)    =        Fi(i,j)
       invG(i+N,j+N)  =-conjg( Gi(i,j) )
    enddo
    !
    call inv(invG)
    !
    Wi = z - nSi - invG(1:N,1:N)     ![w+mu]1 - Hloc - Sigma - [G^-1]
    Ti =   - aSi - invG(1:N,N+1:2*N) !               - Self  - [F^-1]
  end subroutine dmft_delta_superc


  !                                     !
  ![G0]^-1 = [ [Gloc]^-1 + [Sigma] ]^-1 ! 
  !                                     !
  subroutine dmft_weiss_normal(Gi,Si,Wi)
    complex(8),dimension(:,:)                       :: Gi
    complex(8),dimension(size(Gi,1),size(Gi,2))     :: Si
    complex(8),dimension(size(Gi,1),size(Gi,2))     :: Wi
    complex(8),dimension(size(Gi,1),size(Gi,2))     :: invG
    !
    invG  = Gi
    call inv(invG)
    !
    Wi = invG + Si
    call inv(Wi)
    !
  end subroutine dmft_weiss_normal

  subroutine dmft_weiss_superc(Gi,Fi,nSi,aSi,Wi,Ti)
    complex(8),dimension(:,:)                       :: Gi
    complex(8),dimension(size(Gi,1),size(Gi,1))     :: Fi
    complex(8),dimension(size(Gi,1),size(Gi,1))     :: nSi
    complex(8),dimension(size(Gi,1),size(Gi,1))     :: aSi
    complex(8),dimension(size(Gi,1),size(Gi,1))     :: Wi
    complex(8),dimension(size(Gi,1),size(Gi,1))     :: Ti
    integer                                         :: i,j,N
    complex(8),dimension(2*size(Gi,1),2*size(Gi,1)) :: S,G
    !
    N=size(Gi,1)
    !
    do concurrent(i=1:N,j=1:N)
       S(i,j)      =        nSi(i,j)
       S(i,j+N)    =        aSi(i,j)
       S(i+N,j)    =        aSi(i,j)
       S(i+N,j+N)  =-conjg( nSi(i,j) )
       !
       G(i,j)      =        Gi(i,j)
       G(i,j+N)    =        Fi(i,j)
       G(i+N,j)    =        Fi(i,j)
       G(i+N,j+N)  =-conjg( Gi(i,j) )
    enddo
    !
    call inv(G)
    !
    !G0^-1 = G^-1 + S
    G = G + S
    !
    call inv(G)
    !
    do concurrent(i=1:N,j=1:N)
       Wi(i,j) = G(i,j)
       Ti(i,j) = G(i,j+N)
    enddo
  end subroutine dmft_weiss_superc


  !####################################################################
  ! SET MPI
  !####################################################################
  subroutine set_gf_mpi()
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
  end subroutine set_gf_mpi




  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################

  !--------------------------------------------------------------------!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !--------------------------------------------------------------------!
  function select_block_Nlso(ilat,Matrix) result(Vblock)
    integer                                 :: ilat
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    complex(8),dimension(Nso,Nso)           :: Vblock
    integer                                 :: i,j
    Vblock=zero
    i = 1+(ilat-1)*Nso
    j =       ilat*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ilat,Matrix) result(Vblock)
    integer                                          :: ilat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Vblock(is,js) = Matrix(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function select_block_nnn


  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
  function blocks_to_matrix(Vblocks) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ilat
    Matrix=zero
    do ilat=1,Nlat
       i = 1 + (ilat-1)*Nso
       j = ilat*Nso
       Matrix(i:j,i:j) =  Vblocks(ilat,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ilat
    Vblocks=zero
    do ilat=1,Nlat
       i = 1 + (ilat-1)*Nso
       j = ilat*Nso
       Vblocks(ilat,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks



  !####################################################################
  !####################################################################


  !DBLE
  !MATRIX --> RANK3_7
  function d_matrix_TO_rank3(Hmat) result(Hrank)
    real(8),dimension(Nlat*Nso,Nlat*Nso)                                    :: Hmat
    real(8),dimension(Nlat,Nso,Nso)                                         :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hrank(ilat,is,js) = Hmat(i,j)
    enddo
  end function d_matrix_TO_rank3

  function d_matrix_TO_rank4(Hmat) result(Hrank)
    real(8),dimension(Nspin*Norb,Nspin*Norb)                                :: Hmat
    real(8),dimension(Nspin,Nspin,Norb,Norb)                                :: Hrank
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank4

  function d_matrix_TO_rank5(Hmat) result(Hrank)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                           :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hrank(ilat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank5

  function d_matrix_TO_rank6(Hmat) result(Hrank)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                      :: Hrank
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank6

  function d_matrix_TO_rank7(Hmat) result(Hrank)
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)          :: Hmat
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)                :: Hrank
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank7

  function d_matrix_TO_rank3L(Hmat) result(Hrank)
    real(8),dimension(Nlat*Nso,Nlat*Nso,Lfreq)                              :: Hmat
    real(8),dimension(Nlat,Nso,Nso,Lfreq)                                   :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hrank(ilat,is,js,:) = Hmat(i,j,:)
    enddo
  end function d_matrix_TO_rank3L

  function d_matrix_TO_rank4L(Hmat) result(Hrank)
    real(8),dimension(Nspin*Norb,Nspin*Norb,Lfreq)                          :: Hmat
    real(8),dimension(Nspin,Nspin,Norb,Norb,Lfreq)                          :: Hrank
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank4L

  function d_matrix_TO_rank5L(Hmat) result(Hrank)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                     :: Hrank
    integer                                                                 :: iorb,ispin,ilat,is
    integer                                                                 :: jorb,jspin,jlat,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hrank(ilat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank5L

  function d_matrix_TO_rank6L(Hmat) result(Hrank)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                :: Hrank
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank6L

#if __GFORTRAN__ &&  __GNUC__ > 8
  function d_matrix_TO_rank7L(Hmat) result(Hrank)
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,Lfreq)    :: Hmat
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)          :: Hrank
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank7L
#endif

  !##################################################################


  !RANK4_7 --> MATRIX:
  function d_rank3_TO_matrix(Hrank) result(Hmat)
    real(8),dimension(Nlat*Nso,Nlat*Nso)                                    :: Hmat
    real(8),dimension(Nlat,Nso,Nso)                                         :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j) = Hrank(ilat,is,js)
    enddo
  end function d_rank3_TO_matrix

  function d_rank4_TO_matrix(Hrank) result(Hmat)
    real(8),dimension(Nspin,Nspin,Norb,Norb)                                :: Hrank
    real(8),dimension(Nspin*Norb,Nspin*Norb)                                :: Hmat
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js) = Hrank(ispin,jspin,iorb,jorb)
    enddo
  end function d_rank4_TO_matrix


  function d_rank5_TO_matrix(Hrank) result(Hmat)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                           :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    integer                                                                 :: iorb,ispin,ilat,is
    integer                                                                 :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb, jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hmat(is,js) = Hrank(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function d_rank5_TO_matrix

  function d_rank6_TO_matrix(Hrank) result(Hmat)
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                      :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function d_rank6_TO_matrix

  function d_rank7_TO_matrix(Hrank) result(Hmat)
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)                :: Hrank
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)          :: Hmat
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function d_rank7_TO_matrix


  function d_rank3_TO_matrixL(Hrank) result(Hmat)
    real(8),dimension(Nlat*Nso,Nlat*Nso,Lfreq)                              :: Hmat
    real(8),dimension(Nlat,Nso,Nso,Lfreq)                                   :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j,:) = Hrank(ilat,is,js,:)
    enddo
  end function d_rank3_TO_matrixL


  function d_rank4_TO_matrixL(Hrank) result(Hmat)
    real(8),dimension(Nspin,Nspin,Norb,Norb,Lfreq)                          :: Hrank
    real(8),dimension(Nspin*Norb,Nspin*Norb,Lfreq)                          :: Hmat
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js,:) = Hrank(ispin,jspin,iorb,jorb,:)
    enddo
  end function d_rank4_TO_matrixL

  function d_rank5_TO_matrixL(Hrank) result(Hmat)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                     :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    integer                                                                 :: iorb,ispin,ilat,is
    integer                                                                 :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hmat(is,js,:) = Hrank(ilat,ispin,jspin,iorb,jorb,:)
    enddo
  end function d_rank5_TO_matrixL

  function d_rank6_TO_matrixL(Hrank) result(Hmat)
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js,:) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:)
    enddo
  end function d_rank6_TO_matrixL

#if __GFORTRAN__ &&  __GNUC__ > 8
  function d_rank7_TO_matrixL(Hrank) result(Hmat)
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)          :: Hrank
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,Lfreq)    :: Hmat
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js,:) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) 
    enddo
  end function d_rank7_TO_matrixL
#endif



  !##################################################################
  !##################################################################
  !##################################################################







  !COMPLEX
  !MATRIX --> RANK3_7
  function c_matrix_TO_rank3(Hmat) result(Hrank)
    complex(8),dimension(Nlat*Nso,Nlat*Nso)                                    :: Hmat
    complex(8),dimension(Nlat,Nso,Nso)                                         :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hrank(ilat,is,js) = Hmat(i,j)
    enddo
  end function c_matrix_TO_rank3

  function c_matrix_TO_rank4(Hmat) result(Hrank)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                                :: Hmat
    complex(8),dimension(Nspin,Nspin,Norb,Norb)                                :: Hrank
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank4

  function c_matrix_TO_rank5(Hmat) result(Hrank)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                           :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hrank(ilat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank5

  function c_matrix_TO_rank6(Hmat) result(Hrank)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                      :: Hrank
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank6

  function c_matrix_TO_rank7(Hmat) result(Hrank)
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb) :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)       :: Hrank
    integer                                                           :: iineq
    integer                                                           :: ilat,jlat
    integer                                                           :: iorb,jorb
    integer                                                           :: ispin,jspin
    integer                                                           :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank7

  function c_matrix_TO_rank3L(Hmat) result(Hrank)
    complex(8),dimension(Nlat*Nso,Nlat*Nso,Lfreq)                              :: Hmat
    complex(8),dimension(Nlat,Nso,Nso,Lfreq)                                   :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hrank(ilat,is,js,:) = Hmat(i,j,:)
    enddo
  end function c_matrix_TO_rank3L

  function c_matrix_TO_rank4L(Hmat) result(Hrank)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lfreq)                          :: Hmat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lfreq)                          :: Hrank
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank4L

  function c_matrix_TO_rank5L(Hmat) result(Hrank)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                     :: Hrank
    integer                                                                 :: iorb,ispin,ilat,is
    integer                                                                 :: jorb,jspin,jlat,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hrank(ilat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank5L

  function c_matrix_TO_rank6L(Hmat) result(Hrank)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                :: Hrank
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank6L

#if __GFORTRAN__ &&  __GNUC__ > 8
  function c_matrix_TO_rank7L(Hmat) result(Hrank)
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,Lfreq)    :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)          :: Hrank
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank7L
#endif


  !##################################################################


  !RANK4_7 --> MATRIX:
  function c_rank3_TO_matrix(Hrank) result(Hmat)
    complex(8),dimension(Nlat*Nso,Nlat*Nso)                                    :: Hmat
    complex(8),dimension(Nlat,Nso,Nso)                                         :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j) = Hrank(ilat,is,js)
    enddo
  end function c_rank3_TO_matrix

  function c_rank4_TO_matrix(Hrank) result(Hmat)
    complex(8),dimension(Nspin,Nspin,Norb,Norb)                                :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                                :: Hmat
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js) = Hrank(ispin,jspin,iorb,jorb)
    enddo
  end function c_rank4_TO_matrix


  function c_rank5_TO_matrix(Hrank) result(Hmat)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                           :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    integer                                                                 :: iorb,ispin,ilat,is
    integer                                                                 :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb, jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hmat(is,js) = Hrank(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function c_rank5_TO_matrix

  function c_rank6_TO_matrix(Hrank) result(Hmat)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                      :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)                      :: Hmat
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function c_rank6_TO_matrix

  function c_rank7_TO_matrix(Hrank) result(Hmat)
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)                :: Hrank
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)          :: Hmat
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function c_rank7_TO_matrix


  function c_rank3_TO_matrixL(Hrank) result(Hmat)
    complex(8),dimension(Nlat*Nso,Nlat*Nso,Lfreq)                              :: Hmat
    complex(8),dimension(Nlat,Nso,Nso,Lfreq)                                   :: Hrank
    integer                                                                 :: ilat
    integer                                                                 :: i,is
    integer                                                                 :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j,:) = Hrank(ilat,is,js,:)
    enddo
  end function c_rank3_TO_matrixL


  function c_rank4_TO_matrixL(Hrank) result(Hmat)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lfreq)                          :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lfreq)                          :: Hmat
    integer                                                                 :: iorb,ispin,is
    integer                                                                 :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js,:) = Hrank(ispin,jspin,iorb,jorb,:)
    enddo
  end function c_rank4_TO_matrixL

  function c_rank5_TO_matrixL(Hrank) result(Hmat)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                     :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    integer                                                                 :: iorb,ispin,ilat,is
    integer                                                                 :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hmat(is,js,:) = Hrank(ilat,ispin,jspin,iorb,jorb,:)
    enddo
  end function c_rank5_TO_matrixL

  function c_rank6_TO_matrixL(Hrank) result(Hmat)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)                :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lfreq)                :: Hmat
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js,:) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:)
    enddo
  end function c_rank6_TO_matrixL

#if __GFORTRAN__ &&  __GNUC__ > 8
  function c_rank7_TO_matrixL(Hrank) result(Hmat)
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq)          :: Hrank
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,Lfreq)    :: Hmat
    integer                                                                 :: iineq
    integer                                                                 :: ilat,jlat
    integer                                                                 :: iorb,jorb
    integer                                                                 :: ispin,jspin
    integer                                                                 :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js,:) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) 
    enddo
  end function c_rank7_TO_matrixL
#endif




















end module SC_COMMON
