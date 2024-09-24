module GF_COMMON
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,str,txtfy,splot,sread,file_gzip,file_gunzip,file_targz,file_untargz,to_lower
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

  ![Nlat,Nso,Nso] <-> [N,N]
  interface reshape_matrix_to_rank3
     module procedure :: d_matrix_TO_rank3
     module procedure :: c_matrix_TO_rank3
     module procedure :: d_matrix_TO_rank3L
     module procedure :: c_matrix_TO_rank3L
  end interface reshape_matrix_to_rank3
  !
  interface reshape_rank3_to_matrix
     module procedure :: d_rank3_TO_matrix
     module procedure :: c_rank3_TO_matrix
     module procedure :: d_rank3_TO_matrixL
     module procedure :: c_rank3_TO_matrixL
  end interface reshape_rank3_to_matrix


  ![Nspin,Nspin,Norb,Norb] <-> [N,N]
  interface reshape_matrix_to_rank4
     module procedure :: d_matrix_TO_rank4
     module procedure :: c_matrix_TO_rank4
     module procedure :: d_matrix_TO_rank4L
     module procedure :: c_matrix_TO_rank4L
  end interface reshape_matrix_to_rank4
  !
  interface reshape_rank4_to_matrix
     module procedure :: d_rank4_TO_matrix
     module procedure :: c_rank4_TO_matrix
     module procedure :: d_rank4_TO_matrixL
     module procedure :: c_rank4_TO_matrixL
  end interface reshape_rank4_to_matrix


  ![Nlat,Nspin,Nspin,Norb,Norb] <-> [N,N]
  interface reshape_matrix_to_rank5
     module procedure :: d_matrix_TO_rank5
     module procedure :: c_matrix_TO_rank5
     module procedure :: d_matrix_TO_rank5L
     module procedure :: c_matrix_TO_rank5L
  end interface reshape_matrix_to_rank5
  !
  interface reshape_rank5_to_matrix
     module procedure :: d_rank5_TO_matrix
     module procedure :: c_rank5_TO_matrix
     module procedure :: d_rank5_TO_matrixL
     module procedure :: c_rank5_TO_matrixL
  end interface reshape_rank5_to_matrix

  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb] <-> [N,N]
  interface reshape_matrix_to_rank6
     module procedure :: d_matrix_TO_rank6
     module procedure :: c_matrix_TO_rank6
     module procedure :: d_matrix_TO_rank6L
     module procedure :: c_matrix_TO_rank6L
  end interface reshape_matrix_to_rank6
  !
  interface reshape_rank6_to_matrix
     module procedure :: d_rank6_TO_matrix
     module procedure :: c_rank6_TO_matrix
     module procedure :: d_rank6_TO_matrixL
     module procedure :: c_rank6_TO_matrixL
  end interface reshape_rank6_to_matrix


  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb] <-> [N,N]
  interface reshape_matrix_to_rank7
     module procedure :: d_matrix_TO_rank7
     module procedure :: c_matrix_TO_rank7
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: d_matrix_TO_rank7L
     module procedure :: c_matrix_TO_rank7L
#endif
  end interface reshape_matrix_to_rank7

  interface reshape_rank7_to_matrix
     module procedure :: d_rank7_TO_matrix
     module procedure :: c_rank7_TO_matrix
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: d_rank7_TO_matrixL
     module procedure :: c_rank7_TO_matrixL
#endif
  end interface reshape_rank7_to_matrix



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


  integer                             :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Nineq,Nilso,Ntot
  integer                             :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js,iineq,iel,jel
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
  real(8),dimension(:),allocatable    :: wio
  logical                             :: pushed=.false.

  character(len=128)                  :: suffix
  character(len=128)                  :: gf_suffix='.dat'
  character(len=8)                    :: w_suffix

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
    if(allocated(wio))deallocate(wio)
    allocate(wio(Lfreq))
    select case(to_lower(axis(1:1)))
    case default;
       stop "gf_push_zeta error: axis undefined. axis=[matsubara,realaxis]"
    case("m")
       wfreq = zeta
       wio   = dimag(zeta)
    case("r")
       wfreq = zeta
       wio   = dreal(zeta)
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
          wio   = pi/beta*(2*arange(1,Lfreq)-1)
          wfreq = dcmplx(0d0,wio)
       case("r")
          call get_ctrl_var(wini,"WINI")
          call get_ctrl_var(wfin,"WFIN")
          call get_ctrl_var(eps,"EPS")
          wio   = linspace(wini,wfin,Lfreq)
          wfreq = dcmplx(wio,eps)
       end select
    endif
    select case(to_lower(str(axis(1:1))))
    case default;stop "dmft_gfio_build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
    case("m");w_suffix="_iw"
    case("r") ;w_suffix="_realw"
    end select
    return
  end subroutine build_frequency_array

  
  subroutine set_gf_suffix(string)
    character(len=*) :: string
    gf_suffix=reg(string)
  end subroutine set_gf_suffix




  !####################################################################
  !####################################################################
  !
  !        GET LOCAL SUPERCONDUCTING Z(2,2,:,:)
  !
  !####################################################################
  !####################################################################
  function sc_zi_rank2(i,S,axis) result(zi)
    integer                                            :: i,L
    integer                                            :: N
    complex(8),dimension(2,Ntot,Ntot,Lfreq),intent(in) :: S
    character(len=*),intent(in)                        :: axis
    complex(8),dimension(2,2,Ntot,Ntot)                :: zi
    !
    if(i<0 .OR. i>Lfreq)stop  "sc_zi_rank2 error: i<0 OR i>L "
    !
    select case(axis(1:1))
    case default; stop "sc_zi_rank2 error: axis != {m,r}"
    case("m","M")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) -        S(1,:,:,i)
       zi(1,2,:,:) =                          -        S(2,:,:,i)
       zi(2,1,:,:) =                          -        S(2,:,:,i)
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + conjg( S(1,:,:,i))
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - S(1,:,:,i)
       zi(1,2,:,:) =                                         - S(2,:,:,i)
       zi(2,1,:,:) =                                         - S(2,:,:,i)
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + conjg( S(1,:,:,Lfreq+1-i) )
    end select
  end function sc_zi_rank2

  function sc_zi_rank3(i,S,axis) result(zi)
    integer                                               :: i,L
    integer                                               :: N
    complex(8),dimension(2,Nlat,Nso,Nso,Lfreq),intent(in) :: S
    character(len=*),intent(in)                           :: axis
    complex(8),dimension(2,2,Ntot,Ntot)                   :: zi
    !
    if(i<0 .OR. i>Lfreq)stop  "sc_zi_rank3 error: i<0 OR i>L "
    !
    select case(axis(1:1))
    case default; stop "sc_zi_rank3 error: axis != {m,r}"
    case("m","M")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - reshape_rank3_to_matrix(S(1,:,:,:,i),Nlat,Nso)
       zi(1,2,:,:) =                          - reshape_rank3_to_matrix(S(2,:,:,:,i),Nlat,Nso)
       zi(2,1,:,:) =                          - reshape_rank3_to_matrix(S(2,:,:,:,i),Nlat,Nso)
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + reshape_rank3_to_matrix(conjg(S(1,:,:,:,i)),Nlat,Nso)
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - reshape_rank3_to_matrix(S(1,:,:,:,i),Nlat,Nso)
       zi(1,2,:,:) =                                         - reshape_rank3_to_matrix(S(2,:,:,:,i),Nlat,Nso)
       zi(2,1,:,:) =                                         - reshape_rank3_to_matrix(S(2,:,:,:,i),Nlat,Nso)
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + reshape_rank3_to_matrix(conjg(S(1,:,:,:,Lfreq+1-i)),Nlat,Nso)
    end select
  end function sc_zi_rank3

  function sc_zi_rank4(i,S,axis) result(zi)
    integer                                                        :: i,L
    integer                                                        :: N
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Lfreq),intent(in) :: S
    character(len=*),intent(in)                                    :: axis
    complex(8),dimension(2,2,Ntot,Ntot)                            :: zi
    !
    if(i<0 .OR. i>Lfreq)stop  "sc_zi_rank4 error: i<0 OR i>L "
    !
    select case(axis(1:1))
    case default; stop "sc_zi_rank4 error: axis != {m,r}"
    case("m","M")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - reshape_rank4_to_matrix(S(1,:,:,:,:,i),Nspin,Norb)
       zi(1,2,:,:) =                          - reshape_rank4_to_matrix(S(2,:,:,:,:,i),Nspin,Norb)
       zi(2,1,:,:) =                          - reshape_rank4_to_matrix(S(2,:,:,:,:,i),Nspin,Norb)
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + reshape_rank4_to_matrix(conjg(S(1,:,:,:,:,i)),Nspin,Norb)
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - reshape_rank4_to_matrix(S(1,:,:,:,:,i),Nspin,Norb)
       zi(1,2,:,:) =                                         - reshape_rank4_to_matrix(S(2,:,:,:,:,i),Nspin,Norb)
       zi(2,1,:,:) =                                         - reshape_rank4_to_matrix(S(2,:,:,:,:,i),Nspin,Norb)
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + reshape_rank4_to_matrix(conjg(S(1,:,:,:,:,Lfreq+1-i)),Nspin,Norb)
    end select
  end function sc_zi_rank4


  function sc_zi_rank5(i,S,axis) result(zi)
    integer                                                             :: i,L
    integer                                                             :: N
    complex(8),dimension(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq),intent(in) :: S
    character(len=*),intent(in)                                         :: axis
    complex(8),dimension(2,2,Ntot,Ntot)                                 :: zi
    !
    if(i<0 .OR. i>Lfreq)stop  "sc_zi_rank5 error: i<0 OR i>L "
    !
    select case(axis(1:1))
    case default; stop "sc_zi_rank5 error: axis != {m,r}"
    case("m","M")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - reshape_rank5_to_matrix(S(1,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(1,2,:,:) =                          - reshape_rank5_to_matrix(S(2,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,1,:,:) =                          - reshape_rank5_to_matrix(S(2,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + reshape_rank5_to_matrix(conjg(S(1,:,:,:,:,:,i)),Nlat,Nspin,Norb)
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - reshape_rank5_to_matrix(S(1,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(1,2,:,:) =                                         - reshape_rank5_to_matrix(S(2,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,1,:,:) =                                         - reshape_rank5_to_matrix(S(2,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + reshape_rank5_to_matrix(conjg(S(1,:,:,:,:,:,Lfreq+1-i)),Nlat,Nspin,Norb)
    end select
  end function sc_zi_rank5

#if __GFORTRAN__ &&  __GNUC__ > 8
  function sc_zi_rank6(i,S,axis) result(zi)
    integer                                                                  :: i,L
    integer                                                                  :: N
    complex(8),dimension(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq),intent(in) :: S
    character(len=*),intent(in)                                              :: axis
    complex(8),dimension(2,2,Ntot,Ntot)                                      :: zi
    !
    if(i<0 .OR. i>Lfreq)stop  "sc_zi_rank5 error: i<0 OR i>L "
    !
    select case(axis(1:1))
    case default; stop "sc_zi_rank5 error: axis != {m,r}"
    case("m","M")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - reshape_rank6_to_matrix(S(1,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(1,2,:,:) =                          - reshape_rank6_to_matrix(S(2,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,1,:,:) =                          - reshape_rank6_to_matrix(S(2,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + reshape_rank6_to_matrix(conjg(S(1,:,:,:,:,:,:,i)),Nlat,Nspin,Norb)
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - reshape_rank6_to_matrix(S(1,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(1,2,:,:) =                                         - reshape_rank6_to_matrix(S(2,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,1,:,:) =                                         - reshape_rank6_to_matrix(S(2,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + reshape_rank6_to_matrix(conjg(S(1,:,:,:,:,:,:,Lfreq+1-i)),Nlat,Nspin,Norb)
    end select
  end function sc_zi_rank6
#endif





  !####################################################################
  !####################################################################
  !
  !    INVERT G_k(w) for rank2-->8
  !
  !####################################################################
  !####################################################################
  !
  !####################################################################
  !    NORMAL - HK
  !####################################################################
  ![N,N]
  function invert_gki_normal_rank2(z,Hk,Si) result(Gki)
    complex(8)                                             :: z
    complex(8),dimension(:,:),intent(in)                   :: Hk
    complex(8),dimension(size(Hk,1),size(Hk,1)),intent(in) :: Si
    complex(8),dimension(size(Hk,1),size(Hk,1))            :: Gki
    integer                                                :: N
    !
    N   = size(Hk,1)
    Gki = z*eye(N) - Si - Hk
    call inv(Gki)
  end function invert_gki_normal_rank2

  ![Nlat,Nso,Nso]
  function invert_gki_normal_rank3(z,Hk,Si) result(Gki)
    complex(8)                                    :: z
    complex(8),dimension(:,:),intent(in)          :: Hk
    complex(8),dimension(Nlat,Nso,Nso),intent(in) :: Si
    complex(8),dimension(Nlat,Nso,Nso)            :: Gki
    complex(8),dimension(size(Hk,1),size(Hk,1))   :: G
    integer                                       :: N
    !
    N   = size(Hk,1)
    G   = z*eye(N) - reshape_rank3_to_matrix(Si,Nlat,Nso) - Hk
    call inv(G)
    Gki = reshape_matrix_to_rank3(G,Nlat,Nso)
  end function invert_gki_normal_rank3

  ![Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_rank4(z,Hk,Si) result(Gki)
    complex(8)                                                :: z
    complex(8),dimension(:,:),intent(in)                      :: Hk
    complex(8),dimension(Nspin,Nspin,Norb,Norb),intent(in)    :: Si    
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: G
    integer                                                   :: N
    !
    N    = size(Hk,1)
    G   = z*eye(N) - reshape_rank4_to_matrix(Si,Nspin,Norb) - Hk
    call inv(G)
    Gki = reshape_matrix_to_rank4(G,Nspin,Norb)
  end function invert_gki_normal_rank4

  ![Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_rank5(z,Hk,Si) result(Gki)
    complex(8)                                                     :: z
    complex(8),dimension(:,:),intent(in)                           :: Hk
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb),intent(in)    :: Si    
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(size(Hk,1),size(Hk,1))      :: G
    integer                                                        :: N
    !
    N    = size(Hk,1)
    G   = z*eye(N) - reshape_rank5_to_matrix(Si,Nlat,Nspin,Norb) - Hk
    call inv(G)
    Gki = reshape_matrix_to_rank5(G,Nlat,Nspin,Norb)
  end function invert_gki_normal_rank5

  ![Nlat,Nspin,Nspin,Norb,Norb] -> ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_rank5_6(z,Hk,Si) result(Gki)
    complex(8)                                                          :: z
    complex(8),dimension(:,:),intent(in)                                :: Hk
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb),intent(in)         :: Si    
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(size(Hk,1),size(Hk,1))           :: G
    integer                                                             :: N
    !
    N    = size(Hk,1)
    G   = z*eye(N) - reshape_rank5_to_matrix(Si,Nlat,Nspin,Norb) - Hk
    call inv(G)
    Gki = reshape_matrix_to_rank6(G,Nlat,Nspin,Norb)
  end function invert_gki_normal_rank5_6

  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_rank6(z,Hk,Si) result(Gki)
    complex(8)                                                          :: z
    complex(8),dimension(:,:),intent(in)                                :: Hk
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb),intent(in)    :: Si
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(size(Hk,1),size(Hk,1))           :: G
    integer                                                             :: N
    !
    N    = size(Hk,1)
    G   = z*eye(N) - reshape_rank6_to_matrix(Si,Nlat,Nspin,Norb) - Hk
    call inv(G)
    Gki = reshape_matrix_to_rank6(G,Nlat,Nspin,Norb)
  end function invert_gki_normal_rank6

  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_rank7(z,Hk,Si) result(Gki)
    complex(8)                                                                :: z
    complex(8),dimension(:,:),intent(in)                                      :: Hk
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb),intent(in)    :: Si
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(size(Hk,1),size(Hk,1))                 :: G
    integer                                                                   :: N
    !
    N    = size(Hk,1)
    G   = z*eye(N) - reshape_rank7_to_matrix(Si,Nineq,Nlat,Nspin,Norb) - Hk
    call inv(G)
    Gki = reshape_matrix_to_rank7(G,Nineq,Nlat,Nspin,Norb)
  end function invert_gki_normal_rank7



  !####################################################################
  !    NORMAL - TRIDIAG H
  !####################################################################
  ![N,N]
  function invert_gki_normal_tridiag_rank2(z,Hk,Si,Nsites,Ncell) result(Gki)
    integer                                                    :: Nsites,Ncell
    complex(8)                                                 :: z
    complex(8),dimension(Nsites*Ncell,Nsites*Ncell),intent(in) :: Hk
    complex(8),dimension(Nsites*Ncell,Nsites*Ncell),intent(in) :: Si
    complex(8),dimension(Nsites*Ncell,Nsites*Ncell)            :: Gki
    complex(8),dimension(Nsites,  Ncell,Ncell)                 :: D
    complex(8),dimension(Nsites-1,Ncell,Ncell)                 :: Sub
    complex(8),dimension(Nsites-1,Ncell,Ncell)                 :: Over
    complex(8),dimension(Nsites,Ncell,Ncell)                   :: Gmatrix(Nsites,Ncell,Ncell)
    integer                                                    :: N,i,isite,ic,jc,io,jo
    !
    N   = size(Hk,1)    
    call get_tridiag(Nsites,Ncell,(z*eye(N) - Si - Hk),Sub,D,Over)
    call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
    do concurrent(isite=1:Nsites,ic=1:Ncell,jc=1:Ncell)
       io = ic + (isite-1)*Ncell
       jo = jc + (isite-1)*Ncell
       Gki(io,jo) = Gmatrix(isite,ic,jc)
    enddo
  end function invert_gki_normal_tridiag_rank2

  ![Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_tridiag_rank5(z,Hk,Si,Nsites,Ncell) result(Gki)
    integer                                                        :: Nsites,Ncell
    complex(8)                                                     :: z
    complex(8),dimension(Nsites*Ncell,Nsites*Ncell),intent(in)     :: Hk
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb),intent(in)    :: Si    
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(Nsites,  Ncell,Ncell)                     :: D
    complex(8),dimension(Nsites-1,Ncell,Ncell)                     :: Sub
    complex(8),dimension(Nsites-1,Ncell,Ncell)                     :: Over
    complex(8),dimension(Nsites,Ncell,Ncell)                       :: Gmatrix(Nsites,Ncell,Ncell)
    complex(8),dimension(size(Hk,1),Nsites*Ncell)    :: G
    integer                                                        :: N,i,isite,ic,jc,io,jo
    !
    N   = size(Hk,1)
    call get_tridiag(Nsites,Ncell,( z*eye(N) - reshape_rank5_to_matrix(Si,Nlat,Nspin,Norb) - Hk ),Sub,D,Over)
    call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
    do concurrent(isite=1:Nsites,ic=1:Ncell,jc=1:Ncell)
       io = ic + (isite-1)*Ncell
       jo = jc + (isite-1)*Ncell
       G(io,jo) = Gmatrix(isite,ic,jc)
    enddo
    Gki = reshape_matrix_to_rank5(G,Nlat,Nspin,Norb)
  end function invert_gki_normal_tridiag_rank5

  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_normal_tridiag_rank7(z,Hk,Si,Nsites,Ncell) result(Gki)
    integer                                                                :: Nsites,Ncell
    complex(8)                                                             :: z
    complex(8),dimension(Nsites*Ncell,Nsites*Ncell),intent(in)             :: Hk
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb),intent(in) :: Si    
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)            :: Gki
    complex(8),dimension(Nsites,  Ncell,Ncell)                             :: D
    complex(8),dimension(Nsites-1,Ncell,Ncell)                             :: Sub
    complex(8),dimension(Nsites-1,Ncell,Ncell)                             :: Over
    complex(8),dimension(Nsites,Ncell,Ncell)                               :: Gmatrix(Nsites,Ncell,Ncell)
    complex(8),dimension(size(Hk,1),Nsites*Ncell)                          :: G
    integer                                                                :: N,i,isite,ic,jc,io,jo
    !
    N    = size(Hk,1)
    call get_tridiag(Nsites,Ncell,( z*eye(N) - reshape_rank7_to_matrix(Si,Nineq,Nlat,Nspin,Norb) - Hk ),Sub,D,Over)
    call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
    do concurrent(isite=1:Nsites,ic=1:Ncell,jc=1:Ncell)
       io = ic + (isite-1)*Ncell
       jo = jc + (isite-1)*Ncell
       G(io,jo) = Gmatrix(isite,ic,jc)
    enddo
    Gki = reshape_matrix_to_rank7(G,Nineq,Nlat,Nspin,Norb)
  end function invert_gki_normal_tridiag_rank7







  !####################################################################
  !    SUPERC - Hk
  !####################################################################
  ![2][N,N]
  function invert_gki_superc_rank2(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in) :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)   :: Hk
    complex(8),dimension(2,Ntot,Ntot)              :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)            :: Gmatrix
    integer                                        :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "invert_gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:) = Gmatrix(1:Ntot,1:Ntot)
    Gki(2,:,:) = Gmatrix(1:Ntot,Ntot+1:2*Ntot)
  end function invert_gki_superc_rank2

  ![2][Nlat,Nso,Nso]
  function invert_gki_superc_rank3(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in) :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)   :: Hk
    complex(8),dimension(2,Nlat,Nso,Nso)           :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)            :: Gmatrix
    integer                                        :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "invert_gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:) = reshape_matrix_to_rank3(Gmatrix(1:N,1:N),Nlat,Nso)
    Gki(2,:,:,:) = reshape_matrix_to_rank3(Gmatrix(1:N,N+1:2*N),Nlat,Nso)
  end function invert_gki_superc_rank3

  ![2][Nspin,Nspin,Norb,Norb]
  function invert_gki_superc_rank4(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in) :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)   :: Hk
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb)  :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)            :: Gmatrix
    integer                                        :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "invert_gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:,:) = reshape_matrix_to_rank4(Gmatrix(1:Ntot,1:Ntot),Nspin,Norb)
    Gki(2,:,:,:,:) = reshape_matrix_to_rank4(Gmatrix(1:Ntot,Ntot+1:2*Ntot),Nspin,Norb)
  end function invert_gki_superc_rank4

  ![2][Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_superc_rank5(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in)     :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)       :: Hk
    complex(8),dimension(2,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)                :: Gmatrix
    integer                                            :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "invert_gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:,:,:) = reshape_matrix_to_rank5(Gmatrix(1:Ntot,1:Ntot),Nlat,Nspin,Norb)
    Gki(2,:,:,:,:,:) = reshape_matrix_to_rank5(Gmatrix(1:Ntot,Ntot+1:2*Ntot),Nlat,Nspin,Norb)
  end function invert_gki_superc_rank5

  ![2][Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  function invert_gki_superc_rank6(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in)          :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)            :: Hk
    complex(8),dimension(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)                     :: Gmatrix
    integer                                                 :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "invert_gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:,:,:,:) = reshape_matrix_to_rank6(Gmatrix(1:Ntot,1:Ntot),Nlat,Nspin,Norb)
    Gki(2,:,:,:,:,:,:) = reshape_matrix_to_rank6(Gmatrix(1:Ntot,Ntot+1:2*Ntot),Nlat,Nspin,Norb)
  end function invert_gki_superc_rank6

  !
  ! INVERT_GK_SUPERC(_MPI)
  !
  subroutine invert_gk_superc(zeta,Hk,Gkout)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: zeta    ![2][2][N][N][Lfreq]
    complex(8),dimension(:,:,:),intent(in)          :: Hk      ![2][N][N]
    complex(8),dimension(:,:,:,:)     :: Gkout   ![2][N][N][Lfreq]
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
    complex(8),dimension(:,:,:,:)                              :: Gkout   ![2][N][N][Lfreq]
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

  ! #ifdef _MPI
  !   function MPI_Get_size(comm) result(size)
  !     integer :: comm
  !     integer :: size,ierr
  !     call MPI_Comm_size(comm,size,ierr)
  !   end function MPI_Get_size

  !   function MPI_Get_rank(comm) result(rank)
  !     integer :: comm
  !     integer :: rank,ierr
  !     call MPI_Comm_rank(comm,rank,ierr)
  !   end function MPI_Get_rank

  !   function MPI_Get_master(comm) result(master)
  !     integer :: comm
  !     logical :: master
  !     integer :: rank,ierr
  !     call MPI_Comm_rank(comm,rank,ierr)
  !     master=.false.
  !     if(rank==0)master=.true.
  !   end function MPI_Get_master
  ! #endif








  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !DBLE
  !MATRIX --> RANK3_7
  function d_matrix_TO_rank3(Hmat,Nlat,Nso) result(Hrank)
    integer                                                          :: Nso,Nlat
    real(8),dimension(Nlat*Nso,Nlat*Nso)                             :: Hmat
    real(8),dimension(Nlat,Nso,Nso)                                  :: Hrank
    integer                                                          :: ilat
    integer                                                          :: i,is
    integer                                                          :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hrank(ilat,is,js) = Hmat(i,j)
    enddo
  end function d_matrix_TO_rank3

  function d_matrix_TO_rank4(Hmat,Nspin,Norb) result(Hrank)
    integer                                                          :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb)                         :: Hmat
    real(8),dimension(Nspin,Nspin,Norb,Norb)                         :: Hrank
    integer                                                          :: iorb,ispin,is
    integer                                                          :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank4

  function d_matrix_TO_rank5(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                                          :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                    :: Hrank
    integer                                                          :: ilat
    integer                                                          :: iorb,ispin,is
    integer                                                          :: jorb,jspin,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hrank(ilat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank5

  function d_matrix_TO_rank6(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                                          :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)               :: Hrank
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank6

  function d_matrix_TO_rank7(Hmat,Nineq,Nlat,Nspin,Norb) result(Hrank)
    integer                                                          :: Nineq,Nlat,Nspin,Norb
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)   :: Hmat
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: Hrank
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function d_matrix_TO_rank7

  function d_matrix_TO_rank3L(Hmat,Nlat,Nso,L) result(Hrank)
    integer                                                          :: Nso,Nlat,L
    real(8),dimension(Nlat*Nso,Nlat*Nso,L)                           :: Hmat
    real(8),dimension(Nlat,Nso,Nso,L)                                :: Hrank
    integer                                                          :: ilat
    integer                                                          :: i,is
    integer                                                          :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hrank(ilat,is,js,:) = Hmat(i,j,:)
    enddo
  end function d_matrix_TO_rank3L

  function d_matrix_TO_rank4L(Hmat,Nspin,Norb,L) result(Hrank)
    integer                                                          :: Nspin,Norb,L
    real(8),dimension(Nspin*Norb,Nspin*Norb,L)                       :: Hmat
    real(8),dimension(Nspin,Nspin,Norb,Norb,L)                       :: Hrank
    integer                                                          :: iorb,ispin,is
    integer                                                          :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank4L

  function d_matrix_TO_rank5L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                          :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)                  :: Hrank
    integer                                                          :: iorb,ispin,ilat,is
    integer                                                          :: jorb,jspin,jlat,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hrank(ilat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank5L

  function d_matrix_TO_rank6L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                          :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)             :: Hrank
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank6L

#if __GFORTRAN__ &&  __GNUC__ > 8
  function d_matrix_TO_rank7L(Hmat,Nineq,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function d_matrix_TO_rank7L
#endif



  !##################################################################
  !##################################################################
  !##################################################################



  !RANK4_7 --> MATRIX:
  function d_rank3_TO_matrix(Hrank,Nlat,Nso) result(Hmat)
    integer                                                          :: Nso,Nlat
    real(8),dimension(Nlat*Nso,Nlat*Nso)                             :: Hmat
    real(8),dimension(Nlat,Nso,Nso)                                  :: Hrank
    integer                                                          :: ilat
    integer                                                          :: i,is
    integer                                                          :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j) = Hrank(ilat,is,js)
    enddo
  end function d_rank3_TO_matrix

  function d_rank4_TO_matrix(Hrank,Nspin,Norb) result(Hmat)
    integer                                                          :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb)                         :: Hrank
    real(8),dimension(Nspin*Norb,Nspin*Norb)                         :: Hmat
    integer                                                          :: iorb,ispin,is
    integer                                                          :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js) = Hrank(ispin,jspin,iorb,jorb)
    enddo
  end function d_rank4_TO_matrix


  function d_rank5_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                                          :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                    :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    integer                                                          :: iorb,ispin,ilat,is
    integer                                                          :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb, jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hmat(is,js) = Hrank(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function d_rank5_TO_matrix

  function d_rank6_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                                          :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)               :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function d_rank6_TO_matrix

  function d_rank7_TO_matrix(Hrank,Nineq,Nlat,Nspin,Norb) result(Hmat)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: Hrank
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)   :: Hmat
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function d_rank7_TO_matrix


  function d_rank3_TO_matrixL(Hrank,Nlat,Nso,L) result(Hmat)
    integer                                :: Nso,Nlat,L
    real(8),dimension(Nlat*Nso,Nlat*Nso,L) :: Hmat
    real(8),dimension(Nlat,Nso,Nso,L)      :: Hrank
    integer                                :: ilat
    integer                                :: i,is
    integer                                :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j,:) = Hrank(ilat,is,js,:)
    enddo
  end function d_rank3_TO_matrixL


  function d_rank4_TO_matrixL(Hrank,Nspin,Norb,L) result(Hmat)
    integer                                                          :: Nspin,Norb,L
    real(8),dimension(Nspin,Nspin,Norb,Norb,L)                       :: Hrank
    real(8),dimension(Nspin*Norb,Nspin*Norb,L)                       :: Hmat
    integer                                                          :: iorb,ispin,is
    integer                                                          :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js,:) = Hrank(ispin,jspin,iorb,jorb,:)
    enddo
  end function d_rank4_TO_matrixL

  function d_rank5_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                          :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)                  :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    integer                                                          :: iorb,ispin,ilat,is
    integer                                                          :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hmat(is,js,:) = Hrank(ilat,ispin,jspin,iorb,jorb,:)
    enddo
  end function d_rank5_TO_matrixL

  function d_rank6_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                          :: Nlat,Nspin,Norb,L
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)             :: Hrank
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js,:) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:)
    enddo
  end function d_rank6_TO_matrixL

#if __GFORTRAN__ &&  __GNUC__ > 8
  function d_rank7_TO_matrixL(Hrank,Nineq,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                          :: Nineq,Nlat,Nspin,Norb,L
    real(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    real(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    integer                                                          :: iineq
    integer                                                          :: ilat,jlat
    integer                                                          :: iorb,jorb
    integer                                                          :: ispin,jspin
    integer                                                          :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js,:) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) 
    enddo
  end function d_rank7_TO_matrixL
#endif






  !COMPLEX
  !MATRIX --> RANK3_7
  function c_matrix_TO_rank3(Hmat,Nlat,Nso) result(Hrank)
    integer                                                             :: Nso,Nlat
    complex(8),dimension(Nlat*Nso,Nlat*Nso)                             :: Hmat
    complex(8),dimension(Nlat,Nso,Nso)                                  :: Hrank
    integer                                                             :: ilat
    integer                                                             :: i,is
    integer                                                             :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nspin*Norb
       j = js + (ilat-1)*Nspin*Norb
       Hrank(ilat,is,js) = Hmat(i,j)
    enddo
  end function c_matrix_TO_rank3

  function c_matrix_TO_rank4(Hmat,Nspin,Norb) result(Hrank)
    integer                                                             :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                         :: Hmat
    complex(8),dimension(Nspin,Nspin,Norb,Norb)                         :: Hrank
    integer                                                             :: iorb,ispin,is
    integer                                                             :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank4

  function c_matrix_TO_rank5(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                                             :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                    :: Hrank
    integer                                                             :: ilat
    integer                                                             :: iorb,ispin,is
    integer                                                             :: jorb,jspin,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hrank(ilat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank5

  function c_matrix_TO_rank6(Hmat,Nlat,Nspin,Norb) result(Hrank)
    integer                                                             :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)               :: Hrank
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank6

  function c_matrix_TO_rank7(Hmat,Nineq,Nlat,Nspin,Norb) result(Hrank)
    integer                                                             :: Nineq,Nlat,Nspin,Norb
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)   :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: Hrank
    integer                                                             :: iineq
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb) = Hmat(is,js)
    enddo
  end function c_matrix_TO_rank7

  function c_matrix_TO_rank3L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                             :: Nspin,Norb,Nlat,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb,L)                  :: Hrank
    integer                                                             :: ilat
    integer                                                             :: i,is
    integer                                                             :: j,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,is=1:Nspin*Norb,js=1:Nspin*Norb)
       i = is + (ilat-1)*Nspin*Norb
       j = js + (ilat-1)*Nspin*Norb
       Hrank(ilat,is,js,:) = Hmat(i,j,:)
    enddo
  end function c_matrix_TO_rank3L

  function c_matrix_TO_rank4L(Hmat,Nspin,Norb,L) result(Hrank)
    integer                                                             :: Nspin,Norb,L
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L)                       :: Hmat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L)                       :: Hrank
    integer                                                             :: iorb,ispin,is
    integer                                                             :: jorb,jspin,js
    Hrank=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hrank(ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank4L

  function c_matrix_TO_rank5L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                             :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)                  :: Hrank
    integer                                                             :: iorb,ispin,ilat,is
    integer                                                             :: jorb,jspin,jlat,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hrank(ilat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank5L

  function c_matrix_TO_rank6L(Hmat,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                             :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)             :: Hrank
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hrank=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank6L

#if __GFORTRAN__ &&  __GNUC__ > 8
  function c_matrix_TO_rank7L(Hmat,Nineq,Nlat,Nspin,Norb,L) result(Hrank)
    integer                                                             :: Nineq,Nlat,Nspin,Norb,L
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    integer                                                             :: iineq
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hrank=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) = Hmat(is,js,:)
    enddo
  end function c_matrix_TO_rank7L
#endif



  !##################################################################
  !##################################################################
  !##################################################################



  !RANK4_7 --> MATRIX:
  function c_rank3_TO_matrix(Hrank,Nlat,Nso) result(Hmat)
    integer                                                             :: Nso,Nlat
    complex(8),dimension(Nlat*Nso,Nlat*Nso)               :: Hmat
    complex(8),dimension(Nlat,Nso,Nso)                    :: Hrank
    integer                                                             :: ilat
    integer                                                             :: i,is
    integer                                                             :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j) = Hrank(ilat,is,js)
    enddo
  end function c_rank3_TO_matrix

  function c_rank4_TO_matrix(Hrank,Nspin,Norb) result(Hmat)
    integer                                                             :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb)                         :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                         :: Hmat
    integer                                                             :: iorb,ispin,is
    integer                                                             :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js) = Hrank(ispin,jspin,iorb,jorb)
    enddo
  end function c_rank4_TO_matrix


  function c_rank5_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                                             :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)                    :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    integer                                                             :: iorb,ispin,ilat,is
    integer                                                             :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb, jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin 
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin 
       Hmat(is,js) = Hrank(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function c_rank5_TO_matrix

  function c_rank6_TO_matrix(Hrank,Nlat,Nspin,Norb) result(Hmat)
    integer                                                             :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)               :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)               :: Hmat
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function c_rank6_TO_matrix

  function c_rank7_TO_matrix(Hrank,Nineq,Nlat,Nspin,Norb) result(Hmat)
    integer                                                             :: Nineq,Nlat,Nspin,Norb
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: Hrank
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)   :: Hmat
    integer                                                             :: iineq
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb)
    enddo
  end function c_rank7_TO_matrix


  function c_rank3_TO_matrixL(Hrank,Nlat,Nso,L) result(Hmat)
    integer                                   :: Nso,Nlat,L
    complex(8),dimension(Nlat*Nso,Nlat*Nso,L) :: Hmat
    complex(8),dimension(Nlat,Nso,Nso,L)      :: Hrank
    integer                                   :: ilat
    integer                                   :: i,is
    integer                                   :: j,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,is=1:Nso,js=1:Nso)
       i = is + (ilat-1)*Nso
       j = js + (ilat-1)*Nso
       Hmat(i,j,:) = Hrank(ilat,is,js,:)
    enddo
  end function c_rank3_TO_matrixL


  function c_rank4_TO_matrixL(Hrank,Nspin,Norb,L) result(Hmat)
    integer                                                             :: Nspin,Norb,L
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L)                       :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L)                       :: Hmat
    integer                                                             :: iorb,ispin,is
    integer                                                             :: jorb,jspin,js
    Hmat=zero
    do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb
       js = jorb + (jspin-1)*Norb
       Hmat(is,js,:) = Hrank(ispin,jspin,iorb,jorb,:)
    enddo
  end function c_rank4_TO_matrixL

  function c_rank5_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                             :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)                  :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    integer                                                             :: iorb,ispin,ilat,is
    integer                                                             :: jorb,jspin,jlat,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Hmat(is,js,:) = Hrank(ilat,ispin,jspin,iorb,jorb,:)
    enddo
  end function c_rank5_TO_matrixL

  function c_rank6_TO_matrixL(Hrank,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                             :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)             :: Hrank
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L)             :: Hmat
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hmat=zero
    do concurrent(ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
       js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
       Hmat(is,js,:) = Hrank(ilat,jlat,ispin,jspin,iorb,jorb,:)
    enddo
  end function c_rank6_TO_matrixL

#if __GFORTRAN__ &&  __GNUC__ > 8
  function c_rank7_TO_matrixL(Hrank,Nineq,Nlat,Nspin,Norb,L) result(Hmat)
    integer                                                             :: Nineq,Nlat,Nspin,Norb,L
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,L)       :: Hrank
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb,L) :: Hmat
    integer                                                             :: iineq
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: ispin,jspin
    integer                                                             :: is,js
    Hmat=zero
    do concurrent(iineq=1:Nineq,ilat=1:Nlat,jlat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ilat-1)*Norb + (ispin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       js = jorb + (jlat-1)*Norb + (jspin-1)*Nlat*Norb + (iineq-1)*Nlat*Nspin*Norb
       Hmat(is,js,:) = Hrank(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:) 
    enddo
  end function c_rank7_TO_matrixL
#endif












  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
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
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
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
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
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
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
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











!   ! INVERT_GK_NORMAL(_MPI)
!   !SERIAL (OR PARALLEL ON K):
!   subroutine invert_gk_normal(zeta,Hk,Gkout)
!     complex(8),dimension(:,:,:),intent(in)    :: zeta    ![N,N,Lfreq]
!     complex(8),dimension(:,:),intent(in)      :: Hk      ![N,N]
!     complex(8),dimension(:,:,:),intent(inout) :: Gkout   ![N,N,Lfreq]
!     complex(8)                                :: Gmatrix(size(Hk,1),size(Hk,2))
!     integer                                   :: Ntot,Lfreq
!     integer                                   :: i
!     !
!     Ntot  = size(Hk,1)
!     Lfreq = size(zeta,3)
!     do i=1,Lfreq
!        Gmatrix  = zeta(:,:,i) - Hk
!        call inv(Gmatrix)
!        Gkout(:,:,i) = Gmatrix
!     enddo
!   end subroutine invert_gk_normal

!   !PARALLEL ON FREQ:
!   subroutine invert_gk_normal_mpi(zeta,Hk,Gkout)
!     complex(8),dimension(:,:,:),intent(in)    :: zeta    ![N,N,Lfreq]
!     complex(8),dimension(:,:),intent(in)      :: Hk      ![N,N]
!     complex(8),dimension(:,:,:),intent(inout) :: Gkout   ![N,N,Lfreq]
!     complex(8)                                :: Gktmp(size(Hk,1),size(Hk,2),size(zeta,3))   !as Gkout
!     complex(8)                                :: Gmatrix(size(Hk,1),size(Hk,2))
!     integer                                   :: Ntot,Lfreq
!     integer                                   :: i
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
!     Ntot  = size(Hk,1)
!     Lfreq = size(zeta,3)
!     !
!     Gktmp=zero
!     do i=1+mpi_rank,Lfreq,mpi_size
!        Gmatrix  = zeta(:,:,i) - Hk
!        call inv(Gmatrix)
!        Gktmp(:,:,i) = Gmatrix
!     enddo
!     !
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
!   end subroutine invert_gk_normal_mpi




!   !
!   ! INVERT_GK_NORMAL_TRIDIAG(_MPI)
!   !SERIAL (OR PARALLEL ON K)
!   subroutine invert_gk_normal_tridiag(zeta,Hk,Gkout,Nsites,Ncell)
!     complex(8),dimension(:,:,:),intent(in)     :: zeta    ![N,N,Lfreq]
!     complex(8),dimension(:,:),intent(in)       :: Hk      ![N,N]
!     complex(8),dimension(:,:,:),intent(inout)  :: Gkout   ![N,N,Lfreq]
!     integer                                    :: Nsites,Ncell
!     complex(8),dimension(Nsites,  Ncell,Ncell) :: D
!     complex(8),dimension(Nsites-1,Ncell,Ncell) :: Sub
!     complex(8),dimension(Nsites-1,Ncell,Ncell) :: Over
!     complex(8),dimension(Nsites,Ncell,Ncell)   :: Gmatrix(Nsites,Ncell,Ncell)
!     integer                                    :: Ntot,Lfreq
!     integer                                    :: i,isite,ic,jc,io,jo
!     !
!     Ntot  = size(Hk,1)
!     Lfreq = size(zeta,3)
!     !
!     Gkout=zero
!     do i=1,Lfreq      
!        call get_tridiag(Nsites,Ncell,(zeta(:,:,i)-Hk),Sub,D,Over)
!        call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
!        do isite=1,Nsites
!           do ic=1,Ncell
!              do jc=1,Ncell
!                 io = ic + (isite-1)*Ncell
!                 jo = jc + (isite-1)*Ncell
!                 Gkout(io,jo,i) = Gmatrix(isite,ic,jc)
!              enddo
!           enddo
!        enddo
!     enddo
!   end subroutine invert_gk_normal_tridiag

!   !PARALLEL ON FREQ:
!   subroutine invert_gk_normal_tridiag_mpi(zeta,Hk,Gkout,Nsites,Ncell)
!     complex(8),dimension(:,:,:),intent(in)     :: zeta    ![N,N,Lfreq]
!     complex(8),dimension(:,:),intent(in)       :: Hk      ![N,N]
!     complex(8),dimension(:,:,:),intent(inout)  :: Gkout   ![N,N,Lfreq]
!     integer                                    :: Nsites,Ncell
!     complex(8),dimension(Nsites,  Ncell,Ncell) :: D
!     complex(8),dimension(Nsites-1,Ncell,Ncell) :: Sub
!     complex(8),dimension(Nsites-1,Ncell,Ncell) :: Over
!     complex(8)                                 :: Gktmp(size(Hk,1),size(Hk,2),size(zeta,3))   !as Gkout
!     complex(8),dimension(Nsites,Ncell,Ncell)   :: Gmatrix(Nsites,Ncell,Ncell)
!     integer                                    :: Ntot,Lfreq
!     integer                                    :: i,isite,ic,jc,io,jo
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
!     !
!     Ntot  = size(Hk,1)
!     Lfreq = size(zeta,3)
!     !
!     Gktmp=zero
!     do i=1+mpi_rank,Lfreq,mpi_size
!        call get_tridiag(Nsites,Ncell,(zeta(:,:,i)-Hk),Sub,D,Over)
!        call inv_tridiag(Nsites,Ncell,-Sub,D,-Over,Gmatrix)
!        do isite=1,Nsites
!           do ic=1,Ncell
!              do jc=1,Ncell
!                 io = ic + (isite-1)*Ncell
!                 jo = jc + (isite-1)*Ncell
!                 Gktmp(io,jo,i) = Gmatrix(isite,ic,jc)
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
!   end subroutine invert_gk_normal_tridiag_mpi
