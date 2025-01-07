module GF_COMMON
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS
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




  interface sc_zi
     module procedure :: sc_zi_rank2
     module procedure :: sc_zi_rank3
     module procedure :: sc_zi_rank4
     module procedure :: sc_zi_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: sc_zi_rank6
#endif
  end interface sc_zi



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
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - from_rank3(S(1,:,:,:,i))
       zi(1,2,:,:) =                          - from_rank3(S(2,:,:,:,i))
       zi(2,1,:,:) =                          - from_rank3(S(2,:,:,:,i))
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + from_rank3(conjg(S(1,:,:,:,i)))
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - from_rank3(S(1,:,:,:,i))
       zi(1,2,:,:) =                                         - from_rank3(S(2,:,:,:,i))
       zi(2,1,:,:) =                                         - from_rank3(S(2,:,:,:,i))
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + from_rank3(conjg(S(1,:,:,:,Lfreq+1-i)))
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
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - from_rank4(S(1,:,:,:,:,i))
       zi(1,2,:,:) =                          - from_rank4(S(2,:,:,:,:,i))
       zi(2,1,:,:) =                          - from_rank4(S(2,:,:,:,:,i))
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + from_rank4(conjg(S(1,:,:,:,:,i)))
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - from_rank4(S(1,:,:,:,:,i))
       zi(1,2,:,:) =                                         - from_rank4(S(2,:,:,:,:,i))
       zi(2,1,:,:) =                                         - from_rank4(S(2,:,:,:,:,i))
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + from_rank4(conjg(S(1,:,:,:,:,Lfreq+1-i)))
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
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - from_rank5(S(1,:,:,:,:,:,i))
       zi(1,2,:,:) =                          - from_rank5(S(2,:,:,:,:,:,i))
       zi(2,1,:,:) =                          - from_rank5(S(2,:,:,:,:,:,i))
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + from_rank5(conjg(S(1,:,:,:,:,:,i)))
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - from_rank5(S(1,:,:,:,:,:,i))
       zi(1,2,:,:) =                                         - from_rank5(S(2,:,:,:,:,:,i))
       zi(2,1,:,:) =                                         - from_rank5(S(2,:,:,:,:,:,i))
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + from_rank5(conjg(S(1,:,:,:,:,:,Lfreq+1-i)))
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
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot) - from_rank6(S(1,:,:,:,:,:,:,i))
       zi(1,2,:,:) =                          - from_rank6(S(2,:,:,:,:,:,:,i))
       zi(2,1,:,:) =                          - from_rank6(S(2,:,:,:,:,:,:,i))
       zi(2,2,:,:) = (wfreq(i)-xmu)*eye(Ntot) + from_rank6(conjg(S(1,:,:,:,:,:,:,i)))
    case("r","R")
       zi(1,1,:,:) = (wfreq(i)+xmu)*eye(Ntot)                - from_rank6(S(1,:,:,:,:,:,:,i))
       zi(1,2,:,:) =                                         - from_rank6(S(2,:,:,:,:,:,:,i))
       zi(2,1,:,:) =                                         - from_rank6(S(2,:,:,:,:,:,:,i))
       zi(2,2,:,:) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Ntot)  + from_rank6(conjg(S(1,:,:,:,:,:,:,Lfreq+1-i)))
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
  function gki_normal(z,Hk,Si) result(Gki)
    complex(8)                                             :: z
    complex(8),dimension(:,:),intent(in)                   :: Hk
    complex(8),dimension(size(Hk,1),size(Hk,1)),intent(in) :: Si
    complex(8),dimension(size(Hk,1),size(Hk,1))            :: Gki
    integer                                                :: N
    !
    N   = size(Hk,1)
    Gki = z*eye(N) - Si - Hk
    call inv(Gki)
  end function gki_normal


  !####################################################################
  !    NORMAL - TRIDIAG H
  !####################################################################
  ![N,N]
  function gki_tridiag(z,Hk,Si,Nsites,Ncell) result(Gki)
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
  end function gki_tridiag






  !####################################################################
  !    SUPERC - Hk
  !####################################################################
  ![2][N,N]
  function gki_superc_rank2(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in) :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)   :: Hk
    complex(8),dimension(2,Ntot,Ntot)              :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)            :: Gmatrix
    integer                                        :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:) = Gmatrix(1:Ntot,1:Ntot)
    Gki(2,:,:) = Gmatrix(1:Ntot,Ntot+1:2*Ntot)
  end function gki_superc_rank2

  ![2][Nlat,Nso,Nso]
  function gki_superc_rank3(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in) :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)   :: Hk
    complex(8),dimension(2,Nlat,Nso,Nso)           :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)            :: Gmatrix
    integer                                        :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:) = to_rank3(Gmatrix(1:N,1:N))
    Gki(2,:,:,:) = to_rank3(Gmatrix(1:N,N+1:2*N))
  end function gki_superc_rank3

  ![2][Nspin,Nspin,Norb,Norb]
  function gki_superc_rank4(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in) :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)   :: Hk
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb)  :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)            :: Gmatrix
    integer                                        :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:,:) = to_rank4(Gmatrix(1:Ntot,1:Ntot))
    Gki(2,:,:,:,:) = to_rank4(Gmatrix(1:Ntot,Ntot+1:2*Ntot))
  end function gki_superc_rank4

  ![2][Nlat,Nspin,Nspin,Norb,Norb]
  function gki_superc_rank5(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in)     :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)       :: Hk
    complex(8),dimension(2,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)                :: Gmatrix
    integer                                            :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:,:,:) = to_rank5(Gmatrix(1:Ntot,1:Ntot))
    Gki(2,:,:,:,:,:) = to_rank5(Gmatrix(1:Ntot,Ntot+1:2*Ntot))
  end function gki_superc_rank5

  ![2][Nlat,Nlat,Nspin,Nspin,Norb,Norb]
  function gki_superc_rank6(Z,Hk) result(Gki)
    complex(8),dimension(2,2,Ntot,Ntot),intent(in)          :: Z
    complex(8),dimension(2,Ntot,Ntot),intent(in)            :: Hk
    complex(8),dimension(2,Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Gki
    complex(8),dimension(2*Ntot,2*Ntot)                     :: Gmatrix
    integer                                                 :: N
    !
    N = size(Hk,2);if(N/=Ntot)stop "gki_superc_rank_2 error: size(Hk,1)!=Ntot"
    Gmatrix(1:N,1:N)         = Z(1,1,:,:) - Hk(1,:,:)
    Gmatrix(1:N,N+1:2*N)     = Z(1,2,:,:) 
    Gmatrix(N+1:2*N,1:N)     = Z(2,1,:,:)
    Gmatrix(N+1:2*N,N+1:2*N) = Z(2,2,:,:) - Hk(2,:,:)
    call inv(Gmatrix)
    Gki(1,:,:,:,:,:,:) = to_rank6(Gmatrix(1:Ntot,1:Ntot))
    Gki(2,:,:,:,:,:,:) = to_rank6(Gmatrix(1:Ntot,Ntot+1:2*Ntot))
  end function gki_superc_rank6












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
    complex(8),dimension(Nineq*Nlat*Nspin*Norb,Nineq*Nlat*Nspin*Norb)          :: Hmat
    complex(8),dimension(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb)                :: Hrank
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







end module GF_COMMON

