module SC_COMMON
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only:reg,txtfy
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


  character(len=128)                  :: suffix
  character(len=128)                  :: weiss_suffix=".dat"
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
    call get_ctrl_var(xmu,"XMU")
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    integer                                       :: Nspin,Norb,L
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Hrank
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Hmat
    integer                                       :: iorb,ispin,is
    integer                                       :: jorb,jspin,js
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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
    Hmat=zero
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








end module SC_COMMON
