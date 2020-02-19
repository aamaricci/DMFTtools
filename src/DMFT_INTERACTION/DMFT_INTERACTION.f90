module DMFT_INTERACTION
  USE SF_TIMER
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_IOTOOLS,   only:reg,txtfy
  USE SF_ARRAYS,    only:arange
  USE SF_LINALG,    only:eye,inv,diag,kronecker_product
  USE SF_MISC,      only:assert_shape
  USE SF_IOTOOLS,   only:free_unit,reg,txtfy
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  interface dmft_interaction_setKanamori
     module procedure :: dmft_interaction_setKanamori_serial
#ifdef _MPI
     module procedure :: dmft_interaction_setKanamori_mpi
#endif
  end interface dmft_interaction_setKanamori


  interface dmft_interaction_rotate
     module procedure :: dmft_interaction_rotate_serial
#ifdef _MPI
     module procedure :: dmft_interaction_rotate_mpi
#endif
  end interface dmft_interaction_rotate


  public :: dmft_interaction_setKanamori
  public :: dmft_interaction_rotate

  public :: dmft_interaction_print
  public :: dmft_interaction_read




  !##################################################################
  !##################################################################
  !##################################################################



  interface tensor2matrix_reshape
     module procedure  :: d_tensor2matrix_reshape
     module procedure  :: c_tensor2matrix_reshape
  end interface tensor2matrix_reshape

  interface matrix2tensor_reshape
     module procedure d_matrix2tensor_reshape
     module procedure c_matrix2tensor_reshape
  end interface matrix2tensor_reshape


  integer                          :: mpi_ierr
  integer                          :: mpi_rank
  integer                          :: mpi_size
  logical                          :: mpi_master



contains


  !--------------------------------------------------------------------!
  !PURPOSE: Set/rotate the rank-4 tensor to the standard Kanamori interaction
  !--------------------------------------------------------------------!

  include "dmft_interaction_serial.f90"
#ifdef _MPI
  include "dmft_interaction_mpi.f90"
#endif



  !--------------------------------------------------------------------!
  !                          INTERACTION I/O
  !--------------------------------------------------------------------!
  subroutine dmft_interaction_read(Utensor,file)
    implicit none
    real(8),dimension(:,:,:,:,:),intent(inout)  :: Utensor
    character(len=*),intent(in)                 :: file
    !aux
    real(8),dimension(:,:,:),allocatable        :: Umatrix
    integer                                     :: Norb
    integer                                     :: LOGfile
    integer                                     :: io,jo
    integer                                     :: unit
    !
    call get_ctrl_var(Norb,"NORB")
    !call get_ctrl_var(LOGfile,"LOGFILE") This is not working for some reason
    LOGfile=6
    !
    call assert_shape(Utensor,[Norb,Norb,Norb,Norb,2],"dmft_interaction_read","Utensor")
    allocate(Umatrix(Norb*Norb,Norb*Norb,2));Umatrix=0d0
    !
    write(LOGfile,*)
    write(LOGfile,"(A)")"Read interaction tensor from file: "//reg(file)
    open(free_unit(unit),file=reg(file))
    !
    read(unit,*)
    do io=1,Norb*Norb
      read(unit,"(90F21.12)") (Umatrix(io,jo,1),jo=1,Norb*Norb)
    enddo
    read(unit,*)
    read(unit,*)
    do io=1,Norb*Norb
      read(unit,"(90F21.12)") (Umatrix(io,jo,2),jo=1,Norb*Norb)
    enddo
    close(unit)
    !
    Utensor=matrix2tensor_reshape(Umatrix,Norb)
    !
    deallocate(Umatrix)
    !
  end subroutine dmft_interaction_read



  subroutine dmft_interaction_print(Utensor,file)
    implicit none
    real(8),dimension(:,:,:,:,:),intent(in)     :: Utensor
    character(len=*),optional                   :: file
    !aux
    real(8),dimension(:,:,:),allocatable        :: Umatrix
    integer                                     :: Norb
    integer                                     :: LOGfile
    integer                                     :: io,jo
    integer                                     :: unit
    !
    call get_ctrl_var(Norb,"NORB")
    !call get_ctrl_var(LOGfile,"LOGFILE") This is not working for some reason
    LOGfile=6
    !
    call assert_shape(Utensor,[Norb,Norb,Norb,Norb,2],"dmft_interaction_print","Utensor")
    allocate(Umatrix(Norb*Norb,Norb*Norb,2));Umatrix=0d0
    Umatrix=tensor2matrix_reshape(Utensor,Norb)
    !
    write(LOGfile,*)
    write(LOGfile,"(A)")"Interaction tensor - Block (up,up)(dw,dw):"
    call print_interaction_states(LOGfile,Norb)
    do io=1,Norb*Norb
       write(LOGfile,"(90F8.3)") (Umatrix(io,jo,1),jo=1,Norb*Norb)
    enddo
    write(LOGfile,*)
    write(LOGfile,"(A)")"Interaction tensor - Block (up,up)(up,up):"
    call print_interaction_states(LOGfile,Norb)
    do io=1,Norb*Norb
       write(LOGfile,"(90F8.3)") (Umatrix(io,jo,2),jo=1,Norb*Norb)
    enddo
    write(LOGfile,*)
    !
    if(present(file))then
      !
       open(free_unit(unit),file=reg(file),action='write')
       write(LOGfile,"(A)")"Print interaction tensor on file: "//reg(file)
       !
       call print_interaction_states(unit,Norb)
       do io=1,Norb*Norb
          write(unit,"(90F21.12)") (Umatrix(io,jo,1),jo=1,Norb*Norb)
       enddo
       write(unit,*)
       call print_interaction_states(unit,Norb)
       do io=1,Norb*Norb
          write(unit,"(90F21.12)") (Umatrix(io,jo,2),jo=1,Norb*Norb)
       enddo
       close(unit)
       !
    endif
    !
    deallocate(Umatrix)
    !
  end subroutine dmft_interaction_print



  !--------------------------------------------------------------------!
  !                         AUXILIARY STUFF
  !--------------------------------------------------------------------!
  subroutine print_interaction_states(unit,Norb)
    implicit none
    integer,intent(in)                         :: Norb
    integer,intent(in)                         :: unit
    character(len=32)                          :: fmt1,fmt2
    character(len=32)                          :: fmt
    integer                                    :: LOGfile
    !
    !call get_ctrl_var(LOGfile,"LOGFILE") This is not working for some reason
    LOGfile=6
    !
    if(unit==LOGfile)then
      if(Norb==3) write(unit,"(90a8)") "(11)","(22)","(33)", &
                                       "(12)","(21)",        &
                                       "(13)","(31)",        &
                                       "(23)","(32)"
      if(Norb==2) write(unit,"(90a8)") "(11)","(22)", &
                                       "(12)","(21)"
    else
      if(Norb==3) write(unit,"(90a21)")"(11)","(22)","(33)", &
                                       "(12)","(21)",        &
                                       "(13)","(31)",        &
                                       "(23)","(32)"
      if(Norb==2) write(unit,"(90a21)")"(11)","(22)", &
                                       "(12)","(21)"
    endif
  end subroutine print_interaction_states


  function c_tensor2matrix_reshape(tensor,Norb) result(matrix)
    implicit none
    integer                                     :: Norb
    complex(8),dimension(Norb,Norb,Norb,Norb,2) :: tensor
    complex(8),dimension(Norb*Norb,Norb*Norb,2) :: matrix
    integer                                     :: iorb,jorb,korb,lorb
    integer                                     :: io,jo,spin
    matrix = zero
    do iorb=1,Norb
      do jorb=1,Norb
          do korb=1,Norb
             do lorb=1,Norb
                !
                call getUstride(io,iorb,jorb,Norb)
                call getUstride(jo,korb,lorb,Norb)
                !
                do spin=1,2
                   matrix(io,jo,spin) = tensor(iorb,jorb,korb,lorb,spin)
                enddo
                !
             enddo
          enddo
       enddo
    enddo
  end function c_tensor2matrix_reshape

  function d_tensor2matrix_reshape(tensor,Norb) result(matrix)
    implicit none
    integer                                     :: Norb
    real(8),dimension(Norb,Norb,Norb,Norb,2)    :: tensor
    real(8),dimension(Norb*Norb,Norb*Norb,2)    :: matrix
    integer                                     :: iorb,jorb,korb,lorb
    integer                                     :: io,jo,spin
    matrix = 0d0
    do iorb=1,Norb
      do jorb=1,Norb
          do korb=1,Norb
             do lorb=1,Norb
                !
                call getUstride(io,iorb,jorb,Norb)
                call getUstride(jo,korb,lorb,Norb)
                !
                do spin=1,2
                   matrix(io,jo,spin) = tensor(iorb,jorb,korb,lorb,spin)
                enddo
                !
             enddo
          enddo
       enddo
    enddo
  end function d_tensor2matrix_reshape

  function c_matrix2tensor_reshape(matrix,Norb) result(tensor)
    implicit none
    integer                                     :: Norb
    complex(8),dimension(Norb,Norb,Norb,Norb,2) :: tensor
    complex(8),dimension(Norb*Norb,Norb*Norb,2) :: matrix
    integer                                     :: iorb,jorb,korb,lorb
    integer                                     :: io,jo,spin
    tensor = zero
    do io=1,Norb*Norb
       do jo=1,Norb*Norb
          !
          call getMstride(io,iorb,jorb,Norb)
          call getMstride(jo,korb,lorb,Norb)
          !
          do spin=1,2
             tensor(iorb,jorb,korb,lorb,spin) = matrix(io,jo,spin)
          enddo
          !
       enddo
    enddo
  end function c_matrix2tensor_reshape

  function d_matrix2tensor_reshape(matrix,Norb) result(tensor)
    implicit none
    integer                                     :: Norb
    real(8),dimension(Norb,Norb,Norb,Norb,2)    :: tensor
    real(8),dimension(Norb*Norb,Norb*Norb,2)    :: matrix
    integer                                     :: iorb,jorb,korb,lorb
    integer                                     :: io,jo,spin
    tensor = 0d0
    do io=1,Norb*Norb
       do jo=1,Norb*Norb
          !
          call getMstride(io,iorb,jorb,Norb)
          call getMstride(jo,korb,lorb,Norb)
          !
          do spin=1,2
             tensor(iorb,jorb,korb,lorb,spin) = matrix(io,jo,spin)
          enddo
          !
       enddo
    enddo
  end function d_matrix2tensor_reshape


  subroutine getUstride(ndx,iorb,jorb,Norb)
    implicit none
    integer,intent(in)                          :: iorb,jorb,Norb
    integer,intent(out)                         :: ndx
    integer,dimension(Norb,Norb)                :: stride
    integer                                     :: count
    integer                                     :: io,jo
    !
    count=0
    do io=1,Norb
       count=count+1
       stride(io,io)=count
    enddo
    do io=1,Norb
       do jo=io+1,Norb
          count=count+1
          stride(io,jo)=count
          count=count+1
          stride(jo,io)=count
       enddo
    enddo
    !
    ndx = stride(iorb,jorb)
    !
  end subroutine getUstride

  subroutine getMstride(ndx,iorb,jorb,Norb)
    implicit none
    integer,intent(in)                          :: ndx,Norb
    integer,intent(out)                         :: iorb,jorb
    integer,dimension(Norb,Norb)                :: stride
    integer                                     :: count
    integer                                     :: io,jo
    !
    count=0
    do io=1,Norb
       count=count+1
       stride(io,io)=count
    enddo
    do io=1,Norb
       do jo=io+1,Norb
          count=count+1
          stride(io,jo)=count
          count=count+1
          stride(jo,io)=count
       enddo
    enddo
    !
    strideloop:do iorb=1,Norb
       do jorb=1,Norb
          if(stride(iorb,jorb)==ndx) exit strideloop
       enddo
    enddo strideloop
    !
  end subroutine getMstride



end module DMFT_INTERACTION
