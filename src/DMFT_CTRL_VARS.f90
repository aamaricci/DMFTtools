module DMFT_CTRL_VARS
  implicit none
  private

  type ctrl_var
     integer,pointer             :: i
     real(8),pointer             :: d
     logical,pointer             :: l
     character(len=:),pointer    :: ch
  end type ctrl_var

  type ctrl_node
     type(ctrl_var),dimension(:),allocatable :: var
     character(len=3)                        :: type
     character(len=100)                      :: name
     type(ctrl_node),pointer                 :: next !link to next box
  end type ctrl_node

  type ctrl_list
     logical                 :: status=.false.
     integer                 :: size 
     type(ctrl_node),pointer :: root
  end type ctrl_list


  interface add_ctrl_var
     module procedure i_add_ctrl_var
     module procedure d_add_ctrl_var
     module procedure l_add_ctrl_var
     module procedure iv_add_ctrl_var
     module procedure dv_add_ctrl_var
     module procedure lv_add_ctrl_var
     module procedure ch_add_ctrl_var
  end interface add_ctrl_var


  interface get_ctrl_var
     module procedure i_get_ctrl_variable
     module procedure d_get_ctrl_variable
     module procedure l_get_ctrl_variable
     module procedure iv_get_ctrl_variable
     module procedure dv_get_ctrl_variable
     module procedure lv_get_ctrl_variable
     module procedure ch_get_ctrl_variable
  end interface get_ctrl_var



  public :: ctrl_list

  public :: init_ctrl_list
  public :: destroy_ctrl_list
  public :: size_ctrl_list
  public :: print_ctrl_list

  public :: add_ctrl_var
  public :: get_ctrl_var


  type(ctrl_list)              :: dmft_ctrl_list

  character(len=:),allocatable :: p_buffer
  
  character(len=:),allocatable :: name_

  !LOCAL VERSION OF TXTFY//STR
  interface txtfy
     module procedure i_to_ch,r_to_ch,c_to_ch,l_to_ch
  end interface txtfy



contains  





  !+------------------------------------------------------------------+
  !PURPOSE: init the input list
  !+------------------------------------------------------------------+
  subroutine init_ctrl_list(list)
    type(ctrl_list),optional :: list
    if(present(list))then
       allocate(list%root)    
       list%size=0
       list%status=.true.
       list%root%next=>null()
    else
       allocate(dmft_ctrl_list%root)    
       dmft_ctrl_list%size=0
       dmft_ctrl_list%status=.true.
       dmft_ctrl_list%root%next=>null()
    endif
  end subroutine init_ctrl_list





  !+------------------------------------------------------------------+
  !PURPOSE: delete the list
  !+------------------------------------------------------------------+
  subroutine destroy_ctrl_list(list)
    type(ctrl_list),optional :: list
    type(ctrl_node),pointer  :: p,c
    if(present(list))then
       do
          p => list%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       list%status=.false.
    else
       do
          p => dmft_ctrl_list%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       dmft_ctrl_list%status=.false.
    endif
  end subroutine destroy_ctrl_list




  !+------------------------------------------------------------------+
  !PURPOSE: get list size
  !+------------------------------------------------------------------+
  function size_ctrl_list(list) result(size)
    type(ctrl_list),optional :: list
    integer                   :: size
    size=dmft_ctrl_list%size
    if(present(list))size=list%size
  end function size_ctrl_list






  !+------------------------------------------------------------------+
  !PURPOSE: !Append input data to the list:
  !+------------------------------------------------------------------+
  !========================SCALAR==================================
  subroutine i_add_ctrl_var(variable,name,default,list)
    integer,target           :: variable
    integer,optional         :: default
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    type(ctrl_node),pointer  :: p,c
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(1))
    p%next%var(1)%i=>variable
    p%next%name= name_
    p%next%type='i'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine i_add_ctrl_var

  subroutine d_add_ctrl_var(variable,name,default,list)
    real(8),target           :: variable
    real(8),optional         :: default
    character(len=*)         :: name
    type(ctrl_list),optional  :: list
    type(ctrl_node),pointer :: p,c
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(1))
    p%next%var(1)%d=>variable
    p%next%name= name_
    p%next%type='d'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine d_add_ctrl_var

  subroutine l_add_ctrl_var(variable,name,default,list)
    logical,target           :: variable
    logical,optional         :: default
    character(len=*)         :: name
    type(ctrl_list),optional  :: list
    type(ctrl_node),pointer :: p,c
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(1))
    p%next%var(1)%l=>variable
    p%next%name= name_
    p%next%type='l'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine l_add_ctrl_var


  !========================VECTOR==================================
  subroutine iv_add_ctrl_var(variable,name,default,list)
    integer,dimension(:),target   :: variable
    integer,dimension(:),optional :: default
    character(len=*)              :: name
    type(ctrl_list),optional      :: list
    type(ctrl_node),pointer       :: p,c
    integer                       :: i
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(size(variable)))
    do i=1,size(variable)
       p%next%var(i)%i=>variable(i)
    enddo
    p%next%name= name_
    p%next%type='i'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine iv_add_ctrl_var

  subroutine dv_add_ctrl_var(variable,name,default,list)
    real(8),dimension(:),target   :: variable
    real(8),dimension(:),optional :: default
    character(len=*)              :: name
    type(ctrl_list),optional      :: list
    type(ctrl_node),pointer       :: p,c
    integer                       :: i
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(size(variable)))
    do i=1,size(variable)
       p%next%var(i)%d=>variable(i)
    enddo
    p%next%name= name_
    p%next%type='d'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine dv_add_ctrl_var

  subroutine lv_add_ctrl_var(variable,name,default,list)
    logical,dimension(:),target   :: variable
    logical,dimension(:),optional :: default
    character(len=*)              :: name
    type(ctrl_list),optional      :: list
    type(ctrl_node),pointer       :: p,c
    integer                       :: i
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(size(variable)))
    do i=1,size(variable)
       p%next%var(i)%l=>variable(i)
    enddo
    p%next%name= name_
    p%next%type='l'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine lv_add_ctrl_var


  !========================STRING==================================
  subroutine ch_add_ctrl_var(variable,name,default,list)
    character(len=*),target   :: variable
    character(len=*),optional :: default
    character(len=*)          :: name
    type(ctrl_list),optional  :: list
    type(ctrl_node),pointer   :: p,c
    name_=name;call upper_case(name_)
    if(present(default))variable=default
    if(present(list))then
       if(.not.list%status)call init_ctrl_list(list)
       p => list%root
    else
       if(.not.dmft_ctrl_list%status)call init_ctrl_list()
       p => dmft_ctrl_list%root
    endif
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list or beginning of the list
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    allocate(p%next%var(1))
    nullify(p%next%var(1)%ch)
    p%next%var(1)%ch=> variable
    p%next%name= name_
    p%next%type='ch'
    !
    if(present(list))then
       list%size=list%size+1
    else
       dmft_ctrl_list%size=dmft_ctrl_list%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
    p=>null()
    c=>null()
  end subroutine ch_add_ctrl_var














  !+------------------------------------------------------------------+
  !PURPOSE:   !Get input variable from the list:
  !+------------------------------------------------------------------+
  !========================0-dimension==================================
  subroutine i_get_ctrl_variable(variable,name,list)
    integer                  :: variable
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             variable = c%var(1)%i
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine i_get_ctrl_variable

  subroutine d_get_ctrl_variable(variable,name,list)
    real(8)                  :: variable
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             variable = c%var(1)%d
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine d_get_ctrl_variable

  subroutine l_get_ctrl_variable(variable,name,list)
    logical                  :: variable
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             variable = c%var(1)%l
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine l_get_ctrl_variable


  !========================1-dimension==================================
  subroutine iv_get_ctrl_variable(variable,name,list)
    integer,dimension(:)     :: variable
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             if(size(variable)/=size(c%var))write(*,"(A)")"get_ctrl_variable warning: variable has wrong dimensions"
             do i=1,size(variable)
                variable(i) = c%var(i)%i
             enddo
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine iv_get_ctrl_variable

  subroutine dv_get_ctrl_variable(variable,name,list)
    real(8),dimension(:)      :: variable
    character(len=*)          :: name
    type(ctrl_list),optional :: list
    integer                   :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             if(size(variable)/=size(c%var))write(*,"(A)")"get_ctrl_variable warning: variable has wrong dimensions"
             do i=1,size(variable)
                variable(i) = c%var(i)%d
             enddo
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine dv_get_ctrl_variable

  subroutine lv_get_ctrl_variable(variable,name,list)
    logical,dimension(:)     :: variable
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             if(size(variable)/=size(c%var))write(*,"(A)")"get_ctrl_variable warning: variable has wrong dimensions"
             do i=1,size(variable)
                variable(i) = c%var(i)%l
             enddo
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine lv_get_ctrl_variable


  !========================STRING==================================
  subroutine ch_get_ctrl_variable(variable,name,list)
    character(len=*)         :: variable
    character(len=*)         :: name
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size_
    type(ctrl_node),pointer  :: c
    character(len=len(name)) :: name_
    name_=name;call upper_case(name_)
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size_=dmft_ctrl_list%size
    if(present(list))size_=list%size
    if(size_>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          if(trim(c%name)==trim(name_))then
             variable = c%var(1)%ch
             c=>null()
             return
          endif
          c => c%next
       enddo
       write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
    else
       write(*,"(A)")"get ctrl variable: empty. Stopping"
       stop
    endif
    c => null()
  end subroutine ch_get_ctrl_variable








  !+------------------------------------------------------------------+
  !PURPOSE: print the list to file
  !+------------------------------------------------------------------+
  subroutine print_ctrl_list(list)
    type(ctrl_list),optional :: list
    integer                  :: i,counter,size
    type(ctrl_node),pointer  :: c
    if(present(list))then
       c => list%root%next
    else
       c => dmft_ctrl_list%root%next
    endif
    counter = 0 
    size=dmft_ctrl_list%size
    if(present(list))size=list%size
    if(size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          call print_ctrl_node(c)
          c => c%next
       enddo
    else
       write(*,*)"print ctrl list: empty"
       return
    endif
    c => null()
  end subroutine print_ctrl_list
  !---------------------------------------------------------------------
  subroutine print_ctrl_node(c)
    type(ctrl_node) :: c
    integer         :: i
    !
    call s_blank_delete(c%name)
    select case(c%type)
    case('ch')
       p_buffer=trim(c%name)//"="//trim(adjustl(trim(c%var(1)%ch)))

    case('i')
       if(size(c%var)==1)then   !scalar
          p_buffer=trim(c%name)//"="//txtfy(c%var(1)%i)
       else                     !vector
          p_buffer=trim(c%name)//"="
          do i=1,size(c%var)-1
             p_buffer=trim(p_buffer)//trim(txtfy(c%var(i)%i))//","
          end do
          p_buffer=trim(p_buffer)//trim(txtfy(c%var(size(c%var))%i))
       endif

    case('d')
       if(size(c%var)==1)then   !scalar
          p_buffer=trim(c%name)//"="//txtfy(c%var(1)%d)
       else                     !vector
          p_buffer=trim(c%name)//"="
          do i=1,size(c%var)-1
             p_buffer=trim(p_buffer)//trim(txtfy(c%var(i)%d))//","
          end do
          p_buffer=trim(p_buffer)//trim(txtfy(c%var(size(c%var))%d))
       endif

    case('l')
       if(size(c%var)==1)then   !scalar
          p_buffer=trim(c%name)//"="//txtfy(c%var(1)%l)
       else                     !vector
          p_buffer=trim(c%name)//"="
          do i=1,size(c%var)-1
             p_buffer=trim(p_buffer)//trim(txtfy(c%var(i)%l))//","
          end do
          p_buffer=trim(p_buffer)//trim(txtfy(c%var(size(c%var))%l))
       endif
    end select
    !
    call s_blank_delete(p_buffer)
    !
    write(*,"(1x,A)")trim(p_buffer)
    p_buffer=""
  end subroutine print_ctrl_node











  !+------------------------------------------------------------------+
  !PURPOSE: ANCILLARY routines
  !+------------------------------------------------------------------+
  !Auxiliary routines:
  subroutine upper_case(s)
    character              ch
    integer   ( kind = 4 ) i
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len_trim ( s )
    do i = 1, s_length
       ch = s(i:i)
       call ch_cap ( ch )
       s(i:i) = ch
    end do
  end subroutine upper_case

  subroutine lower_case(s)
    integer   ( kind = 4 ) i
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len_trim ( s )
    do i = 1, s_length
       call ch_low ( s(i:i) )
    end do
  end subroutine lower_case

  subroutine ch_cap(ch)
    character              ch
    integer   ( kind = 4 ) itemp
    itemp = iachar ( ch )
    if ( 97 <= itemp .and. itemp <= 122 ) then
       ch = achar ( itemp - 32 )
    end if
  end subroutine ch_cap

  subroutine ch_low ( ch )
    character ch
    integer ( kind = 4 ) i
    i = iachar ( ch )
    if ( 65 <= i .and. i <= 90 ) then
       ch = achar ( i + 32 )
    end if
  end subroutine ch_low


  subroutine s_blank_delete ( s )
    !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
    !    All TAB characters are also removed.
    !    Input/output, character ( len = * ) S, the string to be transformed.
    implicit none
    character              ch
    integer   ( kind = 4 ) get
    integer   ( kind = 4 ) put
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    character, parameter :: tab = achar ( 9 )
    put = 0
    s_length = len_trim ( s )
    do get = 1, s_length
       ch = s(get:get)
       if ( ch /= ' ' .and. ch /= tab ) then
          put = put + 1
          s(put:put) = ch
       end if
    end do
    s(put+1:s_length) = ' '
    return
  end subroutine s_blank_delete


  function i_to_ch(i4) result(string)
    character(len=32) :: string
    integer           :: i4
    call i4_to_s_left(i4,string)
  end function i_to_ch

  function r_to_ch(r8) result(string)
    character(len=32) :: string
    character(len=16) :: string_
    real(8)           :: r8
    call r8_to_s_left(r8,string_)
    string=adjustl(string_)
  end function r_to_ch

  function c_to_ch(c) result(string)
    character(len=32+3) :: string
    character(len=16) :: sre,sim
    complex(8)        :: c
    real(8)           :: re,im
    re=real(c,8);im=aimag(c)
    call r8_to_s_left(re,sre)
    call r8_to_s_left(im,sim)
    string="("//trim(sre)//","//trim(sim)//")"
  end function c_to_ch

  function l_to_ch(bool) result(string)
    logical :: bool
    character(len=1) :: string
    string='F'
    if(bool)string='T'
  end function l_to_ch

  subroutine i4_to_s_left ( i4, s )
    character :: c
    integer   :: i
    integer   :: i4
    integer   :: idig
    integer   :: ihi
    integer   :: ilo
    integer   :: ipos
    integer   :: ival
    character(len=*) ::  s
    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !  Make a copy of the integer.
    ival = i4
    !  Handle the negative sign.
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !  The absolute value of the integer goes into S(ILO:IHI).
    ipos = ihi
    !  Find the last digit of IVAL, strip it off, and stick it into the string.
    do
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
       if ( ival == 0 ) then
          exit
       end if
    end do
    !  Shift the string to the left.
    s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
    s(ilo+ihi-ipos:ihi) = ' '
  end subroutine i4_to_s_left

  subroutine r8_to_s_left ( r8, s )
    integer :: i
    real(8) :: r8
    integer :: s_length
    character(len=*) ::  s
    character(len=16) :: s2
    s_length = len ( s )
    if ( s_length < 16 ) then
       do i = 1, s_length
          s(i:i) = '*'
       end do
    else if ( r8 == 0.0D+00 ) then
       s(1:16) = '     0.d0     '
    else
       if(abs(r8) < 1.d0)then
          write ( s2, '(ES16.9)' ) r8
       else
          write ( s2, '(F16.9)' ) r8
       endif
       s(1:16) = s2
    end if
    !  Shift the string left.
    s = adjustl ( s )
  end subroutine r8_to_s_left


  subroutine digit_to_ch(digit,ch)
    character :: ch
    integer   :: digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
  end subroutine digit_to_ch






end module DMFT_CTRL_VARS
