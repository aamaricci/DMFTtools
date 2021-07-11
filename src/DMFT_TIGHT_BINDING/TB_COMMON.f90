module TB_COMMON
  USE SF_CONSTANTS, only: pi,pi2,xi,one,zero
  USE SF_ARRAYS, only: linspace
  USE SF_IOTOOLS
  USE SF_LINALG, only: eigh,det,eye,zeros,eig,diag,operator(.x.),diagonal,zeros,ones
  USE SF_COLORS
  USE SF_SPECIAL, only: fermi
  USE SF_TIMER, only:start_timer,stop_timer,eta
  USE SF_MISC, only: assert_shape,sort_array
  USE SF_OPTIMIZE,only: fmin_cgminimize,fzero
  USE DMFT_CTRL_VARS
  USE DMFT_GFIO
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  interface add_to
     module procedure :: add_to_A1
     module procedure :: add_to_A2
     module procedure :: add_to_A3
  end interface add_to


  interface TB_reshape_array
     module procedure :: tb_reorder_vec_d,tb_reorder_vec_c
     module procedure :: tb_reorder_mat_d,tb_reorder_mat_c
     module procedure :: tb_reorder_hk_d,tb_reorder_hk_c
  end interface TB_reshape_array
  interface TB_reorder_array
     module procedure :: tb_reorder_vec_d,tb_reorder_vec_c
     module procedure :: tb_reorder_mat_d,tb_reorder_mat_c
     module procedure :: tb_reorder_hk_d,tb_reorder_hk_c
  end interface TB_reorder_array
  interface TB_reshape_hk
     module procedure :: tb_reorder_hk_d,tb_reorder_hk_c
  end interface TB_reshape_hk
  interface TB_reorder_hk
     module procedure :: tb_reorder_hk_d,tb_reorder_hk_c
  end interface TB_reorder_hk

  interface tb_findloc
     module procedure :: findloc_char
     module procedure :: findloc_int
  end interface tb_findloc


  !Some special points in the BZ:
  !we do everything in 3d.
  real(8),dimension(3),parameter :: kpoint_gamma=[0,0,0]*pi
  real(8),dimension(3),parameter :: kpoint_x1=[1,0,0]*pi
  real(8),dimension(3),parameter :: kpoint_x2=[0,1,0]*pi
  real(8),dimension(3),parameter :: kpoint_x3=[0,0,1]*pi
  real(8),dimension(3),parameter :: kpoint_m1=[1,1,0]*pi
  real(8),dimension(3),parameter :: kpoint_m2=[0,1,1]*pi
  real(8),dimension(3),parameter :: kpoint_m3=[1,0,1]*pi
  real(8),dimension(3),parameter :: kpoint_r=[1,1,1]*pi


  real(8),dimension(3),save      :: ei_x=[1d0,0d0,0d0]
  real(8),dimension(3),save      :: ei_y=[0d0,1d0,0d0]
  real(8),dimension(3),save      :: ei_z=[0d0,0d0,1d0]

  real(8),dimension(3),save      :: bk_x=[1d0,0d0,0d0]*pi2
  real(8),dimension(3),save      :: bk_y=[0d0,1d0,0d0]*pi2
  real(8),dimension(3),save      :: bk_z=[0d0,0d0,1d0]*pi2
  real(8),dimension(3),save      :: BZ_origin=[0d0,0d0,0d0]

  logical,save                   :: io_eivec=.false.
  logical,save                   :: io_bkvec=.false.
  logical,save                   :: set_eivec=.false.
  logical,save                   :: set_bkvec=.false.

  integer                        :: mpi_ierr
  integer                        :: mpi_rank
  integer                        :: mpi_size
  logical                        :: mpi_master

  private :: tb_findloc
  private :: add_to_A1,add_to_A2,add_to_A3
  private :: tb_reorder_vec_d,tb_reorder_vec_c
  private :: tb_reorder_mat_d,tb_reorder_mat_c
  private :: findloc_char
  private :: findloc_int
  private :: indx_reorder

contains


  function indices2i(ivec,Nvec) result(istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end function indices2i

  function i2indices(istate,Nvec) result(ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end function i2indices






  function tb_reorder_vec_d(Huser,Nin,OrderIn,OrderOut) result(Hss)
    real(8),dimension(:)                   :: Huser ![Nlat*Nspin*Norb]
    integer,dimension(3)                   :: Nin   !In sequence of Nlat,Nspin,Norb as integers
    character(len=*),dimension(3)          :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
    !
    real(8),dimension(size(Huser))         :: Hss
    integer,dimension(3)                   :: Ivec,Jvec
    integer,dimension(3)                   :: IndexOut
    integer,dimension(3)                   :: Nout
    integer                                :: iss,iuser,i,Nlso
    !
    character(len=*),dimension(3),optional :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
    character(len=5),dimension(3)          :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
    OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
    if(present(OrderOut))then
       do i=1,2
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser)
    !
    !Construct an index array InderOut corresponding to the Out ordering.
    !This is a permutation of the In ordering taken as [1,2,3].
    !For each entry in OrderIn we look for the position of the
    !corresponding entry in OrderOut using tb_findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       IndexOut(i:i)=tb_findloc(OrderIn,OrderOut_(i))
    enddo
    if(any(IndexOut==0))then
       print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iuser = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          !
          Hss(iss) = Huser(iuser)
       enddo
    else
       Hss = Huser
    endif
    return
  end function tb_reorder_vec_d

  function tb_reorder_vec_c(Huser,Nin,OrderIn,OrderOut) result(Hss)
    complex(8),dimension(:)                :: Huser ![Nlat*Nspin*Norb]
    integer,dimension(3)                   :: Nin   !In sequence of Nlat,Nspin,Norb as integers
    character(len=*),dimension(3)          :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
    !
    complex(8),dimension(size(Huser))      :: Hss
    integer,dimension(3)                   :: Ivec,Jvec
    integer,dimension(3)                   :: IndexOut
    integer,dimension(3)                   :: Nout
    integer                                :: iss,iuser,i,Nlso
    !
    character(len=*),dimension(3),optional :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
    character(len=5),dimension(3)          :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
    !
    OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
    if(present(OrderOut))then
       do i=1,2
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser)
    !
    !Construct an index array InderOut corresponding to the Out ordering.
    !This is a permutation of the In ordering taken as [1,2,3].
    !For each entry in OrderIn we look for the position of the
    !corresponding entry in OrderOut using tb_findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       IndexOut(i:i)=tb_findloc(OrderIn,OrderOut_(i))
    enddo
    if(any(IndexOut==0))then
       print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iuser = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          !
          Hss(iss) = Huser(iuser)
       enddo
    else
       Hss = Huser
    endif
    return
  end function tb_reorder_vec_c

  function tb_reorder_mat_d(Huser,Nin,OrderIn,OrderOut) result(Hss)
    real(8),dimension(:,:)                         :: Huser
    integer,dimension(3)                           :: Nin   !In sequence of Nlat,Nspin,Norb as integers
    character(len=*),dimension(3)                  :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
    !
    real(8),dimension(size(Huser,1),size(Huser,2)) :: Hss
    integer,dimension(3)                           :: Ivec,Jvec
    integer,dimension(3)                           :: IndexOut
    integer,dimension(3)                           :: Nout
    integer                                        :: iss,jss,iuser,juser,i,Nlso
    character(len=*),dimension(3),optional         :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
    character(len=5),dimension(3)                  :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
    !
    OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
    if(present(OrderOut))then
       do i=1,2
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser,1)
    call assert_shape(Huser,[Nlso,Nlso],"tb_reorder_mat_c","Huser")
    !
    !Construct an index array InderOut corresponding to the Out ordering.
    !This is a permutation of the In ordering taken as [1,2,3].
    !For each entry in OrderIn we look for the position of the
    !corresponding entry in OrderOut using tb_findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       IndexOut(i:i)=tb_findloc(OrderIn,OrderOut_(i))
    enddo
    if(any(IndexOut==0))then
       print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iuser = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do jss=1,Nlso
             Ivec  = i2indices(jss,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             juser = indices2i(Jvec,Nout)
             !
             Hss(iss,jss) = Huser(iuser,juser)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function tb_reorder_mat_d

  function tb_reorder_mat_c(Huser,Nin,OrderIn,OrderOut) result(Hss)
    complex(8),dimension(:,:)                         :: Huser
    integer,dimension(3)                              :: Nin   !In sequence of Nlat,Nspin,Norb as integers
    character(len=*),dimension(3)                     :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
    !
    complex(8),dimension(size(Huser,1),size(Huser,2)) :: Hss
    integer,dimension(3)                              :: Ivec,Jvec
    integer,dimension(3)                              :: IndexOut
    integer,dimension(3)                              :: Nout
    integer                                           :: iss,jss,iuser,juser,i,Nlso
    character(len=*),dimension(3),optional            :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
    character(len=5),dimension(3)                     :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
    !
    OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
    if(present(OrderOut))then
       do i=1,2
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !

    Nlso = size(Huser,1)
    call assert_shape(Huser,[Nlso,Nlso],"tb_reorder_mat_c","Huser")
    !
    !Construct an index array InderOut corresponding to the Out ordering.
    !This is a permutation of the In ordering taken as [1,2,3].
    !For each entry in OrderIn we look for the position of the
    !corresponding entry in OrderOut using tb_findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       IndexOut(i:i)=tb_findloc(OrderIn,OrderOut_(i))
    enddo
    if(any(IndexOut==0))then
       print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iuser = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do jss=1,Nlso
             Ivec  = i2indices(jss,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             juser = indices2i(Jvec,Nout)
             !
             Hss(iss,jss) = Huser(iuser,juser)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function tb_reorder_mat_c




  function tb_reorder_hk_d(Huser,Nin,OrderIn,OrderOut) result(Hss)
    real(8),dimension(:,:,:)                                     :: Huser
    integer,dimension(3)                                         :: Nin   !In sequence of Nlat,Nspin,Norb as integers
    character(len=*),dimension(3)                                :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
    !
    real(8),dimension(size(Huser,1),size(Huser,2),size(Huser,3)) :: Hss
    integer,dimension(3)                                         :: Ivec,Jvec
    integer,dimension(3)                                         :: IndexOut
    integer,dimension(3)                                         :: Nout
    integer                                                      :: iss,jss,iuser,juser,i,Nlso,Nk
    !
    character(len=*),dimension(3),optional                       :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
    character(len=5),dimension(3)                                :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
    !
    OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
    if(present(OrderOut))then
       do i=1,2
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Nlso,Nlso,Nk],"tb_reorder_hk_d","Huser")
    !
    !Construct an index array InderOut corresponding to the Out ordering.
    !This is a permutation of the In ordering taken as [1,2,3].
    !For each entry in OrderIn we look for the position of the
    !corresponding entry in OrderOut using tb_findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       IndexOut(i:i)=tb_findloc(OrderIn,OrderOut_(i))
    enddo
    if(any(IndexOut==0))then
       print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iuser = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do jss=1,Nlso
             Ivec  = i2indices(jss,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             juser = indices2i(Jvec,Nout)
             !
             Hss(iss,jss,:) = Huser(iuser,juser,:)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function tb_reorder_hk_d

  function tb_reorder_hk_c(Huser,Nin,OrderIn,OrderOut) result(Hss)
    complex(8),dimension(:,:,:)                                     :: Huser
    integer,dimension(3)                                            :: Nin   !In sequence of Nlat,Nspin,Norb as integers
    character(len=*),dimension(3)                                   :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
    !
    complex(8),dimension(size(Huser,1),size(Huser,2),size(Huser,3)) :: Hss
    integer,dimension(3)                                            :: Ivec,Jvec
    integer,dimension(3)                                            :: IndexOut
    integer,dimension(3)                                            :: Nout
    integer                                                         :: iss,jss,iuser,juser,i,Nlso,Nk
    character(len=*),dimension(3),optional                          :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
    character(len=5),dimension(3)                                   :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
    !
    OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
    if(present(OrderOut))then
       do i=1,2
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !   
    Nlso = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Nlso,Nlso,Nk],"tb_reorder_hk_d","Huser")
    !
    !Construct an index array InderOut corresponding to the Out ordering.
    !This is a permutation of the In ordering taken as [1,2,3].
    !For each entry in OrderIn we look for the position of the
    !corresponding entry in OrderOut using tb_findloc.
    !If 0 entries exist, corresponding components are not found. stop. 
    do i=1,3     
       IndexOut(i:i)=tb_findloc(OrderIn,OrderOut_(i))
    enddo
    if(any(IndexOut==0))then
       print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iss=1,Nlso
          Ivec  = i2indices(iss,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iuser = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do jss=1,Nlso
             Ivec  = i2indices(jss,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             juser = indices2i(Jvec,Nout)
             !
             Hss(iss,jss,:) = Huser(iuser,juser,:)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function tb_reorder_hk_c





  subroutine TB_slo2lso_model(Hk,Nlat,Nspin,Norb)
    complex(8),dimension(:,:,:),intent(inout)             :: Hk
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Htmp
    integer                                               :: Nlat,Nspin,Norb,Nlso,Nk,ik
    Nlso = size(Hk,1)
    Nk   = size(Hk,3)
    call assert_shape(Hk,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk],"TB_slo2lso","Hk")
    do ik=1,Nk
       Htmp = slo2lso(Hk(:,:,ik),Nlat,Nspin,Norb)
       Hk(:,:,ik) = Htmp
    enddo
  end subroutine TB_slo2lso_model


  function slo2lso(Hslo,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb),intent(in) :: Hslo
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)            :: Hlso
    integer                                                          :: Nlat,Nspin,Norb
    integer                                                          :: iorb,ispin,ilat
    integer                                                          :: jorb,jspin,jlat
    integer                                                          :: Islo,Jslo
    integer                                                          :: Ilso,Jlso
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      !
                      Islo = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      jslo = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      !
                      Ilso = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      Jlso = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      !
                      Hlso(Ilso,Jlso) = Hslo(Islo,Jslo)
                      !
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function slo2lso

  function lso2slo(Hlso,Nlat,Nspin,Norb) result(Hslo)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb),intent(in) :: Hlso
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)            :: Hslo
    integer                                                          :: Nlat,Nspin,Norb
    integer                                                          :: iorb,ispin,ilat
    integer                                                          :: jorb,jspin,jlat
    integer                                                          :: Islo,Jslo
    integer                                                          :: Ilso,Jlso
    Hslo=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      !
                      Islo = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      jslo = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      !
                      Ilso = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      Jlso = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      !
                      Hslo(Islo,Jslo) = Hlso(Ilso,Jlso)
                      !
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2slo




  function indx_reorder(Ain,Index)  result(Aout)
    integer,dimension(:)         :: Ain
    integer,dimension(size(Ain)) :: Index
    integer,dimension(size(Ain)) :: Aout
    integer                        :: i
    forall(i=1:size(Ain))Aout(Index(i)) = Ain(i)
  end function indx_reorder


  subroutine add_to_A1(vec,val)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),intent(in)                             :: val  
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_A1



  subroutine add_to_A2(vec,val)
    real(8),dimension(:,:),allocatable,intent(inout) :: vec
    real(8),intent(in),dimension(size(vec,2))        :: val  
    real(8),dimension(:,:),allocatable               :: tmp
    integer                                          :: n,ndim
    !
    ndim = size(vec,2)
    !
    if (allocated(vec)) then
       n = size(vec,1)
       allocate(tmp(n+1,ndim))
       tmp(1:n,:) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n,ndim))
    end if
    !
    !Put val as last entry:
    vec(n,:) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_A2


  subroutine add_to_A3(vec,val)
    real(8),dimension(:,:,:),allocatable,intent(inout)    :: vec
    real(8),dimension(size(vec,2),size(vec,3)),intent(in) :: val  
    real(8),dimension(:,:,:),allocatable                  :: tmp
    integer                                               :: n,n2,n3
    !
    n2 = size(vec,2)
    n3 = size(vec,3)
    !
    if (allocated(vec)) then
       n = size(vec,1)
       allocate(tmp(n+1,n2,n3))
       tmp(1:n,:,:) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n,n2,n3))
    end if
    !
    !Put val as last entry:
    vec(n,:,:) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_A3









  function findloc_char(array,val) result(pos)
    character(len=*),dimension(:) :: array
    character(len=*)              :: val
    integer                       :: pos,i
    pos=0
    do i=1,size(array)
       if(array(i)==val)then
          pos = i
          exit
       endif
    enddo
    return
  end function findloc_char

  function findloc_int(array,val) result(pos)
    integer,dimension(:) :: array
    integer              :: val
    integer              :: pos,i
    pos=0
    do i=1,size(array)
       if(array(i)==val)then
          pos = i
          exit
       endif
    enddo
    return
  end function findloc_int




END MODULE TB_COMMON











