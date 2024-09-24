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


  interface TB_reshuffle
     module procedure :: tb_reorder_vec_i
     module procedure :: tb_reorder_vec_d
     module procedure :: tb_reorder_vec_c
     module procedure :: tb_reorder_mat_i
     module procedure :: tb_reorder_mat_d
     module procedure :: tb_reorder_mat_c
     module procedure :: tb_reorder_hk_i
     module procedure :: tb_reorder_hk_d
     module procedure :: tb_reorder_hk_c
  end interface TB_reshuffle

  interface TB_reorder
     module procedure :: legacy_reorder_vec_d
     module procedure :: legacy_reorder_vec_c
     module procedure :: legacy_reorder_mat_d
     module procedure :: legacy_reorder_mat_c
  end interface TB_reorder


  interface TB_reshuffle_hk
     module procedure :: tb_reorder_hk_i
     module procedure :: tb_reorder_hk_d
     module procedure :: tb_reorder_hk_c
  end interface TB_reshuffle_hk

  interface TB_reorder_hk
     module procedure :: legacy_reorder_hk_d
     module procedure :: legacy_reorder_hk_c
  end interface TB_reorder_hk

  interface tb_findloc
     module procedure :: findloc_char
     module procedure :: findloc_int
  end interface tb_findloc


  interface cross2d
     module procedure :: cross_2d_d
     module procedure :: cross_2d_c
  end interface cross2d


  interface cross3d
     module procedure :: cross_3d_d
     module procedure :: cross_3d_c
  end interface cross3d


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


  real(8),dimension(3),save      :: ei_1=[1d0,0d0,0d0]
  real(8),dimension(3),save      :: ei_2=[0d0,1d0,0d0]
  real(8),dimension(3),save      :: ei_3=[0d0,0d0,1d0]

  real(8),dimension(3),save      :: bk_1=[1d0,0d0,0d0]*pi2
  real(8),dimension(3),save      :: bk_2=[0d0,1d0,0d0]*pi2
  real(8),dimension(3),save      :: bk_3=[0d0,0d0,1d0]*pi2
  !
  real(8),dimension(3),save      :: BZ_origin=[0d0,0d0,0d0]

  ! logical,save                   :: io_eivec=.false.
  ! logical,save                   :: io_bkvec=.false.
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




  function cross_2d_d(a,b) result(c)
    real(8),dimension(2) :: a,b
    real(8)              :: c
    c = a(1)*b(2) - a(2)*b(1)
  end function cross_2d_d
  function cross_2d_c(a,b) result(c)
    complex(8),dimension(2) :: a,b
    complex(8)              :: c
    c = a(1)*b(2) - a(2)*b(1)
  end function cross_2d_c
  !
  function cross_3d_d(a,b) result(c)
    real(8),dimension(3) :: a,b
    real(8),dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_3d_d
  function cross_3d_c(a,b) result(c)
    complex(8),dimension(3) :: a,b
    complex(8),dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_3d_c





  subroutine tb_reorder_vec_i(Huser,Nvec,IndexOut)
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    integer,dimension(:)               :: Huser   
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)               :: Nvec
    integer,dimension(:)               :: IndexOut!pi[1,..,D]
    !
    integer,dimension(size(Huser))     :: Hss
    integer,dimension(:,:),allocatable :: Nin,Nout
    integer,dimension(:),allocatable   :: Ivec,Jvec
    integer,dimension(:),allocatable   :: shift
    integer                            :: Ntot,D,M,DM
    integer                            :: iss
    integer                            :: iuser
    integer                            :: i
    integer                            :: im
    integer                            :: id
    !
    Ntot = size(Huser)
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          Hss(iss) = Huser(iuser)
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_vec_i

  subroutine tb_reorder_vec_d(Huser,Nvec,IndexOut)
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    real(8),dimension(:)               :: Huser   
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)               :: Nvec
    integer,dimension(:)               :: IndexOut!pi[1,..,D]
    !
    real(8),dimension(size(Huser))     :: Hss
    integer,dimension(:,:),allocatable :: Nin,Nout
    integer,dimension(:),allocatable   :: Ivec,Jvec
    integer,dimension(:),allocatable   :: shift
    integer                            :: Ntot,D,M,DM
    integer                            :: iss
    integer                            :: iuser
    integer                            :: i
    integer                            :: im
    integer                            :: id
    !
    Ntot = size(Huser)
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          Hss(iss) = Huser(iuser)
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_vec_d

  subroutine tb_reorder_vec_c(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    complex(8),dimension(:)            :: Huser   
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)               :: Nvec
    integer,dimension(:)               :: IndexOut!pi[1,..,D]
    !
    complex(8),dimension(size(Huser))  :: Hss
    integer,dimension(:,:),allocatable :: Nin,Nout
    integer,dimension(:),allocatable   :: Ivec,Jvec
    integer,dimension(:),allocatable   :: shift
    integer                            :: Ntot,D,M,DM
    integer                            :: iss
    integer                            :: iuser
    integer                            :: i
    integer                            :: im
    integer                            :: id
    !
    Ntot = size(Huser)
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          Hss(iss) = Huser(iuser)
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_vec_c







  subroutine tb_reorder_mat_i(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    integer,dimension(:,:)                         :: Huser 
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)                           :: Nvec
    integer,dimension(:)                           :: IndexOut!pi[1,..,D]
    !
    integer,dimension(size(Huser,1),size(Huser,2)) :: Hss
    integer,dimension(:,:),allocatable             :: Nin,Nout
    integer,dimension(:),allocatable               :: Ivec,Jvec
    integer,dimension(:),allocatable               :: Shift
    integer                                        :: Ntot,D,M,DM
    integer                                        :: iss,jss
    integer                                        :: iuser,juser
    integer                                        :: i,j
    integer                                        :: im,jm
    integer                                        :: id,jd
    !
    Ntot = size(Huser,1)
    call assert_shape(Huser,[Ntot,Ntot],"tb_reorder_mat","Huser")
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          !
          do jm=1,M
             do j=1,product(Nin(jm,:))
                juser = shift(jm)+j
                Ivec  = i2indices(j,Nin(jm,:))            !Map iuser to Ivec by Nin
                Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
                jss   = shift(jm)+indices2i(Jvec,Nout(jm,:))  !Map Jvec to iss by Nout
                !
                Hss(iss,jss) = Huser(iuser,juser)
             enddo
          enddo
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_mat_i

  subroutine tb_reorder_mat_d(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    real(8),dimension(:,:)                         :: Huser 
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)                           :: Nvec
    integer,dimension(:)                           :: IndexOut!pi[1,..,D]
    !
    real(8),dimension(size(Huser,1),size(Huser,2)) :: Hss
    integer,dimension(:,:),allocatable             :: Nin,Nout
    integer,dimension(:),allocatable               :: Ivec,Jvec
    integer,dimension(:),allocatable               :: Shift
    integer                                        :: Ntot,D,M,DM
    integer                                        :: iss,jss
    integer                                        :: iuser,juser
    integer                                        :: i,j
    integer                                        :: im,jm
    integer                                        :: id,jd
    !
    Ntot = size(Huser,1)
    call assert_shape(Huser,[Ntot,Ntot],"tb_reorder_mat","Huser")
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          !
          do jm=1,M
             do j=1,product(Nin(jm,:))
                juser = shift(jm)+j
                Ivec  = i2indices(j,Nin(jm,:))            !Map iuser to Ivec by Nin
                Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
                jss   = shift(jm)+indices2i(Jvec,Nout(jm,:))  !Map Jvec to iss by Nout
                !
                Hss(iss,jss) = Huser(iuser,juser)
             enddo
          enddo
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_mat_d

  subroutine tb_reorder_mat_c(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    complex(8),dimension(:,:)                         :: Huser 
    complex(8),dimension(size(Huser,1),size(Huser,2)) :: Hss
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)                           :: Nvec
    integer,dimension(:)                           :: IndexOut!pi[1,..,D]
    !
    integer,dimension(:,:),allocatable             :: Nin,Nout
    integer,dimension(:),allocatable               :: Ivec,Jvec
    integer,dimension(:),allocatable               :: Shift
    integer                                        :: Ntot,D,M,DM
    integer                                        :: iss,jss
    integer                                        :: iuser,juser
    integer                                        :: i,j
    integer                                        :: im,jm
    integer                                        :: id,jd
    !
    Ntot = size(Huser,1)
    call assert_shape(Huser,[Ntot,Ntot],"tb_reorder_mat","Huser")
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          !
          do jm=1,M
             do j=1,product(Nin(jm,:))
                juser = shift(jm)+j
                Ivec  = i2indices(j,Nin(jm,:))            !Map iuser to Ivec by Nin
                Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
                jss   = shift(jm)+indices2i(Jvec,Nout(jm,:))  !Map Jvec to iss by Nout
                !
                Hss(iss,jss) = Huser(iuser,juser)
             enddo
          enddo
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_mat_c





  subroutine tb_reorder_hk_i(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    integer,dimension(:,:,:)                                     :: Huser 
    integer,dimension(size(Huser,1),size(Huser,2),size(Huser,3)) :: Hss
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)                                         :: Nvec
    integer,dimension(:)                                         :: IndexOut!pi[1,..,D]
    !
    integer,dimension(:,:),allocatable                           :: Nin,Nout
    integer,dimension(:),allocatable                             :: Ivec,Jvec
    integer,dimension(:),allocatable                             :: Shift
    integer                                                      :: Ntot,D,M,DM,Nk
    integer                                                      :: iss,jss
    integer                                                      :: iuser,juser
    integer                                                      :: i,j
    integer                                                      :: im,jm
    integer                                                      :: id,jd
    !
    Ntot = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Ntot,Ntot,Nk],"tb_reorder_hk","Huser")
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          !
          do jm=1,M
             do j=1,product(Nin(jm,:))
                juser = shift(jm)+j
                Ivec  = i2indices(j,Nin(jm,:))            !Map iuser to Ivec by Nin
                Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
                jss   = shift(jm)+indices2i(Jvec,Nout(jm,:))  !Map Jvec to iss by Nout
                !
                Hss(iss,jss,:) = Huser(iuser,juser,:)
             enddo
          enddo
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_hk_i

  subroutine tb_reorder_hk_d(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    real(8),dimension(:,:,:)                                     :: Huser 
    real(8),dimension(size(Huser,1),size(Huser,2),size(Huser,3)) :: Hss
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)                                         :: Nvec
    integer,dimension(:)                                         :: IndexOut!pi[1,..,D]
    !
    integer,dimension(:,:),allocatable                           :: Nin,Nout
    integer,dimension(:),allocatable                             :: Ivec,Jvec
    integer,dimension(:),allocatable                             :: Shift
    integer                                                      :: Ntot,D,M,DM,Nk
    integer                                                      :: iss,jss
    integer                                                      :: iuser,juser
    integer                                                      :: i,j
    integer                                                      :: im,jm
    integer                                                      :: id,jd
    !
    Ntot = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Ntot,Ntot,Nk],"tb_reorder_hk","Huser")
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          !
          do jm=1,M
             do j=1,product(Nin(jm,:))
                juser = shift(jm)+j
                Ivec  = i2indices(j,Nin(jm,:))            !Map iuser to Ivec by Nin
                Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
                jss   = shift(jm)+indices2i(Jvec,Nout(jm,:))  !Map Jvec to iss by Nout
                !
                Hss(iss,jss,:) = Huser(iuser,juser,:)
             enddo
          enddo
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_hk_d

  subroutine tb_reorder_hk_c(Huser,Nvec,IndexOut) 
    ![Ntot],!M-elements,D-DegreesOfFreedom(Nspin,Norb,Nlat)
    complex(8),dimension(:,:,:)                                     :: Huser 
    complex(8),dimension(size(Huser,1),size(Huser,2),size(Huser,3)) :: Hss
    ![M*D]=[Nvec_1,...,Nvec_D],Nvec_i=[N_i1,...,N_iM]
    integer,dimension(:)                                         :: Nvec
    integer,dimension(:)                                         :: IndexOut!pi[1,..,D]
    !
    integer,dimension(:,:),allocatable                           :: Nin,Nout
    integer,dimension(:),allocatable                             :: Ivec,Jvec
    integer,dimension(:),allocatable                             :: Shift
    integer                                                      :: Ntot,D,M,DM,Nk
    integer                                                      :: iss,jss
    integer                                                      :: iuser,juser
    integer                                                      :: i,j
    integer                                                      :: im,jm
    integer                                                      :: id,jd
    !
    Ntot = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Ntot,Ntot,Nk],"tb_reorder_hk","Huser")
    DM   = size(Nvec)
    D    = size(IndexOut)
    M    = DM/D
    if(mod(DM,D)/=0)stop "tb_reorder_vec error: size(Nin) is not multiple of size(IndexOut)"
    !
    if(all(IndexOut==(/(i,i=1,D)/) ))return
    !
    allocate(Nin(M,D),Nout(M,D))
    do im=1,M
       do id=1,D
          Nin(im,id) = Nvec(im+(id-1)*M)
       enddo
       Nout(im,:)=indx_reorder(Nin(im,:),IndexOut)
    enddo
    !
    allocate(shift(M))
    shift(1)=0
    do im=2,M
       shift(im)=shift(m-1)+product(Nin(im-1,:))
    enddo
    !
    allocate(Ivec(D),Jvec(D))
    !From IndexOut we can re-order the dimensions array to get the User dimensions array
    do im=1,M
       do i=1,product(Nin(im,:))
          iuser = shift(im)+i
          Ivec  = i2indices(i,Nin(im,:))            !Map iuser to Ivec by Nin
          Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
          iss   = shift(im)+indices2i(Jvec,Nout(im,:))  !Map Jvec to iss by Nout
          !
          do jm=1,M
             do j=1,product(Nin(jm,:))
                juser = shift(jm)+j
                Ivec  = i2indices(j,Nin(jm,:))            !Map iuser to Ivec by Nin
                Jvec  = indx_reorder(Ivec,IndexOut)       !Reorder according to Out ordering Ivec-->Jvec
                jss   = shift(jm)+indices2i(Jvec,Nout(jm,:))  !Map Jvec to iss by Nout
                !
                Hss(iss,jss,:) = Huser(iuser,juser,:)
             enddo
          enddo
       enddo
    enddo
    Huser = Hss
    return
  end subroutine tb_reorder_hk_c












  !VEC
  function legacy_reorder_vec_d(Huser,Nin,OrderIn,OrderOut) result(Hss)
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
       do i=1,3
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
       print*,"Legacy_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iuser=1,Nlso
          Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iss = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          !
          Hss(iss) = Huser(iuser)
       enddo
    else
       Hss = Huser
    endif
    return
  end function legacy_reorder_vec_d

  function legacy_reorder_vec_c(Huser,Nin,OrderIn,OrderOut) result(Hss)
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
       do i=1,3
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
       print*,"Legacy_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iuser=1,Nlso
          Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iss = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          !
          Hss(iss) = Huser(iuser)
       enddo
    else
       Hss = Huser
    endif
    return
  end function legacy_reorder_vec_c






  !MAT
  function legacy_reorder_mat_d(Huser,Nin,OrderIn,OrderOut) result(Hss)
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
       do i=1,3
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser,1)
    call assert_shape(Huser,[Nlso,Nlso],"legacy_reorder_mat_c","Huser")
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
       print*,"Legacy_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iuser=1,Nlso
          Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iss   = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do juser=1,Nlso
             Ivec  = i2indices(juser,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             jss = indices2i(Jvec,Nout)
             !
             Hss(iss,jss) = Huser(iuser,juser)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function legacy_reorder_mat_d

  function legacy_reorder_mat_c(Huser,Nin,OrderIn,OrderOut) result(Hss)
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
       do i=1,3
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser,1)
    call assert_shape(Huser,[Nlso,Nlso],"legacy_reorder_mat_c","Huser")
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
       print*,"Legacy_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iuser=1,Nlso
          Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iss   = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do juser=1,Nlso
             Ivec  = i2indices(juser,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             jss = indices2i(Jvec,Nout)
             !
             Hss(iss,jss) = Huser(iuser,juser)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function legacy_reorder_mat_c





  !HK
  function legacy_reorder_hk_d(Huser,Nin,OrderIn,OrderOut) result(Hss)
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
       do i=1,3
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !
    Nlso = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Nlso,Nlso,Nk],"legacy_reorder_hk_d","Huser")
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
       print*,"Legacy_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iuser=1,Nlso
          Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iss   = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do juser=1,Nlso
             Ivec  = i2indices(juser,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             jss = indices2i(Jvec,Nout)
             !
             Hss(iss,jss,:) = Huser(iuser,juser,:)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function legacy_reorder_hk_d

  function legacy_reorder_hk_c(Huser,Nin,OrderIn,OrderOut) result(Hss)
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
       do i=1,3
          OrderOut_(i) = trim(OrderOut(i))
       enddo
    endif
    !   
    Nlso = size(Huser,1)
    Nk   = size(Huser,3)
    call assert_shape(Huser,[Nlso,Nlso,Nk],"legacy_reorder_hk_d","Huser")
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
       print*,"Legacy_Reorder_vec ERROR: wrong entry in IndexOut at: ",tb_findloc(IndexOut,0)
       stop
    endif
    !
    !From IndexOut we can re-order the dimensions array to get the User dimensions array 
    Nout=indx_reorder(Nin,IndexOut)
    !
    if(any(IndexOut/=[1,2,3]))then
       do iuser=1,Nlso
          Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
          Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
          iss   = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
          do juser=1,Nlso
             Ivec  = i2indices(juser,Nin)
             Jvec  = indx_reorder(Ivec,IndexOut)
             jss = indices2i(Jvec,Nout)
             !
             Hss(iss,jss,:) = Huser(iuser,juser,:)
          enddo
       enddo
    else
       Hss = Huser
    endif
    return
  end function legacy_reorder_hk_c










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
    do i=1,size(Ain)
       Aout(i) = Ain(Index(i))!Aout(Index(i)) = Ain(i)
    enddo
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
    real(8),intent(in),dimension(:)                  :: val
    real(8),dimension(:,:),allocatable               :: tmp
    integer                                          :: n,ndim
    !
    ndim = size(val)
    !
    if (allocated(vec)) then
       if(size(vec,2)/=ndim)stop "add_to_A2 error: size(vec,2)!=ndim"
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
    real(8),dimension(:,:,:),allocatable,intent(inout) :: vec
    real(8),dimension(:,:),intent(in)                  :: val  
    real(8),dimension(:,:,:),allocatable               :: tmp
    integer                                            :: n,n2,n3
    !
    n2 = size(val,1)
    n3 = size(val,2)
    !
    if (allocated(vec)) then
       if(n2/=size(vec,2))stop "add_to_A3 error: size(vec,2)!=n2"
       if(n3/=size(vec,3))stop "add_to_A3 error: size(vec,3)!=n3"
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











