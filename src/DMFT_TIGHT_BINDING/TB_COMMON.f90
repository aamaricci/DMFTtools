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
  USE DMFT_GLOC
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



  !Some special points in the BZ:
  !we do everything in 3d.
  real(8),dimension(3),public,parameter         :: kpoint_gamma=[0,0,0]*pi
  real(8),dimension(3),public,parameter         :: kpoint_x1=[1,0,0]*pi
  real(8),dimension(3),public,parameter         :: kpoint_x2=[0,1,0]*pi
  real(8),dimension(3),public,parameter         :: kpoint_x3=[0,0,1]*pi
  real(8),dimension(3),public,parameter         :: kpoint_m1=[1,1,0]*pi
  real(8),dimension(3),public,parameter         :: kpoint_m2=[0,1,1]*pi
  real(8),dimension(3),public,parameter         :: kpoint_m3=[1,0,1]*pi
  real(8),dimension(3),public,parameter         :: kpoint_r=[1,1,1]*pi


  real(8),dimension(3),save                     :: ei_x=[1d0,0d0,0d0]
  real(8),dimension(3),save                     :: ei_y=[0d0,1d0,0d0]
  real(8),dimension(3),save                     :: ei_z=[0d0,0d0,1d0]

  real(8),dimension(3),save                     :: bk_x=[1d0,0d0,0d0]*pi2
  real(8),dimension(3),save                     :: bk_y=[0d0,1d0,0d0]*pi2
  real(8),dimension(3),save                     :: bk_z=[0d0,0d0,1d0]*pi2
  real(8),dimension(3),save                     :: BZ_origin=[0d0,0d0,0d0]

  logical,save                                  :: io_eivec=.false.
  logical,save                                  :: io_bkvec=.false.
  logical,save                                  :: set_eivec=.false.
  logical,save                                  :: set_bkvec=.false.

  type(ctrl_list)                               :: dos_ctrl_list
  real(8),dimension(2),save                     :: dos_range=[-10d0,10d0]
  integer,save                                  :: dos_Lreal=2048
  real(8),save                                  :: dos_eps=0.01d0
  character(len=128)                            :: dos_file="dos_Hk"
  real(8),dimension(:),allocatable              :: dos_wtk
  complex(8),dimension(:,:,:,:,:,:),allocatable :: dos_Greal ![Nlat,Nspin,Nspin,Norb,Norb,dos_Lreal]

  integer                                       :: mpi_ierr
  integer                                       :: mpi_rank
  integer                                       :: mpi_size
  logical                                       :: mpi_master


contains



  subroutine TB_get_FermiLevel(Hk,filling,Ef,Nspin,verbose)
    complex(8),dimension(:,:,:)                 :: Hk
    real(8)                                     :: filling
    real(8)                                     :: Ef,Dmin,Dmax
    integer,optional                            :: Nspin
    logical,optional                            :: verbose
    complex(8),dimension(:,:),allocatable       :: Uvec
    real(8),dimension(:),allocatable            :: Eval
    complex(8),dimension(:,:,:),allocatable     :: Rk,Rk_tmp
    real(8),dimension(:,:),allocatable          :: Ek,Ek_tmp
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Rho,diagRho
    integer                                     :: Nk,Nso,ik,info,Nspin_
    logical                                     :: verbose_
    !
    verbose_ = .true. ;if(present(verbose))verbose_=verbose
    Nspin_   = 1      ;if(present(Nspin))Nspin_=Nspin
    !
    Nso = size(Hk,1)
    Nk = size(Hk,3)
    call assert_shape(Hk,[Nso,Nso,Nk],"TB_FermiLevel","Hk")
    !
    !MPI setup:
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#ifdef _MPI 
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    endif
#endif
    !
    allocate(Uvec(Nso,Nso))
    allocate(Eval(Nso))
    allocate(Rk_tmp(Nso,Nso,Nk),Rk(Nso,Nso,Nk))
    allocate(Ek_tmp(Nso,Nk),Ek(Nso,Nk))
    Rk_tmp = zero
    Ek_tmp = 0d0
    do ik=1+mpi_rank,Nk,mpi_size
       Uvec = Hk(:,:,ik)
       call eigh(Uvec,Eval)
       Rk_tmp(:,:,ik) = Uvec
       Ek_tmp(:,ik)   = Eval
    enddo
#ifdef _MPI
    if(check_MPI())then
       Rk = zero
       Ek = 0d0
       call AllReduce_MPI(MPI_COMM_WORLD,Ek_tmp,Ek)
       call AllReduce_MPI(MPI_COMM_WORLD,Rk_tmp,Rk)
       deallocate(Ek_tmp,Rk_tmp)
    else
       Rk = Rk_tmp
       Ek = Ek_tmp
    endif
#else
    Rk = Rk_tmp
    Ek = Ek_tmp
#endif
    !
    !
    Dmin = minval(Ek)
    Dmax = maxval(Ek)
    Ef   = Dmin
    call fzero(get_dens,Ef,Dmax,info,rguess=Dmin+0.5d0*(Dmax-Dmin))
    if(info/=1)then
       if(mpi_master)write(*,*)"ERROR TB_get_Fermi: fzero returned info>1 ",info
       stop
    endif
    if(mpi_master)write(*,*)"w90 Fermi Level: ",Ef    
  contains
    function get_dens(ef) result(ndens)
      real(8),intent(in)            :: ef
      real(8)                       :: dens,dens_tmp
      real(8)                       :: ndens
      dens_tmp = 0d0
      do ik=1+mpi_rank,Nk,mpi_size
         diagRho = diag(fermi(Ek(:,ik)-ef,1000d0))
         Uvec    = Rk(:,:,ik)
         Rho     = (Uvec .x. diagRho) .x. (conjg(transpose(Uvec)))
         dens_tmp= dens_tmp + sum(diagonal(Rho))/Nk
      enddo
#ifdef _MPI
      if(check_MPI())then
         dens = 0d0
         call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens)
      else
         dens = dens_tmp
      endif
#else
      dens = dens_tmp
#endif
      dens  = (3-Nspin_)*dens
      ndens = dens-filling
      if(mpi_master.AND.verbose_)write(*,"(A9,3G18.9)")"Ef,N0   =",ef,dens,filling
    end function get_dens
  end subroutine TB_Get_FermiLevel


  !   subroutine TB_get_FermiLevel(Hk,filling,Ef,nspin,verbose)
  !     complex(8),dimension(:,:,:)           :: Hk
  !     real(8)                               :: filling
  !     real(8)                               :: Ef
  !     integer,optional                            :: Nspin
  !     logical,optional                            :: verbose
  !     complex(8),dimension(:,:),allocatable :: Uk
  !     real(8),dimension(:),allocatable      :: Ek
  !     real(8),dimension(:),allocatable      :: Ek_all,Ek_tmp
  !     integer,dimension(:),allocatable      :: Ek_indx
  !     integer                               :: stride,Nk,Nso,ik,indx
  !     Nk = size(Hk,3)
  !     Nso = size(Hk,1)
  !     call assert_shape(Hk,[Nso,Nso,Nk],"TB_FermiLevel","Hk")
  !     allocate(Uk(Nso,Nso),Ek(Nso),Ek_all(Nso*Nk),Ek_indx(Nk*Nso),Ek_tmp(Nk*Nso))
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
  !     Ek_tmp = 0d0
  !     stride = 0
  !     do ik=1+mpi_rank,Nk,mpi_size
  !        stride = (ik-1)*Nso
  !        Uk = Hk(:,:,ik)
  !        call eigh(Uk,Ek)
  !        Ek_tmp(stride+1:stride+Nso) = Ek
  !     enddo
  ! #ifdef _MPI
  !     if(check_MPI())then
  !        Ek_all = 0d0
  !        call AllReduce_MPI(MPI_COMM_WORLD,Ek_tmp,Ek_all)
  !     else
  !        Ek_all = Ek_tmp
  !     endif
  ! #else
  !     Ek_all = Ek_tmp
  ! #endif
  !     call sort_array(Ek_all,Ek_indx)
  !     if(filling==0d0)filling=Nso
  !     indx  = ceiling(filling*Nk/2d0)
  !     Ef    = Ek_all(indx)
  !   end subroutine TB_Get_FermiLevel




  !< Set properties for non=interacting DOS calculation 
  subroutine TB_set_dos_range(range)
    real(8),dimension(2) :: range
    dos_range = range
  end subroutine TB_set_dos_range

  subroutine TB_set_dos_lreal(lreal)
    integer :: lreal
    dos_lreal = lreal
  end subroutine TB_set_dos_lreal

  subroutine TB_set_dos_eps(eps)
    real(8) :: eps
    dos_eps = eps
  end subroutine TB_set_dos_eps

  subroutine TB_set_dos_file(file)
    character(len=*) :: file
    dos_file = reg(file)
  end subroutine TB_set_dos_file





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


END MODULE TB_COMMON











