module TB_COMMON
  USE SF_CONSTANTS, only: pi,pi2,xi,one,zero
  USE SF_IOTOOLS
  USE SF_LINALG, only: eigh,det,eye,zeros,eig,diag,operator(.x.)
  USE SF_COLORS
  USE SF_TIMER, only:start_timer,stop_timer,eta
  USE SF_MISC, only: assert_shape,sort_array
  USE SF_OPTIMIZE,only: fmin_cgminimize
  USE DMFT_CTRL_VARS
  USE DMFT_GLOC
  USE DMFT_GFIO
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


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




  subroutine TB_get_FermiLevel(Hk,filling,Ef)
    complex(8),dimension(:,:,:)           :: Hk
    real(8)                               :: filling
    real(8)                               :: Ef
    complex(8),dimension(:,:),allocatable :: Uk
    real(8),dimension(:),allocatable      :: Ek
    real(8),dimension(:),allocatable      :: Ek_all,Ek_tmp
    integer,dimension(:),allocatable      :: Ek_indx
    integer                               :: stride,Nk,Nso,ik,indx
    Nk = size(Hk,3)
    Nso = size(Hk,1)
    call assert_shape(Hk,[Nso,Nso,Nk],"TB_FermiLevel","Hk")
    allocate(Uk(Nso,Nso),Ek(Nso),Ek_all(Nso*Nk),Ek_indx(Nk*Nso),Ek_tmp(Nk*Nso))
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

    Ek_tmp = 0d0
    stride = 0
    do ik=1+mpi_rank,Nk,mpi_size
       stride = (ik-1)*Nso
       Uk = Hk(:,:,ik)
       call eigh(Uk,Ek)
       Ek_tmp(stride+1:stride+Nso) = Ek
    enddo
#ifdef _MPI
    if(check_MPI())then
       Ek_all = 0d0
       call AllReduce_MPI(MPI_COMM_WORLD,Ek_tmp,Ek_all)
    else
       Ek_all = Ek_tmp
    endif
#else
    Ek_all = Ek_tmp
#endif
    call sort_array(Ek_all,Ek_indx)
    if(filling==0d0)filling=Nso
    indx  = ceiling(filling*Nk/2d0)
    Ef    = Ek_all(indx)
  end subroutine TB_Get_FermiLevel




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





END MODULE TB_COMMON











