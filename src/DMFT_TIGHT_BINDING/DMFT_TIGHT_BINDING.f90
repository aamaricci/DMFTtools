module DMFT_TIGHT_BINDING
  USE SF_CONSTANTS, only: pi,pi2,xi,one,zero
  USE SF_IOTOOLS
  USE SF_LINALG, only: eigh,det,eye,zeros,eig
  USE SF_COLORS
  USE SF_TIMER, only:start_timer,stop_timer,eta
  USE SF_MISC, only: assert_shape,sort_array
  USE SF_OPTIMIZE,only: fmin_cgminimize
  USE DMFT_CTRL_VARS
  USE DMFT_GLOC
  USE DMFT_GFIO
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private


  interface TB_build_model
     module procedure :: build_hk_model_kgrid_d
     module procedure :: build_hk_model_kgrid_c
     module procedure :: build_hk_model_nkvec_d
     module procedure :: build_hk_model_nkvec_c
     !
     module procedure :: build_hkR_model_kgrid_d
     module procedure :: build_hkR_model_kgrid_c
     module procedure :: build_hkR_model_nkvec_d
     module procedure :: build_hkR_model_nkvec_c
     !
     module procedure :: build_hk_path_d
     module procedure :: build_hk_path_c
     !
     module procedure :: build_hkR_path_d
     module procedure :: build_hkR_path_c
     !
     module procedure :: build_Hij_Nrvec
     !
     module procedure :: build_hk_w90
  end interface TB_build_model


  interface TB_solve_model
     module procedure :: solve_Hk_along_BZpath
     module procedure :: solve_w90Hk_along_BZpath !this is redundant
     module procedure :: solve_HkR_along_BZpath
     !< obsolete
     module procedure :: read_Hr_w90_solve_Hk_along_BZpath
  end interface TB_solve_model


  interface TB_w90_setup
     module procedure :: setup_w90
  end interface TB_w90_setup

  interface TB_w90_delete
     module procedure :: delete_w90
  end interface TB_w90_delete

  interface TB_w90_FermiLevel
     module procedure :: FermiLevel_w90
  end interface TB_w90_FermiLevel

  interface TB_write_hk
     module procedure :: write_hk_w90_func
     module procedure :: write_hk_w90_array
     module procedure :: write_hk_w90_path
  end interface TB_write_hk


  interface TB_read_hk
     module procedure :: read_hk_w90_array
     module procedure :: read_hk_w90_path
  end interface TB_read_hk


  interface TB_write_Hloc
     module procedure :: write_Hloc_1
     module procedure :: write_Hloc_2
  end interface TB_write_Hloc


  interface TB_read_Hloc
     module procedure :: read_Hloc_1
     module procedure :: read_Hloc_2
  end interface TB_read_Hloc


  interface TB_build_kgrid
     module procedure :: build_kgrid
     module procedure :: build_kgrid_generic
     module procedure :: kgrid_from_path_grid
     module procedure :: kgrid_from_path_dim
  end interface TB_build_kgrid


  interface TB_build_Rgrid
     module procedure ::   build_Rgrid
  end interface TB_build_Rgrid


  interface TB_print_bk
     module procedure :: print_bk
  end interface TB_print_bk

  interface TB_print_ei
     module procedure :: print_ei
  end interface TB_print_ei



  interface TB_hr_to_hk
     module procedure :: hk_from_w90_hr
     module procedure :: hloct_from_w90_hr
#ifdef _MPI
     module procedure :: hk_from_w90_hr_mpi
     module procedure :: hkt_from_w90_hr_mpi
#endif
  end interface TB_hr_to_hk


  interface TB_dipole
     module procedure :: dipole_t2g_LDA
#ifdef _MPI
     module procedure :: dipole_t2g_LDA_mpi
#endif
  end interface TB_dipole


  abstract interface
     function w90_hk(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: w90_hk
     end function w90_hk

     function d_hk(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       real(8),dimension(N,N)    :: hk_model
     end function d_hk

     function c_hk(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function c_hk
  end interface
  procedure(w90_hk),pointer      :: TB_w90_model=>w90_hk_model


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


  logical,save                                  :: set_eivec=.false.
  logical,save                                  :: set_bkvec=.false.

  type(ctrl_list)                               :: dos_ctrl_list
  real(8),dimension(2),save                     :: dos_range=[-10d0,10d0]
  integer,save                                  :: dos_Lreal=2048
  real(8),save                                  :: dos_eps=0.01d0
  character(len=128)                            :: dos_file="dos_Hk"
  real(8),dimension(:),allocatable              :: dos_wtk
  complex(8),dimension(:,:,:,:,:,:),allocatable :: dos_Greal ![Nlat,Nspin,Nspin,Norb,Norb,dos_Lreal]

  type,public :: w90_structure
     character(len=:),allocatable               :: w90_file
     integer                                    :: Num_wann=0
     integer                                    :: Nrpts=0
     integer                                    :: N15=15
     integer                                    :: Qst=0
     integer                                    :: Rst=0
     integer,allocatable,dimension(:)           :: Ndegen
     integer                                    :: Nlat=0
     integer                                    :: Norb=0   
     integer                                    :: Nspin=0
     integer,allocatable,dimension(:,:)         :: Rvec
     real(8),allocatable,dimension(:,:)         :: Rgrid
     complex(8),allocatable,dimension(:,:,:)    :: Hij
     complex(8),allocatable,dimension(:,:)      :: Hloc
     real(8)                                    :: Efermi
     logical                                    :: iFermi=.false.
     logical                                    :: verbose=.false.       
     logical                                    :: status=.false.
  end type w90_structure



  type(w90_structure) :: TB_w90

  public :: TB_set_ei
  public :: TB_set_bk
  !
  public :: TB_get_ei
  public :: TB_get_bk
  !
  public :: TB_reset_ei
  public :: TB_reset_bk
  !
  public :: TB_build_ei
  public :: TB_build_bk
  !
  public :: TB_print_ei
  public :: TB_print_bk
  !
  public :: TB_reciprocal_basis
  !
  public :: TB_build_kgrid
  public :: TB_build_Rgrid
  !
  public :: TB_write_grid
  !
  public :: TB_build_model
  public :: TB_solve_model
  !
  public :: TB_get_FermiLevel
  !
  public :: TB_w90_setup
  public :: TB_w90_delete
  public :: TB_w90_FermiLevel
  public :: TB_w90_model
  !
  public :: TB_write_hk
  public :: TB_read_hk
  !
  public :: TB_write_hloc
  public :: TB_read_hloc
  !
  public :: TB_build_CoordGrid
  public :: TB_find_IndxCoord
  !
  public :: TB_dipole
  public :: TB_hr_to_hk

  public :: TB_set_dos_range
  public :: TB_set_dos_lreal
  public :: TB_set_dos_eps
  public :: TB_set_dos_file

contains


  !< Setup, reset and build the lattices basis set (direct and reciprocal)
  include "tight_binding_basis.f90"


  !< build the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) Hamiltonian matrix
  ! from the function user defined hk_model procedure.
  include "tight_binding_build_hk_model.f90"
  include "tight_binding_w90_helper.f90"


  !<  solve the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) along a given linear
  ! path in the Brillouin Zone. A GNUPLOT script to plot the bands together with their
  ! character is generated.
  include "tight_binding_solve_hk.f90"



  !< read/write the Hamiltonian matrix H(k) and its local part 
  include "tight_binding_io_Hk.f90"


  !< read/write the local part of the Hamiltonian to a file
  include "tight_binding_io_Hloc.f90"


  !< construct a grid of k-points with minimum information
  include "tight_binding_grid.f90"


  !OLD STUFF:
  !< read the real space hopping matrix from Wannier90 output and create H(k)
  include "w90hr/tight_binding_build_hk_from_w90hr.f90"
#ifdef _MPI
  include "w90hr/tight_binding_build_hk_from_w90hr_mpi.f90"
#endif
  !< read the real space lattice position and compute the space integrals
  !          <t2g| [x,y,z] |t2g> using atomic orbital as a subset.
  include "w90hr/dipole_w90hr.f90"
#ifdef _MPI
  include "w90hr/dipole_w90hr_mpi.f90"
#endif




  subroutine TB_get_FermiLevel(Hk,filling,Ef)
    complex(8),dimension(:,:,:)           :: Hk
    real(8)                               :: filling
    real(8)                               :: Ef
    complex(8),dimension(:,:),allocatable :: Uk
    real(8),dimension(:),allocatable      :: Ek
    real(8),dimension(:),allocatable      :: Ek_all
    integer,dimension(:),allocatable      :: Ek_indx
    integer                               :: stride,Nk,Nso,ik,indx
    Nk = size(Hk,3)
    Nso = size(Hk,1)
    call assert_shape(Hk,[Nso,Nso,Nk],"TB_FermiLevel","Hk")
    allocate(Uk(Nso,Nso),Ek(Nso),Ek_all(Nso*Nk),Ek_indx(Nk*Nso))
    stride = 0
    do ik = 1,Nk 
       Uk = Hk(:,:,ik)
       call eigh(Uk,Ek)
       Ek_all(stride+1:stride+Nso) = Ek
       stride = stride+Nso
    enddo
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





#ifdef _MPI
  function MPI_Get_size(comm) result(size)
    integer :: comm
    integer :: size,ierr
    call MPI_Comm_size(comm,size,ierr)
  end function MPI_Get_size

  function MPI_Get_rank(comm) result(rank)
    integer :: comm
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
  end function MPI_Get_rank

  function MPI_Get_master(comm) result(master)
    integer :: comm
    logical :: master
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    master=.false.
    if(rank==0)master=.true.
  end function MPI_Get_master
#endif


  subroutine herm_check(A)
    complex(8),intent(in) ::   A(:,:)
    integer               ::   row,col,i,j
    row=size(A,1)
    col=size(A,2)
    do i=1,col
       do j=1,row
          if(abs(A(i,j))-abs(A(j,i)) .gt. 1d-12)then
             write(*,'(1A)') "--> NON HERMITIAN MATRIX <--"
             write(*,'(2(1A7,I3))') " row: ",i," col: ",j
             write(*,'(1A)') "  A(i,j)"
             write(*,'(2F22.18)') real(A(i,j)),aimag(A(i,j))
             write(*,'(1A)') "  A(j,i)"
             write(*,'(2F22.18)') real(A(j,i)),aimag(A(j,i))
             stop
          endif
       enddo
    enddo
  end subroutine herm_check


  function slo2lso(Hslo,Nlat,Nspin,Norb) result(Hlso)
    implicit none
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hslo
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    integer                                               :: iorb,ispin,ilat
    integer                                               :: jorb,jspin,jlat
    integer                                               :: iI,jI,iO,jO
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      !
                      iI = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                      jI = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                      !
                      iO = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                      jO = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                      !
                      Hlso(iO,jO) = Hslo(iI,jI)
                      !
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function slo2lso


END MODULE DMFT_TIGHT_BINDING











