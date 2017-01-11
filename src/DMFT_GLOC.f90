module DMFT_GLOC
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,txtfy,splot,store_data
  USE SF_LINALG,    only: eye,inv,inv_sym,inv_tridiag,get_tridiag
  USE SF_ARRAYS,    only: linspace,arange
  USE SF_MISC,      only: assert_shape
  USE DMFT_CTRL_VARS

#ifdef _MPI
  USE MPI
#endif
  implicit none
  private

  interface dmft_gloc_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_main
     module procedure :: dmft_get_gloc_matsubara_normal_1band
     module procedure :: dmft_get_gloc_matsubara_normal_Nband
     module procedure :: dmft_get_gloc_matsubara_normal_lattice_main
     module procedure :: dmft_get_gloc_matsubara_normal_lattice_1band
     module procedure :: dmft_get_gloc_matsubara_normal_lattice_Nband
#ifdef _MPI
     module procedure :: dmft_get_gloc_matsubara_normal_main_mpi
     module procedure :: dmft_get_gloc_matsubara_normal_1band_mpi
     module procedure :: dmft_get_gloc_matsubara_normal_Nband_mpi
     module procedure :: dmft_get_gloc_matsubara_normal_lattice_main_mpi
     module procedure :: dmft_get_gloc_matsubara_normal_lattice_1band_mpi     
     module procedure :: dmft_get_gloc_matsubara_normal_lattice_Nband_mpi
#endif
  end interface dmft_gloc_matsubara




  interface dmft_gloc_matsubara_superc
     module procedure :: dmft_get_gloc_matsubara_superc_main
     module procedure :: dmft_get_gloc_matsubara_superc_1band
     module procedure :: dmft_get_gloc_matsubara_superc_Nband
     module procedure :: dmft_get_gloc_matsubara_superc_lattice_main          
     module procedure :: dmft_get_gloc_matsubara_superc_lattice_1band
     module procedure :: dmft_get_gloc_matsubara_superc_lattice_Nband     
#ifdef _MPI
     module procedure :: dmft_get_gloc_matsubara_superc_main_mpi
     module procedure :: dmft_get_gloc_matsubara_superc_1band_mpi
     module procedure :: dmft_get_gloc_matsubara_superc_Nband_mpi     
     module procedure :: dmft_get_gloc_matsubara_superc_lattice_main_mpi     
     module procedure :: dmft_get_gloc_matsubara_superc_lattice_1band_mpi
     module procedure :: dmft_get_gloc_matsubara_superc_lattice_Nband_mpi
#endif
  end interface dmft_gloc_matsubara_superc






  interface dmft_gloc_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_main
     module procedure :: dmft_get_gloc_realaxis_normal_1band
     module procedure :: dmft_get_gloc_realaxis_normal_Nband
     module procedure :: dmft_get_gloc_realaxis_normal_lattice_main
     module procedure :: dmft_get_gloc_realaxis_normal_lattice_1band     
     module procedure :: dmft_get_gloc_realaxis_normal_lattice_Nband
#ifdef _MPI
     module procedure :: dmft_get_gloc_realaxis_normal_main_mpi
     module procedure :: dmft_get_gloc_realaxis_normal_1band_mpi
     module procedure :: dmft_get_gloc_realaxis_normal_Nband_mpi
     module procedure :: dmft_get_gloc_realaxis_normal_lattice_main_mpi
     module procedure :: dmft_get_gloc_realaxis_normal_lattice_1band_mpi
     module procedure :: dmft_get_gloc_realaxis_normal_lattice_Nband_mpi
#endif
  end interface dmft_gloc_realaxis





  interface dmft_gloc_realaxis_superc
     module procedure :: dmft_get_gloc_realaxis_superc_main
     module procedure :: dmft_get_gloc_realaxis_superc_1band
     module procedure :: dmft_get_gloc_realaxis_superc_Nband
     module procedure :: dmft_get_gloc_realaxis_superc_lattice_main
     module procedure :: dmft_get_gloc_realaxis_superc_lattice_1band     
     module procedure :: dmft_get_gloc_realaxis_superc_lattice_Nband
#ifdef _MPI
     module procedure :: dmft_get_gloc_realaxis_superc_main_mpi
     module procedure :: dmft_get_gloc_realaxis_superc_1band_mpi
     module procedure :: dmft_get_gloc_realaxis_superc_Nband_mpi
     module procedure :: dmft_get_gloc_realaxis_superc_lattice_main_mpi          
     module procedure :: dmft_get_gloc_realaxis_superc_lattice_1band_mpi
     module procedure :: dmft_get_gloc_realaxis_superc_lattice_Nband_mpi
#endif
  end interface dmft_gloc_realaxis_superc









  interface dmft_gij_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_gij_main
#ifdef _MPI
     module procedure :: dmft_get_gloc_matsubara_normal_gij_main_mpi
#endif     
  end interface dmft_gij_matsubara






  interface dmft_gij_matsubara_superc
     module procedure :: dmft_get_gloc_matsubara_superc_gij_main
#ifdef _MPI
     module procedure :: dmft_get_gloc_matsubara_superc_gij_main_mpi
#endif     
  end interface dmft_gij_matsubara_superc






  interface dmft_gij_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_gij_main
#ifdef _MPI
     module procedure :: dmft_get_gloc_realaxis_normal_gij_main_mpi
#endif     
  end interface dmft_gij_realaxis




  interface dmft_gij_realaxis_superc
     module procedure :: dmft_get_gloc_realaxis_superc_gij_main
#ifdef _MPI
     module procedure :: dmft_get_gloc_realaxis_superc_gij_main_mpi
#endif     
  end interface dmft_gij_realaxis_superc




  interface dmft_set_Gamma_matsubara
     module procedure :: dmft_set_Gamma_matsubara_local
     module procedure :: dmft_set_Gamma_matsubara_lattice
  end interface dmft_set_Gamma_matsubara



  interface dmft_set_Gamma_realaxis
     module procedure :: dmft_set_Gamma_realaxis_local
     module procedure :: dmft_set_Gamma_realaxis_lattice
  end interface dmft_set_Gamma_realaxis




  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block
  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape
  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape
  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape
  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape


  real(8),dimension(:),allocatable          :: wm !Matsubara frequencies
  real(8),dimension(:),allocatable          :: wr !Real frequencies
  character(len=128)                        :: suffix
  character(len=128)                        :: gloc_suffix=".dat"
  complex(8),dimension(:,:,:),allocatable   :: lattice_Gamma_mats ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:,:),allocatable :: local_Gamma_mats   ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable   :: lattice_Gamma_real ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:,:),allocatable :: local_Gamma_real   ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  integer                                   :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats
  integer                                   :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
  !
  integer                                   :: mpi_ierr
  integer                                   :: mpi_rank
  integer                                   :: mpi_size
  logical                                   :: mpi_master
  !
  real(8)                                   :: beta
  real(8)                                   :: xmu,eps
  real(8)                                   :: wini,wfin 


  !PUBLIC IN DMFT:
  public :: dmft_gloc_matsubara
  public :: dmft_gloc_matsubara_superc

  public :: dmft_gloc_realaxis
  public :: dmft_gloc_realaxis_superc

  public :: dmft_gij_matsubara
  public :: dmft_gij_matsubara_superc

  public :: dmft_gij_realaxis
  public :: dmft_gij_realaxis_superc

  public :: dmft_set_Gamma_matsubara
  public :: dmft_set_Gamma_realaxis



  !####################################################################
  !####################################################################

contains


  !####################################################################
  !####################################################################



  !--------------------------------------------------------------------!
  ! PURPOSE: evaluate the normal/superconducting local Green's function:
  ! input:
  ! 1.  H(k): [(Nlat*)Nspin*Norb][(Nlat*)Nspin*Norb][Lk]
  ! 2. DMFT  Sigma & Self: ([2][Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! Additional interfaces with different shapes of the Sigma are
  ! appended to each file:
  !  1 BAND
  ! - Sigma shape: ([Nlat])[L] (no Nspin[=1], no Norb[=1])
  ! NN FORM
  ! - Sigma shape: ([Nlat])[Norb][Norb][L] (Nspin=1)
  !--------------------------------------------------------------------!
#define FROUTINE 'dmft_get_gloc_matsubara_normal_main'
  subroutine dmft_get_gloc_matsubara_normal_main(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                            :: iprint
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
    !
    real(8)                                       :: beta
    real(8)                                       :: xmu,eps
    real(8)                                       :: wini,wfin
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gkmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    call start_timer
    Gmats=zero
    do ik=1,Lk
       call invert_gk_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
       Gmats = Gmats + Gkmats*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_matsubara(wm,Gmats,"Gloc",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_main'
  subroutine dmft_get_gloc_matsubara_normal_lattice_main(Hk,Wtk,Gmats,Smats,iprint,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                              :: iprint    !
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(.not.tridiag_)then
       write(*,"(A)")"Direct Inversion algorithm:"
    else
       write(*,"(A)")"Quantum Zipper algorithm:"
    endif
    !
    allocate(Gkmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    call start_timer
    Gmats=zero
    if(.not.tridiag_)then
       do ik=1,Lk
          call invert_gk_normal_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats*Wtk(ik)
          call eta(ik,Lk)
       end do
    else
       do ik=1,Lk
          call invert_gk_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats*Wtk(ik)
          call eta(ik,Lk)
       end do
    endif
    call stop_timer
    call dmft_gloc_print_matsubara_lattice(wm,Gmats,"LG",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_lattice_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_gij_main'
  subroutine dmft_get_gloc_matsubara_normal_gij_main(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats     !      [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer,intent(in)                                :: iprint
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats));zeta_mats=zero
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    call start_timer
    Gmats=zero
    do ik=1,Lk
       call invert_gk_normal_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       Gmats = Gmats + Gkmats*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_matsubara_gij(wm,Gmats,"Gij",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_gij_main
#undef FROUTINE


#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_matsubara_normal_main_mpi'
  subroutine dmft_get_gloc_matsubara_normal_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                       :: MpiComm
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                            :: iprint
    character(len=*),optional                     :: mpi_split  
    character(len=1)                              :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkmats    !as Smats  
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Smats
    complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gkmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)call start_timer
    Gmats=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_matsubara_normal_main_mpi: ! mpi_split_ in ['w','k'] "
    case ('w')
       do ik=1,Lk
          call invert_gk_normal_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
          Gmats = Gmats + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
    case ('k')
       allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_matsubara(wm,Gmats,"Gloc",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_main_mpi'
  subroutine dmft_get_gloc_matsubara_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                              :: iprint    !
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    character(len=*),optional                       :: mpi_split  
    character(len=1)                                :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gkmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gmats=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_matsubara_normal_main_lattice_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_lattice_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gmats = Gmats + Gkmats*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_tridiag_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gmats = Gmats + Gkmats*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
       !
    case ('k')
       allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gtmp = Gtmp + Gkmats*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gtmp = Gtmp + Gkmats*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
       call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_matsubara_lattice(wm,Gmats,"LG",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_lattice_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_gij_main_mpi'
  subroutine dmft_get_gloc_matsubara_normal_gij_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats     !      [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer,intent(in)                                :: iprint
    character(len=*),optional                         :: mpi_split  
    character(len=1)                                  :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gmats=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_matsubara_normal_gij_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_normal_gij_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero     
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp, Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_matsubara_gij(wm,Gmats,"Gij",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_gij_main_mpi
#undef FROUTINE
#endif





  !##################################################################
  !##################################################################
  !##################################################################



#define FROUTINE 'dmft_get_gloc_matsubara_normal_1band'
  subroutine dmft_get_gloc_matsubara_normal_1band(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
    real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                     :: Smats(:)
    complex(8),intent(inout)                  :: Gmats(size(Smats))
    logical,optional                          :: hk_symm(size(Hk,1))
    logical                                   :: hk_symm_(size(Hk,1))
    integer                                   :: iprint
    !
    complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(1,1,1,1,size(Smats)) :: Gmats_
    complex(8),dimension(1,1,1,1,size(Smats)) :: Smats_
    !
    Hk_(1,1,:)        = Hk(:)
    Smats_(1,1,1,1,:) = Smats(:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_normal_main(Hk_,Wtk,Gmats_,Smats_,iprint,hk_symm_)
    Gmats(:) = Gmats_(1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_normal_1band
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_1band'
  subroutine dmft_get_gloc_matsubara_normal_lattice_1band(Hk,Wtk,Gmats,Smats,iprint,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),dimension(:,:),intent(in)                            :: Smats           ![Nlat][Lmats]
    complex(8),dimension(size(Smats,1),size(Smats,2)),intent(inout) :: Gmats
    integer                                                         :: iprint
    logical,optional                                                :: tridiag
    logical                                                         :: tridiag_
    logical,optional                                                :: hk_symm(size(Hk,1))
    logical                                                         :: hk_symm_(size(Hk,1))
    !
    complex(8),dimension(size(Smats,1),1,1,1,1,size(Smats,2))       :: Gmats_
    complex(8),dimension(size(Smats,1),1,1,1,1,size(Smats,2))       :: Smats_
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[size(Hk,1),size(Smats,2)],FROUTINE,"Smats")
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_gloc_matsubara_normal_lattice_main(Hk,Wtk,Gmats_,Smats_,iprint,tridiag_,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_normal_lattice_1band
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_matsubara_normal_1band_mpi'
  subroutine dmft_get_gloc_matsubara_normal_1band_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                   :: MpiComm
    complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
    real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                     :: Smats(:)
    complex(8),intent(inout)                  :: Gmats(size(Smats))
    logical,optional                          :: hk_symm(size(Hk,1))
    logical                                   :: hk_symm_(size(Hk,1))
    integer                                   :: iprint
    character(len=*),optional                 :: mpi_split  
    character(len=1)                          :: mpi_split_ 
    !
    complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(1,1,1,1,size(Smats)) :: Gmats_
    complex(8),dimension(1,1,1,1,size(Smats)) :: Smats_
    !
    Hk_(1,1,:)        = Hk(:)
    Smats_(1,1,1,1,:) = Smats(:)
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_normal_main_mpi(MpiComm,Hk_,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
    Gmats(:) = Gmats_(1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_normal_1band_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_1band_mpi'
  subroutine dmft_get_gloc_matsubara_normal_lattice_1band_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),dimension(:,:),intent(in)                            :: Smats           ![Nlat][Lmats]
    complex(8),dimension(size(Smats,1),size(Smats,2)),intent(inout) :: Gmats
    integer                                                         :: iprint
    logical,optional                                                :: tridiag
    logical                                                         :: tridiag_
    character(len=*),optional                                       :: mpi_split  
    character(len=1)                                                :: mpi_split_ 
    logical,optional                                                :: hk_symm(size(Hk,1))
    logical                                                         :: hk_symm_(size(Hk,1))
    !
    complex(8),dimension(size(Smats,1),1,1,1,1,size(Smats,2))       :: Gmats_
    complex(8),dimension(size(Smats,1),1,1,1,1,size(Smats,2))       :: Smats_
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[size(Hk,1),size(Smats,2)],FROUTINE,"Smats")
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_gloc_matsubara_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,tridiag_,mpi_split_,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_normal_lattice_1band_mpi
#undef FROUTINE
#endif






  !##################################################################
  !##################################################################
  !##################################################################







#define FROUTINE 'dmft_get_gloc_matsubara_normal_Nband'
  subroutine dmft_get_gloc_matsubara_normal_Nband(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                        :: Hk              ![Norb][Norb][Nk]
    real(8)                                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in),dimension(:,:,:)                        :: Smats(:,:,:)    ![Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:)                     :: Gmats
    logical,optional                                              :: hk_symm(size(Hk,3))
    logical                                                       :: hk_symm_(size(Hk,3))
    integer                                                       :: iprint
    !
    complex(8),&
         dimension(1,1,size(Smats,1),size(Smats,2),size(Smats,3)) :: Gmats_
    complex(8),&
         dimension(1,1,size(Smats,1),size(Smats,2),size(Smats,3)) :: Smats_
    !
    integer                                                       :: Nspin,Norb,Nso,Lmats
    !
    Nspin = 1
    Norb  = size(Smats,1)
    Lmats = size(Smats,3)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[Nso,Nso,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nso,Nso,Lmats],FROUTINE,"Gmats")
    !
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_normal_main(Hk,Wtk,Gmats_,Smats_,iprint,hk_symm_)
    Gmats(:,:,:) = Gmats_(1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_normal_Nband
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_Nband'
  subroutine dmft_get_gloc_matsubara_normal_lattice_Nband(Hk,Wtk,Gmats,Smats,iprint,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                                      :: Hk              ![Nlat*Norb][Nlat*Norb][Nk]
    real(8)                                                                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                                       :: Smats(:,:,:,:)  ![Nlat][Norb][Norb][Lmats]
    complex(8),intent(inout)                                                    :: Gmats(:,:,:,:)
    logical,optional                                                            :: hk_symm(size(Hk,3))
    logical                                                                     :: hk_symm_(size(Hk,3))
    integer                                                                     :: iprint
    logical,optional                                                            :: tridiag
    logical                                                                     :: tridiag_
    !
    complex(8),&
         dimension(size(Smats,1),1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Gmats_ ![Nlat][1][1][Norb][Norb][Lmats]
    complex(8),&
         dimension(size(Smats,1),1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Smats_![Nlat][1][1][Norb][Norb][Lmats]
    integer                                                                     :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1    
    Nlat  = size(Smats,1)
    Norb  = size(Smats,2)
    Lmats = size(Smats,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[Nlat,Nso,Nso,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nso,Nso,Lmats],FROUTINE,"Gmats")
    !
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_gloc_matsubara_normal_lattice_main(Hk,Wtk,Gmats_,Smats_,iprint,tridiag_,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_normal_lattice_Nband
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_matsubara_normal_Nband_mpi'
  subroutine dmft_get_gloc_matsubara_normal_Nband_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                                             :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                              :: Hk              ![Norb][Norb][Nk]
    real(8)                                                             :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in),dimension(:,:,:)                              :: Smats(:,:,:)    ![Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:)                           :: Gmats
    logical,optional                                                    :: hk_symm(size(Hk,3))
    logical                                                             :: hk_symm_(size(Hk,3))
    integer                                                             :: iprint
    character(len=*),optional                                           :: mpi_split  
    character(len=1)                                                    :: mpi_split_ 
    !
    complex(8),dimension(1,1,size(Smats,1),size(Smats,2),size(Smats,3)) :: Gmats_
    complex(8),dimension(1,1,size(Smats,1),size(Smats,2),size(Smats,3)) :: Smats_
    integer                                                             :: Nspin,Norb,Nso,Lmats
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Norb  = size(Smats,1)
    Lmats = size(Smats,3)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[Nso,Nso,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nso,Nso,Lmats],FROUTINE,"Gmats")
    !
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    call dmft_get_gloc_matsubara_normal_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
    Gmats(:,:,:) = Gmats_(1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_normal_Nband_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_Nband_mpi'
  subroutine dmft_get_gloc_matsubara_normal_lattice_Nband_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
    integer                                                                     :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                                      :: Hk              ![Nlat*Norb][Nlat*Norb][Nk]
    real(8)                                                                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                                       :: Smats(:,:,:,:)  ![Nlat][Norb][Norb][Lmats]
    complex(8),intent(inout)                                                    :: Gmats(:,:,:,:)
    logical,optional                                                            :: hk_symm(size(Hk,3))
    logical                                                                     :: hk_symm_(size(Hk,3))
    integer                                                                     :: iprint
    logical,optional                                                            :: tridiag
    logical                                                                     :: tridiag_
    character(len=*),optional                                                   :: mpi_split  
    character(len=1)                                                            :: mpi_split_ 
    !
    complex(8),&
         dimension(size(Smats,1),1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Gmats_ ![Nlat][1][1][Norb][Norb][Lmats]
    complex(8),&
         dimension(size(Smats,1),1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Smats_![Nlat][1][1][Norb][Norb][Lmats]
    !
    integer                                                                     :: Nspin,Norb,Nso,Nlso,Lmats
    !
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1    
    Nlat  = size(Smats,1)
    Norb  = size(Smats,2)
    Lmats = size(Smats,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[Nlat,Nso,Nso,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nso,Nso,Lmats],FROUTINE,"Gmats")
    !
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_gloc_matsubara_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,tridiag_,mpi_split_,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_normal_lattice_Nband_mpi
#undef FROUTINE
#endif










#define FROUTINE 'dmft_get_gloc_realaxis_normal_main'
  subroutine dmft_get_gloc_realaxis_normal_main(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                            :: iprint
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                       :: beta
    real(8)                                       :: xmu,eps
    real(8)                                       :: wini,wfin
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(Gkreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    call start_timer
    Greal=zero
    do ik=1,Lk
       call invert_gk_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
       Greal = Greal + Gkreal*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_realaxis(wr,Greal,"Gloc",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_main'
  subroutine dmft_get_gloc_realaxis_normal_lattice_main(Hk,Wtk,Greal,Sreal,iprint,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint    !
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(.not.tridiag_)then
       write(*,"(A)")"Direct Inversion algorithm:"
    else
       write(*,"(A)")"Quantum Zipper algorithm:"
    endif
    !
    allocate(Gkreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    call start_timer
    Greal=zero
    if(.not.tridiag_)then
       do ik=1,Lk
          call invert_gk_normal_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal*Wtk(ik)
          call eta(ik,Lk)
       end do
    else
       do ik=1,Lk
          call invert_gk_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal*Wtk(ik)
          call eta(ik,Lk)
       end do
    endif
    call stop_timer
    call dmft_gloc_print_realaxis_lattice(wr,Greal,"LG",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_lattice_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_gij_main'
  subroutine dmft_get_gloc_realaxis_normal_gij_main(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer,intent(in)                                :: iprint
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(zeta_real(Nlat,Nso,Nso,Lreal));zeta_real=zero
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    call start_timer
    Greal=zero
    do ik=1,Lk
       call invert_gk_normal_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Greal = Greal + Gkreal*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_realaxis_gij(wm,Greal,"Gij",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_gij_main
#undef FROUTINE


#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_realaxis_normal_main_mpi'
  subroutine dmft_get_gloc_realaxis_normal_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                       :: MpiComm
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                            :: iprint
    character(len=*),optional                     :: mpi_split  
    character(len=1)                              :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkreal    !as Sreal  
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal] 
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(Gkreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)call start_timer
    Greal=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_realaxis_normal_main_mpi: ! mpi_split_ in ['w','k'] "
    case ('w')
       do ik=1,Lk
          call invert_gk_normal_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
          Greal = Greal + Gkreal*Wtk(ik)
          call eta(ik,Lk)
       end do
    case ('k')
       allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_realaxis(wr,Greal,"Gloc",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_main_mpi'
  subroutine dmft_get_gloc_realaxis_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,tridiag,mpi_split,hk_symm)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint    !
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    character(len=*),optional                       :: mpi_split  
    character(len=1)                                :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin  
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gkreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Greal=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_realaxis_normal_lattice_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_lattice_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Greal = Greal + Gkreal*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_tridiag_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Greal = Greal + Gkreal*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
       !
    case ('k')
       allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Gtmp = Gtmp + Gkreal*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Gtmp = Gtmp + Gkreal*Wtk(ik)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
       call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_realaxis_lattice(wr,Greal,"LG",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_lattice_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_gij_main_mpi'
  subroutine dmft_get_gloc_realaxis_normal_gij_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer,intent(in)                                :: iprint
    character(len=*),optional                         :: mpi_split  
    character(len=1)                                  :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Greal=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_realaxis_normal_gij_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_normal_gij_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero     
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp, Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_realaxis_gij(wm,Greal,"Gij",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_gij_main_mpi
#undef FROUTINE
#endif





  !##################################################################
  !##################################################################
  !##################################################################


#define FROUTINE 'dmft_get_gloc_realaxis_normal_1band'
  subroutine dmft_get_gloc_realaxis_normal_1band(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
    real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                     :: Sreal(:)
    complex(8),intent(inout)                  :: Greal(size(Sreal))
    logical,optional                          :: hk_symm(size(Hk,1))
    logical                                   :: hk_symm_(size(Hk,1))
    integer                                   :: iprint
    !
    complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(1,1,1,1,size(Sreal)) :: Greal_
    complex(8),dimension(1,1,1,1,size(Sreal)) :: Sreal_
    !
    Hk_(1,1,:)        = Hk(:)
    Sreal_(1,1,1,1,:) = Sreal(:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_normal_main(Hk_,Wtk,Greal_,Sreal_,iprint,hk_symm_)
    Greal(:) = Greal_(1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_normal_1band
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_1band'
  subroutine dmft_get_gloc_realaxis_normal_lattice_1band(Hk,Wtk,Greal,Sreal,iprint,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),dimension(:,:),intent(in)                            :: Sreal           ![Nlat][Lreal]
    complex(8),dimension(size(Sreal,1),size(Sreal,2)),intent(inout) :: Greal
    integer                                                         :: iprint
    logical,optional                                                :: tridiag
    logical                                                         :: tridiag_
    logical,optional                                                :: hk_symm(size(Hk,1))
    logical                                                         :: hk_symm_(size(Hk,1))
    !
    complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Greal_
    complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Sreal_
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[size(Hk,1),size(Sreal,2)],FROUTINE,"Sreal")
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    call dmft_get_gloc_realaxis_normal_lattice_main(Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,hk_symm_)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_normal_lattice_1band
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_realaxis_normal_1band_mpi'
  subroutine dmft_get_gloc_realaxis_normal_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                   :: MpiComm
    complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
    real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                     :: Sreal(:)
    complex(8),intent(inout)                  :: Greal(size(Sreal))
    logical,optional                          :: hk_symm(size(Hk,1))
    logical                                   :: hk_symm_(size(Hk,1))
    integer                                   :: iprint
    character(len=*),optional                 :: mpi_split  
    character(len=1)                          :: mpi_split_ 
    !
    complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(1,1,1,1,size(Sreal)) :: Greal_
    complex(8),dimension(1,1,1,1,size(Sreal)) :: Sreal_
    !
    Hk_(1,1,:)        = Hk(:)
    Sreal_(1,1,1,1,:) = Sreal(:)
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_normal_main_mpi(MpiComm,Hk_,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
    Greal(:) = Greal_(1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_normal_1band_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_1band_mpi'
  subroutine dmft_get_gloc_realaxis_normal_lattice_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,tridiag,mpi_split,hk_symm)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),dimension(:,:),intent(in)                            :: Sreal           ![Nlat][Lreal]
    complex(8),dimension(size(Sreal,1),size(Sreal,2)),intent(inout) :: Greal
    integer                                                         :: iprint
    logical,optional                                                :: tridiag
    logical                                                         :: tridiag_
    character(len=*),optional                                       :: mpi_split  
    character(len=1)                                                :: mpi_split_ 
    logical,optional                                                :: hk_symm(size(Hk,1))
    logical                                                         :: hk_symm_(size(Hk,1))
    !
    complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Greal_
    complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Sreal_
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[size(Hk,1),size(Sreal,2)],FROUTINE,"Sreal")
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    call dmft_get_gloc_realaxis_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,mpi_split_,hk_symm_)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_normal_lattice_1band_mpi
#undef FROUTINE
#endif




  !##################################################################
  !##################################################################
  !##################################################################



#define FROUTINE 'dmft_get_gloc_realaxis_normal_Nband'
  subroutine dmft_get_gloc_realaxis_normal_Nband(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                        :: Hk              ![Norb][Norb][Nk]
    real(8)                                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in),dimension(:,:,:)                        :: Sreal(:,:,:)    ![Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:)                     :: Greal
    logical,optional                                              :: hk_symm(size(Hk,3))
    logical                                                       :: hk_symm_(size(Hk,3))
    integer                                                       :: iprint
    !
    complex(8),&
         dimension(1,1,size(Sreal,1),size(Sreal,2),size(Sreal,3)) :: Greal_
    complex(8),&
         dimension(1,1,size(Sreal,1),size(Sreal,2),size(Sreal,3)) :: Sreal_
    !
    integer                                                       :: Nspin,Norb,Nso,Lreal
    !
    Nspin = 1
    Norb  = size(Sreal,1)
    Lreal = size(Sreal,3)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nso,Nso,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nso,Nso,Lreal],FROUTINE,"Greal")
    !
    Sreal_(1,1,:,:,:) = Sreal(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_normal_main(Hk,Wtk,Greal_,Sreal_,iprint,hk_symm_)
    Greal(:,:,:) = Greal_(1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_normal_Nband
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_Nband'
  subroutine dmft_get_gloc_realaxis_normal_lattice_Nband(Hk,Wtk,Greal,Sreal,iprint,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                                      :: Hk              ![Nlat*Norb][Nlat*Norb][Nk]
    real(8)                                                                     :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                                       :: Sreal(:,:,:,:)  ![Nlat][Norb][Norb][Lreal]
    complex(8),intent(inout)                                                    :: Greal(:,:,:,:)
    logical,optional                                                            :: hk_symm(size(Hk,3))
    logical                                                                     :: hk_symm_(size(Hk,3))
    integer                                                                     :: iprint
    logical,optional                                                            :: tridiag
    logical                                                                     :: tridiag_
    !
    complex(8),&
         dimension(size(Sreal,1),1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Greal_ ![Nlat][1][1][Norb][Norb][Lreal]
    complex(8),&
         dimension(size(Sreal,1),1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Sreal_![Nlat][1][1][Norb][Norb][Lreal]
    integer                                                                     :: Nlat,Nspin,Norb,Nso,Nlso,Lreal
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1    
    Nlat  = size(Sreal,1)
    Norb  = size(Sreal,2)
    Lreal = size(Sreal,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nlat,Nso,Nso,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nso,Nso,Lreal],FROUTINE,"Greal")
    !
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    call dmft_get_gloc_realaxis_normal_lattice_main(Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,hk_symm_)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_normal_lattice_Nband
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_realaxis_normal_Nband_mpi'
  subroutine dmft_get_gloc_realaxis_normal_Nband_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                                          :: Hk              ![Norb][Norb][Nk]
    real(8)                                                                         :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in),dimension(:,:,:)                                          :: Sreal(:,:,:)    ![Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:)                                       :: Greal
    logical,optional                                                                :: hk_symm(size(Hk,3))
    logical                                                                         :: hk_symm_(size(Hk,3))
    integer                                                                         :: iprint
    character(len=*),optional                                                       :: mpi_split  
    character(len=1)                                                                :: mpi_split_ 
    !
    complex(8),dimension(1,1,size(Sreal,1),size(Sreal,2),size(Sreal,3)) :: Greal_
    complex(8),dimension(1,1,size(Sreal,1),size(Sreal,2),size(Sreal,3)) :: Sreal_
    integer                                                                         :: Nspin,Norb,Nso,Lreal
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Norb  = size(Sreal,1)
    Lreal = size(Sreal,3)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nso,Nso,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nso,Nso,Lreal],FROUTINE,"Greal")
    !
    Sreal_(1,1,:,:,:) = Sreal(:,:,:)
    call dmft_get_gloc_realaxis_normal_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
    Greal(:,:,:) = Greal_(1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_normal_Nband_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_Nband_mpi'
  subroutine dmft_get_gloc_realaxis_normal_lattice_Nband_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,tridiag,mpi_split,hk_symm)
    integer                                                                           :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                                            :: Hk              ![Nlat*Norb][Nlat*Norb][Nk]
    real(8)                                                                           :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                                             :: Sreal(:,:,:,:)  ![Nlat][Norb][Norb][Lreal]
    complex(8),intent(inout)                                                          :: Greal(:,:,:,:)
    logical,optional                                                                  :: hk_symm(size(Hk,3))
    logical                                                                           :: hk_symm_(size(Hk,3))
    integer                                                                           :: iprint
    logical,optional                                                                  :: tridiag
    logical                                                                           :: tridiag_
    character(len=*),optional                                                         :: mpi_split  
    character(len=1)                                                                  :: mpi_split_ 
    !
    complex(8),dimension(size(Sreal,1),1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Greal_ ![Nlat][1][1][Norb][Norb][Lreal]
    complex(8),dimension(size(Sreal,1),1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Sreal_![Nlat][1][1][Norb][Norb][Lreal]
    !
    integer                                                                           :: Nspin,Norb,Nso,Nlso,Lreal
    !
    !
    tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1    
    Nlat  = size(Sreal,1)
    Norb  = size(Sreal,2)
    Lreal = size(Sreal,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[Nlat,Nso,Nso,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nso,Nso,Lreal],FROUTINE,"Greal")
    !
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    call dmft_get_gloc_realaxis_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,mpi_split_,hk_symm_)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_normal_lattice_Nband_mpi
#undef FROUTINE
#endif









#define FROUTINE 'dmft_get_gloc_matsubara_superc_main'
  subroutine dmft_get_gloc_matsubara_superc_main(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                              :: iprint
    logical,dimension(size(Hk,3)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    allocate(Gkmats(2,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(2,2,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    call start_timer
    Gmats=zero
    do ik=1,Lk
       call invert_gk_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       Gmats = Gmats + Gkmats*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_matsubara(wm,Gmats(1,:,:,:,:,:),"Gloc",iprint)
    call dmft_gloc_print_matsubara(wm,Gmats(2,:,:,:,:,:),"Floc",iprint)
  end subroutine dmft_get_gloc_matsubara_superc_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_lattice_main'
  subroutine dmft_get_gloc_matsubara_superc_lattice_main(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                                :: iprint
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    !
    real(8)                                           :: beta
    real(8)                                           :: xmu,eps
    real(8)                                           :: wini,wfin
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))    
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
       !G22(iw) = -[G11[iw]]*
       !G21(iw) =   G12[w]
       do i=1,Lmats
          zeta_mats(1,1,ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,2,ilat,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb) )
       enddo
    enddo
    !
    call start_timer
    Gmats=zero
    do ik=1,Lk
       call invert_gk_superc_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
       Gmats = Gmats + Gkmats*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_matsubara_lattice(wm,Gmats(1,:,:,:,:,:,:),"LG",iprint)
    call dmft_gloc_print_matsubara_lattice(wm,Gmats(2,:,:,:,:,:,:),"LF",iprint)
  end subroutine dmft_get_gloc_matsubara_superc_lattice_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_gij_main'
  subroutine dmft_get_gloc_matsubara_superc_gij_main(Hk,Wtk,Gmats,Fmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:)                       :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Fmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Smats           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer                                           :: iprint
    logical,optional                                  :: hk_symm(size(Hk,3))
    logical                                           :: hk_symm_(size(Hk,3))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkmats          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer                                           :: ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    integer                                           :: Lk,Nlso,Lmats,Nlat,Nspin,Norb,Nso
    !
    real(8)                                           :: beta
    real(8)                                           :: xmu,eps
    real(8)                                           :: wini,wfin
    !
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    call assert_shape(Fmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Fmats")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(Fkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Fkmats=zero
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats));zeta_mats=zero
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    zeta_mats=zero
    !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
    !G22(iw) = -[G11[iw]]*
    !G21(iw) =   G12[w]
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg( Smats(1,ilat,ispin,jspin,iorb,jorb,:) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call start_timer
    Gmats=zero
    do ik=1,Lk
       call invert_gk_superc_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
       Gmats = Gmats + Gkmats*Wtk(ik)
       Fmats = Fmats + Fkmats*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_matsubara_gij(wm,Gmats,"Gij",iprint)
    call dmft_gloc_print_matsubara_gij(wm,Fmats,"Fij",iprint)
  end subroutine dmft_get_gloc_matsubara_superc_gij_main
#undef FROUTINE


#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_matsubara_superc_main_mpi'
  subroutine dmft_get_gloc_matsubara_superc_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                              :: iprint
    character(len=*),optional                       :: mpi_split  
    character(len=1)                                :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    allocate(Gkmats(2,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(2,2,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)call start_timer
    Gmats=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_matsubara_superc_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_superc_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lmats))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_matsubara(wm,Gmats(1,:,:,:,:,:),"Gloc",iprint)
    if(mpi_master)call dmft_gloc_print_matsubara(wm,Gmats(2,:,:,:,:,:),"Floc",iprint)
  end subroutine dmft_get_gloc_matsubara_superc_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_lattice_main_mpi'
  subroutine dmft_get_gloc_matsubara_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    integer,intent(in)                                :: iprint
    character(len=*),optional                         :: mpi_split  
    character(len=1)                                  :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)call start_timer
    !
    allocate(Gkmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))    
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
       !G22(iw) = -[G11[iw]]*
       !G21(iw) =   G12[w]
       do i=1,Lmats
          zeta_mats(1,1,ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,2,ilat,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb) )
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Gmats=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_matsubara_superc_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_superc_lattice_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_matsubara_lattice(wm,Gmats(1,:,:,:,:,:,:),"LG",iprint)
    if(mpi_master)call dmft_gloc_print_matsubara_lattice(wm,Gmats(2,:,:,:,:,:,:),"LF",iprint)
  end subroutine dmft_get_gloc_matsubara_superc_lattice_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_gij_main_mpi'
  subroutine dmft_get_gloc_matsubara_superc_gij_main_mpi(MpiComm,Hk,Wtk,Gmats,Fmats,Smats,iprint,mpi_split,hk_symm)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:)                       :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Fmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Smats           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer                                           :: iprint
    character(len=*),optional                         :: mpi_split  
    character(len=1)                                  :: mpi_split_ 
    logical,optional                                  :: hk_symm(size(Hk,3))
    logical                                           :: hk_symm_(size(Hk,3))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkmats,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer                                           :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
    call assert_shape(Fmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Fmats")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(Fkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Fkmats=zero
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats));zeta_mats=zero
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
    !G22(iw) = -[G11[iw]]*
    !G21(iw) =   G12[w]
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg( Smats(1,ilat,ispin,jspin,iorb,jorb,:) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Gmats=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_matsubara_superc_gij_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_superc_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
          Gmats = Gmats + Gkmats*Wtk(ik)
          Fmats = Fmats + Fkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Ftmp=zero
       do ik=1,Lk
          call invert_gk_superc_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
          Gtmp = Gtmp + Gkmats*Wtk(ik)
          Ftmp = Ftmp + Fkmats*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp, Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       call Mpi_AllReduce(Ftmp, Fmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp,Ftmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_matsubara_gij(wm,Gmats,"Gij",iprint)
    if(mpi_master)call dmft_gloc_print_matsubara_gij(wm,Fmats,"Fij",iprint)
  end subroutine dmft_get_gloc_matsubara_superc_gij_main_mpi
#undef FROUTINE
#endif





  !##################################################################
  !##################################################################
  !##################################################################


#define FROUTINE 'dmft_get_gloc_matsubara_superc_1band'
  subroutine dmft_get_gloc_matsubara_superc_1band(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:),intent(in)            :: Hk              ![Nk]
    real(8),intent(in)                            :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                         :: Smats(:,:)
    complex(8),intent(inout)                      :: Gmats(2,size(Smats,2))
    logical,optional                              :: hk_symm(size(Hk,1))
    logical                                       :: hk_symm_(size(Hk,1))
    integer                                       :: iprint
    !
    complex(8),dimension(1,1,size(Hk))            :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(2,1,1,1,1,size(Smats,2)) :: Gmats_
    complex(8),dimension(2,1,1,1,1,size(Smats,2)) :: Smats_
    !
    call assert_shape(Smats,[2,size(Smats,2)],FROUTINE,"Smats")
    !
    Hk_(1,1,:)          = Hk
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_superc_main(Hk_,Wtk,Gmats_,Smats_,iprint,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_superc_1band
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_lattice_1band'
  subroutine dmft_get_gloc_matsubara_superc_lattice_1band(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                   :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8),intent(in)                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                    :: Smats(:,:,:)
    complex(8),intent(inout)                                 :: Gmats(2,size(Hk,1),size(Smats,3))
    logical,optional                                         :: hk_symm(size(Hk,1))
    logical                                                  :: hk_symm_(size(Hk,1))
    integer                                                  :: iprint
    !
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Smats,3)) :: Gmats_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Smats,3)) :: Smats_
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[2,size(Hk,1),size(Smats,3)],FROUTINE,"Smats")
    !
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_superc_lattice_main(Hk,Wtk,Gmats_,Smats_,iprint,hk_symm_)
    Gmats(:,:,:) = Gmats_(:,:,1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_superc_lattice_1band
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_matsubara_superc_1band_mpi'
  subroutine dmft_get_gloc_matsubara_superc_1band_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                       :: MpiComm
    complex(8),dimension(:),intent(in)            :: Hk              ![Nk]
    real(8),intent(in)                            :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                         :: Smats(:,:)
    complex(8),intent(inout)                      :: Gmats(2,size(Smats,2))
    logical,optional                              :: hk_symm(size(Hk,1))
    logical                                       :: hk_symm_(size(Hk,1))
    integer                                       :: iprint
    character(len=*),optional                     :: mpi_split  
    character(len=1)                              :: mpi_split_ 
    !
    complex(8),dimension(1,1,size(Hk))            :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(2,1,1,1,1,size(Smats,2)) :: Gmats_
    complex(8),dimension(2,1,1,1,1,size(Smats,2)) :: Smats_
    !
    call assert_shape(Smats,[2,size(Smats,2)],FROUTINE,"Smats")
    !
    Hk_(1,1,:)          = Hk
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_superc_main_mpi(MpiComm,Hk_,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
    Gmats(:,:) = Gmats_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_superc_1band_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_lattice_1band_mpi'
  subroutine dmft_get_gloc_matsubara_superc_lattice_1band_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                                  :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                   :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8),intent(in)                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                    :: Smats(:,:,:)
    complex(8),intent(inout)                                 :: Gmats(2,size(Hk,1),size(Smats,3))
    logical,optional                                         :: hk_symm(size(Hk,1))
    logical                                                  :: hk_symm_(size(Hk,1))
    integer                                                  :: iprint
    character(len=*),optional                                :: mpi_split  
    character(len=1)                                         :: mpi_split_ 
    !
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Smats,3)) :: Gmats_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Smats,3)) :: Smats_
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[2,size(Hk,1),size(Smats,3)],FROUTINE,"Smats")
    !
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    call dmft_get_gloc_matsubara_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
    Gmats(:,:,:) = Gmats_(:,:,1,1,1,1,:)
  end subroutine dmft_get_gloc_matsubara_superc_lattice_1band_mpi
#undef FROUTINE
#endif




  !##################################################################
  !##################################################################
  !##################################################################

#define FROUTINE 'dmft_get_gloc_matsubara_superc_Nband'
  subroutine dmft_get_gloc_matsubara_superc_Nband(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Norb][Norb][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                           :: Smats(:,:,:,:)  ![2][Norb][Norb][Lmats]
    complex(8),intent(inout)                                        :: Gmats(:,:,:,:)  !as Smats
    logical,optional                                                :: hk_symm(size(Hk,3))
    logical                                                         :: hk_symm_(size(Hk,3))
    integer                                                         :: iprint
    !
    complex(8),&
         dimension(2,1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Gmats_
    complex(8),&
         dimension(2,1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Smats_
    integer                                                         :: Nspin,Norb,Nso,Lmats,Lreal
    !
    Nspin = 1
    Norb  = size(Smats,2)
    Lmats = size(Smats,4)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_matsubara_superc_main(Hk,Wtk,Gmats_,Smats_,iprint,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_superc_Nband
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_lattice_Nband'
  subroutine dmft_get_gloc_matsubara_superc_lattice_Nband(Hk,Wtk,Gmats,Smats,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                                        :: Hk               ![Nlat*Norb][Nlat*Norb][Nk]
    real(8),intent(in)                                                            :: Wtk(size(Hk,3))  ![Nk]
    complex(8),intent(in)                                                         :: Smats(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    complex(8),intent(inout)                                                      :: Gmats(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    logical,optional                                                              :: hk_symm(size(Hk,3))
    logical                                                                       :: hk_symm_(size(Hk,3))
    integer                                                                       :: iprint
    !
    complex(8),&
         dimension(2,size(Smats,2),1,1,size(Smats,3),size(Smats,4),size(Smats,5)) :: Gmats_
    complex(8),&
         dimension(2,size(Smats,2),1,1,size(Smats,3),size(Smats,4),size(Smats,5)) :: Smats_
    integer                                                                       :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Nlat  = size(Smats,2)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Nlat,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call dmft_get_gloc_matsubara_superc_lattice_main(Hk,Wtk,Gmats_,Smats_,iprint,hk_symm_)
    Gmats(:,:,:,:,:) = Gmats_(:,:,1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_superc_lattice_Nband
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_matsubara_superc_Nband_mpi'
  subroutine dmft_get_gloc_matsubara_superc_Nband_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Norb][Norb][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                           :: Smats(:,:,:,:)  ![2][Norb][Norb][Lmats]
    complex(8),intent(inout)                                        :: Gmats(:,:,:,:)  !as Smats
    logical,optional                                                :: hk_symm(size(Hk,3))
    logical                                                         :: hk_symm_(size(Hk,3))
    integer                                                         :: iprint
    character(len=*),optional                                       :: mpi_split  
    character(len=1)                                                :: mpi_split_ 
    !
    complex(8),&
         dimension(2,1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Gmats_
    complex(8),&
         dimension(2,1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Smats_
    integer                                                         :: Nspin,Norb,Nso,Lmats,Lreal
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Norb  = size(Smats,2)
    Lmats = size(Smats,4)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_gloc_matsubara_superc_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
    Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_superc_Nband_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_superc_lattice_Nband_mpi'
  subroutine dmft_get_gloc_matsubara_superc_lattice_Nband_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
    integer                                                                       :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                                        :: Hk               ![Nlat*Norb][Nlat*Norb][Nk]
    real(8),intent(in)                                                            :: Wtk(size(Hk,3))  ![Nk]
    complex(8),intent(in)                                                         :: Smats(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    complex(8),intent(inout)                                                      :: Gmats(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lmats]
    logical,optional                                                              :: hk_symm(size(Hk,3))
    logical                                                                       :: hk_symm_(size(Hk,3))
    integer                                                                       :: iprint
    character(len=*),optional                                                     :: mpi_split  
    character(len=1)                                                              :: mpi_split_ 
    !
    complex(8),&
         dimension(2,size(Smats,2),1,1,size(Smats,3),size(Smats,4),size(Smats,5)) :: Gmats_
    complex(8),&
         dimension(2,size(Smats,2),1,1,size(Smats,3),size(Smats,4),size(Smats,5)) :: Smats_
    integer                                                                       :: Nspin,Norb,Nso,Nlso,Lmats
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Nlat  = size(Smats,2)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Gmats,[2,Nlat,Norb,Norb,Lmats],FROUTINE,"Gmats")
    !
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call dmft_get_gloc_matsubara_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
    Gmats(:,:,:,:,:) = Gmats_(:,:,1,1,:,:,:)
  end subroutine dmft_get_gloc_matsubara_superc_lattice_Nband_mpi
#undef FROUTINE
#endif






#define FROUTINE 'dmft_get_gloc_realaxis_superc_main'
  subroutine dmft_get_gloc_realaxis_superc_main(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint
    logical,dimension(size(Hk,3)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin  
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    allocate(Gkreal(2,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(2,2,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
            nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
            conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    call start_timer
    Greal=zero
    do ik=1,Lk
       call invert_gk_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Greal = Greal + Gkreal*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_realaxis(wr,Greal(1,:,:,:,:,:),"Gloc",iprint)
    call dmft_gloc_print_realaxis(wr,Greal(2,:,:,:,:,:),"Floc",iprint)
  end subroutine dmft_get_gloc_realaxis_superc_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_lattice_main'
  subroutine dmft_get_gloc_realaxis_superc_lattice_main(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                                :: iprint
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                           :: beta
    real(8)                                           :: xmu,eps
    real(8)                                           :: wini,wfin  
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       !SYMMETRIES in real-frequencies   [assuming a real order parameter]
       !G22(w)  = -[G11[-w]]*
       !G21(w)  =   G12[w]   
       do i=1,Lreal
          zeta_real(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
               nn2so_reshape(Sreal(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(1,2,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,1,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
               conjg( nn2so_reshape(Sreal(1,ilat,:,:,:,:,Lreal+1-i),Nspin,Norb) )
       enddo
    enddo
    !
    call start_timer
    Greal=zero
    do ik=1,Lk
       call invert_gk_superc_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
       Greal = Greal + Gkreal*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_realaxis_lattice(wr,Greal(1,:,:,:,:,:,:),"LG",iprint)
    call dmft_gloc_print_realaxis_lattice(wr,Greal(2,:,:,:,:,:,:),"LF",iprint)
  end subroutine dmft_get_gloc_realaxis_superc_lattice_main
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_gij_main'
  subroutine dmft_get_gloc_realaxis_superc_gij_main(Hk,Wtk,Greal,Freal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:)                       :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Greal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Freal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sreal           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer                                           :: iprint
    logical,optional                                  :: hk_symm(size(Hk,3))
    logical                                           :: hk_symm_(size(Hk,3))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkreal          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer                                           :: ik,Lk,Nlso,Nlat,Nspin,Nso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js,Lreal
    real(8)                                           :: beta
    real(8)                                           :: xmu,eps
    real(8)                                           :: wini,wfin  
    !
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    call assert_shape(Freal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Freal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Fkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Fkreal=zero
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal));zeta_real=zero
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    zeta_real=zero
    !SYMMETRIES in real-frequencies   [assuming a real order parameter]
    !G22(w)  = -[G11[-w]]*
    !G21(w)  =   G12[w]    
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu
             zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu )
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - &
                        Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - &
                        Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - &
                        Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) + &
                        conjg( Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call start_timer
    Greal=zero
    do ik=1,Lk
       call invert_gk_superc_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
       Greal = Greal + Gkreal*Wtk(ik)
       Freal = Freal + Fkreal*Wtk(ik)
       call eta(ik,Lk)
    end do
    call stop_timer
    call dmft_gloc_print_realaxis_gij(wm,Greal,"Gij",iprint)
    call dmft_gloc_print_realaxis_gij(wm,Freal,"Fij",iprint)
  end subroutine dmft_get_gloc_realaxis_superc_gij_main
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_realaxis_superc_main_mpi'
  subroutine dmft_get_gloc_realaxis_superc_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                              :: iprint
    character(len=*),optional                       :: mpi_split  
    character(len=1)                                :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    allocate(Gkreal(2,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(2,2,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
            nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
            conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)call start_timer
    Greal=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_realaxis_superc_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_superc_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lreal))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_realaxis(wr,Greal(1,:,:,:,:,:),"Gloc",iprint)
    if(mpi_master)call dmft_gloc_print_realaxis(wr,Greal(2,:,:,:,:,:),"Floc",iprint)
  end subroutine dmft_get_gloc_realaxis_superc_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_lattice_main_mpi'
  subroutine dmft_get_gloc_realaxis_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    integer,intent(in)                                :: iprint
    character(len=*),optional                         :: mpi_split  
    character(len=1)                                  :: mpi_split_ 
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    if(mpi_master)call start_timer
    !
    allocate(Gkreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       !SYMMETRIES in real-frequencies   [assuming a real order parameter]
       !G22(w)  = -[G11[-w]]*
       !G21(w)  =   G12[w]   
       do i=1,Lreal
          zeta_real(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
               nn2so_reshape(Sreal(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(1,2,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,1,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
               conjg( nn2so_reshape(Sreal(1,ilat,:,:,:,:,Lreal+1-i),Nspin,Norb) )
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Greal=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_realaxis_superc_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_superc_lattice_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_realaxis_lattice(wr,Greal(1,:,:,:,:,:,:),"LG",iprint)
    if(mpi_master)call dmft_gloc_print_realaxis_lattice(wr,Greal(2,:,:,:,:,:,:),"LF",iprint)
  end subroutine dmft_get_gloc_realaxis_superc_lattice_main_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_main_gij_mpi'
  subroutine dmft_get_gloc_realaxis_superc_gij_main_mpi(MpiComm,Hk,Wtk,Greal,Freal,Sreal,iprint,mpi_split,hk_symm)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:)                       :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Greal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Freal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sreal           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer                                           :: iprint
    character(len=*),optional                         :: mpi_split  
    character(len=1)                                  :: mpi_split_ 
    logical,optional                                  :: hk_symm(size(Hk,3))
    logical                                           :: hk_symm_(size(Hk,3))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkreal,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer                                           :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
    call assert_shape(Freal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Freal")
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
    !
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Fkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Fkreal=zero
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal));zeta_real=zero
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    !SYMMETRIES in real-frequencies   [assuming a real order parameter]
    !G22(w)  = -[G11[-w]]*
    !G21(w)  =   G12[w]
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu
             zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu )
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - &
                        Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - &
                        Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - &
                        Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) + &
                        conjg( Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Greal=zero
    select case(mpi_split_)
    case default
       stop "dmft_get_gloc_normal_realaxis_gij_main_mpi: ! mpi_split_ in ['w','k'] "
       !
    case ('w')
       do ik=1,Lk
          call invert_gk_superc_gij_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
          Greal = Greal + Gkreal*Wtk(ik)
          Freal = Freal + Fkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       !
    case ('k')
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Ftmp=zero
       do ik=1,Lk
          call invert_gk_superc_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
          Gtmp = Gtmp + Gkreal*Wtk(ik)
          Ftmp = Ftmp + Fkreal*Wtk(ik)
          if(mpi_master)call eta(ik,Lk)
       end do
       call Mpi_AllReduce(Gtmp, Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       call Mpi_AllReduce(Ftmp, Freal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
       deallocate(Gtmp,Ftmp)
       !
    end select
    if(mpi_master)call stop_timer
    if(mpi_master)call dmft_gloc_print_realaxis_gij(wm,Greal,"Gij",iprint)
    if(mpi_master)call dmft_gloc_print_realaxis_gij(wm,Freal,"Fij",iprint)
  end subroutine dmft_get_gloc_realaxis_superc_gij_main_mpi
#undef FROUTINE
#endif




  !##################################################################
  !##################################################################
  !##################################################################


#define FROUTINE 'dmft_get_gloc_realaxis_superc_1band'
  subroutine dmft_get_gloc_realaxis_superc_1band(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:),intent(in)            :: Hk              ![Nk]
    real(8),intent(in)                            :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                         :: Sreal(:,:)
    complex(8),intent(inout)                      :: Greal(2,size(Sreal,2))
    logical,optional                              :: hk_symm(size(Hk,1))
    logical                                       :: hk_symm_(size(Hk,1))
    integer                                       :: iprint
    !
    complex(8),dimension(1,1,size(Hk))            :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Greal_
    complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Sreal_
    !
    call assert_shape(Sreal,[2,size(Sreal,2)],FROUTINE,"Sreal")
    !
    Hk_(1,1,:)          = Hk
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_superc_main(Hk_,Wtk,Greal_,Sreal_,iprint,hk_symm_)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_superc_1band
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_lattice_1band'
  subroutine dmft_get_gloc_realaxis_superc_lattice_1band(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                   :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8),intent(in)                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                    :: Sreal(:,:,:)
    complex(8),intent(inout)                                 :: Greal(2,size(Hk,1),size(Sreal,3))
    logical,optional                                         :: hk_symm(size(Hk,1))
    logical                                                  :: hk_symm_(size(Hk,1))
    integer                                                  :: iprint
    !
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Sreal,3)) :: Greal_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Sreal,3)) :: Sreal_
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,size(Hk,1),size(Sreal,3)],FROUTINE,"Sreal")
    !
    Sreal_(:,:,1,1,1,1,:) = Sreal(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_superc_lattice_main(Hk,Wtk,Greal_,Sreal_,iprint,hk_symm_)
    Greal(:,:,:) = Greal_(:,:,1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_superc_lattice_1band
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_realaxis_superc_1band_mpi'
  subroutine dmft_get_gloc_realaxis_superc_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                       :: MpiComm
    complex(8),dimension(:),intent(in)            :: Hk              ![Nk]
    real(8),intent(in)                            :: Wtk(size(Hk))   ![Nk]
    complex(8),intent(in)                         :: Sreal(:,:)
    complex(8),intent(inout)                      :: Greal(2,size(Sreal,2))
    logical,optional                              :: hk_symm(size(Hk,1))
    logical                                       :: hk_symm_(size(Hk,1))
    integer                                       :: iprint
    character(len=*),optional                     :: mpi_split  
    character(len=1)                              :: mpi_split_ 
    !
    complex(8),dimension(1,1,size(Hk))            :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
    complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Greal_
    complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Sreal_
    !
    call assert_shape(Sreal,[2,size(Sreal,2)],FROUTINE,"Sreal")
    !
    Hk_(1,1,:)          = Hk
    Sreal_(:,1,1,1,1,:) = Sreal(:,:)
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_superc_main_mpi(MpiComm,Hk_,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
    Greal(:,:) = Greal_(:,1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_superc_1band_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_lattice_1band_mpi'
  subroutine dmft_get_gloc_realaxis_superc_lattice_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                                  :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                   :: Hk              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8),intent(in)                                       :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                    :: Sreal(:,:,:)
    complex(8),intent(inout)                                 :: Greal(2,size(Hk,1),size(Sreal,3))
    logical,optional                                         :: hk_symm(size(Hk,1))
    logical                                                  :: hk_symm_(size(Hk,1))
    integer                                                  :: iprint
    character(len=*),optional                                :: mpi_split  
    character(len=1)                                         :: mpi_split_ 
    !
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Sreal,3)) :: Greal_
    complex(8),dimension(2,size(Hk,1),1,1,1,1,size(Sreal,3)) :: Sreal_
    !
    call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,size(Hk,1),size(Sreal,3)],FROUTINE,"Sreal")
    !
    Sreal_(:,:,1,1,1,1,:) = Sreal(:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    call dmft_get_gloc_realaxis_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
    Greal(:,:,:) = Greal_(:,:,1,1,1,1,:)
  end subroutine dmft_get_gloc_realaxis_superc_lattice_1band_mpi
#undef FROUTINE
#endif




  !##################################################################
  !##################################################################
  !##################################################################



#define FROUTINE 'dmft_get_gloc_realaxis_superc_Nband'
  subroutine dmft_get_gloc_realaxis_superc_Nband(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Norb][Norb][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                           :: Sreal(:,:,:,:)  ![2][Norb][Norb][Lreal]
    complex(8),intent(inout)                                        :: Greal(:,:,:,:)  !as Sreal
    logical,optional                                                :: hk_symm(size(Hk,3))
    logical                                                         :: hk_symm_(size(Hk,3))
    integer                                                         :: iprint
    !
    complex(8),&
         dimension(2,1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Greal_
    complex(8),&
         dimension(2,1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Sreal_
    integer                                                         :: Nspin,Norb,Nso,Lreal
    !
    Nspin = 1
    Norb  = size(Sreal,2)
    Lreal = size(Sreal,4)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    call dmft_get_gloc_realaxis_superc_main(Hk,Wtk,Greal_,Sreal_,iprint,hk_symm_)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_superc_Nband
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_lattice_Nband'
  subroutine dmft_get_gloc_realaxis_superc_lattice_Nband(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
    complex(8),dimension(:,:,:),intent(in)                                        :: Hk               ![Nlat*Norb][Nlat*Norb][Nk]
    real(8),intent(in)                                                            :: Wtk(size(Hk,3))  ![Nk]
    complex(8),intent(in)                                                         :: Sreal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lreal]
    complex(8),intent(inout)                                                      :: Greal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lreal]
    logical,optional                                                              :: hk_symm(size(Hk,3))
    logical                                                                       :: hk_symm_(size(Hk,3))
    integer                                                                       :: iprint
    !
    complex(8),&
         dimension(2,size(Sreal,2),1,1,size(Sreal,3),size(Sreal,4),size(Sreal,5)) :: Greal_
    complex(8),&
         dimension(2,size(Sreal,2),1,1,size(Sreal,3),size(Sreal,4),size(Sreal,5)) :: Sreal_
    integer                                                                       :: Nlat,Nspin,Norb,Nso,Nlso,Lreal
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Nlat  = size(Sreal,2)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nlat,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Nlat,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    Sreal_(:,:,1,1,:,:,:) = Sreal(:,:,:,:,:)
    call dmft_get_gloc_realaxis_superc_lattice_main(Hk,Wtk,Greal_,Sreal_,iprint,hk_symm_)
    Greal(:,:,:,:,:) = Greal_(:,:,1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_superc_lattice_Nband
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_gloc_realaxis_superc_Nband_mpi'
  subroutine dmft_get_gloc_realaxis_superc_Nband_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Norb][Norb][Nk]
    real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
    complex(8),intent(in)                                           :: Sreal(:,:,:,:)  ![2][Norb][Norb][Lreal]
    complex(8),intent(inout)                                        :: Greal(:,:,:,:)  !as Sreal
    logical,optional                                                :: hk_symm(size(Hk,3))
    logical                                                         :: hk_symm_(size(Hk,3))
    integer                                                         :: iprint
    character(len=*),optional                                       :: mpi_split  
    character(len=1)                                                :: mpi_split_ 
    !
    complex(8),&
         dimension(2,1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Greal_
    complex(8),&
         dimension(2,1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Sreal_
    integer                                                         :: Nspin,Norb,Nso,Lreal
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Norb  = size(Sreal,2)
    Lreal = size(Sreal,4)
    Nso   = Nspin*Norb
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
    call dmft_get_gloc_realaxis_superc_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
    Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_superc_Nband_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_superc_lattice_Nband_mpi'
  subroutine dmft_get_gloc_realaxis_superc_lattice_Nband_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
    integer                                                                       :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                                        :: Hk               ![Nlat*Norb][Nlat*Norb][Nk]
    real(8),intent(in)                                                            :: Wtk(size(Hk,3))  ![Nk]
    complex(8),intent(in)                                                         :: Sreal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lreal]
    complex(8),intent(inout)                                                      :: Greal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lreal]
    logical,optional                                                              :: hk_symm(size(Hk,3))
    logical                                                                       :: hk_symm_(size(Hk,3))
    integer                                                                       :: iprint
    character(len=*),optional                                                     :: mpi_split  
    character(len=1)                                                              :: mpi_split_ 
    !
    complex(8),&
         dimension(2,size(Sreal,2),1,1,size(Sreal,3),size(Sreal,4),size(Sreal,5)) :: Greal_
    complex(8),&
         dimension(2,size(Sreal,2),1,1,size(Sreal,3),size(Sreal,4),size(Sreal,5)) :: Sreal_
    integer                                                                       :: Nspin,Norb,Nso,Nlso,Lreal
    !
    mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = 1
    Nlat  = size(Sreal,2)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
    call assert_shape(Sreal,[2,Nlat,Norb,Norb,Lreal],FROUTINE,"Sreal")
    call assert_shape(Greal,[2,Nlat,Norb,Norb,Lreal],FROUTINE,"Greal")
    !
    Sreal_(:,:,1,1,:,:,:) = Sreal(:,:,:,:,:)
    call dmft_get_gloc_realaxis_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
    Greal(:,:,:,:,:) = Greal_(:,:,1,1,:,:,:)
  end subroutine dmft_get_gloc_realaxis_superc_lattice_Nband_mpi
#undef FROUTINE
#endif











  !--------------------------------------------------------------------!
  ! PURPOSE: invert the Gk matrix in all cases:
  ! + normal/superc
  ! + serial/parallel
  ! + single site/tridiagonal/lattice/full gij
  !--------------------------------------------------------------------!
  !
  ! INVERT_GK_NORMAL(_MPI)
  !
  !SERIAL (OR PARALLEL ON K):
  subroutine invert_gk_normal(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:),intent(in)        :: zeta    ![Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)          :: Hk      ![Nspin*Norb][Nspin*Norb]
    logical,intent(in)                            :: hk_symm                
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkout   ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gktmp   ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable         :: Gmatrix ![Nspin*Norb][Nspin*Norb]
    integer                                       :: Nspin,Norb,Nso,Lfreq
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    !
    Nspin = size(Gkout,1)
    Norb  = size(Gkout,3)
    Lfreq = size(zeta,3)
    Nso   = Nspin*Norb
    !testing
    call assert_shape(zeta,[Nso,Nso,Lfreq],"invert_gk_normal","zeta")
    call assert_shape(Hk,[Nso,Nso],"invert_gk_normal","Hk")
    call assert_shape(Gkout,[Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal","Gkout")
    !
    allocate(Gktmp(Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nso,Nso))
    Gktmp=zero
    do i=1,Lfreq
       Gmatrix  = zeta(:,:,i) - Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gktmp(ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine invert_gk_normal

  !PARALLEL ON FREQ:
#ifdef _MPI
  subroutine invert_gk_normal_mpi(MpiComm,zeta,Hk,hk_symm,Gkout)
    integer                                       :: MpiComm
    complex(8),dimension(:,:,:),intent(in)        :: zeta    ![Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)          :: Hk      ![Nspin*Norb][Nspin*Norb]
    logical,intent(in)                            :: hk_symm                
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkout   ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gktmp   ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable         :: Gmatrix ![Nspin*Norb][Nspin*Norb]
    integer                                       :: Nspin,Norb,Nso,Lfreq
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nspin = size(Gkout,1)
    Norb  = size(Gkout,3)
    Lfreq = size(zeta,3)
    Nso   = Nspin*Norb
    !testing
    call assert_shape(zeta,[Nso,Nso,Lfreq],"invert_gk_normal_mpi","zeta")
    call assert_shape(Hk,[Nso,Nso],"invert_gk_normal_mpi","Hk")
    call assert_shape(Gkout,[Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_mpi","Gkout")
    !
    allocate(Gktmp(Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nso,Nso))
    Gktmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = zeta(:,:,i) - Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gktmp(ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_normal_mpi
#endif







  !
  ! INVERT_GK_NORMAL_LATTICE(_MPI)
  !
  !SERIAL (OR PARALLEL ON K)
  subroutine invert_gk_normal_lattice(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                              :: hk_symm
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    !
    complex(8),dimension(:,:,:),allocatable         :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing:
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_lattice","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_lattice","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_lattice","Gkout")
    !
    if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
         stop "invert_gk_normal_lattice: lattice_Gamma_mats & lattice_Gamma_real both allocated"
    if(allocated(lattice_gamma_mats))then
       call assert_shape(lattice_gamma_mats,[Nlso,Nlso,Lfreq],"invert_gk_normal_lattice","lattice_gamma_mats")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_mats
    endif
    if(allocated(lattice_gamma_real))then
       call assert_shape(lattice_gamma_real,[Nlso,Nlso,Lfreq],"invert_gk_normal_lattice","lattice_gamma_real")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_real
    endif
    !
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1,Lfreq
       Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nlat,Nso) - Hk
       if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                      Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine invert_gk_normal_lattice

  !PARALLEL ON FREQ:
#ifdef _MPI
  subroutine invert_gk_normal_lattice_mpi(MpiComm,zeta,Hk,hk_symm,Gkout)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                              :: hk_symm
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    !
    complex(8),dimension(:,:,:),allocatable         :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_lattice_mpi","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_lattice_mpi","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_lattice_mpi","Gkout")
    !
    if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
         stop "invert_gk_normal_lattice_mpi: lattice_Gamma_mats & lattice_Gamma_real both allocated"
    if(allocated(lattice_gamma_mats))then
       call assert_shape(lattice_gamma_mats,[Nlso,Nlso,Lfreq],"invert_gk_normal_lattice_mpi","lattice_gamma_mats")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_mats
    endif
    if(allocated(lattice_gamma_real))then
       call assert_shape(lattice_gamma_real,[Nlso,Nlso,Lfreq],"invert_gk_normal_lattice_mpi","lattice_gamma_real")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_real
    endif
    !
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nlat,Nso) - Hk
       if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                      Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_normal_lattice_mpi
#endif












  !
  ! INVERT_GK_NORMAL_TRIDIAG(_MPI)
  !
  !SERIAL (OR PARALLEL ON K)
  subroutine invert_gk_normal_tridiag(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                              :: hk_symm
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    !
    complex(8),dimension(:,:,:,:),allocatable       :: Gembed  ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:),allocatable         :: Diag
    complex(8),dimension(:,:,:),allocatable         :: Sub
    complex(8),dimension(:,:,:),allocatable         :: Over
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:),allocatable         :: Gmatrix ![Nlat][Nspin*Norb][Nspin*Norb]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_tridiag","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_tridiag","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_tridiag","Gkout")
    !
    if(allocated(local_gamma_mats).AND.allocated(local_gamma_real))&
         stop "invert_gk_normal_tridiag: local_Gamma_mats & local_Gamma_real both allocated"
    if(allocated(local_gamma_mats))then
       call assert_shape(local_gamma_mats,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_tridiag","local_gamma_mats")
       allocate(Gembed(Nlat,Nso,Nso,Lfreq));Gembed=local_gamma_mats
    endif
    if(allocated(lattice_gamma_real))then
       call assert_shape(lattice_gamma_real,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_tridiag","local_gamma_real")
       allocate(Gembed(Nlat,Nso,Nso,Lfreq));Gembed=local_gamma_real
    endif
    !
    allocate(Sub(Nlat-1,Nso,Nso))
    allocate(Diag(Nlat,Nso,Nso))
    allocate(Over(Nlat-1,Nso,Nso))
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlat,Nso,Nso))
    Gktmp=zero
    do i=1,Lfreq
       call get_tridiag(Nlat,Nso,Hk,Sub,Diag,Over)
       Diag = zeta(:,:,:,i) - Diag
       if(allocated(Gembed))Diag = Diag - Gembed(:,:,:,i)
       call inv_tridiag(Nlat,Nso,-Sub,Diag,-Over,Gmatrix)
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(ilat,io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine invert_gk_normal_tridiag

  !PARALLEL ON FREQ:
#ifdef _MPI
  subroutine invert_gk_normal_tridiag_mpi(MpiComm,zeta,Hk,hk_symm,Gkout)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:),intent(in)        :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                              :: hk_symm
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    !
    complex(8),dimension(:,:,:),allocatable         :: Diag
    complex(8),dimension(:,:,:),allocatable         :: Sub
    complex(8),dimension(:,:,:),allocatable         :: Over
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:),allocatable         :: Gmatrix ![Nlat][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:),allocatable       :: Gembed  ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]    
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                         :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_tridiag_mpi","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_tridiag_mpi","Hk")
    call assert_shape(Gkout,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_tridiag_mpi","Gkout")
    !
    if(allocated(local_gamma_mats).AND.allocated(local_gamma_real))&
         stop "invert_gk_normal_tridiag_mpi: local_Gamma_mats & local_Gamma_real both allocated"
    if(allocated(local_gamma_mats))then
       call assert_shape(local_gamma_mats,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_tridiag_mpi","local_gamma_mats")
       allocate(Gembed(Nlat,Nso,Nso,Lfreq));Gembed=local_gamma_mats
    endif
    if(allocated(lattice_gamma_real))then
       call assert_shape(lattice_gamma_real,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_tridiag_mpi","local_gamma_real")
       allocate(Gembed(Nlat,Nso,Nso,Lfreq));Gembed=local_gamma_real
    endif
    !
    allocate(Sub(Nlat-1,Nso,Nso))
    allocate(Diag(Nlat,Nso,Nso))
    allocate(Over(Nlat-1,Nso,Nso))
    allocate(Gktmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlat,Nso,Nso))
    Gktmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       call get_tridiag(Nlat,Nso,Hk,Sub,Diag,Over)
       Diag = zeta(:,:,:,i) - Diag
       if(allocated(Gembed))Diag = Diag - Gembed(:,:,:,i)
       call inv_tridiag(Nlat,Nso,-Sub,Diag,-Over,Gmatrix)
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      Gktmp(ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(ilat,io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_normal_tridiag_mpi
#endif








  !
  ! INVERT_GK_NORMAL_GIJ(_MPI)
  !
  !SERIAL (OR PARALLEL ON K)
  subroutine invert_gk_normal_gij(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:),intent(in)          :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable           :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    integer                                           :: Lfreq
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_gij","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_gij","Hk")
    call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_gij","Gkout")
    !
    if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
         stop "invert_gk_normal_lattice: lattice_Gamma_mats & lattice_Gamma_real both allocated"
    if(allocated(lattice_gamma_mats))then
       call assert_shape(lattice_gamma_mats,[Nlso,Nlso,Lfreq],"invert_gk_normal_gij","lattice_gamma_mats")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_mats
    endif
    if(allocated(lattice_gamma_real))then
       call assert_shape(lattice_gamma_real,[Nlso,Nlso,Lfreq],"invert_gk_normal_gij","lattice_gamma_real")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_real
    endif
    !
    allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1,Lfreq
       Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nlat,Nso) - Hk
       if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                         Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine invert_gk_normal_gij

  !PARALLEL ON FREQ:
#ifdef _MPI
  subroutine invert_gk_normal_gij_mpi(MpiComm,zeta,Hk,hk_symm,Gkout)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:,:),intent(in)          :: zeta    ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable           :: Gembed  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    integer                                           :: Lfreq
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[Nlat,Nso,Nso,Lfreq],"invert_gk_normal_gij_mpi","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_normal_gij_mpi","Hk")
    call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_normal_gij_mpi","Gkout")
    !
    if(allocated(lattice_gamma_mats).AND.allocated(lattice_gamma_real))&
         stop "invert_gk_normal_gij_mpi: lattice_Gamma_mats & lattice_Gamma_real both allocated"
    if(allocated(lattice_gamma_mats))then
       call assert_shape(lattice_gamma_mats,[Nlso,Nlso,Lfreq],"invert_gk_normal_gij_mpi","lattice_gamma_mats")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_mats
    endif
    if(allocated(lattice_gamma_real))then
       call assert_shape(lattice_gamma_real,[Nlso,Nlso,Lfreq],"invert_gk_normal_gij","lattice_gamma_real")
       allocate(Gembed(Nlso,Nlso,Lfreq));Gembed=lattice_gamma_real
    endif
    !
    allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(Nlso,Nlso))
    Gktmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = blocks_to_matrix(zeta(:,:,:,i),Nlat,Nso) - Hk
       if(allocated(Gembed))Gmatrix = Gmatrix - Gembed(:,:,i)
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                         Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_normal_gij_mpi
#endif








  !
  ! INVERT_GK_NORMAL(_MPI)
  !
  !SERIAL (OR PARALLEL ON K):
  subroutine invert_gk_superc(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nspin*Norb][Nspin*Norb]
    logical,intent(in)                              :: hk_symm !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix ![2*Nspin*Norb][2*Nspin*Norb]
    integer                                         :: Nspin,Norb,Nso,Lfreq
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
    !
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,5)
    Nso   = Nspin*Norb
    call assert_shape(zeta,[2,2,Nso,Nso,Lfreq],"invert_gk_superc","zeta")
    call assert_shape(Hk,[Nso,Nso],"invert_gk_superc","Hk")
    call assert_shape(Gkout,[2,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc","Gkout")
    !
    allocate(Gktmp(2,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nso,2*Nso))
    Gkout = zero
    Gktmp = zero
    do i=1,Lfreq
       Gmatrix  = zero
       Gmatrix(1:Nso,1:Nso)             = zeta(1,1,:,:,i) - Hk
       Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta(1,2,:,:,i)
       Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta(2,1,:,:,i)
       Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta(2,2,:,:,i) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gktmp(1,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                   Gktmp(2,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nso+jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine invert_gk_superc

  !PARALLEL ON FREQ:
#ifdef _MPI
  subroutine invert_gk_superc_mpi(MpiComm,zeta,Hk,hk_symm,Gkout)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)            :: Hk      ![Nspin*Norb][Nspin*Norb]
    logical,intent(in)                              :: hk_symm !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkout   ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gktmp   ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix ![2*Nspin*Norb][2*Nspin*Norb]
    integer                                         :: Nspin,Norb,Nso,Lfreq
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nspin = size(Gkout,2)
    Norb  = size(Gkout,4)
    Lfreq = size(zeta,5)
    Nso   = Nspin*Norb
    call assert_shape(zeta,[2,2,Nso,Nso,Lfreq],"invert_gk_superc_mpi","zeta")
    call assert_shape(Hk,[Nso,Nso],"invert_gk_superc_mpi","Hk")
    call assert_shape(Gkout,[2,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_mpi","Gkout")
    !
    allocate(Gktmp(2,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nso,2*Nso))
    Gkout = zero
    Gktmp = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = zero
       Gmatrix(1:Nso,1:Nso)             = zeta(1,1,:,:,i) - Hk
       Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta(1,2,:,:,i)
       Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta(2,1,:,:,i)
       Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta(2,2,:,:,i) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Gktmp(1,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                   Gktmp(2,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nso+jo)
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_superc_mpi
#endif




  !
  ! INVERT_GK_NORMAL_LATTICE(_MPI)
  !
  !SERIAL (OR PARALLEL ON K)
  subroutine invert_gk_superc_lattice(zeta,Hk,hk_symm,Gkout)
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm !
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![2*Nlat*Nspin*Norb][2*Nlat*Nspin*Norb]
    integer                                           :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                           :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlat  = size(zeta,3)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[2,2,Nlat,Nso,Nso,Lfreq],"invert_gk_superc_lattice","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_superc_lattice","Hk")
    call assert_shape(Gkout,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_lattice","Gkout")
    !
    allocate(Gktmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nlso,2*Nlso))
    Gkout = zero
    Gktmp = zero
    do i=1,Lfreq
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta(1,1,:,:,:,i),Nlat,Nso) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(1,2,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta(2,1,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(2,2,:,:,:,i),Nlat,Nso) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                      Gktmp(1,ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                      Gktmp(2,ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout = Gktmp
  end subroutine invert_gk_superc_lattice


  !PARALLEL ON FREQ:
#ifdef _MPI
  subroutine invert_gk_superc_lattice_mpi(MpiComm,zeta,Hk,hk_symm,Gkout)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: zeta    ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),intent(in)              :: Hk      ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    logical,intent(in)                                :: hk_symm !
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![2*Nlat*Nspin*Norb][2*Nlat*Nspin*Norb]
    integer                                           :: Nlat,Nspin,Norb,Nso,Nlso,Lfreq
    integer                                           :: i,is,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nlat  = size(zeta,3)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[2,2,Nlat,Nso,Nso,Lfreq],"invert_gk_superc_lattice_mpi","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_superc_lattice_mpi","Hk")
    call assert_shape(Gkout,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_lattice_mpi","Gkout")
    !
    allocate(Gktmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nlso,2*Nlso))
    Gkout = zero
    Gktmp = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta(1,1,:,:,:,i),Nlat,Nso) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(1,2,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta(2,1,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(2,2,:,:,:,i),Nlat,Nso) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                      Gktmp(1,ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                      Gktmp(2,ilat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    Gkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_superc_lattice_mpi
#endif




  subroutine invert_gk_superc_gij(zeta,Hk,hk_symm,Gkout,Fkout)
    complex(8)                                        :: zeta(:,:,:,:,:,:)              ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)                                        :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                                           :: hk_symm                
    !output:
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Fkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![2*Nlat*Nspin*Norb][2*Nlat*Nspin*Norb]
    integer                                           :: Lfreq
    !
    Nlat  = size(zeta,3)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[2,2,Nlat,Nso,Nso,Lfreq],"invert_gk_superc_gij","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_superc_gij","Hk")
    call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_gij","Gkout")
    call assert_shape(Fkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_gij","Fkout")
    !
    allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Fktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nlso,2*Nlso))
    Gkout = zero
    Fkout = zero
    Gktmp = zero
    Fktmp = zero
    do i=1,Lfreq
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta(1,1,:,:,:,i),Nlat,Nso) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(1,2,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta(2,1,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(2,2,:,:,:,i),Nlat,Nso) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                         Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                         Fktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    enddo
    Gkout = Gktmp
    Fkout = Fktmp
  end subroutine invert_gk_superc_gij

#ifdef _MPI
  subroutine invert_gk_superc_gij_mpi(MpiComm,zeta,Hk,hk_symm,Gkout,Fkout)
    integer                                           :: MpiComm
    complex(8)                                        :: zeta(:,:,:,:,:,:)              ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8)                                        :: Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb) 
    logical                                           :: hk_symm                
    !output:
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Fkout   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fktmp   ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:),allocatable             :: Gmatrix ![2*Nlat*Nspin*Norb][2*Nlat*Nspin*Norb]
    integer                                           :: Lfreq
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    Nlat  = size(zeta,3)
    Nspin = size(Gkout,3)
    Norb  = size(Gkout,5)
    Lfreq = size(zeta,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[2,2,Nlat,Nso,Nso,Lfreq],"invert_gk_superc_gij","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gk_superc_gij","Hk")
    call assert_shape(Gkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_gij","Gkout")
    call assert_shape(Fkout,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"invert_gk_superc_gij","Fkout")
    !
    allocate(Gktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Fktmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nlso,2*Nlso))
    Gktmp = zero
    Fktmp = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gmatrix  = zero
       Gmatrix(1:Nlso,1:Nlso)        = blocks_to_matrix(zeta(1,1,:,:,:,i),Nlat,Nso) - Hk
       Gmatrix(1:Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(1,2,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,1:Nlso) = blocks_to_matrix(zeta(2,1,:,:,:,i),Nlat,Nso)
       Gmatrix(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = blocks_to_matrix(zeta(2,2,:,:,:,i),Nlat,Nso) + Hk
       if(hk_symm) then
          call inv_sym(Gmatrix)
       else
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
       end if
       !store the diagonal blocks directly into the tmp output 
       do ilat=1,Nlat
          do jlat=1,Nlat
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                         Gktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,jo)
                         Fktmp(ilat,jlat,ispin,jspin,iorb,jorb,i) = Gmatrix(io,Nlso+jo)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    enddo
    Gkout=zero
    Fkout=zero
    call MPI_AllReduce(Gktmp, Gkout, size(Gkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
    call MPI_AllReduce(Fktmp, Fkout, size(Fkout), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  end subroutine invert_gk_superc_gij_mpi
#endif







  !--------------------------------------------------------------------!
  ! PURPOSE: Set the lattice/local Gamma functions for Embedding
  !--------------------------------------------------------------------!
  subroutine dmft_set_Gamma_matsubara_lattice(Gamma_mats)
    complex(8),dimension(:,:,:) :: Gamma_mats![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    integer                     :: Nlso,Lmats
    Nlso  = size(Gamma_mats,1)
    Lmats = size(Gamma_mats,3)
    !Testing part:
    call assert_shape(Gamma_mats,[Nlso,Nlso,Lmats],"dmft_set_Gamma_matsubara_lattice","Gamma_mats")
    !
    if(allocated(lattice_Gamma_mats))deallocate(lattice_Gamma_mats)
    allocate(lattice_Gamma_mats(Nlso,Nlso,Lmats))
    !
    write(*,"(A)")"Set lattice Gamma_mats"
    lattice_Gamma_mats = Gamma_mats
  end subroutine dmft_set_Gamma_matsubara_lattice

  subroutine dmft_set_Gamma_matsubara_local(Gamma_mats)
    complex(8),dimension(:,:,:,:),optional          :: Gamma_mats![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    integer                                         :: Nlat,Nso,Lmats
    !
    Nlat  = size(Gamma_mats,1)
    Nso   = size(Gamma_mats,2)
    Lmats = size(Gamma_mats,4)
    !Testing part:
    call assert_shape(Gamma_mats,[Nlat,Nso,Nso,Lmats],"dmft_set_Gamma_matsubara_local","Gamma_mats")         
    if(allocated(local_Gamma_mats))deallocate(local_Gamma_mats)
    allocate(local_Gamma_mats(Nlat,Nso,Nso,Lmats))
    !
    write(*,"(A)")"Set local Gamma_mats"
    local_Gamma_mats = Gamma_mats
  end subroutine dmft_set_Gamma_matsubara_local

  subroutine dmft_set_Gamma_realaxis_lattice(Gamma_real)
    complex(8),dimension(:,:,:) :: Gamma_real![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
    integer                     :: Nlso,Lreal
    Nlso  = size(Gamma_real,1)
    Lreal = size(Gamma_real,3)
    !Testing part:
    call assert_shape(Gamma_real,[Nlso,Nlso,Lreal],"dmft_set_Gamma_realaxis_lattice","Gamma_real")
    !
    if(allocated(lattice_Gamma_real))deallocate(lattice_Gamma_real)
    allocate(lattice_Gamma_real(Nlso,Nlso,Lreal))
    !
    write(*,"(A)")"Set lattice Gamma_real"
    lattice_Gamma_real = Gamma_real
  end subroutine dmft_set_Gamma_realaxis_lattice

  subroutine dmft_set_Gamma_realaxis_local(Gamma_real)
    complex(8),dimension(:,:,:,:),optional          :: Gamma_real![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    integer                                         :: Nlat,Nso,Lreal
    !
    Nlat  = size(Gamma_real,1)
    Nso   = size(Gamma_real,2)
    Lreal = size(Gamma_real,4)
    !Testing part:
    call assert_shape(Gamma_real,[Nlat,Nso,Nso,Lreal],"dmft_set_Gamma_realaxis_local","Gamma_real")         
    if(allocated(local_Gamma_real))deallocate(local_Gamma_real)
    allocate(local_Gamma_real(Nlat,Nso,Nso,Lreal))
    !
    write(*,"(A)")"Set local Gamma_real"
    local_Gamma_real = Gamma_real
  end subroutine dmft_set_Gamma_realaxis_local







  !--------------------------------------------------------------------!
  ! PURPOSE: Print Glocal and Gij
  !--------------------------------------------------------------------!
  !SET THE GLOC SUFFIX
  subroutine dmft_gloc_set_suffix(string)
    character(len=*) :: string
    gloc_suffix=reg(string)
  end subroutine dmft_gloc_set_suffix




  !PRINT GLOC SINGLE SITE
  subroutine dmft_gloc_print_matsubara(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats
    real(8),dimension(size(Gmats,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    Nspin = size(Gmats,1)
    Norb  = size(Gmats,3)
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,size(Gmats,5)],"dmft_gloc_print_matsubara",reg(fname)//"_mats")
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Gloc matsubara: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,"(A)") "Gloc matsubara: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gloc_suffix)
             call splot(reg(suffix),w,Gmats(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,"(A)") "Gloc matsubara: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gloc_suffix)
                call splot(reg(suffix),w,Gmats(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default                  !print all off-diagonals
       write(*,"(A)") "Gloc matsubara: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gloc_suffix)
                   call splot(reg(suffix),w,Gmats(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gloc_print_matsubara

  subroutine dmft_gloc_print_realaxis(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Greal
    real(8),dimension(size(Greal,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    Nspin = size(Greal,1)
    Norb  = size(Greal,3)
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,size(Greal,5)],"dmft_gloc_print_realaxis",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Gloc real: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,"(A)") "Gloc real: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gloc_suffix)
             call splot(reg(suffix),w,Greal(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,"(A)") "Gloc real: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gloc_suffix)
                call splot(reg(suffix),w,Greal(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default                  !print all off-diagonals
       write(*,"(A)") "Gloc real: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gloc_suffix)
                   call splot(reg(suffix),w,Greal(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gloc_print_realaxis



  !PRINT GLOC LATTICE
  subroutine dmft_gloc_print_matsubara_lattice(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Gmats
    real(8),dimension(size(Gmats,6))             :: w        
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,2)
    Norb  = size(Gmats,4)
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,size(Gmats,6)],"dmft_gloc_print_matsubara_lattice",reg(fname)//"_mats")
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Gloc matsubara: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,*)"Gloc matsubara: write spin-orbital diagonal elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gloc_suffix)
             call store_data(reg(suffix),Gmats(:,ispin,ispin,iorb,iorb,:),w)
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc matsubara: write spin diagonal and all orbitals elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gloc_suffix)
                call store_data(reg(suffix),Gmats(:,ispin,ispin,iorb,jorb,:),w)
             enddo
          enddo
       enddo
       !
    case(3)                  !print all off-diagonals
       write(*,*)"Gloc matsubara: write all elements. No Split."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gloc_suffix)
                   call store_data(reg(suffix),Gmats(:,ispin,jspin,iorb,jorb,:),w)
                enddo
             enddo
          enddo
       enddo
       !
    case(4)                  !print only diagonal elements
       write(*,*)"Gloc matsubara: write spin-orbital diagonal elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw_indx"//reg(txtfy(ilat,6))//reg(gloc_suffix)
                call splot(reg(suffix),w,Gmats(ilat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
       !
    case(5)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc matsubara: write spin diagonal and all orbitals elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//&
                        "_iw_indx"//reg(txtfy(ilat,6))//reg(gloc_suffix)
                   call splot(reg(suffix),w,Gmats(ilat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    case default
       write(*,*)"Gloc matsubara: write all elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                           "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                           "_iw_indx"//reg(txtfy(ilat,6))//reg(gloc_suffix)
                      call splot(reg(suffix),w,Gmats(ilat,ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gloc_print_matsubara_lattice

  subroutine dmft_gloc_print_realaxis_lattice(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Greal
    real(8),dimension(size(Greal,6))             :: w
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,2)
    Norb  = size(Greal,4)
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,size(Greal,6)],"dmft_gloc_print_realaxis_lattice",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Gloc real: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,*)"Gloc real: write spin-orbital diagonal elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gloc_suffix)
             call store_data(reg(suffix),Greal(:,ispin,ispin,iorb,iorb,:),w)
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc real: write spin diagonal and all orbitals elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gloc_suffix)
                call store_data(reg(suffix),Greal(:,ispin,ispin,iorb,jorb,:),w)
             enddo
          enddo
       enddo
       !
    case(3)
       write(*,*)"Gloc real: write all elements. No Split."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gloc_suffix)
                   call store_data(reg(suffix),Greal(:,ispin,jspin,iorb,jorb,:),w)
                enddo
             enddo
          enddo
       enddo
       !
    case(4)                  !print only diagonal elements
       write(*,*)"Gloc real: write spin-orbital diagonal elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw_indx"//reg(txtfy(ilat,6))//reg(gloc_suffix)
                call splot(reg(suffix),w,Greal(ilat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
       !
    case(5)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc real: write spin diagonal and all orbitals elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//&
                        "_realw_indx"//reg(txtfy(ilat,6))//reg(gloc_suffix)
                   call splot(reg(suffix),w,Greal(ilat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    case default
       write(*,*)"Gloc real: write all elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                           "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                           "_realw_indx"//reg(txtfy(ilat,6))//reg(gloc_suffix)
                      call splot(reg(suffix),w,Greal(ilat,ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gloc_print_realaxis_lattice


  !PRINT GLOC FULL LATTICE
  subroutine dmft_gloc_print_matsubara_gij(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Gmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats/Lreal]
    real(8),dimension(size(Gmats,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,3)
    Norb  = size(Gmats,5)
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Gmats,7)],"dmft_gloc_print_matsubara_gij",reg(fname)//"_mats")
    !
    select case(iprint)
    case (0)
       write(*,*)"Gloc matsubara: not written on file."
    case(1)                  !print only diagonal elements     
       write(*,*)"Gloc matsubara: spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gloc_suffix)
             call store_data(reg(suffix),Gmats(:,:,ispin,ispin,iorb,iorb,:),w)
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc matsubara: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gloc_suffix)
                call store_data(reg(suffix),Gmats(:,:,ispin,ispin,iorb,jorb,:),w)
             enddo
          enddo
       enddo
    case default              !print all off-diagonals
       write(*,*)"Gloc matsubara: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gloc_suffix)
                   call store_data(reg(suffix),Gmats(:,:,ispin,jspin,iorb,jorb,:),w)
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine dmft_gloc_print_matsubara_gij

  subroutine dmft_gloc_print_realaxis_gij(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Greal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal/Lreal]
    real(8),dimension(size(Greal,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    integer                                        :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb,ilat,jlat
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,3)
    Norb  = size(Greal,5)
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Greal,7)],"dmft_gloc_print_realaxis_gij",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(*,*)"Gloc real: not written on file."
    case(1)                  !print only diagonal elements     
       write(*,*)"Gloc real: spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gloc_suffix)
             call store_data(reg(suffix),Greal(:,:,ispin,ispin,iorb,iorb,:),w)
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc real: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gloc_suffix)
                call store_data(reg(suffix),Greal(:,:,ispin,ispin,iorb,jorb,:),w)
             enddo
          enddo
       enddo
    case default              !print all off-diagonals
       write(*,*)"Gloc real: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gloc_suffix)
                   call store_data(reg(suffix),Greal(:,:,ispin,jspin,iorb,jorb,:),w)
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine dmft_gloc_print_realaxis_gij



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





  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
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








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso


end module DMFT_GLOC
