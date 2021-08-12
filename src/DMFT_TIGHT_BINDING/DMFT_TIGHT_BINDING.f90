module DMFT_TIGHT_BINDING
  USE TB_COMMON
  USE TB_BASIS
  USE TB_IO
  USE TB_EFERMI
  USE TB_WANNIER90
  USE TB_BUILD
  USE TB_SOLVE
  USE TB_FSURF
  implicit none
  private


  !TB_BASIS
  interface TB_build_kgrid
     module procedure :: build_kgrid
     module procedure :: build_kgrid_generic
     module procedure :: kgrid_from_path_grid
     module procedure :: kgrid_from_path_dim
     module procedure :: klen_from_path
  end interface TB_build_kgrid
  interface TB_kgrid
     module procedure :: build_kgrid
     module procedure :: build_kgrid_generic
     module procedure :: kgrid_from_path_grid
     module procedure :: kgrid_from_path_dim
     module procedure :: klen_from_path
  end interface TB_kgrid

  interface TB_build_Rgrid
     module procedure ::   build_Rgrid
  end interface TB_build_Rgrid
  interface TB_Rgrid
     module procedure ::   build_Rgrid
  end interface TB_Rgrid


  interface TB_print_bk
     module procedure :: print_bk
  end interface TB_print_bk

  interface TB_print_ei
     module procedure :: print_ei
  end interface TB_print_ei



  !TB_IO
  interface TB_write_hk
     module procedure :: write_hk_func
     module procedure :: write_hk_array
     module procedure :: write_hk_w90_func
     module procedure :: write_hk_w90_array
  end interface TB_write_hk

  interface TB_read_hk
     module procedure :: read_hk_array
     module procedure :: read_hk_w90_array
  end interface TB_read_hk

  interface TB_write_Hloc
     module procedure :: write_Hloc_1
     module procedure :: write_Hloc_2
  end interface TB_write_Hloc

  interface TB_read_Hloc
     module procedure :: read_Hloc_1
     module procedure :: read_Hloc_2
  end interface TB_read_Hloc



  !TB_BUILD
  interface TB_build_model
     module procedure :: build_hk_model_kgrid
     module procedure :: build_hk_model_nkvec
     module procedure :: build_hkR_model_kgrid
     module procedure :: build_hkR_model_nkvec
     module procedure :: build_hk_path
     module procedure :: build_hkR_path
     module procedure :: build_hk_w90
     !
     module procedure :: build_Hij_Nrvec
  end interface TB_build_model


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


  !TB_SOLVE
  interface TB_solve_model
     module procedure :: solve_Hk_along_BZpath
     module procedure :: solve_w90Hk_along_BZpath !this is redundant
     module procedure :: solve_HkR_along_BZpath
     !< obsolete
     module procedure :: read_Hr_w90_solve_Hk_along_BZpath
  end interface TB_solve_model



  !TB_WANNIER
  interface TB_w90_setup
     module procedure :: setup_w90
  end interface TB_w90_setup

  interface TB_fix_w90file
     module procedure :: fix_w90file
  end interface TB_fix_w90file

  interface TB_w90_delete
     module procedure :: delete_w90
  end interface TB_w90_delete

  interface TB_w90_Hloc
     module procedure :: Hloc_w90
  end interface TB_w90_Hloc

  interface TB_w90_Zeta
     module procedure :: Zeta_w90_matrix
     module procedure :: Zeta_w90_vector
  end interface TB_w90_Zeta

  interface TB_w90_Self
     module procedure :: Self0_w90
  end interface TB_w90_Self

  interface TB_w90_FermiLevel
     module procedure :: FermiLevel_w90
  end interface TB_w90_FermiLevel


  !TB_FSURF
  interface TB_fsurface
     module procedure :: TB_fsurf_nkvec
     module procedure :: TB_fsurf_w90_nkvec
  end interface TB_fsurface



  abstract interface
     function w90_hk(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: w90_hk
     end function w90_hk
  end interface
  procedure(w90_hk),pointer      :: TB_w90_model=>w90_hk_model



  public :: w90_structure

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
  public :: TB_build_Kgrid
  public :: TB_build_Rgrid
  public :: TB_Kgrid
  public :: TB_Rgrid
  !
  public :: TB_bk_length
  public :: TB_ei_length
  !
  public :: TB_write_grid
  !
  public :: TB_build_model
  public :: TB_solve_model
  !
  public :: TB_FermiLevel
  !
  public :: TB_fsurface
  !
  public :: TB_w90_setup
  public :: TB_w90_delete
  public :: TB_fix_w90file
  public :: TB_w90_Hloc
  public :: TB_w90_FermiLevel
  public :: TB_w90_Zeta
  public :: TB_w90_Self
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

  public :: TB_slo2lso_model
  public :: TB_reshape_array
  public :: TB_reorder_array
  public :: TB_reshape_hk
  public :: TB_reorder_hk
  !
  public :: kpoint_gamma
  public :: kpoint_x1
  public :: kpoint_x2
  public :: kpoint_x3
  public :: kpoint_m1
  public :: kpoint_m2
  public :: kpoint_m3
  public :: kpoint_r

END MODULE DMFT_TIGHT_BINDING











