module DMFT_TIGHT_BINDING
  USE TB_COMMON
  USE TB_BASIS
  USE TB_IO
  USE TB_WANNIER90
  USE TB_BUILD
  USE TB_SOLVE
  USE TB_FSURF
  implicit none
  private

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
  public :: TB_build_kgrid
  public :: TB_build_Rgrid
  public :: TB_refine_kgrid
  !
  public :: TB_bk_length
  public :: TB_ei_length
  !
  public :: TB_write_grid
  !
  public :: TB_build_model
  public :: TB_solve_model
  !
  public :: TB_get_FermiLevel
  !
  public :: TB_fsurface
  !
  public :: TB_w90_setup
  public :: TB_w90_delete
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

  public :: TB_set_dos_range
  public :: TB_set_dos_lreal
  public :: TB_set_dos_eps
  public :: TB_set_dos_file

  public :: kpoint_gamma
  public :: kpoint_x1
  public :: kpoint_x2
  public :: kpoint_x3
  public :: kpoint_m1
  public :: kpoint_m2
  public :: kpoint_m3
  public :: kpoint_r

END MODULE DMFT_TIGHT_BINDING











