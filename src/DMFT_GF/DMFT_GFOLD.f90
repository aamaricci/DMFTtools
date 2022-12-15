module DMFT_GFOLD
  ! USE GF_COMMON_OLD
  USE GK_MATSUBARA
  USE GK_REALAXIS
  USE GLOC_MATSUBARA
  USE GLOC_REALAXIS
  implicit none
  private

  public :: get_gk_matsubara
  public :: get_gk_realaxis
  public :: get_gloc_matsubara
  public :: get_gloc_realaxis

  public :: dmft_gk_matsubara
  public :: dmft_gk_realaxis
  public :: dmft_gloc_matsubara
  public :: dmft_gloc_realaxis


  public :: get_gij_matsubara
  public :: get_gij_realaxis
  public :: dmft_gij_matsubara
  public :: dmft_gij_realaxis

  ! public :: dmft_set_Gamma_matsubara
  ! public :: dmft_set_Gamma_realaxis

end module DMFT_GFOLD
