module LEGACY_DMFT_GF
  USE GK_MATSUBARA
  USE GK_REALAXIS
  USE GLOC_MATSUBARA
  USE GLOC_REALAXIS
  implicit none
  private

  public :: dmft_gk_matsubara
  public :: dmft_gk_realaxis
  public :: dmft_gloc_matsubara
  public :: dmft_gloc_realaxis


  public :: get_gij_matsubara
  public :: get_gij_realaxis
  public :: dmft_gij_matsubara
  public :: dmft_gij_realaxis

end module LEGACY_DMFT_GF
