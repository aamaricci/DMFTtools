module DMFT_GF
  USE GF_COMMON
  USE GF_GLOC
  USE GF_GK
  implicit none
  private


  !GK:
  public :: get_gk
  public :: get_gk_matsubara
  public :: get_gk_realaxis

  public :: dmft_gk
  public :: dmft_gk_matsubara
  public :: dmft_gk_realaxis

  public :: dmft_get_gk
  public :: dmft_get_gk_matsubara
  public :: dmft_get_gk_realaxis


  !GLOC
  public :: get_gloc
  public :: get_gloc_matsubara
  public :: get_gloc_realaxis

  public :: dmft_gloc
  public :: dmft_gloc_matsubara
  public :: dmft_gloc_realaxis

  public :: dmft_get_gloc
  public :: dmft_get_gloc_matsubara
  public :: dmft_get_gloc_realaxis


  !Gij
  public :: get_gij
  public :: get_gij_matsubara
  public :: get_gij_realaxis

  public :: dmft_gij
  public :: dmft_gij_matsubara
  public :: dmft_gij_realaxis

  public :: dmft_get_gij
  public :: dmft_get_gij_matsubara
  public :: dmft_get_gij_realaxis


  public :: dmft_set_Gamma_matsubara
  public :: dmft_set_Gamma_realaxis

end module DMFT_GF
