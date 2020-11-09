module DMFT_GLOC
  USE GLOC_COMMON
  USE GLOC_MATSUBARA
  USE GLOC_REALAXIS
  implicit none
  private



  !PUBLIC IN DMFT:
  public :: get_gloc_matsubara
  public :: get_gloc_realaxis
  public :: get_gij_matsubara
  public :: get_gij_realaxis

  public :: dmft_gloc_matsubara
  public :: dmft_gloc_realaxis
  public :: dmft_gij_matsubara
  public :: dmft_gij_realaxis

  public :: dmft_set_Gamma_matsubara
  public :: dmft_set_Gamma_realaxis
end module DMFT_GLOC
