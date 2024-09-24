module ARCHIVED_DMFT_GF
  USE GK_MATSUBARA, archived_dmft_gk_matsubara => dmft_gk_matsubara
  USE GK_REALAXIS,  archived_dmft_gk_realaxis => dmft_gk_realaxis
  USE GLOC_MATSUBARA, archived_dmft_gloc_matsubara => dmft_gloc_matsubara
  USE GLOC_REALAXIS, archived_dmft_gloc_realaxis => dmft_gloc_realaxis
  implicit none
  private

  public :: archived_dmft_gk_matsubara
  public :: archived_dmft_gk_realaxis
  public :: archived_dmft_gloc_matsubara
  public :: archived_dmft_gloc_realaxis
  ! public :: archived_dmft_gij_matsubara
  ! public :: archived_dmft_gij_realaxis

end module ARCHIVED_DMFT_GF
