module DMFT_GK
  USE GK_MATSUBARA
  USE GK_REALAXIS
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private

  public :: dmft_gk_matsubara
  public :: dmft_gk_realaxis

contains

end module DMFT_GK
