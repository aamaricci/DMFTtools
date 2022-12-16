module LEGACY_DMFT_WEISS_FIELD
  USE LEGACY_SC_COMMON
  USE LEGACY_SC_WEISS, legacy_dmft_weiss => dmft_weiss 
  USE LEGACY_SC_DELTA, legacy_dmft_delta => dmft_delta
  USE LEGACY_SC_GLOBAL, legacy_dmft_self_consistency => dmft_self_consistency
  implicit none
  private

  public :: legacy_dmft_self_consistency
  public :: legacy_dmft_weiss
  public :: legacy_dmft_delta

end module LEGACY_DMFT_WEISS_FIELD
