module ARCHIVED_DMFT_WEISS_FIELD
  USE ARCHIVED_SC_COMMON
  USE ARCHIVED_SC_WEISS, archived_dmft_weiss => dmft_weiss 
  USE ARCHIVED_SC_DELTA, archived_dmft_delta => dmft_delta
  USE ARCHIVED_SC_GLOBAL, archived_dmft_self_consistency => dmft_self_consistency
  implicit none
  private

  public :: archived_dmft_self_consistency
  public :: archived_dmft_weiss
  public :: archived_dmft_delta

end module ARCHIVED_DMFT_WEISS_FIELD
