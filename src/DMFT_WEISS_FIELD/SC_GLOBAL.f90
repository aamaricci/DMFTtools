module SC_GLOBAL
  USE DMFT_CTRL_VARS
  USE SC_COMMON
  USE SC_WEISS
  USE SC_DELTA
  implicit none
  private


  interface dmft_weiss
     module procedure :: dmft_get_weiss_normal_main
     module procedure :: dmft_get_weiss_normal_rank4
     module procedure :: dmft_get_weiss_normal_rank5
     module procedure :: dmft_get_weiss_normal_rank6
     module procedure :: dmft_get_weiss_normal_rank7
     !
     module procedure :: dmft_get_weiss_superc_main
     module procedure :: dmft_get_weiss_superc_rank4
     module procedure :: dmft_get_weiss_superc_rank5
     module procedure :: dmft_get_weiss_superc_rank6
  end interface dmft_weiss


  interface dmft_delta
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_rank4
     module procedure :: dmft_get_delta_normal_rank5
     module procedure :: dmft_get_delta_normal_rank6
     module procedure :: dmft_get_delta_normal_rank7
     !
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_rank4
     module procedure :: dmft_get_delta_superc_rank5
     module procedure :: dmft_get_delta_superc_rank6
  end interface dmft_delta

  interface dmft_self_consistency
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_rank4
     module procedure :: dmft_get_delta_normal_rank5
     module procedure :: dmft_get_delta_normal_rank6
     module procedure :: dmft_get_delta_normal_rank7
     !
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_rank4
     module procedure :: dmft_get_delta_superc_rank5
     module procedure :: dmft_get_delta_superc_rank6
     !
     module procedure :: dmft_get_weiss_normal_main
     module procedure :: dmft_get_weiss_normal_rank4
     module procedure :: dmft_get_weiss_normal_rank5
     module procedure :: dmft_get_weiss_normal_rank6
     module procedure :: dmft_get_weiss_normal_rank7
     !
     module procedure :: dmft_get_weiss_superc_main
     module procedure :: dmft_get_weiss_superc_rank4
     module procedure :: dmft_get_weiss_superc_rank5
     module procedure :: dmft_get_weiss_superc_rank6
  end interface dmft_self_consistency

  public :: dmft_weiss
  public :: dmft_delta
  public :: dmft_self_consistency



end module SC_GLOBAL
