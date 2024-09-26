module LEGACY_DMFT_WEISS_FIELD
  USE LEGACY_SC_COMMON
  USE LEGACY_SC_WEISS
  USE LEGACY_SC_DELTA
  implicit none
  private


  interface legacy_dmft_weiss
     module procedure :: dmft_get_weiss_normal_main
     module procedure :: dmft_get_weiss_normal_rank4
     module procedure :: dmft_get_weiss_normal_rank5
     module procedure :: dmft_get_weiss_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_weiss_normal_rank7
#endif
     !
     module procedure :: dmft_get_weiss_superc_main
     module procedure :: dmft_get_weiss_superc_rank4
     module procedure :: dmft_get_weiss_superc_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_weiss_superc_rank6
#endif
  end interface legacy_dmft_weiss


  interface legacy_dmft_delta
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_rank4
     module procedure :: dmft_get_delta_normal_rank5
     module procedure :: dmft_get_delta_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_delta_normal_rank7
#endif
     !
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_rank4
     module procedure :: dmft_get_delta_superc_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_delta_superc_rank6
#endif
  end interface legacy_dmft_delta



  interface legacy_dmft_self_consistency
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_rank4
     module procedure :: dmft_get_delta_normal_rank5
     module procedure :: dmft_get_delta_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_delta_normal_rank7
#endif
     !
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_rank4
     module procedure :: dmft_get_delta_superc_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_delta_superc_rank6
#endif
     !
     module procedure :: dmft_get_weiss_normal_main
     module procedure :: dmft_get_weiss_normal_rank4
     module procedure :: dmft_get_weiss_normal_rank5
     module procedure :: dmft_get_weiss_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_weiss_normal_rank7
#endif
     !
     module procedure :: dmft_get_weiss_superc_main
     module procedure :: dmft_get_weiss_superc_rank4
     module procedure :: dmft_get_weiss_superc_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_weiss_superc_rank6
#endif
  end interface LEGACY_dmft_self_consistency



  public :: legacy_dmft_weiss
  public :: legacy_dmft_delta
  public :: legacy_dmft_self_consistency

end module LEGACY_DMFT_WEISS_FIELD
