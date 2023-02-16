module DMFT_GF
  USE GF_COMMON, only:gf_push_zeta
  USE GF_GLOC
  USE GF_GK
  implicit none
  private

  interface dmft_get_gloc
     !push zeta array
     module procedure :: gf_push_zeta
     !main procedures normal
     module procedure :: get_gloc_normal_main
     module procedure :: get_gloc_normal_tridiag
     module procedure :: get_gloc_normal_dos
     !interfaces normal
     module procedure :: get_gloc_normal_hk_rank4
     module procedure :: get_gloc_normal_dos_rank4
     module procedure :: get_gloc_normal_hk_rank5
     module procedure :: get_gloc_normal_tridiag_rank5
     module procedure :: get_gloc_normal_hk_rank5_6
     module procedure :: get_gloc_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: get_gloc_normal_hk_rank7
     module procedure :: get_gloc_normal_tridiag_rank7
#endif
     !main procedures superc
     module procedure :: get_gloc_superc_main
     module procedure :: get_gloc_superc_dos
     !interfaces superc
     module procedure :: get_gloc_superc_hk_rank4
     module procedure :: get_gloc_superc_dos_rank4
     module procedure :: get_gloc_superc_hk_rank5
     module procedure :: get_gloc_superc_hk_rank5_6
  end interface dmft_get_gloc


  interface dmft_get_gk
     !push zeta array
     module procedure :: gf_push_zeta
     !main procedures normal
     module procedure :: get_gk_normal_main
     module procedure :: get_gk_normal_tridiag
     module procedure :: get_gk_normal_dos
     !interfaces normal
     module procedure :: get_gk_normal_hk_rank4
     module procedure :: get_gk_normal_dos_rank4
     module procedure :: get_gk_normal_hk_rank5
     module procedure :: get_gk_normal_tridiag_rank5
     module procedure :: get_gk_normal_hk_rank5_6
     module procedure :: get_gk_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: get_gk_normal_hk_rank7
     module procedure :: get_gk_normal_tridiag_rank7
#endif
     !main procedures superc
     module procedure :: get_gk_superc_main
     module procedure :: get_gk_superc_dos
     !interfaces superc
     module procedure :: get_gk_superc_hk_rank4
     module procedure :: get_gk_superc_dos_rank4
     module procedure :: get_gk_superc_hk_rank5
     module procedure :: get_gk_superc_hk_rank5_6
  end interface dmft_get_gk




  interface get_gloc
     !push zeta array
     module procedure :: gf_push_zeta
     !main procedures normal
     module procedure :: get_gloc_normal_main
     module procedure :: get_gloc_normal_tridiag
     module procedure :: get_gloc_normal_dos
     !interfaces normal
     module procedure :: get_gloc_normal_hk_rank4
     module procedure :: get_gloc_normal_dos_rank4
     module procedure :: get_gloc_normal_hk_rank5
     module procedure :: get_gloc_normal_tridiag_rank5
     module procedure :: get_gloc_normal_hk_rank5_6
     module procedure :: get_gloc_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: get_gloc_normal_hk_rank7
     module procedure :: get_gloc_normal_tridiag_rank7
#endif
     !main procedures superc
     module procedure :: get_gloc_superc_main
     module procedure :: get_gloc_superc_dos
     !interfaces superc
     module procedure :: get_gloc_superc_hk_rank4
     module procedure :: get_gloc_superc_dos_rank4
     module procedure :: get_gloc_superc_hk_rank5
     module procedure :: get_gloc_superc_hk_rank5_6
  end interface get_gloc


  interface get_gk
     !push zeta array
     module procedure :: gf_push_zeta
     !main procedures normal
     module procedure :: get_gk_normal_main
     module procedure :: get_gk_normal_tridiag
     module procedure :: get_gk_normal_dos
     !interfaces normal
     module procedure :: get_gk_normal_hk_rank4
     module procedure :: get_gk_normal_dos_rank4
     module procedure :: get_gk_normal_hk_rank5
     module procedure :: get_gk_normal_tridiag_rank5
     module procedure :: get_gk_normal_hk_rank5_6
     module procedure :: get_gk_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: get_gk_normal_hk_rank7
     module procedure :: get_gk_normal_tridiag_rank7
#endif
     !main procedures superc
     module procedure :: get_gk_superc_main
     module procedure :: get_gk_superc_dos
     !interfaces superc
     module procedure :: get_gk_superc_hk_rank4
     module procedure :: get_gk_superc_dos_rank4
     module procedure :: get_gk_superc_hk_rank5
     module procedure :: get_gk_superc_hk_rank5_6
  end interface get_gk



  !GLOC:
  public :: dmft_get_gloc
  public :: get_gloc



  !GK:
  public :: dmft_get_gk
  public :: get_gk



end module DMFT_GF
