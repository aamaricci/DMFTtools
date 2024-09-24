MODULE DMFT_GF
  USE GF_COMMON!, only:gf_push_zeta
  USE GF_IO
  USE GF_GLOC
  USE GF_GK

  implicit none
  private

  interface dmft_get_gloc
     !push zeta array
     module procedure :: gf_push_zeta
     ! NORMAL 
     ! . rank-2   [N,N,:]
     module procedure :: get_gloc_normal_hk_rank2
     module procedure :: get_gloc_normal_tridiag_rank2
     module procedure :: get_gloc_normal_dos_rank2
     ! . rank-3   [Nlat,Nso,Nso,:]
     module procedure :: get_gloc_normal_hk_rank3
     ! . rank-4   [Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_normal_hk_rank4
     module procedure :: get_gloc_normal_dos_rank4
     ! . rank-5   [Nlat,Nspin,Nspin,Norb,Norb,:]
     ! . rank-5_6 [Nlat,Nspin,Nspin,Norb,Norb,:] ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_normal_hk_rank5
     module procedure :: get_gloc_normal_tridiag_rank5
     module procedure :: get_gloc_normal_dos_rank5
     module procedure :: get_gloc_normal_hk_rank5_6
     ! . rank-6   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-7   [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
     module procedure :: get_gloc_normal_hk_rank7
     module procedure :: get_gloc_normal_tridiag_rank7
#endif
     !>SUPERC
     ! . rank-2   [2,N,N,:]
     module procedure :: get_gloc_superc_hk_rank2
     module procedure :: get_gloc_superc_dos_rank2
     ! . rank-3   [2,Nlat,Nso,Nso,:]
     module procedure :: get_gloc_superc_hk_rank3
     ! . rank-4   [2,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_superc_hk_rank4
     module procedure :: get_gloc_superc_dos_rank4
     ! . rank-5   [2,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_superc_hk_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-5_6 [2,Nlat,Nspin,Nspin,Norb,Norb,:]->[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_superc_hk_rank5_6
#endif
  end interface dmft_get_gloc


  interface dmft_get_gk
     !push zeta array
     module procedure :: gf_push_zeta
     ! NORMAL 
     ! . rank-2   [N,N,:]
     module procedure :: get_gk_normal_hk_rank2
     module procedure :: get_gk_normal_tridiag_rank2
     module procedure :: get_gk_normal_dos_rank2
     ! . rank-3   [Nlat,Nso,Nso,:]
     module procedure :: get_gk_normal_hk_rank3
     ! . rank-4   [Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_normal_hk_rank4
     module procedure :: get_gk_normal_dos_rank4
     ! . rank-5   [Nlat,Nspin,Nspin,Norb,Norb,:]
     ! . rank-5_6 [Nlat,Nspin,Nspin,Norb,Norb,:] ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_normal_hk_rank5
     module procedure :: get_gk_normal_tridiag_rank5
     module procedure :: get_gk_normal_dos_rank5
     module procedure :: get_gk_normal_hk_rank5_6
     ! . rank-6   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-7   [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
     module procedure :: get_gk_normal_hk_rank7
     module procedure :: get_gk_normal_tridiag_rank7
#endif
     !>SUPERC
     ! . rank-2   [2,N,N,:]
     module procedure :: get_gk_superc_hk_rank2
     module procedure :: get_gk_superc_dos_rank2
     ! . rank-3   [2,Nlat,Nso,Nso,:]
     module procedure :: get_gk_superc_hk_rank3
     ! . rank-4   [2,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_superc_hk_rank4
     module procedure :: get_gk_superc_dos_rank4
     ! . rank-5   [2,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_superc_hk_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-5_6 [2,Nlat,Nspin,Nspin,Norb,Norb,:]->[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_superc_hk_rank5_6
#endif
  end interface dmft_get_gk




  interface get_gloc
     !push zeta array
     module procedure :: gf_push_zeta
     ! NORMAL 
     ! . rank-2   [N,N,:]
     module procedure :: get_gloc_normal_hk_rank2
     module procedure :: get_gloc_normal_tridiag_rank2
     module procedure :: get_gloc_normal_dos_rank2
     ! . rank-3   [Nlat,Nso,Nso,:]
     module procedure :: get_gloc_normal_hk_rank3
     ! . rank-4   [Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_normal_hk_rank4
     module procedure :: get_gloc_normal_dos_rank4
     ! . rank-5   [Nlat,Nspin,Nspin,Norb,Norb,:]
     ! . rank-5_6 [Nlat,Nspin,Nspin,Norb,Norb,:] ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_normal_hk_rank5
     module procedure :: get_gloc_normal_tridiag_rank5
     module procedure :: get_gloc_normal_dos_rank5
     module procedure :: get_gloc_normal_hk_rank5_6
     ! . rank-6   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-7   [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
     module procedure :: get_gloc_normal_hk_rank7
     module procedure :: get_gloc_normal_tridiag_rank7
#endif
     !>SUPERC
     ! . rank-2   [2,N,N,:]
     module procedure :: get_gloc_superc_hk_rank2
     module procedure :: get_gloc_superc_dos_rank2
     ! . rank-3   [2,Nlat,Nso,Nso,:]
     module procedure :: get_gloc_superc_hk_rank3
     ! . rank-4   [2,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_superc_hk_rank4
     module procedure :: get_gloc_superc_dos_rank4
     ! . rank-5   [2,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_superc_hk_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-5_6 [2,Nlat,Nspin,Nspin,Norb,Norb,:]->[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gloc_superc_hk_rank5_6
#endif
  end interface get_gloc


  interface get_gk
     !push zeta array
     module procedure :: gf_push_zeta
     ! NORMAL 
     ! . rank-2   [N,N,:]
     module procedure :: get_gk_normal_hk_rank2
     module procedure :: get_gk_normal_tridiag_rank2
     module procedure :: get_gk_normal_dos_rank2
     ! . rank-3   [Nlat,Nso,Nso,:]
     module procedure :: get_gk_normal_hk_rank3
     ! . rank-4   [Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_normal_hk_rank4
     module procedure :: get_gk_normal_dos_rank4
     ! . rank-5   [Nlat,Nspin,Nspin,Norb,Norb,:]
     ! . rank-5_6 [Nlat,Nspin,Nspin,Norb,Norb,:] ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_normal_hk_rank5
     module procedure :: get_gk_normal_tridiag_rank5
     module procedure :: get_gk_normal_dos_rank5
     module procedure :: get_gk_normal_hk_rank5_6
     ! . rank-6   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-7   [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
     module procedure :: get_gk_normal_hk_rank7
     module procedure :: get_gk_normal_tridiag_rank7
#endif
     !>SUPERC
     ! . rank-2   [2,N,N,:]
     module procedure :: get_gk_superc_hk_rank2
     module procedure :: get_gk_superc_dos_rank2
     ! . rank-3   [2,Nlat,Nso,Nso,:]
     module procedure :: get_gk_superc_hk_rank3
     ! . rank-4   [2,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_superc_hk_rank4
     module procedure :: get_gk_superc_dos_rank4
     ! . rank-5   [2,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_superc_hk_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
     ! . rank-5_6 [2,Nlat,Nspin,Nspin,Norb,Norb,:]->[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
     module procedure :: get_gk_superc_hk_rank5_6
#endif
  end interface get_gk



  !GLOC:
  public :: dmft_get_gloc
  public :: get_gloc
  public :: dmft_get_gk
  public :: get_gk
  public :: dmft_write_gf
  public :: write_gf

  public :: reshape_matrix_to_rank3
  public :: reshape_matrix_to_rank4
  public :: reshape_matrix_to_rank5
  public :: reshape_matrix_to_rank6
  public :: reshape_matrix_to_rank7

  public :: reshape_rank3_to_matrix
  public :: reshape_rank4_to_matrix
  public :: reshape_rank5_to_matrix
  public :: reshape_rank6_to_matrix
  public :: reshape_rank7_to_matrix

END MODULE DMFT_GF
