module DMFT_WEISS_FIELD
  USE SC_COMMON, only: gf_push_zeta
  USE SC_WEISS
  USE SC_DELTA
  implicit none
  private

  !##################################################################
  !##################################################################
  !                       NORMAL PHASE
  ! Get the local Weiss Field calG0 using DMFT self-consistency 
  !  given G_loc/F_loc and Sigma/Self in different shapes:
  !
  ! . rank-2 DMFT    [N,N,:] check @1Gb memory per Function
  ! . rank-3 DMFT    [Nlat,Nso,Nso,:]
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:]
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5_6 R-DMFT[Nlat,Nspin,Nspin,Norb,Norb,:]
  !                ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-6 CDMFT   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  ! . rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  !##################################################################
  !##################################################################
  !                         SUPERC
  ! . rank-2 DMFT    [][N,N,:] check @1Gb memory per Function
  ! . rank-3 DMFT    [Nlat,Nso,Nso,:] 
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] 
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-6 CDMFT   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  ! . rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  !##################################################################



  interface dmft_weiss
     !push zeta array
     module procedure :: gf_push_zeta
     module procedure :: dmft_get_weiss_normal_rank2
     module procedure :: dmft_get_weiss_normal_rank3
     module procedure :: dmft_get_weiss_normal_rank4
     module procedure :: dmft_get_weiss_normal_rank5
     module procedure :: dmft_get_weiss_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_weiss_normal_rank7
#endif
     module procedure :: dmft_get_weiss_superc_rank2
     module procedure :: dmft_get_weiss_superc_rank3
     module procedure :: dmft_get_weiss_superc_rank4
     module procedure :: dmft_get_weiss_superc_rank5
     module procedure :: dmft_get_weiss_superc_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8  
     module procedure :: dmft_get_weiss_superc_rank7
#endif
  end interface dmft_weiss


  interface dmft_delta
     !push zeta array
     module procedure :: gf_push_zeta
     module procedure :: dmft_get_delta_normal_rank2
     module procedure :: dmft_get_delta_normal_rank3
     module procedure :: dmft_get_delta_normal_rank4
     module procedure :: dmft_get_delta_normal_rank5
     module procedure :: dmft_get_delta_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_delta_normal_rank7
#endif
     module procedure :: dmft_get_delta_superc_rank2
     module procedure :: dmft_get_delta_superc_rank3
     module procedure :: dmft_get_delta_superc_rank4
     module procedure :: dmft_get_delta_superc_rank5
     module procedure :: dmft_get_delta_superc_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8  
     module procedure :: dmft_get_delta_superc_rank7
#endif
  end interface dmft_delta



  interface dmft_self_consistency
     !push zeta array
     module procedure :: gf_push_zeta
     module procedure :: dmft_get_delta_normal_rank2
     module procedure :: dmft_get_delta_normal_rank3
     module procedure :: dmft_get_delta_normal_rank4
     module procedure :: dmft_get_delta_normal_rank5
     module procedure :: dmft_get_delta_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_delta_normal_rank7
#endif
     module procedure :: dmft_get_delta_superc_rank2
     module procedure :: dmft_get_delta_superc_rank3
     module procedure :: dmft_get_delta_superc_rank4
     module procedure :: dmft_get_delta_superc_rank5
     module procedure :: dmft_get_delta_superc_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8  
     module procedure :: dmft_get_delta_superc_rank7
#endif
     module procedure :: dmft_get_weiss_normal_rank2
     module procedure :: dmft_get_weiss_normal_rank3
     module procedure :: dmft_get_weiss_normal_rank4
     module procedure :: dmft_get_weiss_normal_rank5
     module procedure :: dmft_get_weiss_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_weiss_normal_rank7
#endif
     module procedure :: dmft_get_weiss_superc_rank2
     module procedure :: dmft_get_weiss_superc_rank3
     module procedure :: dmft_get_weiss_superc_rank4
     module procedure :: dmft_get_weiss_superc_rank5
     module procedure :: dmft_get_weiss_superc_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8  
     module procedure :: dmft_get_weiss_superc_rank7
#endif

  end interface dmft_self_consistency



  public :: dmft_weiss
  public :: dmft_delta
  public :: dmft_self_consistency



end module DMFT_WEISS_FIELD
