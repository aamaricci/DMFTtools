module ARCHIVED_SC_GLOBAL
  USE DMFT_CTRL_VARS
  USE ARCHIVED_SC_COMMON
  USE ARCHIVED_SC_WEISS
  USE ARCHIVED_SC_DELTA
  implicit none
  private



  interface dmft_self_consistency
     module procedure :: dmft_sc_normal_main
     module procedure :: dmft_sc_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8     
     module procedure :: dmft_sc_normal_cluster_ineq
#endif
     module procedure :: dmft_sc_normal_ineq
     module procedure :: dmft_sc_normal_bethe
     module procedure :: dmft_sc_normal_bethe_ineq
     module procedure :: dmft_sc_superc_main
     module procedure :: dmft_sc_superc_ineq
  end interface dmft_self_consistency


  public :: dmft_self_consistency


contains

  subroutine dmft_sc_normal_main(Gloc,Smats,Weiss,Hloc,SCtype)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    character(len=*)                              :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Smats,Weiss,Hloc)
    case ("delta")
       call dmft_delta(Gloc,Smats,Weiss,Hloc)
    end select
  end subroutine dmft_sc_normal_main


  subroutine dmft_sc_normal_cluster(Gloc,Smats,Weiss,Hloc,SCtype)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    character(len=*)                                  :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Smats,Weiss,Hloc)
    case ("delta")
       call dmft_delta(Gloc,Smats,Weiss,Hloc)
    end select
  end subroutine dmft_sc_normal_cluster


  subroutine dmft_sc_normal_ineq(Gloc,Smats,Weiss,Hloc,SCtype)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    character(len=*)                                :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Smats,Weiss,Hloc)
    case ("delta")
       call dmft_delta(Gloc,Smats,Weiss,Hloc)
    end select
  end subroutine dmft_sc_normal_ineq


#if __GFORTRAN__ &&  __GNUC__ > 8 
  subroutine dmft_sc_normal_cluster_ineq(Gloc,Smats,Weiss,Hloc,SCtype)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Smats        ! [Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Hloc         ! [Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    character(len=*)                                    :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Smats,Weiss,Hloc)
    case ("delta")
       call dmft_delta(Gloc,Smats,Weiss,Hloc)
    end select
  end subroutine dmft_sc_normal_cluster_ineq
#endif

  subroutine dmft_sc_normal_bethe(Gloc,Weiss,Hloc,Wbands,SCtype)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    real(8),dimension(:),intent(in)               :: Wbands ![Nspin*Norb]
    character(len=*)                              :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Weiss,Hloc,Wbands)
    case ("delta")
       call dmft_delta(Gloc,Weiss,Hloc,Wbands)
    end select
  end subroutine dmft_sc_normal_bethe

  subroutine dmft_sc_normal_bethe_ineq(Gloc,Weiss,Hloc,Wbands,SCtype)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    real(8),dimension(:),intent(in)               :: Wbands         ![Nlat*Nspin*Norb]
    character(len=*)                                :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Weiss,Hloc,Wbands)
    case ("delta")
       call dmft_delta(Gloc,Weiss,Hloc,Wbands)
    end select
  end subroutine dmft_sc_normal_bethe_ineq






  subroutine dmft_sc_superc_main(Gloc,Floc,Smats,SAmats,Weiss,aWeiss,Hloc,SCtype)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc         ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Floc         ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats        !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: SAmats        !
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss        !
    complex(8),dimension(:,:,:,:,:),intent(inout) :: aWeiss        !
    complex(8),dimension(:,:,:,:),intent(in)        :: Hloc         ! [Nspin][Nspin][Norb][Norb]
    character(len=*)                                :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Floc,Smats,SAmats,Weiss,aWeiss,Hloc)
    case ("delta")
       call dmft_delta(Gloc,Floc,Smats,SAmats,Weiss,aWeiss,Hloc)
    end select
  end subroutine dmft_sc_superc_main

  subroutine dmft_sc_superc_ineq(Gloc,Floc,Smats,SAmats,Weiss,aWeiss,Hloc,SCtype)
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Floc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats        !
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: SAmats        ! 
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: aWeiss        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    character(len=*)                                  :: SCtype
    select case(SCtype)
    case default
       call dmft_weiss(Gloc,Floc,Smats,SAmats,Weiss,aWeiss,Hloc)
    case ("delta")
       call dmft_delta(Gloc,Floc,Smats,SAmats,Weiss,aWeiss,Hloc)
    end select
  end subroutine dmft_sc_superc_ineq

end module ARCHIVED_SC_GLOBAL
