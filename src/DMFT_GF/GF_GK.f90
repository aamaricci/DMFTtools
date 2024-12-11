module GF_GK
  USE GF_COMMON
  implicit none
  private


  !        NORMAL 
  ! . rank-2   [N,N,:]
  ! . rank-3   [Nlat,Nso,Nso,:]
  ! . rank-4   [Nspin,Nspin,Norb,Norb,:]
  ! . rank-5   [Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5_6 [Nlat,Nspin,Nspin,Norb,Norb,:] ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-6   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-7   [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  public :: get_gk_normal_hk_rank2
  public :: get_gk_normal_hk_rank3 
  public :: get_gk_normal_hk_rank4
  public :: get_gk_normal_hk_rank5
  public :: get_gk_normal_hk_rank5_6
  public :: get_gk_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gk_normal_hk_rank7
#endif
  !
  public :: get_gk_normal_tridiag_rank2
  public :: get_gk_normal_tridiag_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gk_normal_tridiag_rank7
#endif
  !
  public :: get_gk_normal_dos_rank2  
  public :: get_gk_normal_dos_rank4
  public :: get_gk_normal_dos_rank5
  !
  !        SUPERC
  ! . rank-2   [2,N,N,:]
  ! . rank-3   [2,Nlat,Nso,Nso,:]
  ! . rank-4   [2,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5   [2,Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5_6 [2,Nlat,Nspin,Nspin,Norb,Norb,:]->[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  public :: get_gk_superc_hk_rank2
  public :: get_gk_superc_hk_rank3
  public :: get_gk_superc_hk_rank4
  public :: get_gk_superc_hk_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gk_superc_hk_rank5_6
#endif
  !
  public :: get_gk_superc_dos_rank2
  public :: get_gk_superc_dos_rank4



contains



  !##################################################################
  !##################################################################
  !                       NORMAL PHASE
  ! . rank-2 DMFT    [N,N,:] Hk Matrix
  ! . rank-2 DMFT    [N,N,:] Hk TriDiag
  ! . rank-2 DMFT    [N,N,:] Rho(e) DOS
  !
  ! . rank-3 DMFT    [Nlat,Nso,Nso,:] Hk Matrix
  !
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] Hk Matrix
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] Rho(e) DOS
  !
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:] Hk TriDiag
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:] Rho(e) DOS
  ! . rank-5_6 R-DMFT[Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  !                ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  !
  ! . rank-6 CDMFT   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  !                   Nlat=N_cluster
  !
  ! . rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  ! . rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk TriDiag 
  !                    Nineq=N_unit_cell,
  !                    Nlat=N_cluster
  !
  ! Hk Matrix  :: [N,N,Lk]
  !
  ! Hk TriDiag :: [[Ncell,Ncell] x Nsites]  
  ! Block.shape   = [Ncell,Ncell]
  ! Block.number  = Nsites
  ! N             = Nsites*Ncell
  !
  ! Rho(e) DOS :: [N,Le]
  !##################################################################
  ! G/Sigma Shape: [N,N][:]
  !##################################################################
  subroutine get_gk_normal_hk_rank2(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)    :: Hk
    complex(8),dimension(:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:),intent(inout) :: Gk
    character(len=*)                          :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lfreq = size(Sigma,3)
    Nlso  = Ntot
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [N,N], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_hk_rank2","Hk")
    call assert_shape(Sigma,[Ntot,Nlso,Lfreq],"get_gk_normal_hk_rank2","Sigma")
    call assert_shape(Gk, [Ntot,Nlso,Lfreq],"get_gk_normal_hk_rank2","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,i) = gki_normal(wfreq(i)+xmu,Hk(:,:),Sigma(:,:,i))
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank2


  ! G/Sigma Shape: ![Nlat,Nso,Nso][:]
  !##################################################################
  subroutine get_gk_normal_hk_rank3(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)        :: Hk
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:),intent(inout) :: Gk
    character(len=*)                            :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nso   = size(Sigma,2)
    Lfreq = size(Sigma,4)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nso,Nso], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nso],"get_gk_normal_hk_rank3","Hk")
    call assert_shape(Sigma,[Nlat,Nso,Nso,Lfreq],"get_gk_normal_hk_rank3","Sigma")
    call assert_shape(Gk, [Nlat,Nso,Nso,Lfreq],"get_gk_normal_hk_rank3","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,i) =  to_rank3( gki_normal(wfreq(i)+xmu,Hk(:,:),from_rank3(Sigma(:,:,:,i))) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank3


  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_hk_rank4(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)          :: Hk
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                              :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nspin*Norb],"get_gk_normal_hk_rank4","Hk")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank4","Sigma")
    call assert_shape(Gk, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank4","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,i) = to_rank4( gki_normal(wfreq(i)+xmu,Hk(:,:),from_rank4(Sigma(:,:,:,:,i))) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank4

  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_hk_rank5(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)            :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nspin*Norb],"get_gk_normal_hk_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank5","Sigma")
    call assert_shape(Gk, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank5","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,:,i) = to_rank5( gki_normal(wfreq(i)+xmu,Hk(:,:),from_rank5(Sigma(:,:,:,:,:,i))) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank5

  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:] -> [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_hk_rank5_6(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                  :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nspin*Norb],"get_gk_normal_hk_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank5","Sigma")
    call assert_shape(Gk, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank5","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,:,:,i) = to_rank6( gki_normal(wfreq(i)+xmu,Hk(:,:),from_rank5(Sigma(:,:,:,:,:,i))) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank5_6



  !G/Sigma shape: [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_hk_rank6(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                  :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nspin*Norb],"get_gk_normal_hk_rank6","Hk")
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank6","Sigma")
    call assert_shape(Gk, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank6","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,:,:,i) = to_rank6( gki_normal(wfreq(i)+xmu,Hk(:,:),from_rank6(Sigma(:,:,:,:,:,:,i))) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank6


  !G/Sigma shape: ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine get_gk_normal_hk_rank7(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:),intent(in)                :: Hk
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                    :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nineq*Nlat*Nspin*Norb],"get_gk_normal_hk_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank7","Sigma")
    call assert_shape(Gk, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_hk_rank7","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,:,:,:,i) = to_rank7( gki_normal(wfreq(i)+xmu,Hk(:,:),from_rank7(Sigma(:,:,:,:,:,:,:,i))) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_hk_rank7
#endif








  !##################################################################
  !##################################################################
  !##################################################################




  !G/Sigma shape: [N,N][:]
  !##################################################################
  subroutine get_gk_normal_dos_rank2(Ebands,Dbands,Hloc,Gk,Sigma,axis,diagonal)
    real(8),dimension(:),intent(in)            :: Ebands    ![N]
    real(8),dimension(:),intent(in)            :: Dbands    ![N] /[1]
    real(8),dimension(size(Ebands,1)),intent(in) :: Hloc      ![N]
    complex(8),dimension(:,:,:),intent(in)       :: Sigma     ![N,N][Lmats]
    complex(8),dimension(:,:,:),intent(inout)    :: Gk      !as Sigma
    character(len=*)                             :: axis
    !
    complex(8),dimension(:,:),allocatable        :: zeta
    complex(8),dimension(:,:),allocatable        :: Gdos_tmp
    logical                                      :: dos_diag !1. T / 2. F
    logical,optional                             :: diagonal
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,1)
    Nso   = size(Sigma,1)
    Lfreq = size(Sigma,3)
    !
    if(Nso/=Ntot)stop "get_gk_normal_rank2 error: Ntot != Nso"
    if(mpi_master)write(*,"(A)")"Get Green's function [N,N], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm."
    if(mpi_master)call start_timer
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    if(present(diagonal))then
       dos_diag=diagonal
    else
       dos_diag = .not.(size(Dbands) < size(Ebands))
    endif
    !
    !Testing part:
    call assert_shape(Ebands,[Ntot],"dmft_get_gk_normal_dos","Ebands")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"dmft_get_gk_normal_dos","Sigma")
    call assert_shape(Gk,[Ntot,Ntot,Lfreq],"dmft_get_gk_normal_dos","Gk")
    !
    call build_frequency_array(axis)
    !
    !invert (Z-Hk) for each k-point
    allocate(zeta(Ntot,Ntot))
    !
    Gk = zero
    dosdiag: if(dos_diag)then
       do i=1+mpi_rank, Lfreq, mpi_size
          zeta=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
          do io=1,Ntot
             Gk(io,io,i) = Dbands(io)/( zeta(io,io)-Hloc(io)-Ebands(io) )
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       do i = 1+mpi_rank, Lfreq, mpi_size
          zeta     = (wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
          Gdos_tmp = zeta - diag(Hloc(:)) - diag(Ebands(:))
          call inv(Gdos_tmp)
          Gk(:,:,i) = Dbands(1)*Gdos_tmp
          if(mpi_master)call eta(i,Lfreq)
       end do
    end if dosdiag
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
    !
  end subroutine get_gk_normal_dos_rank2


  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_dos_rank4(Ebands,Dbands,Hloc,Gk,Sigma,axis,diagonal)
    real(8),dimension(:),intent(in)               :: Ebands    ![N]
    real(8),dimension(:),intent(in)               :: Dbands    ![N] /[1]
    real(8),dimension(size(Ebands,1)),intent(in)  :: Hloc      ![N]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                              :: axis
    logical,optional                              :: diagonal
    complex(8),dimension(:,:),allocatable         :: zeta
    complex(8),dimension(:,:),allocatable         :: Gdos_tmp
    logical                                       :: dos_diag !1. T / 2. F
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,1)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Nso   = Nspin*Norb
    if(Nso/=Ntot)stop "get_gk_normal_rank4 error: Ntot != Nspin*Norb"
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm."
    if(mpi_master)call start_timer
    !
    !
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank4","Sigma")
    call assert_shape(Gk, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank4","Gk")
    !
    allocate(Zeta(Ntot,Ntot))
    !
    Gk = zero
    dosdiag:if(dos_diag)then
       do i=1+mpi_rank, Lfreq, mpi_size
          zeta = (wfreq(i)+xmu)*eye(Ntot) - from_rank4(Sigma(:,:,:,:,i))
          do concurrent(ispin=1:Nspin,iorb=1:Norb)
             io = iorb + (ispin-1)*Norb
             Gk(ispin,ispin,iorb,iorb,i) = Dbands(io)/( zeta(io,io)-Hloc(io)-Ebands(io) )
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       do i = 1+mpi_rank, Lfreq, mpi_size
          zeta     = (wfreq(i)+xmu)*eye(Ntot) - from_rank4(Sigma(:,:,:,:,i))
          Gdos_tmp = zeta-diag(Hloc(:))-diag(Ebands(:))
          call inv(Gdos_tmp)
          Gk(:,:,:,:,i) = Dbands(1)*to_rank4(Gdos_tmp)
          if(mpi_master)call eta(i,Lfreq)
       end do
    end if dosdiag
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_dos_rank4


  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_dos_rank5(Ebands,Dbands,Hloc,Gk,Sigma,axis,diagonal)
    real(8),dimension(:),intent(in)               :: Ebands ![N]
    real(8),dimension(:),intent(in)               :: Dbands ![N] /[1]
    real(8),dimension(size(Ebands,1)),intent(in)    :: Hloc   ![N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma  ![Nlat,Nspin,Nspin,Norb,Norb][:]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk   
    character(len=*)                                :: axis
    logical,optional                                :: diagonal
    complex(8),dimension(:,:),allocatable           :: zeta
    complex(8),dimension(:,:),allocatable           :: Gdos_tmp
    logical                                         :: dos_diag !1. T / 2. F
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)    
    Nlso  = Nlat*Nspin*Norb
    if(Nlso/=Ntot)stop "get_gk_normal_rank5 error: Ntot != Nlat*Nspin*Norb"
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm:"
    if(mpi_master)call start_timer
    !
    !
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5","Sigma")
    call assert_shape(Gk, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_normal_rank5","Gk")
    !
    allocate(Zeta(Ntot,Ntot))
    !
    Gk = zero
    dosdiag:if(dos_diag)then
       do i=1+mpi_rank, Lfreq, mpi_size
          zeta = (wfreq(i)+xmu)*eye(Ntot)-from_rank5(Sigma(:,:,:,:,:,i))
          do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
             io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
             Gk(ilat,ispin,ispin,iorb,iorb,i) = &
                  Dbands(io)/( zeta(io,io)-Hloc(io)-Ebands(io) )                 
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       do i = 1+mpi_rank, Lfreq, mpi_size
          zeta = (wfreq(i)+xmu)*eye(Ntot)-from_rank5(Sigma(:,:,:,:,:,i))
          Gdos_tmp = zeta-diag(Hloc(:))-diag(Ebands(:))
          call inv(Gdos_tmp)
          Gk(:,:,:,:,:,i) = Dbands(1)*to_rank5(Gdos_tmp)
          if(mpi_master)call eta(i,Lfreq)
       end do
    end if dosdiag
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_dos_rank5




  !##################################################################
  !##################################################################
  !##################################################################




  ! G/Sigma Shape: [N,N][:]
  !##################################################################
  subroutine get_gk_normal_tridiag_rank2(Hk,Gk,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:),intent(in)      :: Hk
    complex(8),dimension(:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:),intent(inout) :: Gk
    character(len=*)                          :: axis
    integer,intent(in)                        :: Nsites,Ncell
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lfreq = size(Sigma,3)
    Nlso  = Ntot
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [N,N], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm."
    if(mpi_master)call start_timer
    !  
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop "get_gk_normal_tridiag_rank2 erro:  Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Nlso],'get_gk_normal_tridiag_rank2',"Hk")
    call assert_shape(Sigma,[Ntot,Nlso,Lfreq],'get_gk_normal_tridiag_rank2',"Sigma")
    call assert_shape(Gk,[Ntot,Nlso,Lfreq],'get_gk_normal_tridiag_rank2',"Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,i) = gki_tridiag(wfreq(i)+xmu,Hk(:,:),Sigma(:,:,i),Nsites,Ncell)
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_tridiag_rank2


  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_tridiag_rank5(Hk,Gk,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:),intent(in)          :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                :: axis
    integer,intent(in)                              :: Nsites,Ncell
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm:"
    if(mpi_master)call start_timer
    !
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop &
         "get_gk_normal_tridiag_rank5 erro:  Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_tridiag_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gk_normal_tridiag_rank5","Sigma")
    call assert_shape(Gk, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gk_normal_tridiag_rank5","Gk")
    !
    call build_frequency_array(axis)
    !
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,:,i) = to_rank5( gki_tridiag(wfreq(i)+xmu,Hk(:,:), from_rank5(Sigma(:,:,:,:,:,i)),Nsites,Ncell) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_tridiag_rank5


#if __GFORTRAN__ &&  __GNUC__ > 8
  !G/Sigma shape: [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_normal_tridiag_rank7(Hk,Gk,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:),intent(in)                :: Hk
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma  
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                    :: axis
    integer,intent(in)                                  :: Nsites,Ncell
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nlso  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm."
    if(mpi_master)call start_timer
    !
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop &
         "get_gk_normal_tridiag_rank7 erro:  Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Nlso],"get_gk_normal_tridiag_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gk_normal_tridiag_rank7","Sigma")
    call assert_shape(Gk, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gk_normal_tridiag_rank7","Gk")
    !
    call build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    Gk = zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Gk(:,:,:,:,:,:,:,i) = to_rank7( gki_tridiag(wfreq(i)+xmu,Hk(:,:),from_rank7(Sigma(:,:,:,:,:,:,:,i)),Nsites,Ncell) )
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_normal_tridiag_rank7
#endif







  !##################################################################
  !##################################################################
  !                         SUPERC
  ! . rank-2 DMFT    [N,N,:] Hk Matrix
  ! . rank-2 DMFT    [N,N,:] Rho(e) DOS
  !
  ! . rank-3 DMFT    [Nlat,Nso,Nso,:] Hk Matrix
  !
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] Hk Matrix
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] Rho(e) DOS
  !
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  !
  ! . rank-5_6 R-DMFT[Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  !                ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  !
  ! Hk Matrix  :: [N,N,Lk]
  !
  ! Rho(e) DOS :: [N,Le]
  !##################################################################
  ! G/Sigma Shape: [2][N,N][:]
  !##################################################################
  subroutine get_gk_superc_hk_rank2(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:),intent(in)            :: Sigma
    complex(8),dimension(:,:,:,:),intent(inout)         :: Gk
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,2)  
    Lfreq = size(Sigma,4)
    Nlso  = Ntot
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][N,N], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Nlso],'get_gk_superc',"Hk")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gk_superc',"Sigma")
    call assert_shape(Gk,[2,Ntot,Ntot,Lfreq],'get_gk_superc',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Z           = sc_zi(i,Sigma,axis)
       Gk(:,:,:,i) = gki_superc_rank2(Z,Hk(:,:,:))
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_hk_rank2


  !G/Sigma shape: [2][Nlat,Nso,Nso][:]
  !##################################################################
  subroutine get_gk_superc_hk_rank3(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:),intent(in)          :: Sigma
    complex(8),dimension(:,:,:,:,:),intent(inout)       :: Gk
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2) 
    Nlat  = size(Sigma,2)
    Nso   = size(Sigma,3)
    Lfreq = size(Sigma,5)
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [Nlat,Nso,Nso], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Nso],'get_gk_superc',"Hk")
    call assert_shape(Sigma,[2,Nlat,Ntot,Nso,Lfreq],'get_gk_superc',"Sigma")
    call assert_shape(Gk,[2,Nlat,Nso,Nso,Lfreq],'get_gk_superc',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Z             = sc_zi(i,Sigma,axis)
       Gk(:,:,:,:,i) = gki_superc_rank3(Z,Hk(:,:,:))
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_hk_rank3

  !G/Sigma shape: [2][Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_superc_hk_rank4(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)        :: Sigma
    complex(8),dimension(:,:,:,:,:,:),intent(inout)     :: Gk
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Nso],'get_gk_superc',"Hk")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],'get_gk_superc',"Sigma")
    call assert_shape(Gk,[2,Nspin,Nspin,Norb,Norb,Lfreq],'get_gk_superc',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Z               = sc_zi(i,Sigma,axis)
       Gk(:,:,:,:,:,i) = gki_superc_rank4(Z,Hk(:,:,:))
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_hk_rank4

  !G/Sigma shape: [2][Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_superc_hk_rank5(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout)   :: Gk
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)  
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nlat,Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Nlso],'get_gk_superc',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gk_superc',"Sigma")
    call assert_shape(Gk,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gk_superc',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Z                 = sc_zi(i,Sigma,axis)
       Gk(:,:,:,:,:,:,i) = gki_superc_rank5(Z,Hk(:,:,:))
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_hk_rank5


#if __GFORTRAN__ &&  __GNUC__ > 8
  !Sigma shape: [2][Nlat,Nspin,Nspin,Norb,Norb][:]
  !G     shape: [2][Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_superc_hk_rank5_6(Hk,Gk,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gk
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2) 
    Nlat = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nlat,Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Nlso],'get_gk_superc',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gk_superc',"Sigma")
    call assert_shape(Gk,[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gk_superc',"Gk")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       Z                   = sc_zi(i,Sigma,axis)
       Gk(:,:,:,:,:,:,:,i) = gki_superc_rank6(Z,Hk(:,:,:))
       if(mpi_master)call eta(i,Lfreq)
    end do
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_hk_rank5_6
#endif  







  ! G/Sigma Shape: [2][N,N][:]
  !##################################################################
  subroutine get_gk_superc_dos_rank2(Ebands,Dbands,Hloc,Gk,Sigma,axis)
    real(8),dimension(:,:),intent(in)                   :: Ebands    ![2,N]
    real(8),dimension(:),intent(in)                     :: Dbands    ![N]
    real(8),dimension(2,size(Ebands,2)),intent(in)      :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:),intent(in)            :: Sigma     ![2,N,N,Lfreq]
    complex(8),dimension(:,:,:,:),intent(inout)         :: Gk      !as Sigma
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z ![2][2][N,N]
    complex(8),dimension(:,:),allocatable               :: G
    integer                                             :: N,N2
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,2)
    Nso  = size(Sigma,2)
    Lfreq = size(Sigma,4)
    !
    if(Nso/=Ntot)stop "get_gk_superc_dos_rank2 error: Ntot != Nso"
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][N,N], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm. WARNING: only diagonal case.. "
    if(mpi_master)call start_timer
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[2,Nso],'get_gk_superc_dos_rank2',"Ebands")
    call assert_shape(Dbands,[Ntot],'get_gk_superc_dos_rank2',"Dbands")
    call assert_shape(Sigma,[2,Ntot,Nso,Lfreq],'get_gk_superc_dos_rank2',"Sigma")
    call assert_shape(Gk,[2,Ntot,Nso,Lfreq],'get_gk_superc_dos_rank2',"Gk")
    !
    call build_frequency_array(axis)
    !
    !
    !invert (Z-Hk) for each k-point
    N  = Ntot
    N2 = 2*Ntot
    allocate(G(N2,N2))
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk = zero
    do i=1+mpi_rank, Lfreq, mpi_size
       Z  = sc_zi_rank2(i,Sigma,axis)
       G(1:N,1:N)       = Z(1,1,:,:) - diag(Hloc(1,:)) - diag(Ebands(1,:))
       G(1:N,N+1:N2)    = Z(1,2,:,:)        
       G(N+1:N2,1:N)    = Z(2,1,:,:)
       G(N+1:N2,N+1:N2) = Z(2,2,:,:) + diag(Hloc(2,:)) - diag(Ebands(2,:))
       call inv(G)
       do concurrent(io=1:Ntot,jo=1:Ntot)
          Gk(1,io,jo,i) = G(io,  jo)*Dbands(io)
          Gk(2,io,jo,i) = G(io,N+jo)*Dbands(io)
       enddo
    enddo
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_dos_rank2

  !G/Sigma shape: [2][Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gk_superc_dos_rank4(Ebands,Dbands,Hloc,Gk,Sigma,axis)
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2,N]
    real(8),dimension(:),intent(in)                 :: Dbands    ![N]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2,N,N,Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gk      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: Z ![2][2][N,N]
    complex(8),dimension(:,:),allocatable           :: G
    integer                                         :: N,N2
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,2)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !
    if(Nso/=Ntot)stop "get_gk_superc_rank4 error: Ntot != Nso"
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm. WARNING: only diagonal case.."
    if(mpi_master)call start_timer
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[2,Nso],'get_gk_superc_dos_rank4',"Ebands")
    call assert_shape(Dbands,[Ntot],'get_gk_superc_dos_rank4',"Dbands")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_dos_rank4","Sigma")
    call assert_shape(Gk, [2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gk_superc_dos_rank4","Gk")
    !
    call build_frequency_array(axis)
    !
    N  = Ntot
    N2 = 2*Ntot
    allocate(G(N2,N2))
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gk = zero
    do i=1+mpi_rank, Lfreq, mpi_size
       Z  = sc_zi_rank4(i,Sigma,axis)
       G(1:N,1:N)       = Z(1,1,:,:) - diag(Hloc(1,:)) - diag(Ebands(1,:))
       G(1:N,N+1:N2)    = Z(1,2,:,:)        
       G(N+1:N2,1:N)    = Z(2,1,:,:)
       G(N+1:N2,N+1:N2) = Z(2,2,:,:) + diag(Hloc(2,:)) - diag(Ebands(2,:))
       call inv(G)
       do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
          io = iorb + (ispin-1)*Norb
          jo = jorb + (jspin-1)*Norb
          Gk(1,ispin,jspin,iorb,jorb,i) = G(io,  jo)*Dbands(io)
          Gk(2,ispin,jspin,iorb,jorb,i) = G(io,N+jo)*Dbands(io)
       enddo
       call eta(i,Lfreq)
    enddo
    call gf_reduce(Gk)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gk_superc_dos_rank4


end module GF_GK
