module GF_GLOC
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
  public :: get_gloc_normal_hk_rank2
  public :: get_gloc_normal_tridiag_rank2
  public :: get_gloc_normal_dos_rank2
  !
  public :: get_gloc_normal_hk_rank3 
  !
  public :: get_gloc_normal_hk_rank4
  public :: get_gloc_normal_dos_rank4
  !
  public :: get_gloc_normal_hk_rank5
  public :: get_gloc_normal_tridiag_rank5
  public :: get_gloc_normal_dos_rank5
  public :: get_gloc_normal_hk_rank5_6
  !
  public :: get_gloc_normal_hk_rank6
  !
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gloc_normal_hk_rank7
  public :: get_gloc_normal_tridiag_rank7
#endif


  !
  !        SUPERC
  ! . rank-2   [2,N,N,:]
  ! . rank-3   [2,Nlat,Nso,Nso,:]
  ! . rank-4   [2,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5   [2,Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5_6 [2,Nlat,Nspin,Nspin,Norb,Norb,:]->[2,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  public :: get_gloc_superc_hk_rank2
  public :: get_gloc_superc_dos_rank2
  public :: get_gloc_superc_hk_rank3
  public :: get_gloc_superc_hk_rank4
  public :: get_gloc_superc_dos_rank4
  public :: get_gloc_superc_hk_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gloc_superc_hk_rank5_6
#endif
  !



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
  subroutine get_gloc_normal_hk_rank2(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)    :: Hk
    complex(8),dimension(:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:),intent(inout) :: Gloc
    character(len=*)                          :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Lfreq = size(Sigma,3)
    Nlso  = Ntot
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [N,N], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_hk_rank2","Hk")
    call assert_shape(Sigma,[Ntot,Nlso,Lfreq],"get_gloc_normal_hk_rank2","Sigma")
    call assert_shape(Gloc, [Ntot,Nlso,Lfreq],"get_gloc_normal_hk_rank2","Gloc")
    !
    call build_frequency_array(axis)
    !
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,i) = Gloc(:,:,i) + invert_gki_normal_rank2(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq        
             Gloc(:,:,i) = Gloc(:,:,i) + invert_gki_normal_rank2(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank2

  subroutine get_gloc_normal_tridiag_rank2(Hk,Gloc,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:,:),intent(in)    :: Hk
    complex(8),dimension(:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:),intent(inout) :: Gloc
    character(len=*)                          :: axis
    integer,intent(in)                        :: Nsites,Ncell
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Lfreq = size(Sigma,3)
    Nlso  = Ntot
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [N,N], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm."
    if(mpi_master)call start_timer
    !  
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop "get_gloc_normal_tridiag_rank2 erro:  Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Nlso,Lk],'get_gloc_normal_tridiag_rank2',"Hk")
    call assert_shape(Sigma,[Ntot,Nlso,Lfreq],'get_gloc_normal_tridiag_rank2',"Sigma")
    call assert_shape(Gloc,[Ntot,Nlso,Lfreq],'get_gloc_normal_tridiag_rank2',"Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,i) = Gloc(:,:,i) + &
                  invert_gki_normal_tridiag_rank2(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,i),Nsites,Ncell)/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,i) = Gloc(:,:,i) + &
                  invert_gki_normal_tridiag_rank2(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,i),Nsites,Ncell)/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
    !
  end subroutine get_gloc_normal_tridiag_rank2

  subroutine get_gloc_normal_dos_rank2(Ebands,Dbands,Hloc,Gloc,Sigma,axis,diagonal)
    real(8),dimension(:,:),intent(in)            :: Ebands    ![N][Lk]
    real(8),dimension(:,:),intent(in)            :: Dbands    ![N][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in) :: Hloc      ![N]
    complex(8),dimension(:,:,:),intent(in)       :: Sigma     ![N,N][Lmats]
    complex(8),dimension(:,:,:),intent(inout)    :: Gloc      !as Sigma
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
    Lk    = size(Ebands,2)
    Nso   = size(Sigma,1)
    Lfreq = size(Sigma,3)
    !
    if(Nso/=Ntot)stop "get_gloc_normal_rank2 error: Ntot != Nso"
    if(mpi_master)write(*,"(A)")"Get Green's function [N,N], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm."
    if(mpi_master)call start_timer
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    if(present(diagonal))then
       dos_diag=diagonal
    else
       dos_diag = .not.(size(Dbands,1) < size(Ebands,1))
    endif
    !
    !Testing part:
    call assert_shape(Ebands,[Ntot,Lk],"dmft_get_gloc_normal_dos","Ebands")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"dmft_get_gloc_normal_dos","Sigma")
    call assert_shape(Gloc,[Ntot,Ntot,Lfreq],"dmft_get_gloc_normal_dos","Gloc")
    !
    call build_frequency_array(axis)
    !
    !invert (Z-Hk) for each k-point
    allocate(zeta(Ntot,Ntot))
    !
    Gloc = zero
    dosdiag: if(dos_diag)then
       dLfMpi: if(Lfreq>=Lk)then
          do i=1+mpi_rank, Lfreq, mpi_size
             zeta=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
             do ik=1,Lk
                do io=1,Ntot
                   Gloc(io,io,i) = Gloc(io,io,i) + Dbands(io,ik)/( zeta(io,io)-Hloc(io)-Ebands(io,ik) )
                enddo
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                zeta=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
                do io=1,Ntot
                   Gloc(io,io,i) = Gloc(io,io,i) + Dbands(io,ik)/( zeta(io,io)-Hloc(io)-Ebands(io,ik) )
                enddo
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if dLfMpi
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       ndLfMpi: if(Lfreq>=Lk)then
          do i = 1+mpi_rank, Lfreq, mpi_size
             zeta  = (wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
             do ik=1,Lk
                Gdos_tmp = zeta - diag(Hloc(:)) - diag(Ebands(:,ik))
                call inv(Gdos_tmp)
                Gloc(:,:,i) = Gloc(:,:,i) + Dbands(1,ik)*Gdos_tmp
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                zeta     = (wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
                Gdos_tmp = zeta - diag(Hloc(:)) - diag(Ebands(:,ik))
                call inv(Gdos_tmp)
                Gloc(:,:,i) = Gloc(:,:,i) + Dbands(1,ik)*Gdos_tmp
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if ndLfMpi
    end if dosdiag
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
    !
  end subroutine get_gloc_normal_dos_rank2


  ! G/Sigma Shape: ![Nlat,Nso,Nso][:]
  !##################################################################
  subroutine get_gloc_normal_hk_rank3(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)      :: Hk
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:),intent(inout) :: Gloc
    character(len=*)                            :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nso   = size(Sigma,2)
    Lfreq = size(Sigma,4)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nso,Nso], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nso,Lk],"get_gloc_normal_hk_rank3","Hk")
    call assert_shape(Sigma,[Nlat,Nso,Nso,Lfreq],"get_gloc_normal_hk_rank3","Sigma")
    call assert_shape(Gloc, [Nlat,Nso,Nso,Lfreq],"get_gloc_normal_hk_rank3","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,i) = Gloc(:,:,:,i) + &
                  invert_gki_normal_rank3(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,i) = Gloc(:,:,:,i) + &
                  invert_gki_normal_rank3(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank3


  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_normal_hk_rank4(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)        :: Hk
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                              :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nspin*Norb,Lk],"get_gloc_normal_hk_rank4","Hk")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank4","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,:,i) = Gloc(:,:,:,:,i) + &
                  invert_gki_normal_rank4(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,i) = Gloc(:,:,:,:,i) + &
                  invert_gki_normal_rank4(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank4

  subroutine get_gloc_normal_dos_rank4(Ebands,Dbands,Hloc,Gloc,Sigma,axis,diagonal)
    real(8),dimension(:,:),intent(in)             :: Ebands    ![N][Lk]
    real(8),dimension(:,:),intent(in)             :: Dbands    ![N][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)  :: Hloc      ![N]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
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
    Lk    = size(Ebands,2)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Nso   = Nspin*Norb
    if(Nso/=Ntot)stop "get_gloc_normal_rank4 error: Ntot != Nspin*Norb"
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm."
    if(mpi_master)call start_timer
    !
    !
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank4","Gloc")
    !
    allocate(Zeta(Ntot,Ntot))
    !
    Gloc = zero
    dosdiag:if(dos_diag)then
       dLfmpi: if(Lfreq>=Lk)then
          do i=1+mpi_rank, Lfreq, mpi_size
             zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank4_to_matrix(Sigma(:,:,:,:,i),Nspin,Norb)
             do concurrent(ispin=1:Nspin,iorb=1:Norb)
                io = iorb + (ispin-1)*Norb
                do ik=1,Lk
                   Gloc(ispin,ispin,iorb,iorb,i) = Gloc(ispin,ispin,iorb,iorb,i) + &
                        Dbands(io,ik)/( zeta(io,io)-Hloc(io)-Ebands(io,ik) )                 
                enddo
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank4_to_matrix(Sigma(:,:,:,:,i),Nspin,Norb)
                do concurrent(ispin=1:Nspin,iorb=1:Norb)
                   io = iorb + (ispin-1)*Norb
                   Gloc(ispin,ispin,iorb,iorb,i) = Gloc(ispin,ispin,iorb,iorb,i) + &
                        Dbands(io,ik)/( zeta(io,io)-Hloc(io)-Ebands(io,ik) ) 
                enddo
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if DLfmpi
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       ndLfmpi: if(Lfreq>=Lk)then
          do i = 1+mpi_rank, Lfreq, mpi_size
             zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank4_to_matrix(Sigma(:,:,:,:,i),Nspin,Norb)
             do ik=1,Lk
                Gdos_tmp = zeta-diag(Hloc(:))-diag(Ebands(:,ik))
                call inv(Gdos_tmp)
                Gloc(:,:,:,:,i) = Gloc(:,:,:,:,i)+Dbands(1,ik)*reshape_matrix_to_rank4(Gdos_tmp,Nspin,Norb)
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                zeta     = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank4_to_matrix(Sigma(:,:,:,:,i),Nspin,Norb)
                Gdos_tmp = zeta - diag(Hloc(:))-diag(Ebands(:,ik))
                call inv(Gdos_tmp)
                Gloc(:,:,:,:,i) = Gloc(:,:,:,:,i)+Dbands(1,ik)*reshape_matrix_to_rank4(Gdos_tmp,Nspin,Norb)
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if ndLfmpi
    end if dosdiag
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_dos_rank4


  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_normal_hk_rank5(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)          :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nspin*Norb,Lk],"get_gloc_normal_hk_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank5","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + &
                  invert_gki_normal_rank5(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + &
                  invert_gki_normal_rank5(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank5

  subroutine get_gloc_normal_tridiag_rank5(Hk,Gloc,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:,:),intent(in)          :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                :: axis
    integer,intent(in)                              :: Nsites,Ncell
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm."
    if(mpi_master)call start_timer
    !
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop &
         "get_gloc_normal_tridiag_rank5 erro:  Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_tridiag_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gloc_normal_tridiag_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gloc_normal_tridiag_rank5","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + invert_gki_normal_tridiag_rank5(&
                  wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,i),Nsites,Ncell)/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + invert_gki_normal_tridiag_rank5(&
                  wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,i),Nsites,Ncell)/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do

    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_tridiag_rank5

  subroutine get_gloc_normal_dos_rank5(Ebands,Dbands,Hloc,Gloc,Sigma,axis,diagonal)
    real(8),dimension(:,:),intent(in)               :: Ebands ![N][Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands ![N][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)    :: Hloc   ![N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma  ![Nlat,Nspin,Nspin,Norb,Norb][:]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc   
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
    Lk    = size(Ebands,2)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)    
    Nlso  = Nlat*Nspin*Norb
    if(Nlso/=Ntot)stop "get_gloc_normal_rank5 error: Ntot != Nlat*Nspin*Norb"
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm."
    if(mpi_master)call start_timer
    !
    !
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5","Gloc")
    !
    allocate(Zeta(Ntot,Ntot))
    !
    Gloc = zero
    dosdiag:if(dos_diag)then
       dLfmpi: if(Lfreq>=Lk)then
          do i=1+mpi_rank, Lfreq, mpi_size
             zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank5_to_matrix(Sigma(:,:,:,:,:,i),Nlat,Nspin,Norb)
             do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
                io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                do ik=1,Lk
                   Gloc(ilat,ispin,ispin,iorb,iorb,i) = Gloc(ilat,ispin,ispin,iorb,iorb,i) + &
                        Dbands(io,ik)/( zeta(io,io)-Hloc(io)-Ebands(io,ik) )                 
                enddo
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank5_to_matrix(Sigma(:,:,:,:,:,i),Nlat,Nspin,Norb)
                do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   Gloc(ilat,ispin,ispin,iorb,iorb,i) = Gloc(ilat,ispin,ispin,iorb,iorb,i) + &
                        Dbands(io,ik)/( zeta(io,io)-Hloc(io)-Ebands(io,ik) )
                enddo
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if DLfmpi
    else
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       ndLfmpi: if(Lfreq>=Lk)then
          do i = 1+mpi_rank, Lfreq, mpi_size
             zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank5_to_matrix(Sigma(:,:,:,:,:,i),Nlat,Nspin,Norb)
             do ik=1,Lk
                Gdos_tmp = zeta-diag(Hloc(:))-diag(Ebands(:,ik))
                call inv(Gdos_tmp)
                Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + &
                     Dbands(1,ik)*reshape_matrix_to_rank5(Gdos_tmp,Nlat,Nspin,Norb)
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                zeta = (wfreq(i)+xmu)*eye(Ntot)-reshape_rank5_to_matrix(Sigma(:,:,:,:,:,i),Nlat,Nspin,Norb)
                Gdos_tmp = zeta - diag(Hloc(:))-diag(Ebands(:,ik))
                call inv(Gdos_tmp)
                Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + &
                     Dbands(1,ik)*reshape_matrix_to_rank5(Gdos_tmp,Nlat,Nspin,Norb)
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if ndLfmpi
    end if dosdiag
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_dos_rank5



  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:] -> [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_normal_hk_rank5_6(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)            :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                  :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nspin*Norb,Lk],"get_gloc_normal_hk_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank5","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,i) + &
                  invert_gki_normal_rank5_6(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,i) + &
                  invert_gki_normal_rank5_6(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank5_6



  !G/Sigma shape: [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_normal_hk_rank6(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)            :: Hk
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                  :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nlat,Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:
    call assert_shape(Hk,[Ntot,Nlat*Nspin*Norb,Lk],"get_gloc_normal_hk_rank6","Hk")
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank6","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do ik=1,Lk
          do i=1+mpi_rank,Lfreq,mpi_size
             Gloc(:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,i) + &
                  invert_gki_normal_rank6(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,i) + &
                  invert_gki_normal_rank6(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank6


  !G/Sigma shape: ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine get_gloc_normal_hk_rank7(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                    :: axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
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
    call assert_shape(Hk,[Ntot,Nineq*Nlat*Nspin*Norb,Lk],"get_gloc_normal_hk_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_hk_rank7","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,:,i) + &
                  invert_gki_normal_rank7(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,:,i) + &
                  invert_gki_normal_rank7(wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,:,:,i))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_hk_rank7


  subroutine get_gloc_normal_tridiag_rank7(Hk,Gloc,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:,:),intent(in)              :: Hk
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma  
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                    :: axis
    integer,intent(in)                                  :: Nsites,Ncell
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nlso  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm."
    !
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop &
         "get_gloc_normal_tridiag_rank7 erro:  Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_tridiag_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gloc_normal_tridiag_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],&
         "get_gloc_normal_tridiag_rank7","Gloc")
    !
    call build_frequency_array(axis)
    !
    Gloc = zero
    if(Lfreq >= Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Gloc(:,:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,:,i) + invert_gki_normal_tridiag_rank7(&
                  wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,:,:,i),Nsites,Ncell)/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Gloc(:,:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,:,i) + invert_gki_normal_tridiag_rank7(&
                  wfreq(i)+xmu,Hk(:,:,ik),Sigma(:,:,:,:,:,:,:,i),Nsites,Ncell)/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    endif
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_tridiag_rank7
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
  !
  ! G/Sigma Shape: [2][N,N][:]
  !##################################################################
  subroutine get_gloc_superc_hk_rank2(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:,:),intent(in)    :: Hk
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:),intent(inout) :: Gloc
    character(len=*)                            :: axis
    complex(8),dimension(:,:,:,:),allocatable   :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Lfreq = size(Sigma,4)
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][N,N], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Ntot,Lk],'get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Ntot,Ntot,Lfreq],'get_gloc_superc',"Gloc")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gloc=zero
    if(Lfreq>=Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Z             = sc_zi_rank2(i,Sigma,axis)
             Gloc(:,:,:,i) = Gloc(:,:,:,i) + invert_gki_superc_rank2(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Z             = sc_zi_rank2(i,Sigma,axis)
             Gloc(:,:,:,i) = Gloc(:,:,:,i) + invert_gki_superc_rank2(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    end if
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_hk_rank2


  subroutine get_gloc_superc_dos_rank2(Ebands,Dbands,Hloc,Gloc,Sigma,axis)
    real(8),dimension(:,:,:),intent(in)            :: Ebands    ![2,N,Lk]
    real(8),dimension(:,:),intent(in)              :: Dbands    ![N,Lk]
    real(8),dimension(2,size(Ebands,2)),intent(in) :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:),intent(in)       :: Sigma     ![2,N,N,Lfreq]
    complex(8),dimension(:,:,:,:),intent(inout)    :: Gloc      !as Sigma
    character(len=*)                               :: axis
    complex(8),dimension(:,:,:,:),allocatable      :: Z ![2][2][N,N]
    complex(8),dimension(:,:),allocatable          :: G
    integer                                        :: N,N2
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,2)
    Lk    = size(Ebands,3)
    Nso   = size(Sigma,2)
    Lfreq = size(Sigma,4)
    !
    if(Nso/=Ntot)stop "get_gloc_superc_dos_rank2 error: Ntot != Nso"
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][N,N], Axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm. WARNING: only diagonal case.. "

    if(mpi_master)call start_timer
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[2,Ntot,Lk],'get_gloc_superc_dos_rank2',"Ebands")
    call assert_shape(Dbands,[Ntot,Lk],'get_gloc_superc_dos_rank2',"Dbands")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gloc_superc_dos_rank2',"Sigma")
    call assert_shape(Gloc,[2,Ntot,Ntot,Lfreq],'get_gloc_superc_dos_rank2',"Gloc")
    !
    call build_frequency_array(axis)
    !
    !invert (Z-Hk) for each k-point
    N  = Ntot
    N2 = 2*Ntot
    allocate(G(N2,N2))
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gloc = zero
    do ik = 1+mpi_rank, Lk, mpi_size
       do i=1,Lfreq
          Z  = sc_zi_rank2(i,Sigma,axis)
          G(1:N,1:N)       = Z(1,1,:,:) - diag(Hloc(1,:)) - diag(Ebands(1,:,ik))
          G(1:N,N+1:N2)    = Z(1,2,:,:)        
          G(N+1:N2,1:N)    = Z(2,1,:,:)
          G(N+1:N2,N+1:N2) = Z(2,2,:,:) + diag(Hloc(2,:)) - diag(Ebands(2,:,ik))
          call inv(G)
          do concurrent(io=1:Ntot,jo=1:Ntot)
             Gloc(1,io,jo,i) = Gloc(1,io,jo,i) + G(io,  jo)*Dbands(io,ik)
             Gloc(2,io,jo,i) = Gloc(2,io,jo,i) + G(io,N+jo)*Dbands(io,ik)
          enddo
       enddo
       call eta(ik,Lk)
    enddo
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_dos_rank2

  !G/Sigma shape: [2][Nlat,Nso,Nso][:]
  !##################################################################
  subroutine get_gloc_superc_hk_rank3(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:,:),intent(in)      :: Hk
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                              :: axis
    complex(8),dimension(:,:,:,:),allocatable     :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nlat  = size(Sigma,2)
    Nso   = size(Sigma,3)
    Lfreq = size(Sigma,5)
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [Nlat,Nso,Nso], Axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Ntot,Lk],'get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nso,Nso,Lfreq],'get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Nlat,Nso,Nso,Lfreq],'get_gloc_superc',"Gloc")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gloc=zero
    if(Lfreq>=Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Z               = sc_zi_rank3(i,Sigma,axis)
             Gloc(:,:,:,:,i) = Gloc(:,:,:,:,i) + invert_gki_superc_rank3(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Z               = sc_zi_rank3(i,Sigma,axis)
             Gloc(:,:,:,:,i) = Gloc(:,:,:,:,i) + invert_gki_superc_rank3(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    end if
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_hk_rank3

  !G/Sigma shape: [2][Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_superc_hk_rank4(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:,:),intent(in)        :: Hk
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Ntot,Lk],'get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],'get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Nspin,Nspin,Norb,Norb,Lfreq],'get_gloc_superc',"Gloc")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gloc=zero
    if(Lfreq>=Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Z                 = sc_zi_rank4(i,Sigma,axis)
             Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + invert_gki_superc_rank4(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Z                 = sc_zi_rank4(i,Sigma,axis)
             Gloc(:,:,:,:,:,i) = Gloc(:,:,:,:,:,i) + invert_gki_superc_rank4(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    end if
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_hk_rank4



  subroutine get_gloc_superc_dos_rank4(Ebands,Dbands,Hloc,Gloc,Sigma,axis)
    real(8),dimension(:,:,:),intent(in)             :: Ebands    ![2,N,Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands    ![N,Lk]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2,N,N,Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: Z ![2][2][N,N]
    complex(8),dimension(:,:),allocatable           :: G
    integer                                         :: N,N2
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Ntot  = size(Ebands,2)
    Lk    = size(Ebands,3)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !
    if(Nso/=Ntot)stop "get_gloc_superc_rank4 error: Ntot != Nso"
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"DOS integration algorithm. WARNING: only diagonal case.."
    if(mpi_master)call start_timer
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[2,Ntot,Lk],'get_gloc_superc_dos_rank4',"Ebands")
    call assert_shape(Dbands,[Ntot,Lk],'get_gloc_superc_dos_rank4',"Dbands")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_dos_rank4","Sigma")
    call assert_shape(Gloc, [2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_dos_rank4","Gloc")
    !
    call build_frequency_array(axis)
    !
    !invert (Z-Hk) for each k-point
    N  = Ntot
    N2 = 2*Ntot
    allocate(G(N2,N2))
    !
    Gloc = zero
    do ik = 1+mpi_rank, Lk, mpi_size
       do i=1,Lfreq
          Z  = sc_zi_rank4(i,Sigma,axis)
          G(1:N,1:N)       = Z(1,1,:,:) - diag(Hloc(1,:)) - diag(Ebands(1,:,ik))
          G(1:N,N+1:N2)    = Z(1,2,:,:)        
          G(N+1:N2,1:N)    = Z(2,1,:,:)
          G(N+1:N2,N+1:N2) = Z(2,2,:,:) + diag(Hloc(2,:)) - diag(Ebands(2,:,ik))
          call inv(G)
          do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
             io = iorb + (ispin-1)*Norb
             jo = jorb + (jspin-1)*Norb
             Gloc(1,ispin,jspin,iorb,jorb,i) = Gloc(1,ispin,jspin,iorb,jorb,i) + G(io,  jo)*Dbands(io,ik)
             Gloc(2,ispin,jspin,iorb,jorb,i) = Gloc(2,ispin,jspin,iorb,jorb,i) + G(io,N+jo)*Dbands(io,ik)
          enddo
       enddo
       call eta(ik,Lk)
    enddo
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_dos_rank4

  !G/Sigma shape: [2][Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_superc_hk_rank5(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:,:),allocatable         :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nlat = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nlat,Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Ntot,Lk],'get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gloc_superc',"Gloc")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gloc=zero
    if(Lfreq>=Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Z                 = sc_zi_rank5(i,Sigma,axis)
             Gloc(:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,i) + invert_gki_superc_rank5(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Z                 = sc_zi_rank5(i,Sigma,axis)
             Gloc(:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,i) + invert_gki_superc_rank5(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    end if
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_hk_rank5


#if __GFORTRAN__ &&  __GNUC__ > 8
  !Sigma shape: [2][Nlat,Nspin,Nspin,Norb,Norb][:]
  !G     shape: [2][Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine get_gloc_superc_hk_rank5_6(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:,:),intent(in)            :: Hk
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:,:),allocatable           :: Z
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nlat = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    !
    if(mpi_master)write(*,"(A)")"Get Nambu Green's function [2][Nlat,Nspin,Nspin,Norb,Norb], axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    !Testing part:  
    call assert_shape(Hk,[2,Ntot,Ntot,Lk],'get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'get_gloc_superc',"Gloc")
    !
    call build_frequency_array(axis)
    !
    allocate(Z(2,2,Ntot,Ntot))
    !
    Gloc=zero
    if(Lfreq>=Lk)then
       do i=1+mpi_rank,Lfreq,mpi_size
          do ik=1,Lk
             Z                     = sc_zi_rank5(i,Sigma,axis)
             Gloc(:,:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,:,i) + invert_gki_superc_rank6(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(i,Lfreq)
       end do
    else
       do ik=1+mpi_rank,Lk,mpi_size
          do i=1,Lfreq
             Z                     = sc_zi_rank5(i,Sigma,axis)
             Gloc(:,:,:,:,:,:,:,i) = Gloc(:,:,:,:,:,:,:,i) + invert_gki_superc_rank6(Z,Hk(:,:,:,ik))/Lk
          enddo
          if(mpi_master)call eta(ik,Lk)
       end do
    end if
    call gf_reduce(Gloc)
    !
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_hk_rank5_6
#endif  


end module GF_GLOC
