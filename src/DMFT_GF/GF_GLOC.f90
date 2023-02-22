module GF_GLOC
  USE GF_COMMON
  implicit none
  private


  public :: get_gloc_normal_main
  public :: get_gloc_normal_tridiag
  public :: get_gloc_normal_dos
  !
  public :: get_gloc_superc_main
  public :: get_gloc_superc_dos
  !
  public :: get_gloc_normal_hk_rank4
  public :: get_gloc_normal_dos_rank4
  public :: get_gloc_normal_hk_rank5
  public :: get_gloc_normal_tridiag_rank5
  public :: get_gloc_normal_hk_rank5_6
  public :: get_gloc_normal_hk_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: get_gloc_normal_hk_rank7
  public :: get_gloc_normal_tridiag_rank7
#endif
  !
  public :: get_gloc_superc_hk_rank4
  public :: get_gloc_superc_dos_rank4
  public :: get_gloc_superc_hk_rank5
  public :: get_gloc_superc_hk_rank5_6







contains



  !##################################################################
  !##################################################################
  !                       NORMAL PHASE
  !##################################################################
  !##################################################################
  subroutine get_gloc_normal_main(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:),intent(in)    :: Hk        ![N,N][Lk]
    complex(8),dimension(:,:,:),intent(in)    :: Sigma     ![N,N][Lfreq]
    complex(8),dimension(:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                          :: axis
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable   :: Gk    !as Sigma  
    complex(8),dimension(:,:,:),allocatable   :: Gtmp  !as Sigma
    complex(8),dimension(:,:,:),allocatable   :: csi   ![N,N][Lfreq]
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    ! !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Lfreq = size(Sigma,3)
    !Testing part:
    call assert_shape(Hk,[Ntot,Ntot,Lk],"get_gloc_normal_main","Hk")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"get_gloc_normal_main","Sigma")
    call assert_shape(Gloc, [Ntot,Ntot,Lfreq],"get_gloc_normal_main","Gloc")
    !
    call build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gk(Ntot,Ntot,Lfreq))
    allocate(csi(Ntot,Ntot,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    if(Lfreq >= Lk)then
       Gloc = zero
       do ik=1,Lk
          call invert_gk_normal_mpi(csi,Hk(:,:,ik),Gk)      
          Gloc = Gloc + Gk/Lk
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Ntot,Ntot,Lfreq));Gtmp=zero
       Gloc = zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal(csi,Hk(:,:,ik),Gk)
          Gtmp = Gtmp + Gk/Lk
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          call Mpi_AllReduce(Gtmp,Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gloc=Gtmp
       endif
#else
       Gloc=Gtmp
#endif
    endif
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_main



  !TRIDIAG Hk
  !Hk has a blocks tridiagonal form with blocks of size [Ncell,Ncell]. Ntot=Nsites*Ncell,
  !with: Nsites the number of sites in the super-cell
  !      Ncell the dimension of the unit cell 
  subroutine get_gloc_normal_tridiag(Hk,Gloc,Sigma,axis,Nsites,Ncell)
    complex(8),dimension(:,:,:),intent(in)          :: Hk      ![N,N,Nk]
    complex(8),dimension(:,:,:),intent(in)          :: Sigma   ![N,N,Lfreq]
    complex(8),dimension(:,:,:),intent(inout)       :: Gloc    !as Sigma
    character(len=*)                                :: axis
    integer,intent(in)                              :: Nsites,Ncell
    !allocatable arrays
    complex(8),dimension(:,:,:),allocatable         :: Gk
    complex(8),dimension(:,:,:),allocatable         :: Gtmp
    complex(8),dimension(:,:,:),allocatable         :: csi
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Lfreq = size(Sigma,3)
    !Testing part:
    if(Nsites*Ncell/=Ntot)stop "get_gloc_normal_tridiag ERROR: passed Nsites*Ncell != size(Hk,1)==Ntot"
    call assert_shape(Hk,[Ntot,Ntot,Lk],'get_gloc_normal_tridiag',"Hk")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],'get_gloc_normal_tridiag',"Sigma")
    call assert_shape(Gloc,[Ntot,Ntot,Lfreq],'get_gloc_normal_tridiag',"Gloc")
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)write(*,"(A)")"Block Tridiagonal Gaussian elimination algorithm:"
    !
    allocate(Gk(Ntot,Ntot,Lfreq))
    allocate(csi(Ntot,Ntot,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i) = (wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
    enddo
    !
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_tridiag_mpi(csi,Hk(:,:,ik),Gk,Nsites,Ncell)
          Gloc = Gloc + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Ntot,Ntot,Lfreq));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_tridiag(csi,Hk(:,:,ik),Gk,Nsites,Ncell)
          Gtmp = Gtmp + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gloc=zero
          call Mpi_AllReduce(Gtmp,Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gloc=Gtmp
       endif
#else
       Gloc=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_tridiag



  !DOS case: Hk--> Es,DOSs
  subroutine get_gloc_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,axis)
    real(8),dimension(:,:),intent(in)            :: Ebands    ![N][Lk]
    real(8),dimension(:,:),intent(in)            :: Dbands    ![N][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in) :: Hloc      ![N]
    complex(8),dimension(:,:,:),intent(in)       :: Sigma     ![N,N][Lmats]
    complex(8),dimension(:,:,:),intent(inout)    :: Gloc      !as Sigma
    character(len=*)                             :: axis
    !
    complex(8)                                   :: gktmp
    complex(8),dimension(:,:,:),allocatable      :: csi ![N,N][Lmats]
    complex(8),dimension(:,:,:),allocatable      :: Gtmp      !as Sigma
    !
    !New
    complex(8),dimension(:,:),allocatable        :: Gdos_tmp ![N][N]
    logical                                      :: dos_diag !1. T / 2. F
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    Ntot  = size(Ebands,1)
    Lk    = size(Ebands,2)
    Lfreq = size(Sigma,3)
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    dos_diag = .not.(size(Dbands,1) < size(Ebands,1))
    !
    !Testing part:
    call assert_shape(Ebands,[Ntot,Lk],"dmft_get_gloc_normal_dos","Ebands")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"dmft_get_gloc_normal_dos","Sigma")
    call assert_shape(Gloc,[Ntot,Ntot,Lfreq],"dmft_get_gloc_normal_dos","Gloc")
    !
    call build_frequency_array(axis)
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Ntot,Ntot,Lfreq))
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Ntot) - Sigma(:,:,i)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    allocate(Gtmp(Ntot,Ntot,Lfreq));Gtmp=zero
    !
    ! diagonal case
    if(dos_diag)then
       !
       if(Lfreq>=Lk)then
          do i=1+mpi_rank, Lfreq, mpi_size
             do io=1,Ntot
                do ik=1,Lk
                   gktmp = Dbands(io,ik)/( csi(io,io,i)-Hloc(io)-Ebands(io,ik) )
                   Gtmp(io,io,i) = Gtmp(io,io,i) + gktmp
                enddo
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do io=1,Ntot
                do i=1,Lfreq
                   gktmp = Dbands(io,ik)/( csi(io,io,i)-Hloc(io)-Ebands(io,ik) )
                   Gtmp(io,io,i) = Gtmp(io,io,i) + gktmp
                enddo
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if
       !
    else
       !
       allocate(Gdos_tmp(Ntot,Ntot)) ;Gdos_tmp=zero
       if(Lfreq>=Lk)then
          do i = 1+mpi_rank, Lfreq, mpi_size                                 !MPI loop over Matsubara frequencies
             do ik=1,Lk                                                      !for all e-value (here named ik)
                Gdos_tmp = csi(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = csi - Hloc - H(e)
                call inv(Gdos_tmp)
                Gtmp(:,:,i) = Gtmp(:,:,i) + Dbands(1,ik)*Gdos_tmp
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                Gdos_tmp = csi(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = csi - Hloc - H(e)
                call inv(Gdos_tmp)
                Gtmp(:,:,i) = Gtmp(:,:,i) + Dbands(1,ik)*Gdos_tmp
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if
    end if
#ifdef _MPI    
    if(check_MPI())then
       Gloc=zero
       call Mpi_AllReduce(Gtmp,Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gloc=Gtmp
    endif
#else
    Gloc=Gtmp
#endif
    if(mpi_master)call stop_timer
  end subroutine get_gloc_normal_dos





  !##################################################################
  !##################################################################
  !                          SUPERC PHASE
  !##################################################################
  !##################################################################
  subroutine get_gloc_superc_main(Hk,Gloc,Sigma,axis)
    complex(8),dimension(:,:,:,:),intent(in)                    :: Hk       ![2][N,N,Nk]
    complex(8),dimension(:,:,:,:),intent(in)                    :: Sigma    ![2][N,N,Lfreq]
    complex(8),dimension(:,:,:,:),intent(inout)                 :: Gloc     !as Sigma
    character(len=*)                                            :: axis
    !allocatable arrays
    complex(8),dimension(2,size(Hk,2),size(Hk,2),size(Sigma,4)) :: Gk      !as Sigma
    complex(8),dimension(2,size(Hk,2),size(Hk,2),size(Sigma,4)) :: Gtmp    !as Sigma
    complex(8),dimension(:,:,:,:,:),allocatable                 :: csi     ![2][2][N][N][Lfreq]
    !
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Lfreq = size(Sigma,4)
    !
    call assert_shape(Hk,[2,Ntot,Ntot,Lk],'get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Ntot,Ntot,Lfreq],'get_gloc_superc',"Gloc")
    !
    call build_frequency_array(axis)
    !        
    allocate(csi(2,2,Ntot,Ntot,Lfreq))
    !
    select case(axis)
    case default; stop "dmft_get_gloc_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Ntot) -        Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                          -        Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                          -        Sigma(2,:,:,i)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Ntot) + conjg( Sigma(1,:,:,i))
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i) + xmu)*eye(Ntot)                - Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                                           - Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                                           - Sigma(2,:,:,i)
          csi(2,2,:,:,i) = -conjg(wfreq(Lfreq+1-i) + xmu)*eye(Ntot)  + conjg( Sigma(1,:,:,Lfreq+1-i) )
       enddo
    end select
    !
    !
    if(mpi_master)write(*,"(A)")"Get local Nambu Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    if(Lfreq>=Lk)then
       Gloc=zero
       do ik=1,Lk
          call invert_gk_superc_mpi(csi,Hk(:,:,:,ik),Gk)
          Gloc = Gloc + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc(csi,Hk(:,:,:,ik),Gk)
          Gtmp = Gtmp + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gloc=zero
          call Mpi_AllReduce(Gtmp,Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gloc=Gtmp
       endif
#else
       Gloc=Gtmp
#endif
    end if
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_main




  subroutine get_gloc_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,axis)
    real(8),dimension(:,:,:),intent(in)             :: Ebands    ![2,N,Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands    ![N,Lk]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:),intent(in)        :: Sigma     ![2,N,N,Lfreq]
    complex(8),dimension(:,:,:,:),intent(inout)     :: Gloc      !as Sigma
    character(len=*)                                :: axis
    ! arrays
    complex(8),dimension(:,:,:,:,:),allocatable     :: csi ![2][2][N,N,Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix 
    complex(8),dimension(:,:,:,:),allocatable       :: Gtmp    !as Sigma
    complex(8),dimension(:,:),allocatable           :: Gdos_tmp ![N][N]
    logical                                         :: dos_diag !1. T / 2. F
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !
    Ntot  = size(Ebands,2)
    Lk    = size(Ebands,3)
    Lfreq = size(Sigma,4)
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[2,Ntot,Lk],'get_gloc_superc_dos',"Ebands")
    call assert_shape(Dbands,[Ntot,Lk],'get_gloc_superc_dos',"Dbands")
    call assert_shape(Sigma,[2,Ntot,Ntot,Lfreq],'get_gloc_superc_dos',"Sigma")
    call assert_shape(Gloc,[2,Ntot,Ntot,Lfreq],'get_gloc_superc_dos',"Gloc")
    !
    call build_frequency_array(axis)
    !
    !
    allocate(csi(2,2,Ntot,Ntot,Lfreq));csi = zero
    select case(axis)
    case default; stop "get_gloc_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Ntot) - diag(Hloc(1,:)) -        Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                                            -        Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                                            -        Sigma(2,:,:,i)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Ntot) + diag(Hloc(2,:)) + conjg( Sigma(1,:,:,i) )
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Ntot)                 - diag(Hloc(1,:)) -        Sigma(1,:,:,i)
          csi(1,2,:,:,i) =                                                            -        Sigma(2,:,:,i)
          csi(2,1,:,:,i) =                                                            -        Sigma(2,:,:,i)
          csi(2,2,:,:,i) = -conjg( wfreq(Lfreq+1-i)+xmu )*eye(Ntot) + diag(Hloc(2,:)) + conjg( Sigma(1,:,:,Lfreq+1-i) )
       enddo
    end select
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"WARNING: This method is limited to diagonal case.. "
    if(mpi_master)write(*,"(A)")"Get local Nambu Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    allocate(Gtmp(2,Ntot,Ntot,Lfreq))
    allocate(Gmatrix(2*Ntot,2*Ntot))
    Gtmp=zero
    !
    !
    do ik = 1+mpi_rank, Lk, mpi_size
       do i=1,Lfreq
          Gmatrix  = zero
          Gmatrix(1:Ntot,1:Ntot)               = csi(1,1,:,:,i) - diag(Ebands(1,:,ik))
          Gmatrix(1:Ntot,Ntot+1:2*Ntot)        = csi(1,2,:,:,i)
          Gmatrix(Ntot+1:2*Ntot,1:Ntot)        = csi(2,1,:,:,i)
          Gmatrix(Ntot+1:2*Ntot,Ntot+1:2*Ntot) = csi(2,2,:,:,i) - diag(Ebands(2,:,ik))
          call inv(Gmatrix)
          do io=1,Ntot
             do jo=1,Ntot
                Gtmp(1,io,jo,i) = Gtmp(1,io,jo,i) + Gmatrix(io,jo)*Dbands(io,ik)
                Gtmp(2,io,jo,i) = Gtmp(2,io,jo,i) + Gmatrix(io,Ntot+jo)*Dbands(io,ik)
             end do
          enddo
       enddo
       call eta(ik,Lk)
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Gloc=zero
       call Mpi_AllReduce(Gtmp,Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gloc=Gtmp
    endif
#else
    Gloc=Gtmp
#endif
    if(mpi_master)call stop_timer
  end subroutine get_gloc_superc_dos











  !##################################################################
  !                        AUX INTERFACES:
  ! 1. rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] 
  ! 2. rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_unit_cell
  ! 3. rank-6 CDMFT   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  ! 4. rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nineq=N_unit_cell,Nlat=N_cluster
  !##################################################################
  !##################################################################
  !##################################################################
  !                         NORMAL
  !##################################################################
  !##################################################################
  subroutine get_gloc_normal_hk_rank4(Hk,Gloc,Sigma,axis) !N=Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:,:,:),allocatable       :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Hk,[Ntot,Nso,Lk],"get_gloc_normal_rank4","Hk")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank4","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    call get_gloc_normal_main(Hk,GF,SF,axis)
    Gloc = reshape_matrix_to_rank4(GF,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_hk_rank4

  subroutine get_gloc_normal_dos_rank4(Ebands,Dbands,Hloc,Gloc,Sigma,axis)!N=Nspin*Norb
    real(8),dimension(:,:),intent(in)             :: Ebands    ![N][Lk]
    real(8),dimension(:,:),intent(in)             :: Dbands    ![N][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)  :: Hloc      ![N]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:,:,:),allocatable       :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    Ntot  = size(Ebands,1)
    Lk    = size(Ebands,2)
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank4","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    call get_gloc_normal_dos(Ebands,Dbands,Hloc,GF,SF,axis)
    Gloc = reshape_matrix_to_rank4(GF,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_dos_rank4

  subroutine get_gloc_normal_hk_rank5(Hk,Gloc,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:),allocatable         :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gloc_normal_main(Hk,GF,SF,axis)
    Gloc = reshape_matrix_to_rank5(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_hk_rank5

  subroutine get_gloc_normal_tridiag_rank5(Hk,Gloc,Sigma,axis,Nsites,Ncell) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                :: axis
    integer,intent(in)                              :: Nsites,Ncell
    complex(8),dimension(:,:,:),allocatable         :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_rank5","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gloc_normal_tridiag(Hk,GF,SF,axis,Nsites,Ncell)
    Gloc = reshape_matrix_to_rank5(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_tridiag_rank5

  subroutine get_gloc_normal_hk_rank5_6(Hk,Gloc,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc      ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:),allocatable           :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_rank5_6","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5_6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank5_6","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gloc_normal_main(Hk,GF,SF,axis)
    Gloc = reshape_matrix_to_rank6(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_hk_rank5_6


  subroutine get_gloc_normal_hk_rank6(Hk,Gloc,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)            :: Hk     ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc   !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:),allocatable           :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin]"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_rank6","Hk")
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank6","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank6_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    call get_gloc_normal_main(Hk,GF,SF,axis)
    Gloc = reshape_matrix_to_rank6(GF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_hk_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine get_gloc_normal_hk_rank7(Hk,Gloc,Sigma,axis) !N=Nineq*Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)              :: Hk     ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc   !as Sigma
    character(len=*)                                    :: axis
    complex(8),dimension(:,:,:),allocatable             :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nlso  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin,Nineq]"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank7","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank7_to_matrix(Sigma,Nineq,Nlat,Nspin,Norb,Lfreq)
    call get_gloc_normal_main(Hk,GF,SF,axis)
    Gloc = reshape_matrix_to_rank7(GF,Nineq,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_hk_rank7

  subroutine get_gloc_normal_tridiag_rank7(Hk,Gloc,Sigma,axis,Nsites,Ncell) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:),intent(in)              :: Hk     ![N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc   !as Sigma
    character(len=*)                                    :: axis
    integer,intent(in)                                  :: Nsites,Ncell
    complex(8),dimension(:,:,:),allocatable             :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,1)
    Lk    = size(Hk,3)
    Nineq = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Nlso  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin,Nineq]"
    call assert_shape(Hk,[Ntot,Nlso,Lk],"get_gloc_normal_rank7","Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_normal_rank7","Gloc")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq))
    SF=zero
    GF=zero
    !
    SF   = reshape_rank7_to_matrix(Sigma,Nineq,Nlat,Nspin,Norb,Lfreq)
    call get_gloc_normal_tridiag(Hk,GF,SF,axis,Nsites,Ncell)
    Gloc = reshape_matrix_to_rank7(GF,Nineq,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_normal_tridiag_rank7
#endif







  !##################################################################
  !##################################################################
  !                         SUPERC
  !##################################################################
  !##################################################################

  subroutine get_gloc_superc_hk_rank4(Hk,Gloc,Sigma,axis) !N=Nspin*Norb
    complex(8),dimension(:,:,:,:),intent(in)        :: Hk        ![2][N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: SF,GF ![2][N,N,Lfreq]
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Hk,[2,Ntot,Nso,Lk],"get_gloc_superc_rank4","Hk")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank4","Sigma")
    call assert_shape(Gloc, [2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank4","Gloc")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank4_to_matrix(Sigma(1,:,:,:,:,:),Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank4_to_matrix(Sigma(2,:,:,:,:,:),Nspin,Norb,Lfreq)
    call get_gloc_superc_main(Hk,GF,SF,axis)
    Gloc(1,:,:,:,:,:) = reshape_matrix_to_rank4(GF(1,:,:,:),Nspin,Norb,Lfreq)
    Gloc(2,:,:,:,:,:) = reshape_matrix_to_rank4(GF(2,:,:,:),Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_superc_hk_rank4

  subroutine get_gloc_superc_dos_rank4(Ebands,Dbands,Hloc,Gloc,Sigma,axis)!N=Nspin*Norb
    real(8),dimension(:,:,:),intent(in)             :: Ebands   ![2][N,Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands    ![N][Lk] /[1][Lk]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2,N]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:,:,:,:),allocatable       :: SF,GF ![2][N,N,Lfreq]
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    Ntot  = size(Ebands,2)
    Lk    = size(Ebands,3)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank4","Sigma")
    call assert_shape(Gloc, [2,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank4","Gloc")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank4_to_matrix(Sigma(1,:,:,:,:,:),Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank4_to_matrix(Sigma(2,:,:,:,:,:),Nspin,Norb,Lfreq)
    call get_gloc_superc_dos(Ebands,Dbands,Hloc,GF,SF,axis)
    Gloc(1,:,:,:,:,:) = reshape_matrix_to_rank4(GF(1,:,:,:),Nspin,Norb,Lfreq)
    Gloc(2,:,:,:,:,:) = reshape_matrix_to_rank4(GF(2,:,:,:),Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_superc_dos_rank4


  subroutine get_gloc_superc_hk_rank5(Hk,Gloc,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF ![2][N,N,Lfreq]
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Hk,[2,Ntot,Nso,Lk],"get_gloc_superc_rank5","Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank5","Sigma")
    call assert_shape(Gloc, [2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank5","Gloc")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank5_to_matrix(Sigma(1,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank5_to_matrix(Sigma(2,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    call get_gloc_superc_main(Hk,GF,SF,axis)
    Gloc(1,:,:,:,:,:,:) = reshape_matrix_to_rank5(GF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Gloc(2,:,:,:,:,:,:) = reshape_matrix_to_rank5(GF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_superc_hk_rank5


  subroutine get_gloc_superc_hk_rank5_6(Hk,Gloc,Floc,Sigma,axis) !N=Nlat*Nspin*Norb
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2,N,N,Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc      ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Floc      ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    character(len=*)                                  :: axis
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Ntot  = size(Hk,2)
    Lk    = size(Hk,4)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Nlso  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Hk,[2,Ntot,Nlso,Lk],"get_gloc_superc_rank5_6","Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank5_6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank5_6","Gloc")
    call assert_shape(Floc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"get_gloc_superc_rank5_6","Floc")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq))
    SF=zero; GF=zero
    !
    SF(1,:,:,:) = reshape_rank5_to_matrix(Sigma(1,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank5_to_matrix(Sigma(2,:,:,:,:,:,:),Nlat,Nspin,Norb,Lfreq)
    call get_gloc_superc_main(Hk,GF,SF,axis)
    Gloc = reshape_matrix_to_rank6(GF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Floc = reshape_matrix_to_rank6(GF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF)
  end subroutine get_gloc_superc_hk_rank5_6



end module GF_GLOC
