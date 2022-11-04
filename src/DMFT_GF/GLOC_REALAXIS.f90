module GLOC_REALAXIS
  USE DMFT_CTRL_VARS  
  USE GF_COMMON
  implicit none
  private


  interface get_gloc_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_main
     module procedure :: dmft_get_gloc_realaxis_normal_cluster
     module procedure :: dmft_get_gloc_realaxis_normal_dos
     module procedure :: dmft_get_gloc_realaxis_normal_ineq
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_realaxis_normal_cluster_ineq
#endif
     module procedure :: dmft_get_gloc_realaxis_superc_main
     module procedure :: dmft_get_gloc_realaxis_superc_dos
     module procedure :: dmft_get_gloc_realaxis_superc_ineq
  end interface get_gloc_realaxis

  interface get_gij_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_gij
     module procedure :: dmft_get_gloc_realaxis_superc_gij
  end interface get_gij_realaxis


  interface dmft_gloc_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_main
     module procedure :: dmft_get_gloc_realaxis_normal_cluster
     module procedure :: dmft_get_gloc_realaxis_normal_dos
     module procedure :: dmft_get_gloc_realaxis_normal_ineq
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_realaxis_normal_cluster_ineq
#endif
     module procedure :: dmft_get_gloc_realaxis_superc_main
     module procedure :: dmft_get_gloc_realaxis_superc_dos
     module procedure :: dmft_get_gloc_realaxis_superc_ineq
  end interface dmft_gloc_realaxis


  interface dmft_gij_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_gij
     module procedure :: dmft_get_gloc_realaxis_superc_gij
  end interface dmft_gij_realaxis


  !PUBLIC IN DMFT:
  public :: get_gloc_realaxis
  public :: get_gij_realaxis

  public :: dmft_gloc_realaxis
  public :: dmft_gij_realaxis


contains



  subroutine dmft_get_gloc_realaxis_normal_main(hk,Greal,Sreal,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkreal    !as Sreal  
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal] 
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],'dmft_get_gloc_realaxis_normal_main_mpi',"Hk")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_main_mpi',"Sreal")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_main_mpi',"Greal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(Gkreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_mpi(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
          Greal = Greal + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_normal_main


  subroutine dmft_get_gloc_realaxis_normal_cluster(hk,Greal,Sreal,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
    complex(8),dimension(:,:,:),allocatable           :: zeta_real ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,3)
    Nlso   = Nlat*Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_cluster_mpi',"Hk")
    call assert_shape(Sreal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_cluster_mpi',"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_cluster_mpi',"Greal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nlso,Nlso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Sreal(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_cluster_mpi(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
          Greal = Greal + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_cluster(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif

    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_normal_cluster


  subroutine dmft_get_gloc_realaxis_normal_dos(Ebands,Dbands,Hloc,Greal,Sreal)
    real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)                           :: Dbands    !1. [Nspin*Norb][Lk] / 2. [1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)                  :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout)               :: Greal     !as Sreal
    !
    complex(8)                                                  :: gktmp
    complex(8),dimension(:,:,:),allocatable                     :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),allocatable                 :: Gtmp      !as Sreal
    !
    real(8)                                                     :: beta
    real(8)                                                     :: xmu,eps
    real(8)                                                     :: wini,wfin
    !
    complex(8),dimension(:,:),allocatable                         :: Gdos_tmp ![Nspin*Norb][Nspin*Norb]
    logical                                                       :: dos_diag !1. T / 2. F
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
    !Retrieve parameters:
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb
    !
    !case F => one DOS, H(e)=diag(Ebands), non-diag
    !case T => more DOS, Ebands, diag
    dos_diag = .not.( size(Dbands,1) < size(Ebands,1) )
    !
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_realaxis_normal_dos_main_mpi',"Ebands")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos_main_mpi',"Sreal")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos_main_mpi',"Greal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
    if(mpi_master)call start_timer
    Greal=zero
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
    !
    ! diagonal case
    if(dos_diag)then
       if(Lmats>=Lk)then
          do i = 1+mpi_rank, Lreal, mpi_size
             do ispin=1,Nspin
                do iorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   do ik=1,Lk
                      gktmp = Dbands(io,ik)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io,ik) )
                      Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
                   enddo
                enddo
             enddo
             if(mpi_master)call eta(i,Lreal)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do ispin=1,Nspin
                do iorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   do i=1,Lreal
                      gktmp = Dbands(io,ik)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io,ik) )
                      Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
                   enddo
                enddo
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if
       ! non diagonal case
    else
       allocate(Gdos_tmp(Nso,Nso)) ;Gdos_tmp=zero
       if(Lmats>=Lk)then
          do i = 1+mpi_rank, Lreal, mpi_size                                 !MPI loop over real frequencies
             do ik=1,Lk                                                      !for all e-value (here named ik)
                Gdos_tmp = zeta_real(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = zeta_real - Hloc - H(e)
                call inv(Gdos_tmp)
                Gtmp(:,:,:,:,i) = Gtmp(:,:,:,:,i) + Dbands(1,ik)*so2nn_reshape(Gdos_tmp,Nspin,Norb)
             enddo
             if(mpi_master)call eta(i,Lreal)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lreal
                Gdos_tmp = zeta_real(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = zeta_real - Hloc - H(e)
                call inv(Gdos_tmp)
                Gtmp(:,:,:,:,i) = Gtmp(:,:,:,:,i) + Dbands(1,ik)*so2nn_reshape(Gdos_tmp,Nspin,Norb)
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if
    end if
#ifdef _MPI    
    if(check_MPI())then
       Greal=zero
       call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
    else
       Greal=Gtmp
    endif
#else
    Greal=Gtmp
#endif
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_normal_dos


  subroutine dmft_get_gloc_realaxis_normal_ineq(hk,Greal,Sreal,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin  
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_ineq_main_mpi',"Sreal")
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_ineq_main_mpi',"Greal")
    !
    if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gkreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_ineq_mpi(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Greal = Greal + Gkreal/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_tridiag_mpi(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Greal = Greal + Gkreal/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
    else
       allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_ineq(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Gtmp = Gtmp + Gkreal/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
             Gtmp = Gtmp + Gkreal/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_normal_ineq


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gloc_realaxis_normal_cluster_ineq(Hk,Greal,Sreal,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)               :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)     :: Sreal     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout)  :: Greal     !as Sreal
    logical,optional                                     :: tridiag
    logical                                              :: tridiag_
    logical,dimension(size(Hk,3)),optional               :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                        :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable    :: Gkmats    !as Sreal
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable    :: Gtmp      !as Sreal
    complex(8),dimension(:,:,:,:),allocatable            :: zeta_mats ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nineq = size(Sreal,1)
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,4)
    Norb  = size(Sreal,6)
    Lreal = size(Sreal,8)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    Nilso = Nineq*Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nilso,Nilso,Lk],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Sreal")
    call assert_shape(Greal,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Greal")
    !
    if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gkmats(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_mats(Nineq,Nlso,Nlso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do iineq=1,Nineq
       do i=1,Lreal
          zeta_mats(iineq,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nlso)     - nnn2lso_cluster_reshape(Sreal(iineq,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_cluster_ineq_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Greal = Greal + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_cluster_ineq_tridiag_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Greal = Greal + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
    else
       allocate(Gtmp(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_cluster_ineq(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gtmp = Gtmp + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_cluster_ineq_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gtmp = Gtmp + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_normal_cluster_ineq
#endif




  subroutine dmft_get_gloc_realaxis_normal_gij(hk,Greal,Sreal,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_gij_main_mpi',"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_gij_main_mpi',"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_gij_main_mpi',"Greal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (no print)"
    !
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(zeta_real(Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       do i=1,Lreal
          zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_gij_mpi(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero     
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_normal_gij









  subroutine dmft_get_gloc_realaxis_superc_main(hk,Greal,Sreal,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)        :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nso,Nso,Lk],'dmft_get_gloc_realaxis_superc_main_mpi',"Hk")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_mpi',"Sreal")
    call assert_shape(Greal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_mpi',"Greal")
    !
    allocate(Gkreal(2,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(2,2,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
            nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
            conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Realax8is Superc Green's function (no print)"
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_mpi(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lreal))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_superc_main


  subroutine dmft_get_gloc_realaxis_superc_dos(Ebands,Dbands,Hloc,Greal,Sreal)
    real(8),dimension(:,:),intent(in)                                     :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in)           :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                          :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)                          :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout)                       :: Greal     !as Sreal
    ! arrays
    complex(8),dimension(2,2,size(Ebands,1),size(Ebands,1),size(Sreal,6)) :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(2*size(Ebands,1),2*size(Ebands,1))               :: Gmatrix
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable                         :: Gtmp    !as Sreal
    !
    real(8)                                                               :: beta
    real(8)                                                               :: xmu,eps
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")

    !
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_realaxis_superc_dos',"Ebands")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main',"Sreal")
    call assert_shape(Greal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main',"Greal")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso) - diag(Hloc)   - &
            nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
       zeta_real(1,2,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,1,:,:,i) =                                                  - &
            nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
       zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + diag(Hloc) + &
            conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Realaxis Superc Green's function (no print)"
    if(mpi_master)call start_timer
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lreal))
    Gtmp=zero
    Greal=zero
    do ik = 1+mpi_rank, Lk, mpi_size
       !
       do i=1,Lreal
          Gmatrix  = zero
          Gmatrix(1:Nso,1:Nso)             = zeta_real(1,1,:,:,i) - diag(Ebands(:,ik))
          Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta_real(1,2,:,:,i)
          Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta_real(2,1,:,:,i)
          Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta_real(2,2,:,:,i) + diag(Ebands(:,ik))
          !
          call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
          !
          do ispin=1,Nspin
             do iorb=1,Norb
                io = iorb + (ispin-1)*Norb
                Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + Gmatrix(io,io)*Dbands(io,ik)
                Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + Gmatrix(io,Nso+io)*Dbands(io,ik)
             enddo
          enddo
          !
       enddo
       !
       call eta(ik,Lk)
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Greal=zero
       call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
    else
       Greal=Gtmp
    endif
#else
    Greal=Gtmp
#endif
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_superc_dos


  subroutine dmft_get_gloc_realaxis_superc_ineq(hk,Greal,Sreal,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
    logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_superc_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_ineq_main_mpi',"Sreal")
    call assert_shape(Greal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_ineq_main_mpi',"Greal")
    !
    if(mpi_master)write(*,"(A)")"Get local Realaxis Superc Green's function (no print)"
    !
    allocate(Gkreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do ilat=1,Nlat
       !SYMMETRIES in real-frequencies   [assuming a real order parameter]
       !G22(w)  = -[G11[-w]]*
       !G21(w)  =   G12[w]   
       do i=1,Lreal
          zeta_real(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
               nn2so_reshape(Sreal(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(1,2,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,1,ilat,:,:,i) =                                                  - &
               nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_real(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
               conjg( nn2so_reshape(Sreal(1,ilat,:,:,:,:,Lreal+1-i),Nspin,Norb) )
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Greal=zero
    if(Lreal>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_ineq_mpi(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
          Greal = Greal + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc_ineq(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
          Gtmp = Gtmp + Gkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          call Mpi_AllReduce(Gtmp,Greal,size(Greal),MPI_Double_Complex,MPI_Sum,MPI_COMM_WORLD,MPI_ierr)
       else
          Greal=Gtmp
       endif
#else
       Greal=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_superc_ineq


  subroutine dmft_get_gloc_realaxis_superc_gij(hk,Greal,Freal,Sreal,hk_symm)
    complex(8),dimension(:,:,:,:)                       :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Greal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Freal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sreal           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    logical,optional                                  :: hk_symm(size(Hk,4))
    logical                                           :: hk_symm_(size(Hk,4))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkreal,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    integer                                           :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
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
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,2)
    Nspin = size(Sreal,3)
    Norb  = size(Sreal,5)
    Lreal = size(Sreal,7)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Greal")
    call assert_shape(Freal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Freal")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (no print)"
    !
    allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Fkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Fkreal=zero
    allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal));zeta_real=zero
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    !SYMMETRIES in real-frequencies   [assuming a real order parameter]
    !G22(w)  = -[G11[-w]]*
    !G21(w)  =   G12[w]
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu
             zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu )
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - &
                        Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - &
                        Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - &
                        Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) + &
                        conjg( Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Greal=zero
    Freal=zero
    if(Lreal>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_gij_mpi(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
          Greal = Greal + Gkreal/dble(Lk)
          Freal = Freal + Fkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
       allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Ftmp=zero
       do ik=1,Lk
          call invert_gk_superc_gij(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
          Gtmp = Gtmp + Gkreal/dble(Lk)
          Ftmp = Ftmp + Fkreal/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Greal=zero
          Freal=zero
          call Mpi_AllReduce(Gtmp, Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
          call Mpi_AllReduce(Ftmp, Freal, size(Greal), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Greal=Gtmp
          Freal=Ftmp
       endif
#else
       Greal=Gtmp
       Freal=Ftmp
#endif
       deallocate(Gtmp,Ftmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_realaxis_superc_gij













end module GLOC_REALAXIS
