module GLOC_MATSUBARA
  USE DMFT_CTRL_VARS  
  USE GF_COMMON
  implicit none
  private


  interface get_gloc_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_main
     module procedure :: dmft_get_gloc_matsubara_normal_cluster
     module procedure :: dmft_get_gloc_matsubara_normal_dos
     module procedure :: dmft_get_gloc_matsubara_normal_ineq
     !
     module procedure :: dmft_get_gloc_matsubara_superc_main
     module procedure :: dmft_get_gloc_matsubara_superc_dos
     module procedure :: dmft_get_gloc_matsubara_superc_ineq
  end interface get_gloc_matsubara


  interface dmft_gloc_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_main
     module procedure :: dmft_get_gloc_matsubara_normal_cluster
     module procedure :: dmft_get_gloc_matsubara_normal_dos
     module procedure :: dmft_get_gloc_matsubara_normal_ineq
     !
     module procedure :: dmft_get_gloc_matsubara_superc_main
     module procedure :: dmft_get_gloc_matsubara_superc_dos
     module procedure :: dmft_get_gloc_matsubara_superc_ineq
  end interface dmft_gloc_matsubara




  interface get_gij_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_gij
     module procedure :: dmft_get_gloc_matsubara_superc_gij
  end interface get_gij_matsubara


  interface dmft_gij_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_gij
     module procedure :: dmft_get_gloc_matsubara_superc_gij
  end interface dmft_gij_matsubara



  !PUBLIC IN DMFT:
  public :: get_gloc_matsubara
  public :: dmft_gloc_matsubara

  public :: get_gij_matsubara
  public :: dmft_gij_matsubara






contains


  subroutine dmft_get_gloc_matsubara_normal_main(Hk,Gmats,Smats,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gkmats    !as Smats  
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Smats
    complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
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
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_get_gloc_matsubara_normal_main_mpi","Hk")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_main_mpi","Smats")
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_main_mpi","Gmats")
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gkmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (no print)"
    if(mpi_master)call start_timer
    Gmats=zero
    if(Lmats>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
          Gmats = Gmats + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
       endif
#else
       Gmats=Gtmp
#endif
    endif
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_normal_main


  subroutine dmft_get_gloc_matsubara_normal_cluster(Hk,Gmats,Smats,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp      !as Smats
    complex(8),dimension(:,:,:),allocatable           :: zeta_mats ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
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
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,3)
    Nlso   = Nlat*Nspin*Norb
    !
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],"dmft_get_gloc_matsubara_normal_cluster_mpi","Hk")
    call assert_shape(Smats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_cluster_mpi","Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_cluster_mpi","Gmats")
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nlso,Nlso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Smats(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (no print)"
    if(mpi_master)call start_timer
    Gmats=zero
    if(Lmats>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_cluster_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
          Gmats = Gmats + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_cluster(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
       endif
#else
       Gmats=Gtmp
#endif
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_normal_cluster


  subroutine dmft_get_gloc_matsubara_normal_dos(Ebands,Dbands,Hloc,Gmats,Smats)
  real(8),dimension(:,:),intent(in)                             :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(:,:),intent(in)                             :: Dbands    !1. [Nspin*Norb][Lk] / 2. [1][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                  :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout)                 :: Gmats     !as Smats
  !
  complex(8)                                                    :: gktmp
  complex(8),dimension(:,:,:),allocatable                       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),allocatable                   :: Gtmp      !as Smats
  !
  real(8)                                                       :: beta
  real(8)                                                       :: xmu,eps
  real(8)                                                       :: wini,wfin
  !
  !New
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
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Lk    = size(Ebands,2)
  Nso   = Nspin*Norb
  !
  !case F  => one DOS, H(e)=diag(Ebands), non-diag
  !case T => more DOS, Ebands, diag
  dos_diag = .not.( size(Dbands,1) < size(Ebands,1) )
  !
  !Testing part:
  call assert_shape(Ebands,[Nso,Lk],"dmft_get_gloc_matsubara_normal_dos_main_mpi","Ebands")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_dos_main_mpi","Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_dos_main_mpi","Gmats")
  !
  !Allocate and setup the Matsubara freq.
  allocate(zeta_mats(Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !
  do i=1,Lmats
     zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (no print)"
  if(mpi_master)call start_timer
  allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
  !
  !
  ! diagonal case
  if(dos_diag)then
     if(Lmats>=Lk)then
        do i = 1+mpi_rank, Lmats, mpi_size
           do ispin=1,Nspin
              do iorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 do ik=1,Lk
                    gktmp = Dbands(io,ik)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io,ik) )
                    Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
                 enddo
              enddo
           enddo
           if(mpi_master)call eta(i,Lmats)
        end do
     else
        do ik = 1+mpi_rank, Lk, mpi_size
           do ispin=1,Nspin
              do iorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 do i=1,Lmats
                    gktmp = Dbands(io,ik)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io,ik) )
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
        do i = 1+mpi_rank, Lmats, mpi_size                                 !MPI loop over Matsubara frequencies
           do ik=1,Lk                                                      !for all e-value (here named ik)
              Gdos_tmp = zeta_mats(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = zeta_mats - Hloc - H(e)
              call inv(Gdos_tmp)
              Gtmp(:,:,:,:,i) = Gtmp(:,:,:,:,i) + Dbands(1,ik)*so2nn_reshape(Gdos_tmp,Nspin,Norb)
           enddo
           if(mpi_master)call eta(i,Lmats)
        end do
     else
        do ik = 1+mpi_rank, Lk, mpi_size
           do i=1,Lmats
              Gdos_tmp = zeta_mats(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = zeta_mats - Hloc - H(e)
              call inv(Gdos_tmp)
              Gtmp(:,:,:,:,i) = Gtmp(:,:,:,:,i) + Dbands(1,ik)*so2nn_reshape(Gdos_tmp,Nspin,Norb)
           enddo
           if(mpi_master)call eta(ik,Lk)
        end do
     end if
  end if
#ifdef _MPI    
  if(check_MPI())then
     Gmats=zero
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
  else
     Gmats=Gtmp
  endif
#else
  Gmats=Gtmp
#endif
  if(mpi_master)call stop_timer
end subroutine dmft_get_gloc_matsubara_normal_dos



  subroutine dmft_get_gloc_matsubara_normal_ineq(Hk,Gmats,Smats,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
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
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Smats")
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Gmats")
    !
    if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (no print)"
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gkmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gmats=zero
    if(Lmats>=Lk)then
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_ineq_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gmats = Gmats + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_tridiag_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gmats = Gmats + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
    else
       allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_ineq(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gtmp = Gtmp + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
             Gtmp = Gtmp + Gkmats/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
       endif
#else
       Gmats=Gtmp
#endif
       deallocate(Gtmp)
       !
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_normal_ineq








  subroutine dmft_get_gloc_matsubara_normal_gij(Hk,Gmats,Smats,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats     !      [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
    complex(8),dimension(:,:,:,:),allocatable         :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
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
    !
    Nlat  = size(Smats,1)
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_normal_gij_main_mpi',"Hk")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_gij_main_mpi',"Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_gij_main_mpi',"Gmats")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (no print)"
    !
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       do i=1,Lmats
          zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gmats=zero
    if(Lmats>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_gij_mpi(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero     
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
       endif
#else
       Gmats=Gtmp
#endif
       deallocate(Gtmp)
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_normal_gij








  subroutine dmft_get_gloc_matsubara_superc_main(Hk,Gmats,Smats,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
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
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nso,Nso,Lk],'dmft_get_gloc_matsubara_superc_main_mpi',"Hk")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main_mpi',"Smats")
    call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main_mpi',"Gmats")
    !
    allocate(zeta_mats(2,2,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
    if(mpi_master)call start_timer
    Gmats=zero
    if(Lmats>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_mpi(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lmats))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
       endif
#else
       Gmats=Gtmp
#endif
       deallocate(Gtmp)
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_superc_main


  subroutine dmft_get_gloc_matsubara_superc_dos(Ebands,Dbands,Hloc,Gmats,Smats)
    real(8),dimension(:,:),intent(in)                                     :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in)           :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                          :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)                          :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout)                       :: Gmats     !as Smats
    ! arrays
    complex(8),dimension(2,2,size(Ebands,1),size(Ebands,1),size(Smats,6)) :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(2*size(Ebands,1),2*size(Ebands,1))               :: Gmatrix
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable                         :: Gtmp    !as Smats
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
    !
    Nspin = size(Smats,2)
    Norb  = size(Smats,4)
    Lmats = size(Smats,6)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_matsubara_superc_dos',"Ebands")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Smats")
    call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Gmats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    zeta_mats=zero
    do i=1,Lmats
       zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) - diag(Hloc) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
       zeta_mats(1,2,:,:,i) =                                      -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,1,:,:,i) =                                      -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
       zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + diag(Hloc) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
    if(mpi_master)call start_timer
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lmats))
    Gtmp=zero
    Gmats=zero
    do ik = 1+mpi_rank, Lk, mpi_size
       !
       do i=1,Lmats
          Gmatrix  = zero
          Gmatrix(1:Nso,1:Nso)             = zeta_mats(1,1,:,:,i) - diag(Ebands(:,ik))
          Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta_mats(1,2,:,:,i)
          Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta_mats(2,1,:,:,i)
          Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta_mats(2,2,:,:,i) + diag(Ebands(:,ik))
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
       Gmats=zero
       call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    else
       Gmats=Gtmp
    endif
#else
    Gmats=Gtmp
#endif
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_superc_dos


  subroutine dmft_get_gloc_matsubara_superc_ineq(Hk,Gmats,Smats,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
    logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
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
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_superc_ineq_main_mpi',"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_ineq_main_mpi',"Smats")
    call assert_shape(Gmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_ineq_main_mpi',"Gmats")
    !
    if(mpi_master)write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
    !
    allocate(Gkmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))    
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    do ilat=1,Nlat
       !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
       !G22(iw) = -[G11[iw]]*
       !G21(iw) =   G12[w]
       do i=1,Lmats
          zeta_mats(1,1,ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
          zeta_mats(2,2,ilat,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb) )
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Gmats=zero
    if(Lmats>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_ineq_mpi(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
          Gmats = Gmats + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc_ineq(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
          Gtmp = Gtmp + Gkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
       endif
#else
       Gmats=Gtmp
#endif
       deallocate(Gtmp)
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_superc_ineq



  subroutine dmft_get_gloc_matsubara_superc_gij(Hk,Gmats,Fmats,Smats,hk_symm)
    complex(8),dimension(:,:,:,:)                       :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Fmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Smats           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    logical,optional                                  :: hk_symm(size(Hk,4))
    logical                                           :: hk_symm_(size(Hk,4))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkmats,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer                                           :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
    !
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
    !
    Nlat  = size(Smats,2)
    Nspin = size(Smats,3)
    Norb  = size(Smats,5)
    Lmats = size(Smats,7)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Hk")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Smats")
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Gmats")
    call assert_shape(Fmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Fmats")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function (no print)"
    !
    allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
    allocate(Fkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Fkmats=zero
    allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats));zeta_mats=zero
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
    !G22(iw) = -[G11[iw]]*
    !G21(iw) =   G12[w]
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu
             zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg( Smats(1,ilat,ispin,jspin,iorb,jorb,:) )
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(mpi_master)call start_timer
    Gmats=zero
    Fmats=zero
    if(Lmats>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_gij(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
          Gmats = Gmats + Gkmats/dble(Lk)
          Fmats = Fmats + Fkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
       allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Ftmp=zero
       do ik=1,Lk
          call invert_gk_superc_gij(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
          Gtmp = Gtmp + Gkmats/dble(Lk)
          Ftmp = Ftmp + Fkmats/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gmats=zero
          Fmats=zero
          call Mpi_AllReduce(Gtmp, Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
          call Mpi_AllReduce(Ftmp, Fmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gmats=Gtmp
          Fmats=Ftmp
       endif
#else
       Gmats=Gtmp
       Fmats=Ftmp
#endif
       deallocate(Gtmp,Ftmp)
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_matsubara_superc_gij















end module GLOC_MATSUBARA
