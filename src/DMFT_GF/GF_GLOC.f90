module GF_GLOC
  USE GF_COMMON
  implicit none
  private


  !##################################################################
  interface get_gloc
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gloc_normal_main
     module procedure :: dmft_get_gloc_normal_dos
     module procedure :: dmft_get_gloc_normal_ineq
     module procedure :: dmft_get_gloc_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_superc_main
     module procedure :: dmft_get_gloc_superc_dos
     module procedure :: dmft_get_gloc_superc_ineq
  end interface get_gloc

  interface get_gloc_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_main
     module procedure :: dmft_get_gloc_matsubara_normal_dos
     module procedure :: dmft_get_gloc_matsubara_normal_ineq
     module procedure :: dmft_get_gloc_matsubara_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_matsubara_superc_main
     module procedure :: dmft_get_gloc_matsubara_superc_dos
     module procedure :: dmft_get_gloc_matsubara_superc_ineq
  end interface get_gloc_matsubara

  interface get_gloc_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_main
     module procedure :: dmft_get_gloc_realaxis_normal_dos
     module procedure :: dmft_get_gloc_realaxis_normal_ineq
     module procedure :: dmft_get_gloc_realaxis_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_realaxis_normal_cluster_ineq
#endif
     module procedure :: dmft_get_gloc_realaxis_superc_main
     module procedure :: dmft_get_gloc_realaxis_superc_dos
     module procedure :: dmft_get_gloc_realaxis_superc_ineq
  end interface get_gloc_realaxis
  !##################################################################






  !##################################################################
  interface dmft_gloc
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gloc_normal_main
     module procedure :: dmft_get_gloc_normal_dos
     module procedure :: dmft_get_gloc_normal_ineq
     module procedure :: dmft_get_gloc_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_superc_main
     module procedure :: dmft_get_gloc_superc_dos
     module procedure :: dmft_get_gloc_superc_ineq
  end interface dmft_gloc

  interface dmft_gloc_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_main
     module procedure :: dmft_get_gloc_matsubara_normal_dos
     module procedure :: dmft_get_gloc_matsubara_normal_ineq
     module procedure :: dmft_get_gloc_matsubara_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_matsubara_superc_main
     module procedure :: dmft_get_gloc_matsubara_superc_dos
     module procedure :: dmft_get_gloc_matsubara_superc_ineq
  end interface dmft_gloc_matsubara

  interface dmft_gloc_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_main
     module procedure :: dmft_get_gloc_realaxis_normal_dos
     module procedure :: dmft_get_gloc_realaxis_normal_ineq
     module procedure :: dmft_get_gloc_realaxis_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_realaxis_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_realaxis_superc_main
     module procedure :: dmft_get_gloc_realaxis_superc_dos
     module procedure :: dmft_get_gloc_realaxis_superc_ineq
  end interface dmft_gloc_realaxis
  !##################################################################






  !##################################################################
  interface dmft_get_gloc
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gloc_normal_main
     module procedure :: dmft_get_gloc_normal_dos
     module procedure :: dmft_get_gloc_normal_ineq
     module procedure :: dmft_get_gloc_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_superc_main
     module procedure :: dmft_get_gloc_superc_dos
     module procedure :: dmft_get_gloc_superc_ineq
  end interface dmft_get_gloc

  interface dmft_get_gloc_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_main
     module procedure :: dmft_get_gloc_matsubara_normal_dos
     module procedure :: dmft_get_gloc_matsubara_normal_ineq
     module procedure :: dmft_get_gloc_matsubara_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_matsubara_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_matsubara_superc_main
     module procedure :: dmft_get_gloc_matsubara_superc_dos
     module procedure :: dmft_get_gloc_matsubara_superc_ineq
  end interface dmft_get_gloc_matsubara

  interface dmft_get_gloc_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_main
     module procedure :: dmft_get_gloc_realaxis_normal_dos
     module procedure :: dmft_get_gloc_realaxis_normal_ineq
     module procedure :: dmft_get_gloc_realaxis_normal_cluster
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_get_gloc_realaxis_normal_cluster_ineq
#endif
     !
     module procedure :: dmft_get_gloc_realaxis_superc_main
     module procedure :: dmft_get_gloc_realaxis_superc_dos
     module procedure :: dmft_get_gloc_realaxis_superc_ineq
  end interface dmft_get_gloc_realaxis
  !##################################################################




  !##################################################################
  interface get_gij
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gloc_normal_gij
     module procedure :: dmft_get_gloc_superc_gij
  end interface get_gij

  interface get_gij_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_gij
     module procedure :: dmft_get_gloc_matsubara_superc_gij
  end interface get_gij_matsubara

  interface get_gij_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_gij
     module procedure :: dmft_get_gloc_realaxis_superc_gij
  end interface get_gij_realaxis
  !##################################################################




  !##################################################################
  interface dmft_gij
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gloc_normal_gij
     module procedure :: dmft_get_gloc_superc_gij
  end interface dmft_gij

  interface dmft_gij_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_gij
     module procedure :: dmft_get_gloc_matsubara_superc_gij
  end interface dmft_gij_matsubara

  interface dmft_gij_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_gij
     module procedure :: dmft_get_gloc_realaxis_superc_gij
  end interface dmft_gij_realaxis
  !##################################################################




  !##################################################################
  interface dmft_get_gij
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_get_gloc_normal_gij
     module procedure :: dmft_get_gloc_superc_gij
  end interface dmft_get_gij

  interface dmft_get_gij_matsubara
     module procedure :: dmft_get_gloc_matsubara_normal_gij
     module procedure :: dmft_get_gloc_matsubara_superc_gij
  end interface dmft_get_gij_matsubara

  interface dmft_get_gij_realaxis
     module procedure :: dmft_get_gloc_realaxis_normal_gij
     module procedure :: dmft_get_gloc_realaxis_superc_gij
  end interface dmft_get_gij_realaxis
  !##################################################################








  !PUBLIC:
  public :: get_gloc
  public :: get_gloc_matsubara
  public :: get_gloc_realaxis


  public :: dmft_gloc
  public :: dmft_gloc_matsubara
  public :: dmft_gloc_realaxis

  public :: dmft_get_gloc
  public :: dmft_get_gloc_matsubara
  public :: dmft_get_gloc_realaxis


  public :: get_gij
  public :: get_gij_matsubara
  public :: get_gij_realaxis


  public :: dmft_gij
  public :: dmft_gij_matsubara
  public :: dmft_gij_realaxis


  public :: dmft_get_gij
  public :: dmft_get_gij_matsubara
  public :: dmft_get_gij_realaxis





contains





  !##################################################################
  !##################################################################
  !                          NORMAL
  !##################################################################
  !##################################################################
  subroutine dmft_get_gloc_normal_main(Hk,Gloc,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc         !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:),optional              :: zeta
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gk    !as Sigma  
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sigma
    complex(8),dimension(:,:,:),allocatable       :: csi ![Nspin*Norb][Nspin*Norb][Lmats]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)    
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nso,Nso,Lk],"dmft_get_gloc_normal_main","Hk")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gloc_normal_main","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gloc_normal_main","Gloc")
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gk(Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(csi(Nso,Nso,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Nso) - nn2so_reshape(Sigma(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    !
    if(Lfreq >= Lk)then
       Gloc = zero
       do ik=1,Lk
          call invert_gk_normal_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)      
          Gloc = Gloc + Gk/Lk
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
       Gloc = zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
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
  end subroutine dmft_get_gloc_normal_main






  subroutine dmft_get_gloc_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,axis,zeta)
    real(8),dimension(:,:),intent(in)             :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)             :: Dbands    ![Nspin*Norb][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                              :: axis
    complex(8),dimension(:),optional              :: zeta
    !
    complex(8)                                    :: gktmp
    complex(8),dimension(:,:,:),allocatable       :: csi ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sigma
    !
    !New
    complex(8),dimension(:,:),allocatable         :: Gdos_tmp ![Nspin*Norb][Nspin*Norb]
    logical                                       :: dos_diag !1. T / 2. F
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lfreq = size(Sigma,5)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    dos_diag = .not.(size(Dbands,1) < size(Ebands,1))
    !
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],"dmft_get_gloc_normal_dos","Ebands")
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gloc_normal_dos","Sigma")
    call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gloc_normal_dos","Gloc")
    !
    !Allocate and setup the Matsubara freq.
    allocate(csi(Nso,Nso,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Nso) - nn2so_reshape(Sigma(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
    !
    !
    ! diagonal case
    if(dos_diag)then
       !
       if(Lfreq>=Lk)then
          do i = 1+mpi_rank, Lfreq, mpi_size
             do ispin=1,Nspin
                do iorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   do ik=1,Lk
                      gktmp = Dbands(io,ik)/( csi(io,io,i)-Hloc(io)-Ebands(io,ik) )
                      Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
                   enddo
                enddo
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do ispin=1,Nspin
                do iorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   do i=1,Lfreq
                      gktmp = Dbands(io,ik)/( csi(io,io,i)-Hloc(io)-Ebands(io,ik) )
                      Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
                   enddo
                enddo
             enddo
             if(mpi_master)call eta(ik,Lk)
          end do
       end if
       !
    else
       !
       allocate(Gdos_tmp(Nso,Nso)) ;Gdos_tmp=zero
       if(Lfreq>=Lk)then
          do i = 1+mpi_rank, Lfreq, mpi_size                                 !MPI loop over Matsubara frequencies
             do ik=1,Lk                                                      !for all e-value (here named ik)
                Gdos_tmp = csi(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = csi - Hloc - H(e)
                call inv(Gdos_tmp)
                Gtmp(:,:,:,:,i) = Gtmp(:,:,:,:,i) + Dbands(1,ik)*so2nn_reshape(Gdos_tmp,Nspin,Norb)
             enddo
             if(mpi_master)call eta(i,Lfreq)
          end do
       else
          do ik = 1+mpi_rank, Lk, mpi_size
             do i=1,Lfreq
                Gdos_tmp = csi(:,:,i)-diag(Hloc(:))-diag(Ebands(:,ik)) !G(e,w) = csi - Hloc - H(e)
                call inv(Gdos_tmp)
                Gtmp(:,:,:,:,i) = Gtmp(:,:,:,:,i) + Dbands(1,ik)*so2nn_reshape(Gdos_tmp,Nspin,Norb)
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
  end subroutine dmft_get_gloc_normal_dos




  subroutine dmft_get_gloc_normal_ineq(Hk,Gloc,Sigma,axis,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gk    !as Sigma
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
    complex(8),dimension(:,:,:,:),allocatable       :: csi ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_normal_ineq',"Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_normal_ineq',"Sigma")
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_normal_ineq',"Gloc")
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gk(Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(csi(Nlat,Nso,Nso,Lfreq))
    !
    do ilat=1,Nlat
       do i=1,Lfreq
          csi(ilat,:,:,i) = (wfreq(i)+xmu)*eye(Nso)     - nn2so_reshape(Sigma(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_ineq_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gloc = Gloc + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_tridiag_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gloc = Gloc + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
    else
       allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_ineq(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gtmp = Gtmp + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_tridiag(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gtmp = Gtmp + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
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
  end subroutine dmft_get_gloc_normal_ineq


  subroutine dmft_get_gloc_normal_cluster(Hk,Gloc,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma ![Nlat:Nlat:Nspin:Nspin:NorbNorb:Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gk    !as Sigma
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp  !as Sigma
    complex(8),dimension(:,:,:),allocatable           :: csi  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Lk    = size(Hk,3)
    Nlso   = Nlat*Nspin*Norb
    !
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],"dmft_get_gloc_normal_cluster","Hk")
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gloc_normal_cluster","Sigma")
    call assert_shape(Gloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_gloc_normal_cluster","Gloc")
    !
    !Allocate and setup the Matsubara freq.
    allocate(Gk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(csi(Nlso,Nlso,Lfreq))
    !
    do i=1,Lfreq
       csi(:,:,i)=(wfreq(i)+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Sigma(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_cluster_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)      
          Gloc = Gloc + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_cluster(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
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
  end subroutine dmft_get_gloc_normal_cluster



#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gloc_normal_cluster_ineq(Hk,Gloc,Sigma,axis,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)              :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    character(len=*)                                    :: axis
    complex(8),dimension(:),optional                    :: zeta
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,dimension(size(Hk,3)),optional              :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                       :: hk_symm_  ![Lk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gk    !as Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),allocatable   :: Gtmp      !as Sigma
    complex(8),dimension(:,:,:,:),allocatable           :: csi ![Nineq][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nineq  = size(Sigma,1)
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,4)
    Norb  = size(Sigma,6)
    Lfreq = size(Sigma,8)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    Nilso = Nineq*Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nilso,Nilso,Lk],'dmft_get_gloc_normal_ineq',"Hk")
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_normal_ineq',"Sigma")
    call assert_shape(Gloc,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_normal_ineq',"Gloc")
    !
    if(mpi_master)write(*,"(A)")"Get local Green's function, axis:"//str(axis)
    if(mpi_master)then
       if(.not.tridiag_)then
          write(*,"(A)")"Direct Inversion algorithm:"
       else
          write(*,"(A)")"Quantum Zipper algorithm:"
       endif
    endif
    !
    allocate(Gk(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(csi(Nineq,Nlso,Nlso,Lfreq))
    !
    do iineq=1,Nineq
       do i=1,Lfreq
          csi(iineq,:,:,i) = (wfreq(i)+xmu)*eye(Nlso)     - nnn2lso_cluster_reshape(Sigma(iineq,:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       if(.not.tridiag_)then
          do ik=1,Lk
             call invert_gk_normal_cluster_ineq_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gloc = Gloc + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1,Lk
             call invert_gk_normal_cluster_ineq_tridiag_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gloc = Gloc + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
    else
       allocate(Gtmp(Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
       if(.not.tridiag_)then
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_cluster_ineq(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gtmp = Gtmp + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       else
          do ik=1+mpi_rank,Lk,mpi_size
             call invert_gk_normal_cluster_ineq_tridiag(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
             Gtmp = Gtmp + Gk/dble(Lk)
             if(mpi_master)call eta(ik,Lk)
          end do
       endif
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
  end subroutine dmft_get_gloc_normal_cluster_ineq
#endif





  subroutine dmft_get_gloc_normal_gij(Hk,Gloc,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     !      [Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    character(len=*)                                  :: axis
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gk    !as Sigma
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
    complex(8),dimension(:,:,:,:),allocatable         :: csi ![Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    Nlat  = size(Sigma,1)
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Lk    = size(Hk,3)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_normal_gij',"Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_normal_gij',"Sigma")
    call assert_shape(Gloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_normal_gij',"Gloc")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Green's function, axis:"//str(axis)
    !
    allocate(Gk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gk=zero
    allocate(csi(Nlat,Nso,Nso,Lfreq))
    !
    do ilat=1,Nlat
       do i=1,Lfreq
          csi(ilat,:,:,i) = (wfreq(i)+xmu)*eye(Nso)     - nn2so_reshape(Sigma(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       do ik=1,Lk
          call invert_gk_normal_gij_mpi(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
          Gloc = Gloc + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero     
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_normal_gij(csi,Hk(:,:,ik),hk_symm_(ik),Gk)
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
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_normal_gij













  !##################################################################
  !##################################################################
  !                          SUPERC
  !##################################################################
  !##################################################################

  subroutine dmft_get_gloc_superc_main(Hk,Gloc,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)        :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:),optional                :: zeta
    logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gk    !as Sigma
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
    complex(8),dimension(:,:,:,:,:),allocatable     :: csi ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nso,Nso,Lk],'dmft_get_gloc_superc',"Hk")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc',"Sigma")
    call assert_shape(Gloc,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc',"Gloc")
    !
    allocate(csi(2,2,Nso,Nso,Lfreq))
    !
    select case(axis)
    case default; stop "dmft_get_gloc_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso) -        nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) =                         -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) =                         -        nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb) )
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i) + xmu)*eye(Nso)                 - nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) =                                           - nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) =                                           - nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = -conjg(wfreq(Lfreq+1-i) + xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,Lfreq+1-i),Nspin,Norb) )
       enddo
    end select
    !
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"Get local Nambu Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_mpi(csi,Hk(:,:,:,ik),hk_symm_(ik),Gk)
          Gloc = Gloc + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lfreq))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc(csi,Hk(:,:,:,ik),hk_symm_(ik),Gk)
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
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_superc_main



  subroutine dmft_get_gloc_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,axis,zeta)
    real(8),dimension(:,:,:),intent(in)             :: Ebands    ![2][Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands    ![Nspin*Norb][Lk]/[1][Lk]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    character(len=*)                                :: axis
    complex(8),dimension(:),optional                :: zeta
    ! arrays
    complex(8),dimension(:,:,:,:,:),allocatable     :: csi ![2][2][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:),allocatable           :: Gmatrix
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
    !
    complex(8),dimension(:,:),allocatable           :: Gdos_tmp ![Nspin*Norb][Nspin*Norb]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    !
    Nspin = size(Sigma,2)
    Norb  = size(Sigma,4)
    Lfreq = size(Sigma,6)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb
    !
    !case F  => 1  DOS, H(e)=diag(Ebands), non-diagonal case
    !case T  => >1 DOS, Ebands, diagonal case
    ! dos_diag = .not.(size(Dbands,1) < size(Ebands,2))
    !
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_superc_dos',"Ebands")
    call assert_shape(Sigma,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_dos',"Sigma")
    call assert_shape(Gloc,[2,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_dos',"Gloc")
    !
    !
    allocate(csi(2,2,Nso,Nso,Lfreq));csi = zero
    select case(axis)
    case default; stop "dmft_get_gloc_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso) - diag(Hloc(1,:)) -&
               nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) =                                           -&
               nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) =                                           -&
               nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = (wfreq(i)-xmu)*eye(Nso) + diag(Hloc(2,:)) +&
               conjg( nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb) )
       enddo
    case("realaxis","real","Realaxis","Real")
       do i=1,Lfreq
          csi(1,1,:,:,i) = (wfreq(i)+xmu)*eye(Nso)                 - diag(Hloc(1,:)) - nn2so_reshape(Sigma(1,:,:,:,:,i),Nspin,Norb)
          csi(1,2,:,:,i) = -nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,1,:,:,i) = -nn2so_reshape(Sigma(2,:,:,:,:,i),Nspin,Norb)
          csi(2,2,:,:,i) = -conjg( wfreq(Lfreq+1-i)+xmu )*eye(Nso) + diag(Hloc(2,:)) + conjg( nn2so_reshape(Sigma(1,:,:,:,:,Lfreq+1-i),Nspin,Norb) )
       enddo
    end select
    !
    !invert (Z-Hk) for each k-point
    if(mpi_master)write(*,"(A)")"WARNING: This method is limited to diagonal case.. "
    if(mpi_master)write(*,"(A)")"Get local Nambu Green's function, axis:"//str(axis)
    if(mpi_master)call start_timer
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(Gmatrix(2*Nso,2*Nso))
    Gtmp=zero
    !
    ! if(dos_diag)then  
    !
    do ik = 1+mpi_rank, Lk, mpi_size
       do i=1,Lfreq
          Gmatrix  = zero
          Gmatrix(1:Nso,1:Nso)             = csi(1,1,:,:,i) - diag(Ebands(1,:,ik))
          Gmatrix(1:Nso,Nso+1:2*Nso)       = csi(1,2,:,:,i)
          Gmatrix(Nso+1:2*Nso,1:Nso)       = csi(2,1,:,:,i)
          Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = csi(2,2,:,:,i) - diag(Ebands(2,:,ik))
          call inv(Gmatrix)
          do ispin=1,Nspin
             do iorb=1,Norb
                io = iorb + (ispin-1)*Norb
                Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + Gmatrix(io,io)*Dbands(io,ik)
                Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + Gmatrix(io,Nso+io)*Dbands(io,ik)
             enddo
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
    !   else
    !      Gtmp=zero
    !      do ik = 1+mpi_rank, Lk, mpi_size
    !         do i=1,Lfreq
    !            Gmatrix  = zero
    !            Gmatrix(1:Nso,1:Nso)             = csi(1,1,:,:,i) - diag(Ebands(1,:,ik))
    !            Gmatrix(1:Nso,Nso+1:2*Nso)       = csi(1,2,:,:,i)
    !            Gmatrix(Nso+1:2*Nso,1:Nso)       = csi(2,1,:,:,i)
    !            Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = csi(2,2,:,:,i) - diag(Ebands(2,:,ik))
    !            call inv(Gmatrix)
    !            do ispin=1,Nspin
    !               do iorb=1,Norb
    !                  io = iorb + (ispin-1)*Norb
    !                  Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + Gmatrix(io,io)*Dbands(1,ik)
    !                  Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + Gmatrix(io,Nso+io)*Dbands(1,ik)
    !               enddo
    !            enddo
    !         enddo
    !         call eta(ik,Lk)
    !      enddo
    ! #ifdef _MPI    
    !      if(check_MPI())then
    !         Gloc=zero
    !         call Mpi_AllReduce(Gtmp,Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
    !      else
    !         Gloc=Gtmp
    !      endif
    ! #else
    !      Gloc=Gtmp
    ! #endif
    !      if(mpi_master)call stop_timer
    !   endif
  end subroutine dmft_get_gloc_superc_dos



  subroutine dmft_get_gloc_superc_ineq(Hk,Gloc,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    character(len=*)                                  :: axis
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gk    !as Sigma
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sigma
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: csi ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_superc_ineq',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_ineq',"Sigma")
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_ineq',"Gloc")
    !
    if(mpi_master)write(*,"(A)")"Get local Nambu Green's function, axis:"//str(axis)
    !
    allocate(Gk(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
    allocate(csi(2,2,Nlat,Nso,Nso,Lfreq))
    !
    select case(axis)
    case default; stop "dmft_get_gloc_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       do ilat=1,Nlat
          !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
          !G22(iw) = -[G11[iw]]*
          !G21(iw) =   G12[w]
          do i=1,Lfreq
             csi(1,1,ilat,:,:,i) = (wfreq(i)+xmu)*eye(Nso) -        nn2so_reshape(Sigma(1,ilat,:,:,:,:,i),Nspin,Norb)
             csi(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,2,ilat,:,:,i) = (wfreq(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,ilat,:,:,:,:,i),Nspin,Norb) )
          enddo
       enddo
    case("realaxis","real","Realaxis","Real")
       do ilat=1,Nlat
          !SYMMETRIES in real-frequencies   [assuming a real order parameter]
          !G22(w)  = -[G11[-w]]*
          !G21(w)  =   G12[w]   
          do i=1,Lfreq
             csi(1,1,ilat,:,:,i) = (wfreq(i)+xmu)*eye(Nso)  - nn2so_reshape(Sigma(1,ilat,:,:,:,:,i),Nspin,Norb)
             csi(1,2,ilat,:,:,i) =                          - nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,1,ilat,:,:,i) =                          - nn2so_reshape(Sigma(2,ilat,:,:,:,:,i),Nspin,Norb)
             csi(2,2,ilat,:,:,i) = -conjg(wfreq(Lfreq+1-i)+xmu)*eye(Nso) + conjg( nn2so_reshape(Sigma(1,ilat,:,:,:,:,Lfreq+1-i),Nspin,Norb) )
          enddo
       enddo
    end select
    !
    if(mpi_master)call start_timer
    Gloc=zero
    if(Lfreq>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_ineq_mpi(csi,Hk(:,:,:,ik),hk_symm_(ik),Gk)
          Gloc = Gloc + Gk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq))
       do ik=1+mpi_rank,Lk,mpi_size
          call invert_gk_superc_ineq(csi,Hk(:,:,:,ik),hk_symm_(ik),Gk)
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
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_superc_ineq




  subroutine dmft_get_gloc_superc_gij(Hk,Gloc,Floc,Sigma,axis,zeta,hk_symm)
    complex(8),dimension(:,:,:,:)                     :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gloc           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Floc           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sigma           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    character(len=*)                                  :: axis
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm(size(Hk,4))
    logical                                           :: hk_symm_(size(Hk,4))
    !
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: csi       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gk,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fk,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
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
    if(present(zeta))then
       call build_frequency_array(axis,zeta)
    else
       call build_frequency_array(axis)
    endif
    !
    !
    Nlat  = size(Sigma,2)
    Nspin = size(Sigma,3)
    Norb  = size(Sigma,5)
    Lfreq = size(Sigma,7)
    Lk    = size(Hk,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_superc_gij',"Hk")
    call assert_shape(Sigma,[2,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_gij',"Sigma")
    call assert_shape(Gloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_gij',"Gloc")
    call assert_shape(Floc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],'dmft_get_gloc_superc_gij',"Floc")
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(mpi_master)write(*,"(A)")"Get full Nambu Green's functions, axis:"//str(axis)
    !
    allocate(Gk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gk=zero
    allocate(Fk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Fk=zero
    allocate(csi(2,2,Nlat,Nso,Nso,Lfreq));csi=zero
    !
    select case(axis)
    case default; stop "dmft_get_gloc_superc_main error: specified axis not valid. axis={matsubara,realaxis}"
    case("matsubara","mats","Matsubara","Mats")
       !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
       !G22(iw) = -[G11[iw]]*
       !G21(iw) =   G12[w]
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                io = iorb + (ispin-1)*Norb
                js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                csi(1,1,ilat,io,io,:) = wfreq(:) + xmu
                csi(2,2,ilat,io,io,:) = wfreq(:) - xmu
             enddo
          enddo
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      csi(1,1,ilat,io,jo,:) = csi(1,1,ilat,io,jo,:) - Sigma(1,ilat,ispin,jspin,iorb,jorb,:)
                      csi(1,2,ilat,io,jo,:) = csi(1,2,ilat,io,jo,:) - Sigma(2,ilat,ispin,jspin,iorb,jorb,:)
                      csi(2,1,ilat,io,jo,:) = csi(2,1,ilat,io,jo,:) - Sigma(2,ilat,ispin,jspin,iorb,jorb,:)
                      csi(2,2,ilat,io,jo,:) = csi(2,2,ilat,io,jo,:) + conjg( Sigma(1,ilat,ispin,jspin,iorb,jorb,:) )
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
       !
    case("realaxis","real","Realaxis","Real")
       !SYMMETRIES in real-frequencies   [assuming a real order parameter]
       !G22(w)  = -[G11[-w]]*
       !G21(w)  =   G12[w]
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                io = iorb + (ispin-1)*Norb
                js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                csi(1,1,ilat,io,io,:) =         wfreq(:)          + xmu
                csi(2,2,ilat,io,io,:) = -conjg( wfreq(Lfreq:1:-1) + xmu )
             enddo
          enddo
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      csi(1,1,ilat,io,jo,:) = csi(1,1,ilat,io,jo,:) - Sigma(1,ilat,ispin,jspin,iorb,jorb,:)
                      csi(1,2,ilat,io,jo,:) = csi(1,2,ilat,io,jo,:) - Sigma(2,ilat,ispin,jspin,iorb,jorb,:)
                      csi(2,1,ilat,io,jo,:) = csi(2,1,ilat,io,jo,:) - Sigma(2,ilat,ispin,jspin,iorb,jorb,:)
                      csi(2,2,ilat,io,jo,:) = csi(2,2,ilat,io,jo,:) + conjg( Sigma(1,ilat,ispin,jspin,iorb,jorb,Lfreq:1:-1) )
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
       !
    end select
    !
    !
    if(mpi_master)call start_timer
    Gloc=zero
    Floc=zero
    if(Lfreq>=Lk)then
       do ik=1,Lk
          call invert_gk_superc_gij(csi,Hk(:,:,:,ik),hk_symm_(ik),Gk,Fk)
          Gloc = Gloc + Gk/dble(Lk)
          Floc = Floc + Fk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
    else
       allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Gtmp=zero
       allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq));Ftmp=zero
       do ik=1,Lk
          call invert_gk_superc_gij(csi,Hk(:,:,:,ik),hk_symm_(ik),Gk,Fk)
          Gtmp = Gtmp + Gk/dble(Lk)
          Ftmp = Ftmp + Fk/dble(Lk)
          if(mpi_master)call eta(ik,Lk)
       end do
#ifdef _MPI    
       if(check_MPI())then
          Gloc=zero
          Floc=zero
          call Mpi_AllReduce(Gtmp, Gloc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
          call Mpi_AllReduce(Ftmp, Floc, size(Gloc), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, MPI_ierr)
       else
          Gloc=Gtmp
          Floc=Ftmp
       endif
#else
       Gloc=Gtmp
       Floc=Ftmp
#endif
       deallocate(Gtmp,Ftmp)
    end if
    if(mpi_master)call stop_timer
  end subroutine dmft_get_gloc_superc_gij



















  !##################################################################
  !##################################################################
  !                         INTERFACES
  !##################################################################
  !##################################################################

  subroutine dmft_get_gloc_matsubara_normal_main(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc         !as Sigma
    complex(8),dimension(:),optional              :: zeta
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !
    hk_symm_  =.false. ; if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_main(Hk,Gloc,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_main(Hk,Gloc,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_normal_main

  subroutine dmft_get_gloc_realaxis_normal_main(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc         !as Sigma
    complex(8),dimension(:),optional              :: zeta
    logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
    !
    hk_symm_  =.false. ; if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_main(Hk,Gloc,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_main(Hk,Gloc,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_normal_main






  !##################################################################
  !##################################################################




  subroutine dmft_get_gloc_matsubara_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,zeta)
    real(8),dimension(:,:),intent(in)             :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)             :: Dbands    ![Nspin*Norb][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    complex(8),dimension(:),optional              :: zeta
    if(present(zeta))then
       call dmft_get_gloc_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"matsubara")
    else
       call dmft_get_gloc_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"matsubara",zeta)
    endif
  end subroutine dmft_get_gloc_matsubara_normal_dos

  subroutine dmft_get_gloc_realaxis_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,zeta)
    real(8),dimension(:,:),intent(in)             :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)             :: Dbands    ![Nspin*Norb][Lk] /[1][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    complex(8),dimension(:),optional              :: zeta
    if(present(zeta))then
       call dmft_get_gloc_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"realaxis")
    else
       call dmft_get_gloc_normal_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"realaxis",zeta)
    endif
  end subroutine dmft_get_gloc_realaxis_normal_dos




  !##################################################################
  !##################################################################


  subroutine dmft_get_gloc_matsubara_normal_ineq(Hk,Gloc,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_ineq(Hk,Gloc,Sigma,"matsubara",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_ineq(Hk,Gloc,Sigma,"matsubara",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_normal_ineq

  subroutine dmft_get_gloc_realaxis_normal_ineq(Hk,Gloc,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    logical,dimension(size(Hk,3)),optional          :: hk_symm
    logical,dimension((size(Hk,3)))                 :: hk_symm_
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_ineq(Hk,Gloc,Sigma,"realaxis",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_ineq(Hk,Gloc,Sigma,"realaxis",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_normal_ineq





  !##################################################################
  !##################################################################
  subroutine dmft_get_gloc_matsubara_normal_cluster(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma ![Nlat:Nlat:Nspin:Nspin:NorbNorb:Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Lk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_cluster(Hk,Gloc,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_cluster(Hk,Gloc,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_normal_cluster

  subroutine dmft_get_gloc_realaxis_normal_cluster(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma ![Nlat:Nlat:Nspin:Nspin:NorbNorb:Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc      !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                     :: hk_symm_  ![Lk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_cluster(Hk,Gloc,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_cluster(Hk,Gloc,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_normal_cluster


  !##################################################################
  !##################################################################


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_gloc_matsubara_normal_cluster_ineq(Hk,Gloc,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)              :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                    :: zeta
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,dimension(size(Hk,3)),optional              :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                       :: hk_symm_  ![Lk]
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_cluster_ineq(Hk,Gloc,Sigma,"matsubara",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_cluster_ineq(Hk,Gloc,Sigma,"matsubara",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_normal_cluster_ineq

  subroutine dmft_get_gloc_realaxis_normal_cluster_ineq(Hk,Gloc,Sigma,zeta,tridiag,hk_symm)
    complex(8),dimension(:,:,:),intent(in)              :: Hk        ![Nineq*Nlat*Nspin*Norb][Nineq*Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma     ![Nineq][Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                    :: zeta
    logical,optional                                    :: tridiag
    logical                                             :: tridiag_
    logical,dimension(size(Hk,3)),optional              :: hk_symm   ![Lk]
    logical,dimension(size(Hk,3))                       :: hk_symm_  ![Lk]
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_cluster_ineq(Hk,Gloc,Sigma,"realaxis",tridiag=tridiag_,hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_cluster_ineq(Hk,Gloc,Sigma,"realaxis",zeta,tridiag=tridiag_,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_normal_cluster_ineq
#endif





  !##################################################################
  !##################################################################


  subroutine dmft_get_gloc_matsubara_normal_gij(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     !      [Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_gij(Hk,Gloc,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_gij(Hk,Gloc,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_normal_gij

  subroutine dmft_get_gloc_realaxis_normal_gij(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sigma     !      [Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,3)),optional            :: hk_symm
    logical,dimension((size(Hk,3)))                   :: hk_symm_
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_normal_gij(Hk,Gloc,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_normal_gij(Hk,Gloc,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_normal_gij


  !##################################################################
  !##################################################################


  subroutine dmft_get_gloc_matsubara_superc_main(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)        :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma

    complex(8),dimension(:),optional                :: zeta
    logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_superc_main(Hk,Gloc,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_superc_main(Hk,Gloc,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_superc_main

  subroutine dmft_get_gloc_realaxis_superc_main(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)        :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_superc_main(Hk,Gloc,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_superc_main(Hk,Gloc,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_superc_main



  !##################################################################
  !##################################################################


  subroutine dmft_get_gloc_matsubara_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,zeta)
    real(8),dimension(:,:,:),intent(in)             :: Ebands    ![2][Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands    ![Nspin*Norb][Lk]/[1][Lk]
    real(8),dimension(2,size(Ebands,1)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    !
    if(present(zeta))then
       call dmft_get_gloc_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"matsubara")
    else
       call dmft_get_gloc_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"matsubara",zeta)
    endif
  end subroutine dmft_get_gloc_matsubara_superc_dos

  subroutine dmft_get_gloc_realaxis_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,zeta)
    real(8),dimension(:,:,:),intent(in)             :: Ebands    ![2][Nspin*Norb][Lk]
    real(8),dimension(:,:),intent(in)               :: Dbands    ![Nspin*Norb][Lk]/[1][Lk]
    real(8),dimension(2,size(Ebands,1)),intent(in)  :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                :: zeta
    !
    if(present(zeta))then
       call dmft_get_gloc_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"realaxis")
    else
       call dmft_get_gloc_superc_dos(Ebands,Dbands,Hloc,Gloc,Sigma,"realaxis",zeta)
    endif
  end subroutine dmft_get_gloc_realaxis_superc_dos



  !##################################################################
  !##################################################################


  subroutine dmft_get_gloc_matsubara_superc_ineq(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_superc_ineq(Hk,Gloc,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_superc_ineq(Hk,Gloc,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_matsubara_superc_ineq

  subroutine dmft_get_gloc_realaxis_superc_ineq(Hk,Gloc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gloc     !as Sigma
    complex(8),dimension(:),optional                  :: zeta
    logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
    logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call dmft_get_gloc_superc_ineq(Hk,Gloc,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call dmft_get_gloc_superc_ineq(Hk,Gloc,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
  end subroutine dmft_get_gloc_realaxis_superc_ineq



  !##################################################################
  !##################################################################


  subroutine dmft_get_gloc_matsubara_superc_gij(Hk,Gloc,Floc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:,:)                     :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gloc           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Floc           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sigma           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm(size(Hk,4))
    logical                                           :: hk_symm_(size(Hk,4))
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call  dmft_get_gloc_superc_gij(Hk,Gloc,Floc,Sigma,"matsubara",hk_symm=hk_symm_)
    else
       call  dmft_get_gloc_superc_gij(Hk,Gloc,Floc,Sigma,"matsubara",zeta,hk_symm=hk_symm_)
    endif
    !
  end subroutine dmft_get_gloc_matsubara_superc_gij

  subroutine dmft_get_gloc_realaxis_superc_gij(Hk,Gloc,Floc,Sigma,zeta,hk_symm)
    complex(8),dimension(:,:,:,:)                     :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
    real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gloc           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Floc           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sigma           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lfreq]
    complex(8),dimension(:),optional                  :: zeta
    logical,optional                                  :: hk_symm(size(Hk,4))
    logical                                           :: hk_symm_(size(Hk,4))
    !
    hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
    !
    if(present(zeta))then
       call  dmft_get_gloc_superc_gij(Hk,Gloc,Floc,Sigma,"realaxis",hk_symm=hk_symm_)
    else
       call  dmft_get_gloc_superc_gij(Hk,Gloc,Floc,Sigma,"realaxis",zeta,hk_symm=hk_symm_)
    endif
    !
  end subroutine dmft_get_gloc_realaxis_superc_gij


















end module GF_GLOC
