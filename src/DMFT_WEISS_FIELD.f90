module DMFT_WEISS_FIELD
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only:reg,txtfy,splot,store_data
  USE SF_ARRAYS,    only:arange
  USE SF_LINALG,    only:eye,inv
  USE SF_MISC,      only:assert_shape
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private


  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block


  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape




  interface dmft_weiss
     module procedure :: dmft_get_weiss_normal_main
     module procedure :: dmft_get_weiss_normal_1b
     module procedure :: dmft_get_weiss_normal_Mb
     module procedure :: dmft_get_weiss_normal_lattice_main
     module procedure :: dmft_get_weiss_normal_lattice_1b
     module procedure :: dmft_get_weiss_normal_lattice_mb
#ifdef _MPI
     module procedure :: dmft_get_weiss_normal_lattice_main_mpi
     module procedure :: dmft_get_weiss_normal_lattice_1b_mpi
     module procedure :: dmft_get_weiss_normal_lattice_mb_mpi
#endif
  end interface dmft_weiss


  interface dmft_weiss_superc
     module procedure :: dmft_get_weiss_superc_main
     module procedure :: dmft_get_weiss_superc_1b
     module procedure :: dmft_get_weiss_superc_Mb
     module procedure :: dmft_get_weiss_superc_lattice_main
     module procedure :: dmft_get_weiss_superc_lattice_1b
     module procedure :: dmft_get_weiss_superc_lattice_Mb
#ifdef _MPI
     module procedure :: dmft_get_weiss_superc_lattice_main_mpi
     module procedure :: dmft_get_weiss_superc_lattice_1b_mpi
     module procedure :: dmft_get_weiss_superc_lattice_Mb_mpi
#endif
  end interface dmft_weiss_superc


  interface dmft_delta
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_1b
     module procedure :: dmft_get_delta_normal_Mb
     module procedure :: dmft_get_delta_normal_lattice_main
     module procedure :: dmft_get_delta_normal_lattice_1b
     module procedure :: dmft_get_delta_normal_lattice_mb
#ifdef _MPI
     module procedure :: dmft_get_delta_normal_lattice_main_mpi
     module procedure :: dmft_get_delta_normal_lattice_1b_mpi
     module procedure :: dmft_get_delta_normal_lattice_mb_mpi
#endif
  end interface dmft_delta

  interface dmft_delta_superc
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_1b
     module procedure :: dmft_get_delta_superc_Mb
     module procedure :: dmft_get_delta_superc_lattice_main
     module procedure :: dmft_get_delta_superc_lattice_1b
     module procedure :: dmft_get_delta_superc_lattice_Mb
#ifdef _MPI
     module procedure :: dmft_get_delta_superc_lattice_main_mpi
     module procedure :: dmft_get_delta_superc_lattice_1b_mpi
     module procedure :: dmft_get_delta_superc_lattice_Mb_mpi
#endif
  end interface dmft_delta_superc


  public :: dmft_weiss
  public :: dmft_delta

  public :: dmft_weiss_superc
  public :: dmft_delta_superc


  real(8),dimension(:),allocatable :: wm !Matsubara frequencies
  character(len=128)               :: suffix
  character(len=128)               :: weiss_suffix=".dat"
  ! complex,parameter                :: zero=dcmplx(0d0,0d0)
  ! complex,parameter                :: xi=dcmplx(0d0,1d0)
  ! complex,parameter                :: one=dcmplx(1d0,0d0)
  ! real(8),parameter                :: pi=3.141592653589793238462643383279502884197d0
  integer                          :: Nlat,Nspin,Norb
  integer                          :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
  !
  integer                          :: mpi_ierr
  integer                          :: mpi_rank
  integer                          :: mpi_size
  logical                          :: mpi_master
  !
  real(8)                          :: beta
  real(8)                          :: xmu


contains



  !--------------------------------------------------------------------!
  !PURPOSE: Get the local Weiss Field calG0 or using self-consistency 
  ! equations and given G_loc and Sigma.
  ! INPUT:
  ! 1. GLOC  : ([2][Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 2. Sigma : ([2][Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 4. Hloc  : local part of the non-interacting Hamiltonian
  ! 5. Iprint: printing flag
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  !
  !                    WEISS FIELD - NORMAL
  !
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
#define FROUTINE 'dmft_get_weiss_normal_main'
  subroutine dmft_get_weiss_normal_main(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable       :: zeta_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: Smats_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: invGloc_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: calG0_site![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                       :: Nspin,Norb,Nso,Lmats
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lmats = size(Gloc,5)
    Nso   = Nspin*Norb
    call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats into a [Norb*Nspin]^2 matrix and create the zeta_site
    do i=1,Lmats
       zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc,Nspin,Norb) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
       invGloc_site(:,:,i) = nn2so_reshape(Gloc(:,:,:,:,i),Nspin,Norb)
       Smats_site(:,:,i)   = nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
    do i=1,Lmats
       call inv(invGloc_site(:,:,i))
    enddo
    !
    ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
    calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
    do i=1,Lmats
       call inv(calG0_site(:,:,i))
    enddo
    !
    !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
    !output structure of [Nspsin,Nspin,Norb,Norb] matrix
    do i=1,Lmats
       Weiss(:,:,:,:,i) = so2nn_reshape(calG0_site(:,:,i),Nspin,Norb)
    enddo
    !
    call dmft_weiss_print_single(wm,Weiss,"WeissG0",iprint)
    !
  end subroutine dmft_get_weiss_normal_main
#undef FROUTINE

#define FROUTINE 'dmft_get_weiss_normal_lattice_main'
  subroutine dmft_get_weiss_normal_lattice_main(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint       
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Weiss_tmp    ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Weiss_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Weiss     = zero
    do ilat=1,Nlat
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       do i=1,Lmats
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc(ilat,:,:,:,:),Nspin,Norb) - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
          invGloc_site(:,:,i) = nn2so_reshape(Gloc(ilat,:,:,:,:,i),Nspin,Norb)
          Smats_site(:,:,i)   = nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
       calG0_site(:,:,:) = invGloc_site(:,:,:) + Smats_site(:,:,:)
       do i=1,Lmats
          call inv(calG0_site(:,:,i))
       enddo
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Weiss(ilat,ispin,jspin,iorb,jorb,1:Lmats) = calG0_site(io,jo,1:Lmats)
                enddo
             enddo
          enddo
       enddo
    end do
    !
    call dmft_weiss_print_lattice(wm,Weiss,"LWeissG0",iprint)
    !
  end subroutine dmft_get_weiss_normal_lattice_main
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_weiss_normal_lattice_main_mpi'
  subroutine dmft_get_weiss_normal_lattice_main_mpi(MpiComm,Gloc,Smats,Weiss,Hloc,iprint)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint       
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Weiss_tmp    ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Weiss_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Weiss_tmp = zero
    Weiss     = zero
    MPIloop: do ilat=1+mpi_rank,Nlat,mpi_size
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       do i=1,Lmats
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc(ilat,:,:,:,:),Nspin,Norb) - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
          invGloc_site(:,:,i) = nn2so_reshape(Gloc(ilat,:,:,:,:,i),Nspin,Norb)
          Smats_site(:,:,i)   = nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
       calG0_site(:,:,:) = invGloc_site(:,:,:) + Smats_site(:,:,:)
       do i=1,Lmats
          call inv(calG0_site(:,:,i))
       enddo
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Weiss_tmp(ilat,ispin,jspin,iorb,jorb,1:Lmats) = calG0_site(io,jo,1:Lmats)
                enddo
             enddo
          enddo
       enddo
    end do MPIloop
    !
    call Mpi_AllReduce(Weiss_tmp, Weiss, size(Weiss), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
    !
    if(mpi_master)call dmft_weiss_print_lattice(wm,Weiss,"LWeissG0",iprint)
    !
  end subroutine dmft_get_weiss_normal_lattice_main_mpi
#undef FROUTINE
#endif


  !##################################################################
  !##################################################################
  !##################################################################


  subroutine dmft_get_weiss_normal_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:),intent(in)             :: Gloc
    complex(8),dimension(size(Gloc)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc)),intent(inout) :: Weiss
    complex(8),dimension(1,1,1,1)                  :: Hloc
    integer                                        :: iprint
    !aux
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Gloc_
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Smats_
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Weiss_
    Gloc_(1,1,1,1,:) = Gloc(:)
    Smats_(1,1,1,1,:) = Smats(:)
    call dmft_get_weiss_normal_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:) = Weiss_(1,1,1,1,:)
  end subroutine dmft_get_weiss_normal_1b

  subroutine dmft_get_weiss_normal_lattice_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:),intent(in)                          :: Gloc
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(inout) :: Weiss
    complex(8),dimension(size(Gloc,1),1,1,1,1)                    :: Hloc
    integer                                                       :: iprint
    !aux
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Weiss_
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_weiss_normal_lattice_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine dmft_get_weiss_normal_lattice_1b

#ifdef _MPI
  subroutine dmft_get_weiss_normal_lattice_1b_mpi(MpiComm,Gloc,Smats,Weiss,Hloc,iprint)
    integer                                                       :: MpiComm
    complex(8),dimension(:,:),intent(in)                          :: Gloc
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(inout) :: Weiss
    complex(8),dimension(size(Gloc,1),1,1,1,1)                    :: Hloc
    integer                                                       :: iprint
    !aux
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Weiss_
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_weiss_normal_lattice_main_mpi(MpiComm,Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine dmft_get_weiss_normal_lattice_1b_mpi
#endif


#define FROUTINE 'dmft_get_weiss_normal_Mb'
  subroutine dmft_get_weiss_normal_Mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:),intent(in)      :: Gloc  ![Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:),intent(in)    :: Hloc  ![Nspin][Nspin][Norb][Norb]
    integer                                     :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:),allocatable :: Gloc_ ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:),allocatable :: Weiss_!
    integer                                     :: Nspin,Norb,Lfreq
    !
    Nspin = 1
    Norb  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    !
    call assert_shape(Gloc,[Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Weiss,[Norb,Norb,Lfreq],FROUTINE,"Weiss")
    call assert_shape(Hloc,[1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(1,1,Norb,Norb,Lfreq))
    allocate(Smats_(1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(1,1,Norb,Norb,Lfreq))
    Gloc_(1,1,:,:,:) = Gloc(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    call dmft_get_weiss_normal_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:) = Weiss_(1,1,:,:,:)
  end subroutine dmft_get_weiss_normal_mb
#undef FROUTINE

#define FROUTINE 'dmft_get_weiss_normal_lattice_mb'
  subroutine dmft_get_weiss_normal_lattice_mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc  ![Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Weiss_!
    integer                                       :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,1)
    Nspin = 1
    Norb  = size(Gloc,2)
    Lfreq = size(Gloc,4)
    call assert_shape(Gloc,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Weiss,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(Nlat,1,1,Norb,Norb,Lfreq))
    !
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_weiss_normal_lattice_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine dmft_get_weiss_normal_lattice_mb
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_weiss_normal_lattice_mb_mpi'
  subroutine dmft_get_weiss_normal_lattice_mb_mpi(MpiComm,Gloc,Smats,Weiss,Hloc,iprint)
    integer                                       :: MpiComm
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc  ![Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Weiss_!
    integer                                       :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,1)
    Nspin = 1
    Norb  = size(Gloc,2)
    Lfreq = size(Gloc,4)
    call assert_shape(Gloc,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Weiss,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(Nlat,1,1,Norb,Norb,Lfreq))
    !
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_weiss_normal_lattice_main_mpi(MpiComm,Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine dmft_get_weiss_normal_lattice_mb_mpi
#undef FROUTINE
#endif




  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  !
  !                    WEISS FIELD - SUPERC
  !
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
#define FROUTINE 'dmft_get_weiss_superc_main'
  subroutine dmft_get_weiss_superc_main(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        !
    complex(8),dimension(:,:,:,:),intent(in)        :: Hloc         ! [Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    integer                                         :: Nspin,Norb,Nso,Nso2,Lmats
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    call assert_shape(Gloc,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Weiss,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
    zeta_site=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          zeta_site(io,io,:)         = xi*wm(:) + xmu 
          zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu 
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Hloc(ispin,jspin,iorb,jorb) - Smats(1,ispin,jspin,iorb,jorb,:)
                zeta_site(io,jo+Nso,:)       =                                                          - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo,:)       =                                                          - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + Hloc(ispin,jspin,iorb,jorb) + conjg(Smats(1,ispin,jspin,iorb,jorb,:))
                !
                invGloc_site(io,jo,:)        = Gloc(1,ispin,jspin,iorb,jorb,:)
                invGloc_site(io,jo+Nso,:)    = Gloc(2,ispin,jspin,iorb,jorb,:)
                invGloc_site(io+Nso,jo,:)    = Gloc(2,ispin,jspin,iorb,jorb,:)
                invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ispin,jspin,iorb,jorb,:))
                !
                Smats_site(io,jo,:)          = Smats(1,ispin,jspin,iorb,jorb,:)
                Smats_site(io,jo+Nso,:)      = Smats(2,ispin,jspin,iorb,jorb,:)
                Smats_site(io+Nso,jo,:)      = Smats(2,ispin,jspin,iorb,jorb,:)
                Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ispin,jspin,iorb,jorb,:))
             enddo
          enddo
       enddo
    enddo
    !
    !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
    do i=1,Lmats
       call inv(invGloc_site(:,:,i))
    enddo
    !
    ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
    calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
    do i=1,Lmats
       call inv(calG0_site(:,:,i))
    enddo
    !
    !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
    !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Weiss(1,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                Weiss(2,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
             enddo
          enddo
       enddo
    enddo
    !
    call dmft_weiss_print_single(wm,Weiss(1,:,:,:,:,:),"WeissG0",iprint)
    call dmft_weiss_print_single(wm,Weiss(2,:,:,:,:,:),"WeissF0",iprint)
    !
  end subroutine dmft_get_weiss_superc_main
#undef FROUTINE

#define FROUTINE 'dmft_get_weiss_superc_lattice_main'
  subroutine dmft_get_weiss_superc_lattice_main(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats        ! 
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                           :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable           :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: Smats_site   !
    complex(8),dimension(:,:,:),allocatable           :: invGloc_site !
    complex(8),dimension(:,:,:),allocatable           :: calG0_site   !
    integer                                           :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
    integer                                           :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lmats = size(Gloc,7)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Weiss,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Weiss       = zero
    do ilat=1,Nlat
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,ispin,iorb,iorb)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,ispin,iorb,iorb)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   !
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(1,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
       calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
       do i=1,Lmats
          call inv(calG0_site(:,:,i))
       enddo
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Weiss(1,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   Weiss(2,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    end do
    !
    call dmft_weiss_print_lattice(wm,Weiss(1,:,:,:,:,:,:),"LWeissG0",iprint)
    call dmft_weiss_print_lattice(wm,Weiss(2,:,:,:,:,:,:),"LWeissF0",iprint)
    !
  end subroutine dmft_get_weiss_superc_lattice_main
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_weiss_superc_lattice_main_mpi'
  subroutine dmft_get_weiss_superc_lattice_main_mpi(MpiComm,Gloc,Smats,Weiss,Hloc,iprint)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats        ! 
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                           :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Weiss_tmp        ! 
    complex(8),dimension(:,:,:),allocatable           :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: Smats_site   !
    complex(8),dimension(:,:,:),allocatable           :: invGloc_site !
    complex(8),dimension(:,:,:),allocatable           :: calG0_site   !
    integer                                           :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
    integer                                           :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lmats = size(Gloc,7)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Weiss,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Weiss_tmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Weiss_tmp   = zero
    Weiss       = zero
    MPIloop: do ilat=1+mpi_rank,Nlat,mpi_size
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,ispin,iorb,iorb)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,ispin,iorb,iorb)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   !
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(1,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
       calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
       do i=1,Lmats
          call inv(calG0_site(:,:,i))
       enddo
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Weiss_tmp(1,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   Weiss_tmp(2,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    end do MPIloop
    !
    call MPI_Allreduce(Weiss_tmp, Weiss, size(Weiss), MPI_DOUBLE_COMPLEX, MPI_SUM, MpiComm, mpi_ierr)
    !
    if(mpi_master)then
       call dmft_weiss_print_lattice(wm,Weiss(1,:,:,:,:,:,:),"LWeissG0",iprint)
       call dmft_weiss_print_lattice(wm,Weiss(2,:,:,:,:,:,:),"LWeissF0",iprint)
    endif
    !
  end subroutine dmft_get_weiss_superc_lattice_main_mpi
#undef FROUTINE
#endif


  !##################################################################
  !##################################################################
  !##################################################################

#define FROUTINE 'dmft_get_weiss_superc_1b'
  subroutine dmft_get_weiss_superc_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:),intent(in)               :: Gloc
    complex(8),dimension(2,size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2)),intent(inout) :: Weiss
    complex(8),dimension(1,1,1,1)                      :: Hloc
    integer                                            :: iprint
    !aux
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Weiss_
    call assert_shape(Gloc,[2,size(Gloc,2)],FROUTINE,"Gloc")
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_weiss_superc_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:) = Weiss_(:,1,1,1,1,:)
  end subroutine dmft_get_weiss_superc_1b
#undef FROUTINE

#define FROUTINE 'dmft_get_weiss_superc_lattice_1b'
  subroutine dmft_get_weiss_superc_lattice_1b(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:),intent(in)                          :: Gloc ![2][Nlat][Lmats]
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(inout) :: Weiss
    complex(8),dimension(size(Gloc,2),1,1,1,1)                      :: Hloc
    integer                                                         :: iprint
    !aux
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Gloc_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Smats_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Weiss_
    call assert_shape(Gloc,[2,size(Gloc,2),size(Gloc,3)],FROUTINE,"Gloc")
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    call dmft_get_weiss_superc_lattice_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:) = Weiss_(:,:,1,1,1,1,:)
  end subroutine dmft_get_weiss_superc_lattice_1b
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_weiss_superc_lattice_1b_mpi'
  subroutine dmft_get_weiss_superc_lattice_1b_mpi(MpiComm,Gloc,Smats,Weiss,Hloc,iprint)
    integer                                                         :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                          :: Gloc ![2][Nlat][Lmats]
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(inout) :: Weiss
    complex(8),dimension(size(Gloc,2),1,1,1,1)                      :: Hloc
    integer                                                         :: iprint
    !aux
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Gloc_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Smats_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Weiss_
    call assert_shape(Gloc,[2,size(Gloc,2),size(Gloc,3)],FROUTINE,"Gloc")
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    call dmft_get_weiss_superc_lattice_main_mpi(MpiComm,Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:) = Weiss_(:,:,1,1,1,1,:)
  end subroutine dmft_get_weiss_superc_lattice_1b_mpi
#undef FROUTINE
#endif

#define FROUTINE 'dmft_get_weiss_superc_Mb'
  subroutine dmft_get_weiss_superc_Mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc   ![2][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats  !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Weiss  !
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc   ![Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_  ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_ !
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Weiss_ !
    integer                                       :: Nspin,Norb,Lfreq
    !
    Nspin=1
    Norb =size(Gloc,2)
    Lfreq=size(Gloc,4)
    call assert_shape(Gloc,[2,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Weiss,[2,Norb,Norb,Lfreq],FROUTINE,"Weiss")
    call assert_shape(Hloc,[1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(2,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(2,1,1,Norb,Norb,Lfreq))
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_weiss_superc_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:) = Weiss_(:,1,1,:,:,:)
  end subroutine dmft_get_weiss_superc_Mb
#undef FROUTINE

#define FROUTINE 'dmft_get_weiss_superc_lattice_Mb'
  subroutine dmft_get_weiss_superc_lattice_Mb(Gloc,Smats,Weiss,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc  ![2][Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Gloc_ ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Weiss_!
    integer                                         :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,2)
    Nspin = 1
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)
    call assert_shape(Gloc,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Weiss,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(2,Nlat,1,1,Norb,Norb,Lfreq))
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call dmft_get_weiss_superc_lattice_main(Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:,:) = Weiss_(:,:,1,1,:,:,:)
  end subroutine dmft_get_weiss_superc_lattice_Mb
#undef FROUTINE


#ifdef _MPI
#define FROUTINE 'dmft_get_weiss_superc_lattice_Mb_mpi'
  subroutine dmft_get_weiss_superc_lattice_Mb_mpi(MpiComm,Gloc,Smats,Weiss,Hloc,iprint)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc  ![2][Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:,:),intent(inout)   :: Weiss !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Gloc_ ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Weiss_!
    integer                                         :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,2)
    Nspin = 1
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)
    call assert_shape(Gloc,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Weiss,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Weiss")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Weiss_(2,Nlat,1,1,Norb,Norb,Lfreq))
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call dmft_get_weiss_superc_lattice_main_mpi(MpiComm,Gloc_,Smats_,Weiss_,Hloc,iprint)
    Weiss(:,:,:,:,:) = Weiss_(:,:,1,1,:,:,:)
  end subroutine dmft_get_weiss_superc_lattice_Mb_mpi
#undef FROUTINE
#endif






















  !--------------------------------------------------------------------!
  !PURPOSE: Get the local Hybridization \Delta functino using 
  ! self-consistency equations and given G_loc and Sigma.
  ! INPUT:
  ! 1. GLOC  : ([2][Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 2. Sigma : ([2][Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 4. Hloc  : local part of the non-interacting Hamiltonian
  ! 5. Iprint: printing flag
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  !
  !                    DELTA FIELD - NORMAL
  !
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
#define FROUTINE 'dmft_get_delta_normal_main'
  subroutine dmft_get_delta_normal_main(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable       :: zeta_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: Smats_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: invGloc_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: calG0_site![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                       :: Nspin,Norb,Nso,Lmats
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lmats = size(Gloc,5)
    Nso   = Nspin*Norb
    call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Delta,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats into a [Norb*Nspin]^2 matrix and create the zeta_site
    do i=1,Lmats
       zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc,Nspin,Norb) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
       invGloc_site(:,:,i) = nn2so_reshape(Gloc(:,:,:,:,i),Nspin,Norb)
       Smats_site(:,:,i)   = nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
    do i=1,Lmats
       call inv(invGloc_site(:,:,i))
    enddo
    !
    ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Gloc]_ilat^-1
    calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
    !
    !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
    !output structure of [Nspsin,Nspin,Norb,Norb] matrix
    do i=1,Lmats
       Delta(:,:,:,:,i) = so2nn_reshape(calG0_site(:,:,i),Nspin,Norb)
    enddo
    !
    call dmft_weiss_print_single(wm,Delta,"DeltaG0",iprint)
    !
  end subroutine dmft_get_delta_normal_main
#undef FROUTINE

#define FROUTINE 'dmft_get_delta_normal_lattice_main'
  subroutine dmft_get_delta_normal_lattice_main(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint       
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Delta_tmp    ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Delta_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta     = zero
    do ilat=1,Nlat
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       do i=1,Lmats
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc(ilat,:,:,:,:),Nspin,Norb) - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
          invGloc_site(:,:,i) = nn2so_reshape(Gloc(ilat,:,:,:,:,i),Nspin,Norb)
          Smats_site(:,:,i)   = nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,:) = zeta_site(:,:,:) - invGloc_site(:,:,:)
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Delta(ilat,ispin,jspin,iorb,jorb,1:Lmats) = calG0_site(io,jo,1:Lmats)
                enddo
             enddo
          enddo
       enddo
    end do
    !
    call dmft_weiss_print_lattice(wm,Delta,"LDeltaG0",iprint)
    !
  end subroutine dmft_get_delta_normal_lattice_main
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_delta_normal_lattice_main_mpi'
  subroutine dmft_get_delta_normal_lattice_main_mpi(MpiComm,Gloc,Smats,Delta,Hloc,iprint)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint       
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Delta_tmp    ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Delta_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso,Nso,Lmats))
    allocate(Smats_site(Nso,Nso,Lmats))
    allocate(invGloc_site(Nso,Nso,Lmats))
    allocate(calG0_site(Nso,Nso,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta_tmp = zero
    Delta     = zero
    MPIloop: do ilat=1+mpi_rank,Nlat,mpi_size
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       do i=1,Lmats
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc(ilat,:,:,:,:),Nspin,Norb) - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
          invGloc_site(:,:,i) = nn2so_reshape(Gloc(ilat,:,:,:,:,i),Nspin,Norb)
          Smats_site(:,:,i)   = nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,:) = zeta_site(:,:,:) - invGloc_site(:,:,:)
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Delta_tmp(ilat,ispin,jspin,iorb,jorb,1:Lmats) = calG0_site(io,jo,1:Lmats)
                enddo
             enddo
          enddo
       enddo
    end do MPIloop
    !
    call Mpi_AllReduce(Delta_tmp, Delta, size(Delta), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
    !
    if(mpi_master)call dmft_weiss_print_lattice(wm,Delta,"LDeltaG0",iprint)
    !
  end subroutine dmft_get_delta_normal_lattice_main_mpi
#undef FROUTINE
#endif


  !##################################################################
  !##################################################################
  !##################################################################


  subroutine dmft_get_delta_normal_1b(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:),intent(in)             :: Gloc
    complex(8),dimension(size(Gloc)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc)),intent(inout) :: Delta
    complex(8),dimension(1,1,1,1)                  :: Hloc
    integer                                        :: iprint
    !aux
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Gloc_
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Smats_
    complex(8),dimension(1,1,1,1,size(Gloc))       :: Delta_
    Gloc_(1,1,1,1,:) = Gloc(:)
    Smats_(1,1,1,1,:) = Smats(:)
    call dmft_get_delta_normal_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:) = Delta_(1,1,1,1,:)
  end subroutine dmft_get_delta_normal_1b

  subroutine dmft_get_delta_normal_lattice_1b(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:),intent(in)                          :: Gloc
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(inout) :: Delta
    complex(8),dimension(size(Gloc,1),1,1,1,1)                    :: Hloc
    integer                                                       :: iprint
    !aux
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Delta_
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_delta_normal_lattice_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:) = Delta_(:,1,1,1,1,:)
  end subroutine dmft_get_delta_normal_lattice_1b

#ifdef _MPI
  subroutine dmft_get_delta_normal_lattice_1b_mpi(MpiComm,Gloc,Smats,Delta,Hloc,iprint)
    integer                                                       :: MpiComm
    complex(8),dimension(:,:),intent(in)                          :: Gloc
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(size(Gloc,1),size(Gloc,2)),intent(inout) :: Delta
    complex(8),dimension(size(Gloc,1),1,1,1,1)                    :: Hloc
    integer                                                       :: iprint
    !aux
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(size(Gloc,1),1,1,1,1,size(Gloc,2))       :: Delta_
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_delta_normal_lattice_main_mpi(MpiComm,Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:) = Delta_(:,1,1,1,1,:)
  end subroutine dmft_get_delta_normal_lattice_1b_mpi
#endif


#define FROUTINE 'dmft_get_delta_normal_Mb'
  subroutine dmft_get_delta_normal_Mb(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:),intent(in)      :: Gloc  ![Norb][Norb][Lmats]
    complex(8),dimension(:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:),intent(inout)   :: Delta !
    complex(8),dimension(:,:,:,:),intent(in)    :: Hloc  ![Nspin][Nspin][Norb][Norb]
    integer                                     :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:),allocatable :: Gloc_ ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:),allocatable :: Delta_!
    integer                                     :: Nspin,Norb,Lfreq
    !
    Nspin = 1
    Norb  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    !
    call assert_shape(Gloc,[Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Delta,[Norb,Norb,Lfreq],FROUTINE,"Delta")
    call assert_shape(Hloc,[1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(1,1,Norb,Norb,Lfreq))
    allocate(Smats_(1,1,Norb,Norb,Lfreq))
    allocate(Delta_(1,1,Norb,Norb,Lfreq))
    Gloc_(1,1,:,:,:) = Gloc(:,:,:)
    Smats_(1,1,:,:,:) = Smats(:,:,:)
    call dmft_get_delta_normal_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:) = Delta_(1,1,:,:,:)
  end subroutine dmft_get_delta_normal_mb
#undef FROUTINE

#define FROUTINE 'dmft_get_delta_normal_lattice_mb'
  subroutine dmft_get_delta_normal_lattice_mb(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc  ![Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Delta !
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Delta_!
    integer                                       :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,1)
    Nspin = 1
    Norb  = size(Gloc,2)
    Lfreq = size(Gloc,4)
    call assert_shape(Gloc,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Delta,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Delta_(Nlat,1,1,Norb,Norb,Lfreq))
    !
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_delta_normal_lattice_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:,:) = Delta_(:,1,1,:,:,:)
  end subroutine dmft_get_delta_normal_lattice_mb
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_delta_normal_lattice_mb_mpi'
  subroutine dmft_get_delta_normal_lattice_mb_mpi(MpiComm,Gloc,Smats,Delta,Hloc,iprint)
    integer                                       :: MpiComm
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc  ![Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Delta !
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_ ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Delta_!
    integer                                       :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,1)
    Nspin = 1
    Norb  = size(Gloc,2)
    Lfreq = size(Gloc,4)
    call assert_shape(Gloc,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Delta,[Nlat,Norb,Norb,Lfreq],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Delta_(Nlat,1,1,Norb,Norb,Lfreq))
    !
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_delta_normal_lattice_main_mpi(MpiComm,Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:,:) = Delta_(:,1,1,:,:,:)
  end subroutine dmft_get_delta_normal_lattice_mb_mpi
#undef FROUTINE
#endif






  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  !
  !                 DELTA FUNCTION - SUPERC
  !
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
#define FROUTINE 'dmft_get_delta_superc_main'
  subroutine dmft_get_delta_superc_main(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta        !
    complex(8),dimension(:,:,:,:),intent(in)        :: Hloc         ! [Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    integer                                         :: Nspin,Norb,Nso,Nso2,Lmats
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    call assert_shape(Gloc,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Delta,[2,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
    zeta_site=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          io = iorb + (ispin-1)*Norb
          zeta_site(io,io,:)         = xi*wm(:) + xmu 
          zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu 
       enddo
    enddo
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Hloc(ispin,jspin,iorb,jorb) - Smats(1,ispin,jspin,iorb,jorb,:)
                zeta_site(io,jo+Nso,:)       =                                                          - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo,:)       =                                                          - Smats(2,ispin,jspin,iorb,jorb,:)
                zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + Hloc(ispin,jspin,iorb,jorb) + conjg(Smats(1,ispin,jspin,iorb,jorb,:))
                !
                invGloc_site(io,jo,:)        = Gloc(1,ispin,jspin,iorb,jorb,:)
                invGloc_site(io,jo+Nso,:)    = Gloc(2,ispin,jspin,iorb,jorb,:)
                invGloc_site(io+Nso,jo,:)    = Gloc(2,ispin,jspin,iorb,jorb,:)
                invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ispin,jspin,iorb,jorb,:))
                !
                Smats_site(io,jo,:)          = Smats(1,ispin,jspin,iorb,jorb,:)
                Smats_site(io,jo+Nso,:)      = Smats(2,ispin,jspin,iorb,jorb,:)
                Smats_site(io+Nso,jo,:)      = Smats(2,ispin,jspin,iorb,jorb,:)
                Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ispin,jspin,iorb,jorb,:))
             enddo
          enddo
       enddo
    enddo
    !
    !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
    do i=1,Lmats
       call inv(invGloc_site(:,:,i))
    enddo
    !
    ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
    calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
    !
    !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
    !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = iorb + (ispin-1)*Norb
                jo = jorb + (jspin-1)*Norb
                Delta(1,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                Delta(2,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
             enddo
          enddo
       enddo
    enddo
    !
    call dmft_weiss_print_single(wm,Delta(1,:,:,:,:,:),"DeltaG0",iprint)
    call dmft_weiss_print_single(wm,Delta(2,:,:,:,:,:),"DeltaF0",iprint)
    !
  end subroutine dmft_get_delta_superc_main
#undef FROUTINE

#define FROUTINE 'dmft_get_delta_superc_lattice_main'
  subroutine dmft_get_delta_superc_lattice_main(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats        ! 
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                           :: iprint
    !aux
    complex(8),dimension(:,:,:),allocatable           :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: Smats_site   !
    complex(8),dimension(:,:,:),allocatable           :: invGloc_site !
    complex(8),dimension(:,:,:),allocatable           :: calG0_site   !
    integer                                           :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
    integer                                           :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lmats = size(Gloc,7)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Delta,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta       = zero
    do ilat=1,Nlat
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,ispin,iorb,iorb)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,ispin,iorb,iorb)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   !
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(1,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Delta(1,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   Delta(2,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    end do
    !
    call dmft_weiss_print_lattice(wm,Delta(1,:,:,:,:,:,:),"LDeltaG0",iprint)
    call dmft_weiss_print_lattice(wm,Delta(2,:,:,:,:,:,:),"LDeltaF0",iprint)
    !
  end subroutine dmft_get_delta_superc_lattice_main
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_delta_superc_lattice_main_mpi'
  subroutine dmft_get_delta_superc_lattice_main_mpi(MpiComm,Gloc,Smats,Delta,Hloc,iprint)
    integer                                           :: MpiComm
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc         ! [2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats        ! 
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    integer                                           :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Delta_tmp        ! 
    complex(8),dimension(:,:,:),allocatable           :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: Smats_site   !
    complex(8),dimension(:,:,:),allocatable           :: invGloc_site !
    complex(8),dimension(:,:,:),allocatable           :: calG0_site   !
    integer                                           :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
    integer                                           :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
    !
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank =  MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    !Testing part:
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lmats = size(Gloc,7)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
    call assert_shape(Delta,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],FROUTINE,"Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Delta_tmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta_tmp   = zero
    Delta       = zero
    MPIloop: do ilat=1+mpi_rank,Nlat,mpi_size
       !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
       zeta_site=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,ispin,iorb,iorb)
             zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,ispin,iorb,iorb)
          enddo
       enddo
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   !
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =-Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(1,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Gloc(2,ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(1,ilat,ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(1,ilat,ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Hloc]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
       !
       !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io = iorb + (ispin-1)*Norb
                   jo = jorb + (jspin-1)*Norb
                   Delta_tmp(1,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   Delta_tmp(2,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    end do MPIloop
    !
    call MPI_Allreduce(Delta_tmp, Delta, size(Delta), MPI_DOUBLE_COMPLEX, MPI_SUM, MpiComm, mpi_ierr)
    !
    if(mpi_master)then
       call dmft_weiss_print_lattice(wm,Delta(1,:,:,:,:,:,:),"LDeltaG0",iprint)
       call dmft_weiss_print_lattice(wm,Delta(2,:,:,:,:,:,:),"LDeltaF0",iprint)
    endif
    !
  end subroutine dmft_get_delta_superc_lattice_main_mpi
#undef FROUTINE
#endif


  !##################################################################
  !##################################################################
  !##################################################################


#define FROUTINE 'dmft_get_delta_superc_1b'
  subroutine dmft_get_delta_superc_1b(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:),intent(in)               :: Gloc
    complex(8),dimension(2,size(Gloc,2)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2)),intent(inout) :: Delta
    complex(8),dimension(1,1,1,1)                      :: Hloc
    integer                                            :: iprint
    !aux
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Gloc_
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Smats_
    complex(8),dimension(2,1,1,1,1,size(Gloc,2))       :: Delta_
    call assert_shape(Gloc,[2,size(Gloc,2)],FROUTINE,"Gloc")
    Gloc_(:,1,1,1,1,:) = Gloc(:,:)
    Smats_(:,1,1,1,1,:) = Smats(:,:)
    call dmft_get_delta_superc_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:) = Delta_(:,1,1,1,1,:)
  end subroutine dmft_get_delta_superc_1b
#undef FROUTINE

#define FROUTINE 'dmft_get_delta_superc_lattice_1b'
  subroutine dmft_get_delta_superc_lattice_1b(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:),intent(in)                          :: Gloc ![2][Nlat][Lmats]
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(inout) :: Delta
    complex(8),dimension(size(Gloc,2),1,1,1,1)                      :: Hloc
    integer                                                         :: iprint
    !aux
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Gloc_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Smats_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Delta_
    call assert_shape(Gloc,[2,size(Gloc,2),size(Gloc,3)],FROUTINE,"Gloc")
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    call dmft_get_delta_superc_lattice_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:) = Delta_(:,:,1,1,1,1,:)
  end subroutine dmft_get_delta_superc_lattice_1b
#undef FROUTINE

#ifdef _MPI
#define FROUTINE 'dmft_get_delta_superc_lattice_1b_mpi'
  subroutine dmft_get_delta_superc_lattice_1b_mpi(MpiComm,Gloc,Smats,Delta,Hloc,iprint)
    integer :: MpiComm
    complex(8),dimension(:,:,:),intent(in)                          :: Gloc ![2][Nlat][Lmats]
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(in)    :: Smats
    complex(8),dimension(2,size(Gloc,2),size(Gloc,3)),intent(inout) :: Delta
    complex(8),dimension(size(Gloc,2),1,1,1,1)                      :: Hloc
    integer                                                         :: iprint
    !aux
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Gloc_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Smats_
    complex(8),dimension(2,size(Gloc,2),1,1,1,1,size(Gloc,3))       :: Delta_
    call assert_shape(Gloc,[2,size(Gloc,2),size(Gloc,3)],FROUTINE,"Gloc")
    Gloc_(:,:,1,1,1,1,:) = Gloc(:,:,:)
    Smats_(:,:,1,1,1,1,:) = Smats(:,:,:)
    call dmft_get_delta_superc_lattice_main_mpi(MpiComm,Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:) = Delta_(:,:,1,1,1,1,:)
  end subroutine dmft_get_delta_superc_lattice_1b_mpi
#undef FROUTINE
#endif

#define FROUTINE 'dmft_get_delta_superc_Mb'
  subroutine dmft_get_delta_superc_Mb(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:),intent(in)      :: Gloc   ![2][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Smats  !
    complex(8),dimension(:,:,:,:),intent(inout)   :: Delta  !
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc   ![Nspin][Nspin][Norb][Norb]
    integer                                       :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gloc_  ![2][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Smats_ !
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Delta_ !
    integer                                       :: Nspin,Norb,Lfreq
    !
    Nspin=1
    Norb =size(Gloc,2)
    Lfreq=size(Gloc,4)
    call assert_shape(Gloc,[2,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Delta,[2,Norb,Norb,Lfreq],FROUTINE,"Delta")
    call assert_shape(Hloc,[1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(2,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,1,1,Norb,Norb,Lfreq))
    allocate(Delta_(2,1,1,Norb,Norb,Lfreq))
    Gloc_(:,1,1,:,:,:) = Gloc(:,:,:,:)
    Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
    call dmft_get_delta_superc_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:,:) = Delta_(:,1,1,:,:,:)
  end subroutine dmft_get_delta_superc_Mb
#undef FROUTINE

#define FROUTINE 'dmft_get_delta_superc_lattice_Mb'
  subroutine dmft_get_delta_superc_lattice_Mb(Gloc,Smats,Delta,Hloc,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc  ![2][Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:,:),intent(inout)   :: Delta !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Gloc_ ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Delta_!
    integer                                         :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,2)
    Nspin = 1
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)
    call assert_shape(Gloc,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Delta,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Delta_(2,Nlat,1,1,Norb,Norb,Lfreq))
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call dmft_get_delta_superc_lattice_main(Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:,:,:) = Delta_(:,:,1,1,:,:,:)
  end subroutine dmft_get_delta_superc_lattice_Mb
#undef FROUTINE


#ifdef _MPI
#define FROUTINE 'dmft_get_delta_superc_lattice_Mb_mpi'
  subroutine dmft_get_delta_superc_lattice_Mb_mpi(MpiComm,Gloc,Smats,Delta,Hloc,iprint)
    integer                                         :: MpiComm
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc  ![2][Nlat][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats !
    complex(8),dimension(:,:,:,:,:),intent(inout)   :: Delta !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc  ![Nlat][Nspin][Nspin][Norb][Norb]
    integer                                         :: iprint
    !aux
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Gloc_ ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Smats_!
    complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Delta_!
    integer                                         :: Nlat,Nspin,Norb,Lfreq
    !
    Nlat  = size(Gloc,2)
    Nspin = 1
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)
    call assert_shape(Gloc,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Gloc")
    call assert_shape(Smats,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Smats")
    call assert_shape(Delta,[2,Nlat,Norb,Norb,Lfreq],FROUTINE,"Delta")
    call assert_shape(Hloc,[Nlat,1,1,Norb,Norb],FROUTINE,"Hloc")
    allocate(Gloc_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Smats_(2,Nlat,1,1,Norb,Norb,Lfreq))
    allocate(Delta_(2,Nlat,1,1,Norb,Norb,Lfreq))
    Gloc_(:,:,1,1,:,:,:) = Gloc(:,:,:,:,:)
    Smats_(:,:,1,1,:,:,:) = Smats(:,:,:,:,:)
    call dmft_get_delta_superc_lattice_main_mpi(MpiComm,Gloc_,Smats_,Delta_,Hloc,iprint)
    Delta(:,:,:,:,:) = Delta_(:,:,1,1,:,:,:)
  end subroutine dmft_get_delta_superc_lattice_Mb_mpi
#undef FROUTINE
#endif













  !##################################################################
  !##################################################################
  !##################################################################

  !##################################################################
  !##################################################################
  !##################################################################

  !##################################################################
  !##################################################################
  !##################################################################



  !SET THE GLOC SUFFIX
  subroutine dmft_weiss_set_suffix(string)
    character(len=*) :: string
    weiss_suffix=reg(string)
  end subroutine dmft_weiss_set_suffix




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: print a local GF according to iprint variable
  !+-----------------------------------------------------------------------------+!
  subroutine dmft_weiss_print_single(w,Green,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Green ![Nspin][Nspin][Norb][Norb][Lfreq]
    real(8),dimension(size(Green,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    integer                                    :: Nspin,Norb,ispin,jspin,iorb,jorb
    !
    Nspin = size(Green,1)
    Norb  = size(Green,3)
    call assert_shape(Green,[Nspin,Nspin,Norb,Norb,size(Green,5)],"print_weiss","Green")
    !
    !
    select case(iprint)
    case (0)
       write(*,*)"Delta//WeissField not written to file."
    case(1)
       write(*,"(A)") "Delta//WeissField: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(weiss_suffix)
             call splot(reg(suffix),wm,Green(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)
       write(*,"(A)") "Delta//WeissField: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(weiss_suffix)
                call splot(reg(suffix),wm,Green(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default
       write(*,"(A)") "Delta//WeissField: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(weiss_suffix)
                   call splot(reg(suffix),wm,Green(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine dmft_weiss_print_single




  !PRINT GLOC LATTICE
  subroutine dmft_weiss_print_lattice(w,Green,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Green
    real(8),dimension(size(Green,6))              :: w        
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    !
    Nlat  = size(Green,1)
    Nspin = size(Green,2)
    Norb  = size(Green,4)
    call assert_shape(Green,[Nlat,Nspin,Nspin,Norb,Norb,size(Green,6)],"dmft_weiss_print_lattice",fname)
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Delta//WeissField: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,*)"Delta//WeissField: write spin-orbital diagonal elements No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(weiss_suffix)
             call store_data(reg(suffix),Green(:,ispin,ispin,iorb,iorb,:),w)
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Delta//WeissField: write spin diagonal and all orbitals elements No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(weiss_suffix)
                call store_data(reg(suffix),Green(:,ispin,ispin,iorb,jorb,:),w)
             enddo
          enddo
       enddo
       !
    case(3)                  !print all off-diagonals
       write(*,*)"Delta//WeissField: write all elements No Split."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(weiss_suffix)
                   call store_data(reg(suffix),Green(:,ispin,jspin,iorb,jorb,:),w)
                enddo
             enddo
          enddo
       enddo
       !
    case(4)                  !print only diagonal elements
       write(*,*)"Delta//WeissField: write spin-orbital diagonal elements Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw_indx"//reg(txtfy(ilat,6))//reg(weiss_suffix)
                call splot(reg(suffix),w,Green(ilat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
       !
    case(5)                  !print spin-diagonal, all orbitals 
       write(*,*)"Delta//WeissField: write spin diagonal and all orbitals elements Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//&
                        "_iw_indx"//reg(txtfy(ilat,6))//reg(weiss_suffix)
                   call splot(reg(suffix),w,Green(ilat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    case default
       write(*,*)"Delta//WeissField: write all elements Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                           "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                           "_iw_indx"//reg(txtfy(ilat,6))//reg(weiss_suffix)
                      call splot(reg(suffix),w,Green(ilat,ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine dmft_weiss_print_lattice




#ifdef _MPI
  function MPI_Get_size(comm) result(size)
    integer :: comm
    integer :: size,ierr
    call MPI_Comm_size(comm,size,ierr)
  end function MPI_Get_size



  function MPI_Get_rank(comm) result(rank)
    integer :: comm
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
  end function MPI_Get_rank



  function MPI_Get_master(comm) result(master)
    integer :: comm
    logical :: master
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    master=.false.
    if(rank==0)master=.true.
  end function MPI_Get_master
#endif




  !####################################################################
  !                    computational ROUTINES
  !####################################################################
  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix,Nlat,Nso) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks






  !--------------------------------------------------------------------!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !--------------------------------------------------------------------!
  function select_block_Nlso(ip,Matrix,Nlat,Nso) result(Vblock)
    integer                                 :: ip
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nso,Nso)           :: Vblock
    integer                                 :: i,j
    Vblock=zero
    i = 1+(ip-1)*Nso
    j =       ip*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ip,Matrix,Nlat,Nspin,Norb) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block_nnn








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso


end module DMFT_WEISS_FIELD
