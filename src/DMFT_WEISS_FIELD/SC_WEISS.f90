module SC_WEISS
  USE DMFT_CTRL_VARS
  USE SC_COMMON
  implicit none
  private


  public :: dmft_get_weiss_normal_main
  public :: dmft_get_weiss_normal_rank4
  public :: dmft_get_weiss_normal_rank5
  public :: dmft_get_weiss_normal_rank6
  public :: dmft_get_weiss_normal_rank7

  public :: dmft_get_weiss_superc_main
  public :: dmft_get_weiss_superc_rank4
  public :: dmft_get_weiss_superc_rank5
  public :: dmft_get_weiss_superc_rank6



contains


  !--------------------------------------------------------------------!
  !PURPOSE: Get the local Weiss Field calG0 or using self-consistency 
  ! equations and given G_loc/F_loc and Sigma/Self
  !--------------------------------------------------------------------!
  subroutine dmft_get_weiss_normal_main(Gloc,Sigma,Weiss,Nsites) !N=Nsites*Nso
    complex(8),dimension(:,:,:),intent(in)    :: Gloc  ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:,:),intent(in)    :: Sigma ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:,:),intent(inout) :: Weiss ! [N,Nsites*Nso][L]
    integer,optional                          :: Nsites
    integer                                   :: Nsites_,Nso
    !aux
    complex(8),dimension(:,:,:),allocatable   :: Weiss_tmp    ![N,N][L]
    complex(8),dimension(:,:,:),allocatable   :: calG0        ![Nsites,Nso,Nso]
    complex(8),dimension(:,:),allocatable     :: Sigma_site   ![Nso,Nso]
    complex(8),dimension(:,:),allocatable     :: invG_site
    complex(8),dimension(:,:),allocatable     :: invG0_site
    !
    Nsites_=1;if(present(Nsites))Nsites_=Nsites
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
    !Testing part:
    Ntot  = size(Gloc,1);
    if(mod(Ntot,Nsites_)/=0)stop "dmft_get_weiss_normal_main: Ntot%Nsites != 0"
    Lfreq = size(Gloc,3)
    Nso   = Ntot/Nsites_
    call assert_shape(Gloc,[Ntot,Nsites_*Nso,Lfreq],"dmft_get_weiss_normal_main","Gloc")
    call assert_shape(Sigma,[Ntot,Nsites_*Nso,Lfreq],"dmft_get_weiss_normal_main","Sigma")
    call assert_shape(Weiss,[Ntot,Nsites_*Nso,Lfreq],"dmft_get_weiss_normal_main","Weiss")
    !
    !Build array to extract diagonal blocks
    allocate(Weiss_tmp(Ntot,Ntot,Lfreq))
    !
    allocate(calG0(Nsites_,Nso,Nso))
    allocate(Sigma_site(Nso,Nso))
    allocate(invG_site(Nso,Nso))
    allocate(invG0_site(Nso,Nso))
    !
    !Work out the frequencies in parallel:
    Weiss_tmp=zero
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       !
       !For a fixed frequency update the local Weiss field
       do ilat=1,Nsites_
          Sigma_site = select_block(ilat,Sigma(:,:,i),Nsites_,Nso)

          invG_site  = select_block(ilat,Gloc(:,:,i),Nsites_,Nso)
          call inv(invG_site)
          !
          ![G0]^-1 = [ [Gloc]^-1 + [Sigma] ]^-1
          invG0_site = invG_site + Sigma_site
          call inv(invG0_site)
          !
          calG0(ilat,:,:) = invG0_site
       enddo
       !Once all sites are obtained we dump back 
       Weiss_tmp(:,:,i) = blocks_to_matrix(calG0,Nsites_,Nso)
    enddo MPIloop
#ifdef _MPI    
    if(check_MPI())then
       Weiss = zero
       call Mpi_AllReduce(Weiss_tmp, Weiss, size(Weiss), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Weiss=Weiss_tmp
    endif
#else
    Weiss=Weiss_tmp
#endif
    !
  end subroutine dmft_get_weiss_normal_main
  


  subroutine dmft_get_weiss_superc_main(Gloc,Floc,Sigma,Self,Weiss,Theta,Nsites) !N=Nsites*Nso
    complex(8),dimension(:,:,:),intent(in)    :: Gloc  ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:,:),intent(in)    :: Floc  ! ..
    complex(8),dimension(:,:,:),intent(in)    :: Sigma 
    complex(8),dimension(:,:,:),intent(in)    :: Self  
    complex(8),dimension(:,:,:),intent(inout) :: Weiss 
    complex(8),dimension(:,:,:),intent(inout) :: Theta 
    integer,optional                          :: Nsites
    integer                                   :: Nsites_,Nso
    !aux
    complex(8),dimension(:,:,:),allocatable   :: Green_site  !(2,Nso,Nso)
    complex(8),dimension(:,:,:),allocatable   :: Sigma_site  !(2,Nso,Nso)
    complex(8),dimension(:,:),allocatable     :: invG_nambu  !(2*Nso,2*Nso)
    complex(8),dimension(:,:),allocatable     :: Sigma_nambu !(2*Nso,2*Nso)
    complex(8),dimension(:,:),allocatable     :: invG0_nambu !(2*Nso,2*Nso)
    complex(8),dimension(:,:,:,:),allocatable :: calG0       !(2,Nsites,Nso,Nso)
    complex(8),dimension(:,:,:),allocatable   :: Weiss_tmp   !(Ntot,Ntot,Lfreq)
    complex(8),dimension(:,:,:),allocatable   :: Theta_tmp   !(Ntot,Ntot,Lfreq)
    !
    Nsites_=1;if(present(Nsites))Nsites_=Nsites
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
    !Testing part:
    Ntot  = size(Gloc,1)
    if(mod(Ntot,Nsites_)/=0)stop "dmft_get_weiss_normal_main: Ntot%Nsites != 0"
    Lfreq = size(Gloc,3)
    Nso   = Ntot/Nsites_
    call assert_shape(Gloc,[Ntot,Ntot,Lfreq],"dmft_get_weiss_superc_main","Gloc")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"dmft_get_weiss_superc_main","Sigma")
    call assert_shape(Weiss,[Ntot,Ntot,Lfreq],"dmft_get_weiss_superc_main","Weiss")
    call assert_shape(Floc,[Ntot,Ntot,Lfreq],"dmft_get_weiss_superc_main","Floc")
    call assert_shape(Self,[Ntot,Ntot,Lfreq],"dmft_get_weiss_superc_main","Self")
    call assert_shape(Theta,[Ntot,Ntot,Lfreq],"dmft_get_weiss_superc_main","Theta")
    !
    !Build array to extract diagonal blocks
    allocate(Green_site(2,Nso,Nso))
    allocate(Sigma_site(2,Nso,Nso))
    allocate(invG_nambu(2*Nso,2*Nso))
    allocate(Sigma_nambu(2*Nso,2*Nso))
    allocate(invG0_nambu(2*Nso,2*Nso))
    allocate(calG0(2,Nsites_,Nso,Nso))
    !
    !Work out the frequencies in parallel:
    allocate(Weiss_tmp(Ntot,Ntot,Lfreq))
    allocate(Theta_tmp(Ntot,Ntot,Lfreq))
    Weiss_tmp=zero
    Theta_tmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       !For a fixed frequency update the local Weiss field
       do ilat=1,Nsites_
          Green_site(1,:,:) = select_block(ilat,Gloc(:,:,i),Nsites_,Nso)
          Green_site(2,:,:) = select_block(ilat,Floc(:,:,i),Nsites_,Nso)
          !
          Sigma_site(1,:,:) = select_block(ilat,Sigma(:,:,i),Nsites_,Nso)
          Sigma_site(2,:,:) = select_block(ilat,Self(:,:,i),Nsites_,Nso)
          do io=1,Nso
             do jo=1,Nso
                Sigma_nambu(io,jo)          =        Sigma_site(1,io,jo)
                Sigma_nambu(io,jo+Nso)      =        Sigma_site(2,io,jo)
                Sigma_nambu(io+Nso,jo)      =        Sigma_site(2,io,jo)
                Sigma_nambu(io+Nso,jo+Nso)  =-conjg( Sigma_site(1,io,jo) )
                !
                invG_nambu(io,jo)           =        Green_site(1,io,jo)
                invG_nambu(io,jo+Nso)       =        Green_site(2,io,jo)
                invG_nambu(io+Nso,jo)       =        Green_site(2,io,jo)
                invG_nambu(io+Nso,jo+Nso)   =-conjg( Green_site(1,io,jo) )
             enddo
          enddo
          !
          call inv(invG_nambu)
          !
          ![calG0]^-1_ilat = [Gloc]_ilat^-1 + [Sigma]_ilat
          invG0_nambu = invG_nambu + Sigma_nambu
          !
          call inv(invG0_nambu)
          do io=1,Nso
             do jo=1,Nso
                calG0(1,ilat,io,jo) = invG0_nambu(io,jo)
                calG0(2,ilat,io,jo) = invG0_nambu(io,jo+Nso)
             enddo
          enddo
       enddo
       !Once all sites are obtained we dump back 
       Weiss_tmp(:,:,i) = blocks_to_matrix(calG0(1,:,:,:),Nsites_,Nso)
       Theta_tmp(:,:,i) = blocks_to_matrix(calG0(2,:,:,:),Nsites_,Nso)
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Weiss = zero
       Theta = zero
       call Mpi_AllReduce(Weiss_tmp, Weiss, size(Weiss), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
       call Mpi_AllReduce(Theta_tmp, Theta, size(Theta), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Weiss=Weiss_tmp
       Theta=Theta_tmp
    endif
#else
    Weiss=Weiss_tmp
    Theta=Theta_tmp
#endif
  end subroutine dmft_get_weiss_superc_main





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
  subroutine dmft_get_weiss_normal_rank4(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gloc      ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Sigma     !..
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss     !..
    !aux
    complex(8),dimension(:,:,:),allocatable       :: SF,GF,WF  ![Nso,Nso][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Ntot  = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank4","Gloc")
    call assert_shape(Weiss, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank4","Weiss")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq))
    !
    SF = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    GF = reshape_rank4_to_matrix(Gloc,Nspin,Norb,Lfreq)
    call dmft_get_weiss_normal_main(GF,SF,WF,Nsites=1)
    Weiss = reshape_matrix_to_rank4(WF,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_normal_rank4


  subroutine dmft_get_weiss_normal_rank5(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc      ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss     !..
    !aux
    complex(8),dimension(:,:,:),allocatable         :: SF,GF,WF  ![N,N][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lfreq = size(Gloc,6)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank5","Gloc")
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank5","Weiss")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq))
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    GF   = reshape_rank5_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    call dmft_get_weiss_normal_main(GF,SF,WF,Nsites=Nlat)
    Weiss = reshape_matrix_to_rank5(WF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_normal_rank5

  subroutine dmft_get_weiss_normal_rank6(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss !..
    !aux
    complex(8),dimension(:,:,:),allocatable           :: SF,GF,WF  ![N,N][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lfreq = size(Gloc,7)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function  [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin]"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank6","Gloc")
    call assert_shape(Weiss,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank6","Weiss")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq))
    !
    SF   = reshape_rank6_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    GF   = reshape_rank6_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    call dmft_get_weiss_normal_main(GF,SF,WF,Nsites=1)
    Weiss = reshape_matrix_to_rank6(WF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_normal_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_weiss_normal_rank7(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma !..
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Weiss !..
    !aux
    complex(8),dimension(:,:,:),allocatable             :: SF,GF,WF  ![N,N][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nineq = size(Gloc,1)
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,4)
    Norb  = size(Gloc,6)
    Lfreq = size(Gloc,8)    
    Ntot  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function  [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin,Nineq]"
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank7","Gloc")
    call assert_shape(Weiss,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank7","Weiss")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq))
    !
    SF   = reshape_rank7_to_matrix(Sigma,Nineq,Nlat,Nspin,Norb,Lfreq)
    GF   = reshape_rank7_to_matrix(Gloc,Nineq,Nlat,Nspin,Norb,Lfreq)
    call dmft_get_weiss_normal_main(GF,SF,WF,Nsites=Nineq)
    Weiss = reshape_matrix_to_rank7(WF,Nineq,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_normal_rank7
#endif


  !##################################################################
  !##################################################################
  !                         SUPERC
  !##################################################################
  !##################################################################
  subroutine dmft_get_weiss_superc_rank4(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc      ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Floc      !..
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Self      !..
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss     !..
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Theta     !..
    !aux
    complex(8),dimension(:,:,:,:),allocatable     :: SF,GF,WF  ![2][Nso,Nso][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Ntot  = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Gloc")
    call assert_shape(Weiss, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Weiss")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Self")
    call assert_shape(Floc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Floc")
    call assert_shape(Theta, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Theta")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq), WF(2,Ntot,Ntot,Lfreq))
    !
    SF(1,:,:,:) = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank4_to_matrix(Self,Nspin,Norb,Lfreq)
    GF(1,:,:,:) = reshape_rank4_to_matrix(Gloc,Nspin,Norb,Lfreq)
    GF(2,:,:,:) = reshape_rank4_to_matrix(Floc,Nspin,Norb,Lfreq)
    call dmft_get_weiss_superc_main(GF(1,:,:,:),GF(2,:,:,:),SF(1,:,:,:),SF(2,:,:,:),WF(1,:,:,:),WF(2,:,:,:),Nsites=1)
    Weiss = reshape_matrix_to_rank4(WF(1,:,:,:),Nspin,Norb,Lfreq)
    Theta = reshape_matrix_to_rank4(WF(2,:,:,:),Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_superc_rank4


  subroutine dmft_get_weiss_superc_rank5(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc      ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Floc      !..
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Self      !..
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss     !..
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Theta     !..
    !aux
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF,WF  ![2][Nso,Nso][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lfreq = size(Gloc,6)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Gloc")
    call assert_shape(Weiss, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Weiss")
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Self")
    call assert_shape(Floc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Floc")
    call assert_shape(Theta, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Theta")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq), WF(2,Ntot,Ntot,Lfreq))
    !
    SF(1,:,:,:) = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank5_to_matrix(Self,Nlat,Nspin,Norb,Lfreq)
    GF(1,:,:,:) = reshape_rank5_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    GF(2,:,:,:) = reshape_rank5_to_matrix(Floc,Nlat,Nspin,Norb,Lfreq)
    call dmft_get_weiss_superc_main(GF(1,:,:,:),GF(2,:,:,:),SF(1,:,:,:),SF(2,:,:,:),WF(1,:,:,:),WF(2,:,:,:),Nsites=Nlat)
    Weiss = reshape_matrix_to_rank5(WF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Theta = reshape_matrix_to_rank5(WF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_superc_rank5


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_weiss_superc_rank6(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc      ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Floc      !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Self      !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss     !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Theta     !..
    !aux
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF,WF  ![2][Nso,Nso][L]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lfreq = size(Gloc,7)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin]"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Gloc")
    call assert_shape(Weiss, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Weiss")
    call assert_shape(Self,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Self")
    call assert_shape(Floc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Floc")
    call assert_shape(Theta, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Theta")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq), WF(2,Ntot,Ntot,Lfreq))
    !
    SF(1,:,:,:) = reshape_rank6_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank6_to_matrix(Self ,Nlat,Nspin,Norb,Lfreq)
    GF(1,:,:,:) = reshape_rank6_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    GF(2,:,:,:) = reshape_rank6_to_matrix(Floc,Nlat,Nspin,Norb,Lfreq)
    call dmft_get_weiss_superc_main(GF(1,:,:,:),GF(2,:,:,:),SF(1,:,:,:),SF(2,:,:,:),WF(1,:,:,:),WF(2,:,:,:),Nsites=1)
    Weiss = reshape_matrix_to_rank6(WF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Theta = reshape_matrix_to_rank6(WF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_weiss_superc_rank6
#endif


end module SC_WEISS







!   subroutine dmft_get_weiss_normal_bethe(Gloc,Weiss,Hloc,Wbands,Nsites,axis) !N=Nsites*Nso
!     complex(8),dimension(:,:,:),intent(in)    :: Gloc  ! [N,Nsites*Nso][L]
!     complex(8),dimension(:,:,:),intent(inout) :: Weiss ! [N,Nsites*Nso][L]
!     complex(8),dimension(:,:),intent(in)      :: Hloc  ! [N,Nsites*Nso]
!     real(8),dimension(:),intent(in)           :: Wbands ![N]
!     integer                                   :: Nsites,Nso
!     character(len=*),optional                 :: axis
!     character(len=1)                          :: axis_
!     !aux
!     complex(8),dimension(:,:,:),allocatable   :: H0
!     complex(8),dimension(:,:),allocatable     :: G_site
!     complex(8),dimension(:,:,:),allocatable   :: Weiss_tmp ! [N,Nsites*Nso][L]
!     complex(8),dimension(:,:,:,:),allocatable :: calG0     ! [Nsites,Nso,Nso][L]
!     complex(8),dimension(size(Gloc,3))        :: invWeiss  ![L]
!     !
!     axis_='m' ; if(present(axis))axis_=axis
!     !
!     !MPI setup:
! #ifdef _MPI    
!     if(check_MPI())then
!        mpi_master= get_master_MPI()
!     else
!        mpi_master=.true.
!     endif
! #else
!     mpi_master=.true.
! #endif
!     !
!     Ntot  = size(Gloc,1);
!     if(mod(Ntot,Nsites)/=0)stop "dmft_get_weiss_normal_main: Ntot%Nsites != 0"
!     Lfreq = size(Gloc,3)
!     Nso   = Ntot/Nsites
!     call assert_shape(Gloc,[Ntot,Nsites*Nso,Lfreq],"dmft_get_weiss_normal_main","Gloc")
!     call assert_shape(Weiss,[Ntot,Nsites*Nso,Lfreq],"dmft_get_weiss_normal_main","Weiss")
!     call assert_shape(Hloc,[Ntot,Nsites*Nso],"dmft_get_weiss_normal_main","Hloc")
!     !
!     call build_frequency_array(axis_)
!     !
!     allocate(H0(Nsites,Nso,Nso))
!     H0 = matrix_to_blocks(Hloc,Nsites,Nso)
!     !
!     allocate(calG0(Nsites,Nso,Nso))
!     allocate(G_site(Nso,Nso))
!     !
!     !\calG0^{-1}_aa = iw + mu - H_0 - d**2/4*Gw
!     !
!     !Work out the frequencies in parallel:
!     allocate(Weiss_tmp(Ntot,Ntot,Lfreq))
!     Weiss_tmp=zero
!     do i=1+mpi_rank,Lfreq,mpi_size
!        !
!        !For a fixed frequency update the local Weiss field
!        do ilat=1,Nsites
!           G_site = select_block(ilat,Gloc(:,:,i),Nsites,Nso)
!           do io=1,Nso
!              invWeiss = wfreq(i) + xmu - H0(ilat,io,io) - 0.25d0*Wbands(io)**2*G_site(io,io)
!              calG0(ilat,io,io) = one/invWeiss
!           enddo
!        enddo
!        !Once all sites are obtained we dump back 
!        Weiss_tmp(:,:,i) = blocks_to_matrix(calG0,Nsites,Nso)
!     enddo
! #ifdef _MPI    
!     if(check_MPI())then
!        Weiss = zero
!        call Mpi_AllReduce(Weiss_tmp, Weiss, size(Weiss), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
!     else
!        Weiss=Weiss_tmp
!     endif
! #else
!     Weiss=Weiss_tmp
! #endif
!   end subroutine dmft_get_weiss_normal_bethe
