module SC_DELTA
  USE DMFT_CTRL_VARS
  USE SC_COMMON
  implicit none
  private


  public :: dmft_get_delta_normal_main
  public :: dmft_get_delta_normal_rank4
  public :: dmft_get_delta_normal_rank5
  public :: dmft_get_delta_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: dmft_get_delta_normal_rank7
#endif
  
  public :: dmft_get_delta_superc_main
  public :: dmft_get_delta_superc_rank4
  public :: dmft_get_delta_superc_rank5
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: dmft_get_delta_superc_rank6
#endif
  
contains




  subroutine dmft_get_delta_normal_main(Gloc,Sigma,Delta,Hloc,Nsites,axis) !N=Nsites*Nso
    complex(8),dimension(:,:,:),intent(in)    :: Gloc  ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:,:),intent(in)    :: Sigma ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:,:),intent(inout) :: Delta ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:),intent(in)      :: Hloc  ! [N,Nsites*Nso]
    integer,optional                          :: Nsites
    integer                                   :: Nsites_,Nso
    character(len=*),optional                 :: axis
    character(len=1)                          :: axis_
    !aux
    complex(8),dimension(:,:,:),allocatable   :: Delta_tmp    ![N,N][L]
    complex(8),dimension(:,:,:),allocatable   :: Delta_site   ![Nsites,Nso,Nso]
    complex(8),dimension(:,:),allocatable     :: Sigma_site   ![Nso,Nso]
    complex(8),dimension(:,:),allocatable     :: invG_site
    complex(8),dimension(:,:),allocatable     :: Hloc_site
    !
    Nsites_=1;if(present(Nsites))Nsites_=Nsites
    axis_='m' ; if(present(axis))axis_=axis
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
    if(mod(Ntot,Nsites_)/=0)stop "dmft_get_delta_normal_main: Ntot%Nsites != 0"
    Lfreq = size(Gloc,3)
    Nso   = Ntot/Nsites_
    call assert_shape(Gloc,[Ntot,Nsites_*Nso,Lfreq],"dmft_get_delta_normal_main","Gloc")
    call assert_shape(Sigma,[Ntot,Nsites_*Nso,Lfreq],"dmft_get_delta_normal_main","Sigma")
    call assert_shape(Delta,[Ntot,Nsites_*Nso,Lfreq],"dmft_get_delta_normal_main","Delta")
    call assert_shape(Hloc,[Ntot,Nsites_*Nso],"dmft_get_delta_normal_main","Hloc")
    !
    call build_frequency_array(axis_)
    !
    !Build array to extract diagonal blocks
    allocate(Delta_tmp(Ntot,Ntot,Lfreq))
    allocate(Delta_site(Nsites_,Nso,Nso))
    allocate(Sigma_site(Nso,Nso))
    allocate(invG_site(Nso,Nso))
    allocate(Hloc_site(Nso,Nso))
    !
    !Work out the frequencies in parallel:
    Delta_tmp=zero
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       !
       !For a fixed frequency update the local Delta field
       do ilat=1,Nsites_
          Sigma_site = select_block(ilat,Sigma(:,:,i),Nsites_,Nso)
          Hloc_site  = select_block(ilat,Hloc(:,:),Nsites_,Nso)
          invG_site  = select_block(ilat,Gloc(:,:,i),Nsites_,Nso)
          call inv(invG_site)
          ! [Delta]_ilat = [z-Hloc-Sigma]_ilat -  [Gloc]^-1_ilat
          Delta_site(ilat,:,:)= (wfreq(i)+xmu)*eye(Nso) - Hloc_site - Sigma_site - invG_site
       enddo
       Delta_tmp(:,:,i) = blocks_to_matrix(Delta_site,Nsites_,Nso)
    enddo MPIloop
#ifdef _MPI    
    if(check_MPI())then
       Delta = zero
       call Mpi_AllReduce(Delta_tmp, Delta, size(Delta), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Delta=Delta_tmp
    endif
#else
    Delta=Delta_tmp
#endif
    !
  end subroutine dmft_get_delta_normal_main





  subroutine dmft_get_delta_superc_main(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,Nsites,axis) !N=Nsites*Nso
    complex(8),dimension(:,:,:),intent(in)    :: Gloc  ! [N,Nsites*Nso][L]
    complex(8),dimension(:,:,:),intent(in)    :: Floc  ! ..
    complex(8),dimension(:,:,:),intent(in)    :: Sigma ! ..
    complex(8),dimension(:,:,:),intent(in)    :: Self  ! ..
    complex(8),dimension(:,:,:),intent(inout) :: Delta ! ..
    complex(8),dimension(:,:,:),intent(inout) :: Theta ! ..
    complex(8),dimension(:,:),intent(in)      :: Hloc  ! [N,Nsites*Nso]
    integer,optional                          :: Nsites
    integer                                   :: Nsites_,Nso
    character(len=*),optional                 :: axis
    character(len=1)                          :: axis_
    !aux
    complex(8),dimension(:,:,:),allocatable   :: Green_site  !(2,Nso,Nso)
    complex(8),dimension(:,:,:),allocatable   :: Sigma_site  !(2,Nso,Nso)
    complex(8),dimension(:,:,:),allocatable   :: zeta_site   !(2,Nso,Nso)
    complex(8),dimension(:,:),allocatable     :: Hloc_site   !(Nso,Nso)
    complex(8),dimension(:,:),allocatable     :: invG_nambu  !(2*Nso,2*Nso)
    complex(8),dimension(:,:,:,:),allocatable :: calG0       !(2,Nsites,Nso,Nso)
    complex(8),dimension(:,:,:),allocatable   :: Delta_tmp   !(Ntot,Ntot,Lfreq)
    complex(8),dimension(:,:,:),allocatable   :: Theta_tmp   !(Ntot,Ntot,Lfreq)
    !
    Nsites_=1;if(present(Nsites))Nsites_=Nsites
    axis_='m' ; if(present(axis))axis_=axis
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
    Ntot  = size(Gloc,2)
    if(mod(Ntot,Nsites_)/=0)stop "dmft_get_delta_normal_main: Ntot%Nsites != 0"
    Lfreq = size(Gloc,3)
    Nso   = Ntot/Nsites_
    call assert_shape(Gloc,[Ntot,Ntot,Lfreq],"dmft_get_delta_superc_main","Gloc")
    call assert_shape(Sigma,[Ntot,Ntot,Lfreq],"dmft_get_delta_superc_main","Sigma")
    call assert_shape(Delta,[Ntot,Ntot,Lfreq],"dmft_get_delta_superc_main","Delta")
    call assert_shape(Floc,[Ntot,Ntot,Lfreq],"dmft_get_delta_superc_main","Floc")
    call assert_shape(Self,[Ntot,Ntot,Lfreq],"dmft_get_delta_superc_main","Self")
    call assert_shape(Theta,[Ntot,Ntot,Lfreq],"dmft_get_delta_superc_main","Theta")
    call assert_shape(Hloc,[Ntot,Ntot],"dmft_get_delta_superc_main","Hloc")
    !
    call build_frequency_array(axis_)
    !
    !Build array to extract diagonal blocks
    allocate(Green_site(2,Nso,Nso))
    allocate(Sigma_site(2,Nso,Nso))
    allocate(zeta_site(2,Nso,Nso))
    allocate(Hloc_site(Nso,Nso))
    allocate(invG_nambu(2*Nso,2*Nso))
    allocate(calG0(2,Nsites_,Nso,Nso))
    !
    !Work out the frequencies in parallel:
    allocate(Delta_tmp(Ntot,Ntot,Lfreq))
    allocate(Theta_tmp(Ntot,Ntot,Lfreq))
    Delta_tmp=zero
    Theta_tmp=zero
    do i=1+mpi_rank,Lfreq,mpi_size
       !For a fixed frequency update the local Delta field
       do ilat=1,Nsites_
          Hloc_site = select_block(ilat,Hloc,Nsites_,Nso)
          !
          Green_site(1,:,:) = select_block(ilat,Gloc(:,:,i),Nsites_,Nso)
          Green_site(2,:,:) = select_block(ilat,Floc(:,:,i),Nsites_,Nso)
          !
          Sigma_site(1,:,:) = select_block(ilat,Sigma(:,:,i),Nsites_,Nso)
          Sigma_site(2,:,:) = select_block(ilat,Self(:,:,i),Nsites_,Nso)
          !
          zeta_site(1,:,:)  = (wfreq(i)+xmu)*eye(Nso) - Hloc_site - Sigma_site(1,:,:)
          zeta_site(2,:,:)  =                                     - Sigma_site(2,:,:)
          do io=1,Nso
             do jo=1,Nso
                invG_nambu(io,jo)           =        Green_site(1,io,jo)
                invG_nambu(io,jo+Nso)       =        Green_site(2,io,jo)
                invG_nambu(io+Nso,jo)       =        Green_site(2,io,jo)
                invG_nambu(io+Nso,jo+Nso)   =-conjg( Green_site(1,io,jo) )
             enddo
          enddo
          !
          call inv(invG_nambu)
          !
          ! [Delta]_ilat = [z-Hloc-Sigma]_ilat -  [Gloc]^-1_ilat
          do io=1,Nso
             do jo=1,Nso
                calG0(1,ilat,io,jo) = zeta_site(1,io,jo) - invG_nambu(io,jo)
                calG0(2,ilat,io,jo) = zeta_site(2,io,jo) - invG_nambu(io,jo+Nso)
             enddo
          enddo
       enddo
       !Once all sites are obtained we dump back 
       Delta_tmp(:,:,i) = blocks_to_matrix(calG0(1,:,:,:),Nsites_,Nso)
       Theta_tmp(:,:,i) = blocks_to_matrix(calG0(2,:,:,:),Nsites_,Nso)
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Delta = zero
       Theta = zero
       call Mpi_AllReduce(Delta_tmp, Delta, size(Delta), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
       call Mpi_AllReduce(Theta_tmp, Theta, size(Theta), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Delta=Delta_tmp
       Theta=Theta_tmp
    endif
#else
    Delta=Delta_tmp
    Theta=Theta_tmp
#endif
  end subroutine dmft_get_delta_superc_main







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
  subroutine dmft_get_delta_normal_rank4(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc      ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta     !..
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc      ![Nspin,Nspin,Norb,Norb]
    character(len=*),optional                     :: axis
    character(len=1)                          :: axis_
    !aux
    complex(8),dimension(:,:,:),allocatable       :: SF,GF,WF  ![Nso,Nso][L]
    complex(8),dimension(:,:),allocatable         :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Ntot  = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank4","Gloc")
    call assert_shape(Delta, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank4","Delta")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    GF = reshape_rank4_to_matrix(Gloc,Nspin,Norb,Lfreq)
    HF = reshape_rank4_to_matrix(Hloc,Nspin,Norb)
    call dmft_get_delta_normal_main(GF,SF,WF,HF,Nsites=1,axis=axis_)
    Delta = reshape_matrix_to_rank4(WF,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_normal_rank4


  subroutine dmft_get_delta_normal_rank5(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc      ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta     !..
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc      ![Nlat,Nspin,Nspin,Norb,Norb]
    character(len=*),optional                       :: axis
    character(len=1)                                :: axis_
    !aux
    complex(8),dimension(:,:,:),allocatable         :: SF,GF,WF  ![N,N][L]
    complex(8),dimension(:,:),allocatable           :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
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
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank5","Gloc")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank5","Delta")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF   = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    GF   = reshape_rank5_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    HF   = reshape_rank5_to_matrix(Hloc,Nlat,Nspin,Norb)
    call dmft_get_delta_normal_main(GF,SF,WF,HF,Nsites=Nlat,axis=axis_)
    Delta = reshape_matrix_to_rank5(WF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_normal_rank5

  subroutine dmft_get_delta_normal_rank6(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta !..
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Hloc  ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    character(len=*),optional                         :: axis
    character(len=1)                                  :: axis_
    !aux
    complex(8),dimension(:,:,:),allocatable           :: SF,GF,WF  ![N,N][L]
    complex(8),dimension(:,:),allocatable             :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
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
    if(mpi_master)write(*,"(A)")"Get Delta function  [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin]"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank6","Gloc")
    call assert_shape(Delta,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank6","Delta")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF   = reshape_rank6_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    GF   = reshape_rank6_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    HF   = reshape_rank6_to_matrix(Hloc,Nlat,Nspin,Norb)
    call dmft_get_delta_normal_main(GF,SF,WF,HF,Nsites=1,axis=axis_)
    Delta = reshape_matrix_to_rank6(WF,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_normal_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_delta_normal_rank7(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma !..
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Delta !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Hloc  ![Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    character(len=*),optional                           :: axis
    character(len=1)                                    :: axis_
    !aux
    complex(8),dimension(:,:,:),allocatable             :: SF,GF,WF  ![N,N][L]
    complex(8),dimension(:,:),allocatable               :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
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
    if(mpi_master)write(*,"(A)")"Get Delta function  [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin,Nineq]"
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank7","Gloc")
    call assert_shape(Delta,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank7","Delta")
    !
    allocate(SF(Ntot,Ntot,Lfreq), GF(Ntot,Ntot,Lfreq), WF(Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF   = reshape_rank7_to_matrix(Sigma,Nineq,Nlat,Nspin,Norb,Lfreq)
    GF   = reshape_rank7_to_matrix(Gloc,Nineq,Nlat,Nspin,Norb,Lfreq)
    HF   = reshape_rank7_to_matrix(Hloc,Nineq,Nlat,Nspin,Norb)
    call dmft_get_delta_normal_main(GF,SF,WF,HF,Nsites=Nineq,axis=axis_)
    Delta = reshape_matrix_to_rank7(WF,Nineq,Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_normal_rank7
#endif


  !##################################################################
  !##################################################################
  !                         SUPERC
  !##################################################################
  !##################################################################
  subroutine dmft_get_delta_superc_rank4(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc      ![Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Floc      !..
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Self      !..
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta     !..
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Theta     !..
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc      ![Nspin,Nspin,Norb,Norb]
    character(len=*),optional                     :: axis
    character(len=1)                              :: axis_
    !aux
    complex(8),dimension(:,:,:,:),allocatable     :: SF,GF,WF  ![2][Nso,Nso][L]
    complex(8),dimension(:,:),allocatable         :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Ntot  = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin]"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Gloc")
    call assert_shape(Delta, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Delta")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Self")
    call assert_shape(Floc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Floc")
    call assert_shape(Theta, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Theta")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq), WF(2,Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF(1,:,:,:) = reshape_rank4_to_matrix(Sigma,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank4_to_matrix(Self,Nspin,Norb,Lfreq)
    GF(1,:,:,:) = reshape_rank4_to_matrix(Gloc,Nspin,Norb,Lfreq)
    GF(2,:,:,:) = reshape_rank4_to_matrix(Floc,Nspin,Norb,Lfreq)
    HF          = reshape_rank4_to_matrix(Hloc,Nspin,Norb)
    call dmft_get_delta_superc_main(GF(1,:,:,:),GF(2,:,:,:),SF(1,:,:,:),SF(2,:,:,:),WF(1,:,:,:),WF(2,:,:,:),&
         HF,Nsites=1,axis=axis_)
    Delta = reshape_matrix_to_rank4(WF(1,:,:,:),Nspin,Norb,Lfreq)
    Theta = reshape_matrix_to_rank4(WF(2,:,:,:),Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_superc_rank4


  subroutine dmft_get_delta_superc_rank5(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc      ![Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Floc      !..
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Self      !..
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta     !..
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Theta     !..
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc     ![Nlat,Nspin,Nspin,Norb,Norb]
    character(len=*),optional                         :: axis
    character(len=1)                                  :: axis_
    !aux
    complex(8),dimension(:,:,:,:),allocatable         :: SF,GF,WF  ![2][Nso,Nso][L]
    complex(8),dimension(:,:),allocatable             :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
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
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nspin,Nlat]"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Gloc")
    call assert_shape(Delta, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Delta")
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Self")
    call assert_shape(Floc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Floc")
    call assert_shape(Theta, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Theta")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq), WF(2,Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF(1,:,:,:) = reshape_rank5_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank5_to_matrix(Self,Nlat,Nspin,Norb,Lfreq)
    GF(1,:,:,:) = reshape_rank5_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    GF(2,:,:,:) = reshape_rank5_to_matrix(Floc,Nlat,Nspin,Norb,Lfreq)
    HF          = reshape_rank5_to_matrix(Hloc,Nlat,Nspin,Norb)
    call dmft_get_delta_superc_main(GF(1,:,:,:),GF(2,:,:,:),SF(1,:,:,:),SF(2,:,:,:),WF(1,:,:,:),WF(2,:,:,:),&
         HF,Nsites=Nlat,axis=axis_)
    Delta = reshape_matrix_to_rank5(WF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Theta = reshape_matrix_to_rank5(WF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_superc_rank5


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_get_delta_superc_rank6(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc      ![Nlat,Nlat,Nspin,Nspin,Norb,Norb][L]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Floc      !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma     !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Self      !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta     !..
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Theta     !..
    complex(8),dimension(:,:,:,:,:,:),intent(in)        :: Hloc       ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    character(len=*),optional                           :: axis
    character(len=1)                                    :: axis_
    !aux
    complex(8),dimension(:,:,:,:),allocatable           :: SF,GF,WF  ![2][Nso,Nso][L]
    complex(8),dimension(:,:),allocatable               :: HF        ![Nso,Nso
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
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
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    if(mpi_master)write(*,"(A)")"The order is (int-->ext):[Norb,Nlat,Nspin]"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Gloc")
    call assert_shape(Delta, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Delta")
    call assert_shape(Self,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Self")
    call assert_shape(Floc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Floc")
    call assert_shape(Theta, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Theta")
    !
    allocate(SF(2,Ntot,Ntot,Lfreq), GF(2,Ntot,Ntot,Lfreq), WF(2,Ntot,Ntot,Lfreq), HF(Ntot,Ntot))
    !
    SF(1,:,:,:) = reshape_rank6_to_matrix(Sigma,Nlat,Nspin,Norb,Lfreq)
    SF(2,:,:,:) = reshape_rank6_to_matrix(Self ,Nlat,Nspin,Norb,Lfreq)
    GF(1,:,:,:) = reshape_rank6_to_matrix(Gloc,Nlat,Nspin,Norb,Lfreq)
    GF(2,:,:,:) = reshape_rank6_to_matrix(Floc,Nlat,Nspin,Norb,Lfreq)
    HF          = reshape_rank6_to_matrix(Hloc,Nlat,Nspin,Norb)
    call dmft_get_delta_superc_main(GF(1,:,:,:),GF(2,:,:,:),SF(1,:,:,:),SF(2,:,:,:),WF(1,:,:,:),WF(2,:,:,:),&         
         HF,Nsites=1,axis=axis_)
    Delta = reshape_matrix_to_rank6(WF(1,:,:,:),Nlat,Nspin,Norb,Lfreq)
    Theta = reshape_matrix_to_rank6(WF(2,:,:,:),Nlat,Nspin,Norb,Lfreq)
    deallocate(SF,GF,WF)
  end subroutine dmft_get_delta_superc_rank6
#endif


end module SC_DELTA
