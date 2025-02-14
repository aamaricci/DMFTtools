module SC_DELTA
  USE DMFT_CTRL_VARS
  USE SC_COMMON
  implicit none
  private


  public :: dmft_get_delta_normal_rank2
  public :: dmft_get_delta_normal_rank3
  public :: dmft_get_delta_normal_rank4
  public :: dmft_get_delta_normal_rank5
  public :: dmft_get_delta_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: dmft_get_delta_normal_rank7
#endif

  public :: dmft_get_delta_superc_rank2
  public :: dmft_get_delta_superc_rank3
  public :: dmft_get_delta_superc_rank4
  public :: dmft_get_delta_superc_rank5
  public :: dmft_get_delta_superc_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8  
  public :: dmft_get_delta_superc_rank7
#endif



  complex(8),dimension(:,:),allocatable   :: Wtmp
  complex(8),dimension(:,:,:),allocatable :: Ttmp


contains



  !##################################################################
  !##################################################################
  !                       NORMAL PHASE
  ! Get the local Delta Field calG0 using DMFT self-consistency 
  !  given G_loc/F_loc and Sigma/Self in different shapes:
  !
  ! . rank-2 DMFT    [N,N,:] check @1Gb memory per Function
  ! . rank-3 DMFT    [Nlat,Nso,Nso,:]
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:]
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-5_6 R-DMFT[Nlat,Nspin,Nspin,Norb,Norb,:]
  !                ->[Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-6 CDMFT   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:]
  !                   Nlat=N_cluster
  ! . rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Hk Matrix
  !
  ! For rank-2 we adopt a somehow general strategy of having N = Nlat*Nso,
  ! with Nlat=1 (default) or Nlat=Nsites user input.
  ! This is a general case when there are inequivalent sites in the unit cell.
  ! >> WARNING: The use of rank-2 for large supercells is highly discouraged.
  !             A quick check shows that for Nlat large enough the memory gets
  !             completely filled by allocating functions F(N,N,L).
  !             e.g. L=1000, N=Nlat*Nso, Nso=4, Nlat=20x20
  !                  N = 400*4 = 1600.
  !                  size(F) = 1600x1600x1000*complex(8) 
  !                          = 2.56x1d9 *(2*8byte) =  40.96Gb
  !
  !##################################################################
  ! G/Sigma Shape: [N,N][:] N = Nlat*Nso
  !##################################################################
  subroutine dmft_get_delta_normal_rank2(Gloc,Sigma,Delta,Hloc,Nsites,axis) !N=Nsites*Nso
    complex(8),dimension(:,:,:),intent(inout) :: Gloc
    complex(8),dimension(:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:),intent(in)      :: Hloc
    integer,optional                          :: Nsites
    character(len=*),optional                 :: axis
    character(len=1)                          :: axis_
    complex(8),dimension(:,:,:),allocatable   :: Delta_site
    !
    Nlat=1;if(present(Nsites))Nlat=Nsites
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !Testing part:
    Ntot  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    if(Ntot*Ntot*Lfreq*16*1d-9 > 1d0)print*, "WARNING dmft_get_delta_normal_rank2: Memory required exceeds 1Gb per Function!"
    Nso   = Ntot/Nlat ; if(mod(Ntot,Nlat)/=0)stop "dmft_get_delta_normal_rank2: Ntot%Nsites != 0"
    call assert_shape(Gloc, [Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_normal_rank2","Gloc")
    call assert_shape(Sigma,[Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_normal_rank2","Sigma")
    call assert_shape(Delta,[Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_normal_rank2","Delta")
    call assert_shape(Hloc, [Ntot,Nlat*Nso],"dmft_get_delta_normal_rank2","Hloc")
    !
    call build_frequency_array(axis_)
    !
    allocate(Delta_site(Nlat,Nso,Nso))
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       do ilat=1,Nlat
          call dmft_delta_normal((wfreq(i)+xmu)*eye(Nso)-get_block(ilat,Hloc),&
               get_block(ilat,Gloc(:,:,i)),   &
               get_block(ilat,Sigma(:,:,i)),  &
               Delta_site(ilat,:,:))
       enddo
       Delta(:,:,i) = from_rank3(Delta_site)
    enddo MPIloop
    call gf_reduce(Delta)
    deallocate(Delta_Site)
    !
  end subroutine dmft_get_delta_normal_rank2

  !G/Sigma shape: [Nlat,Nso,Nso][:]
  !##################################################################
  subroutine dmft_get_delta_normal_rank3(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:),intent(in)      :: Hloc
    character(len=*),optional                   :: axis
    character(len=1)                            :: axis_
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nso   = size(Gloc,2)
    Lfreq = size(Gloc,4)    
    Ntot  = Nlat*Nso
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nso,Nso]:"
    call assert_shape(Sigma,[Nlat,Nso,Nso,Lfreq],"dmft_get_delta_normal_rank3","Sigma")
    call assert_shape(Gloc, [Nlat,Nso,Nso,Lfreq],"dmft_get_delta_normal_rank3","Gloc")
    call assert_shape(Delta,[Nlat,Nso,Nso,Lfreq],"dmft_get_delta_normal_rank3","Delta")
    call assert_shape(Hloc,[Nlat,Nso,Nso],"dmft_get_delta_normal_rank3","Hloc")
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_normal((wfreq(i)+xmu)*eye(Ntot)-from_rank3(Hloc), &
            from_rank3(Gloc(:,:,:,i)), &
            from_rank3(Sigma(:,:,:,i)), Wtmp ) 
       Delta(:,:,:,i) = to_rank3(Wtmp)
    enddo MPIloop
    call gf_reduce(Delta)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_delta_normal_rank3

  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_normal_rank4(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                     :: axis
    character(len=1)                              :: axis_
    !
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Nso   = Nspin*Norb
    Ntot  = Nso
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank4","Gloc")
    call assert_shape(Delta,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank4","Delta")
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_rank4","Hloc")
    !
    Delta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_normal((wfreq(i)+xmu)*eye(Ntot)-from_rank4(Hloc), &
            from_rank4(Gloc(:,:,:,:,i)), &
            from_rank4(Sigma(:,:,:,:,i)), Wtmp)
       Delta(:,:,:,:,i) = to_rank4(Wtmp)
    enddo MPIloop
    call gf_reduce(Delta)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_delta_normal_rank4
  
  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_normal_rank5(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                       :: axis
    character(len=1)                                :: axis_
    !
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lfreq = size(Gloc,6)
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank5","Gloc")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank5","Delta")
    call assert_shape(Hloc, [Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_rank5","Hloc")
    !
    Delta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_normal((wfreq(i)+xmu)*eye(Ntot)-from_rank5(Hloc), &
            from_rank5(Gloc(:,:,:,:,:,i)), &
            from_rank5(Sigma(:,:,:,:,:,i)), Wtmp)
       Delta(:,:,:,:,:,i)= to_rank5(Wtmp)
    enddo MPIloop
    call gf_reduce(Delta)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_delta_normal_rank5

  !G/Sigma shape: [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_normal_rank6(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                         :: axis
    character(len=1)                                  :: axis_
    !
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lfreq = size(Gloc,7)
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function  [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank6","Gloc")
    call assert_shape(Delta,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank6","Delta")
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_rank6","Hloc")
    !
    Delta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_normal((wfreq(i)+xmu)*eye(Ntot)-from_rank6(Hloc), &
            from_rank6(Gloc(:,:,:,:,:,:,i)), &
            from_rank6(Sigma(:,:,:,:,:,:,i)), Wtmp)
       Delta(:,:,:,:,:,:,i) = to_rank6(Wtmp)
    enddo MPIloop
    call gf_reduce(Delta)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_delta_normal_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  !G/Sigma shape: [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_normal_rank7(Gloc,Sigma,Delta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                           :: axis
    character(len=1)                                    :: axis_
    !
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nineq = size(Gloc,1)
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,4)
    Norb  = size(Gloc,6)
    Lfreq = size(Gloc,8)    
    Ntot  = Nineq*Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function  [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank7","Gloc")
    call assert_shape(Delta,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_normal_rank7","Delta")
    call assert_shape(Hloc,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_rank7","Hloc")
    !
    Delta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_normal((wfreq(i)+xmu)*eye(Ntot)-from_rank7(Hloc), &
            from_rank7(Gloc(:,:,:,:,:,:,:,i)), &
            from_rank7(Sigma(:,:,:,:,:,:,:,i)), Wtmp)
       Delta(:,:,:,:,:,:,:,i)= to_rank7(Wtmp)
    enddo MPIloop
    call gf_reduce(Delta)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_delta_normal_rank7
#endif




  !##################################################################
  !##################################################################
  !                         SUPERC
  ! . rank-2 DMFT    [][N,N,:] check @1Gb memory per Function
  ! . rank-3 DMFT    [Nlat,Nso,Nso,:] 
  ! . rank-4 DMFT    [Nspin,Nspin,Norb,Norb,:] 
  ! . rank-5 R-DMFT  [Nlat,Nspin,Nspin,Norb,Norb,:]
  ! . rank-6 CDMFT   [Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  ! . rank-7 R-CDMFT [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,:] Nlat=N_cluster
  !##################################################################
  ! Shape: [N,N][:]
  !##################################################################
  subroutine dmft_get_delta_superc_rank2(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,Nsites,axis) 
    complex(8),dimension(:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:),intent(in)    :: Sigma 
    complex(8),dimension(:,:,:),intent(in)    :: Self  
    complex(8),dimension(:,:,:),intent(inout) :: Delta 
    complex(8),dimension(:,:,:),intent(inout) :: Theta
    complex(8),dimension(:,:),intent(in)      :: Hloc
    integer,optional                          :: Nsites
    character(len=*),optional                 :: axis
    character(len=1)                          :: axis_
    complex(8),dimension(:,:,:,:),allocatable :: calG0
    !
    Nlat=1;if(present(Nsites))Nlat=Nsites
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !Testing part:
    Ntot  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    if(Ntot*Ntot*Lfreq*16*1d-9 > 1d0)print*, "WARNINGL: dmft_get_delta_superc_rank2: Memory required exceeds 1Gb per Function!"
    Nso   = Ntot/Nlat
    if(mod(Ntot,Nlat)/=0)stop "dmft_get_delta_superc_rank2: Ntot%Nsites != 0"
    call assert_shape(Gloc, [Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_superc_rank2","Gloc")
    call assert_shape(Floc, [Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_superc_rank2","Floc")
    call assert_shape(Sigma,[Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_superc_rank2","Sigma")
    call assert_shape(Self, [Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_superc_rank2","Self")
    call assert_shape(Delta,[Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_superc_rank2","Delta")
    call assert_shape(Theta,[Ntot,Nlat*Nso,Lfreq],"dmft_get_delta_superc_rank2","Theta")
    call assert_shape(Hloc, [Ntot,Nlat*Nso],"dmft_get_delta_superc_rank2","Hloc")
    !
    Delta = zero
    Theta = zero
    !
    call build_frequency_array(axis_)
    !
    allocate(calG0(2,Nlat,Nso,Nso))
    calG0 = zero
    !
    MPIloop: do i=1+mpi_rank,Lfreq,mpi_size
       do ilat=1,Nlat        
          call dmft_delta_superc((wfreq(i)+xmu)*eye(Ntot)-get_block(ilat,Hloc), &
               get_block(ilat,Gloc(:,:,i)) , get_block(ilat,Floc(:,:,i)),&
               get_block(ilat,Sigma(:,:,i)), get_block(ilat,Self(:,:,i)),&
               calG0(1,ilat,:,:),calG0(2,ilat,:,:) ) 
       enddo
       Delta(:,:,i) =  from_rank3(calG0(1,:,:,:))
       Theta(:,:,i) =  from_rank3(calG0(2,:,:,:))
    enddo MPIloop
    call gf_reduce(Delta)
    call gf_reduce(Theta)
    !
    deallocate(calG0)
    !
  end subroutine dmft_get_delta_superc_rank2


  !G/Sigma shape: [Nlat,Nso,Nso][:]
  !##################################################################
  subroutine dmft_get_delta_superc_rank3(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis) 
    complex(8),dimension(:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma 
    complex(8),dimension(:,:,:,:),intent(in)    :: Self  
    complex(8),dimension(:,:,:,:),intent(inout) :: Delta 
    complex(8),dimension(:,:,:,:),intent(inout) :: Theta
    complex(8),dimension(:,:,:),intent(in)      :: Hloc
    character(len=*),optional                   :: axis
    character(len=1)                            :: axis_
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nso   = size(Gloc,2)
    Lfreq = size(Gloc,4)    
    Ntot  = Nlat*Nso
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nso,Nso]:"
    call assert_shape(Sigma,[Nlat,Nso,Nso,Lfreq],"dmft_get_delta_superc_rank3","Sigma")
    call assert_shape(Gloc, [Nlat,Nso,Nso,Lfreq],"dmft_get_delta_superc_rank3","Gloc")
    call assert_shape(Delta,[Nlat,Nso,Nso,Lfreq],"dmft_get_delta_superc_rank3","Delta")
    call assert_shape(Self, [Nlat,Nso,Nso,Lfreq],"dmft_get_delta_superc_rank3","Self")
    call assert_shape(Floc, [Nlat,Nso,Nso,Lfreq],"dmft_get_delta_superc_rank3","Floc")
    call assert_shape(Theta,[Nlat,Nso,Nso,Lfreq],"dmft_get_delta_superc_rank3","Theta")
    call assert_shape(Hloc, [Nlat,Nso,Nso],"dmft_get_delta_superc_rank3","Hloc")
    !
    Delta = zero
    Theta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop: do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_superc((wfreq(i)+xmu)*eye(Ntot)-from_rank3(Hloc), &
            from_rank3(Gloc(:,:,:,i)) , from_rank3(Floc(:,:,:,i)), &
            from_rank3(Sigma(:,:,:,i)), from_rank3(Self(:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Delta(:,:,:,i) = to_rank3(Ttmp(1,:,:))
       Theta(:,:,:,i) = to_rank3(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Delta)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_delta_superc_rank3


  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_superc_rank4(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Theta
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                     :: axis
    character(len=1)                              :: axis_
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Nso   = Nspin*Norb
    Ntot  = Nso
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Gloc")
    call assert_shape(Delta,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Delta")
    call assert_shape(Self, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Self")
    call assert_shape(Floc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Floc")
    call assert_shape(Theta,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank4","Theta")
    call assert_shape(Hloc, [Nspin,Nspin,Norb,Norb],"dmft_get_delta_superc_rank4","Hloc")
    !
    Delta = zero
    Theta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_superc((wfreq(i)+xmu)*eye(Ntot)-from_rank4(Hloc), &
            from_rank4(Gloc(:,:,:,:,i)) , from_rank4(Floc(:,:,:,:,i)), &
            from_rank4(Sigma(:,:,:,:,i)), from_rank4(Self(:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Delta(:,:,:,:,i) = to_rank4(Ttmp(1,:,:))
       Theta(:,:,:,:,i) = to_rank4(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Delta)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_delta_superc_rank4


  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_superc_rank5(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Theta
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                       :: axis
    character(len=1)                                :: axis_
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lfreq = size(Gloc,6)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Gloc")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Delta")
    call assert_shape(Self, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Self")
    call assert_shape(Floc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Floc")
    call assert_shape(Theta,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank5","Theta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_superc_rank5","Hloc")
    !
    Delta = zero
    Theta = zero
    !
    call build_frequency_array(axis_)
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_superc((wfreq(i)+xmu)*eye(Ntot)-from_rank5(Hloc), &
            from_rank5(Gloc(:,:,:,:,:,i)) , from_rank5(Floc(:,:,:,:,:,i)), &
            from_rank5(Sigma(:,:,:,:,:,i)), from_rank5(Self(:,:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Delta(:,:,:,:,:,i) = to_rank5(Ttmp(1,:,:))
       Theta(:,:,:,:,:,i) = to_rank5(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Delta)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_delta_superc_rank5

  !G/Sigma shape: [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_superc_rank6(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Theta
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                         :: axis
    character(len=1)                                  :: axis_
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,3)
    Norb  = size(Gloc,5)
    Lfreq = size(Gloc,7)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Gloc")
    call assert_shape(Delta,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Delta")
    call assert_shape(Self, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Self")
    call assert_shape(Floc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Floc")
    call assert_shape(Theta,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank6","Theta")
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_superc_rank6","Hloc")
    !
    Delta = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_superc((wfreq(i)+xmu)*eye(Ntot)-from_rank6(Hloc), &
            from_rank6(Gloc(:,:,:,:,:,:,i)) , from_rank6(Floc(:,:,:,:,:,:,i)), &
            from_rank6(Sigma(:,:,:,:,:,:,i)), from_rank6(Self(:,:,:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Delta(:,:,:,:,:,:,i) = to_rank6(Ttmp(1,:,:))
       Theta(:,:,:,:,:,:,i) = to_rank6(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Delta)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_delta_superc_rank6





#if __GFORTRAN__ &&  __GNUC__ > 8
  !G/Sigma shape: [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_delta_superc_rank7(Gloc,Floc,Sigma,Self,Delta,Theta,Hloc,axis)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Delta
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Theta
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)      :: Hloc
    character(len=*),optional                         :: axis
    character(len=1)                                  :: axis_
    !
    axis_='m' ; if(present(axis))axis_=axis
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nineq = size(Gloc,1)
    Nlat  = size(Gloc,2)
    Nspin = size(Gloc,4)
    Norb  = size(Gloc,6)
    Lfreq = size(Gloc,8)    
    Ntot  = Nlat*Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Delta function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank7","Gloc")
    call assert_shape(Delta,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank7","Delta")
    call assert_shape(Self, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank7","Self")
    call assert_shape(Floc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank7","Floc")
    call assert_shape(Theta,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_delta_superc_rank7","Theta")
    call assert_shape(Hloc,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_superc_rank7","Hloc")
    !
    Delta = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_delta_superc((wfreq(i)+xmu)*eye(Ntot)-from_rank7(Hloc), &
            from_rank7(Gloc(:,:,:,:,:,:,:,i)) , from_rank7(Floc(:,:,:,:,:,:,:,i)), &
            from_rank7(Sigma(:,:,:,:,:,:,:,i)), from_rank7(Self(:,:,:,:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Delta(:,:,:,:,:,:,:,i) = to_rank7(Ttmp(1,:,:))
       Theta(:,:,:,:,:,:,:,i) = to_rank7(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Delta)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_delta_superc_rank7
#endif





end module SC_DELTA
