module SC_WEISS
  USE DMFT_CTRL_VARS
  USE SC_COMMON
  implicit none
  private


  public :: dmft_get_weiss_normal_rank2
  public :: dmft_get_weiss_normal_rank3
  public :: dmft_get_weiss_normal_rank4
  public :: dmft_get_weiss_normal_rank5
  public :: dmft_get_weiss_normal_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
  public :: dmft_get_weiss_normal_rank7
#endif

  public :: dmft_get_weiss_superc_rank2
  public :: dmft_get_weiss_superc_rank3
  public :: dmft_get_weiss_superc_rank4
  public :: dmft_get_weiss_superc_rank5
  public :: dmft_get_weiss_superc_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8  
  public :: dmft_get_weiss_superc_rank7
#endif


  complex(8),dimension(:,:),allocatable   :: Wtmp
  complex(8),dimension(:,:,:),allocatable :: Ttmp
  
contains
  


  !##################################################################
  !##################################################################
  !                       NORMAL PHASE
  ! Get the local Weiss Field calG0 using DMFT self-consistency 
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
  subroutine dmft_get_weiss_normal_rank2(Gloc,Sigma,Weiss,Nsites) !N=Nsites*Nso
    complex(8),dimension(:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:),intent(inout) :: Weiss
    integer,optional                          :: Nsites
    complex(8),dimension(:,:,:),allocatable   :: calG0
    !
    Nlat=1;if(present(Nsites))Nlat=Nsites
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !Testing part:
    Ntot  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    if(Ntot*Ntot*Lfreq*16*1d-9 > 1d0)print*, "WARNINGL: dmft_get_weiss_normal_rank2: Memory required exceeds 1Gb per Function!"
    Nso   = Ntot/Nlat ; if(mod(Ntot,Nlat)/=0)stop "dmft_get_weiss_normal_rank2: Ntot%Nsites != 0"
    call assert_shape(Gloc, [Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_normal_rank2","Gloc")
    call assert_shape(Sigma,[Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_normal_rank2","Sigma")
    call assert_shape(Weiss,[Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_normal_rank2","Weiss")
    !
    Weiss = zero
    !
    allocate(calG0(Nlat,Nso,Nso))
    calG0 = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       do ilat=1,Nlat
          call dmft_weiss_normal( get_block(ilat,Gloc(:,:,i)), get_block(ilat,Sigma(:,:,i)), calG0(ilat,:,:))
       enddo
       Weiss(:,:,i) =  from_rank3(calG0)
    enddo MPIloop
    call gf_reduce(Weiss)
    deallocate(calG0)
    !
  end subroutine dmft_get_weiss_normal_rank2

  !G/Sigma shape: [Nlat,Nso,Nso][:]
  !##################################################################
  subroutine dmft_get_weiss_normal_rank3(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:),intent(inout) :: Weiss
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nso   = size(Gloc,2)
    Lfreq = size(Gloc,4)    
    Ntot  = Nlat*Nso
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nso,Nso]:"
    call assert_shape(Sigma,[Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_normal_rank3","Sigma")
    call assert_shape(Gloc, [Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_normal_rank3","Gloc")
    call assert_shape(Weiss,[Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_normal_rank3","Weiss")
    !
    Weiss = zero
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_normal( from_rank3(Gloc(:,:,:,i)), from_rank3(Sigma(:,:,:,i)), Wtmp)
       Weiss(:,:,:,i) = to_rank3(Wtmp)
    enddo MPIloop
    call gf_reduce(Weiss)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_weiss_normal_rank3

  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_normal_rank4(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss

    !MPI setup:
    call set_gf_mpi()
    !
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Nso   = Nspin*Norb
    Ntot  = Nso
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank4","Gloc")
    call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank4","Weiss")
    !
    Weiss = zero
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_normal(from_rank4(Gloc(:,:,:,:,i)), from_rank4(Sigma(:,:,:,:,i)), Wtmp)
       Weiss(:,:,:,:,i) = to_rank4(Wtmp)
    enddo MPIloop
    call gf_reduce(Weiss)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_weiss_normal_rank4

  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_normal_rank5(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss
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
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank5","Gloc")
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank5","Weiss")
    !
    Weiss = zero
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_normal(from_rank5(Gloc(:,:,:,:,:,i)), from_rank5(Sigma(:,:,:,:,:,i)), Wtmp)
       Weiss(:,:,:,:,:,i) = to_rank5(Wtmp)
    enddo MPIloop
    call gf_reduce(Weiss)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_weiss_normal_rank5

  !G/Sigma shape: [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_normal_rank6(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss
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
    if(mpi_master)write(*,"(A)")"Get Weiss function  [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank6","Gloc")
    call assert_shape(Weiss,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank6","Weiss")
    !
    Weiss = zero
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_normal(from_rank6(Gloc(:,:,:,:,:,:,i)),from_rank6(Sigma(:,:,:,:,:,:,i)),Wtmp)
       Weiss(:,:,:,:,:,:,i) = to_rank6(Wtmp)
       ! Weiss(:,:,:,:,:,:,i) = to_rank6( dmft_weiss_normal( from_rank6(Gloc(:,:,:,:,:,:,i)), from_rank6(Sigma(:,:,:,:,:,:,i)) ) )
    enddo MPIloop
    call gf_reduce(Weiss)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_weiss_normal_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  !G/Sigma shape: [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_normal_rank7(Gloc,Sigma,Weiss)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Weiss
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
    if(mpi_master)write(*,"(A)")"Get Weiss function  [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank7","Gloc")
    call assert_shape(Weiss,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_normal_rank7","Weiss")
    !
    Weiss = zero
    !
    if(allocated(Wtmp))deallocate(Wtmp)
    allocate(Wtmp(Ntot,Ntot))
    Wtmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_normal(from_rank7(Gloc(:,:,:,:,:,:,:,i)), from_rank7(Sigma(:,:,:,:,:,:,:,i)),Wtmp)
       Weiss(:,:,:,:,:,:,:,i) = to_rank7(Wtmp)
    enddo MPIloop
    call gf_reduce(Weiss)
    deallocate(Wtmp)
    !
  end subroutine dmft_get_weiss_normal_rank7
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
  subroutine dmft_get_weiss_superc_rank2(Gloc,Floc,Sigma,Self,Weiss,Theta,Nsites) 
    complex(8),dimension(:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:),intent(in)    :: Sigma 
    complex(8),dimension(:,:,:),intent(in)    :: Self  
    complex(8),dimension(:,:,:),intent(inout) :: Weiss 
    complex(8),dimension(:,:,:),intent(inout) :: Theta 
    integer,optional                          :: Nsites
    complex(8),dimension(:,:,:,:),allocatable :: calG0
    !
    Nlat=1;if(present(Nsites))Nlat=Nsites
    !
    !MPI setup:
    call set_gf_mpi()
    !
    !Testing part:
    Ntot  = size(Gloc,1)
    Lfreq = size(Gloc,3)
    if(Ntot*Ntot*Lfreq*16*1d-9 > 1d0)print*, "WARNINGL: dmft_get_weiss_superc_rank2: Memory required exceeds 1Gb per Function!"
    Nso   = Ntot/Nlat
    if(mod(Ntot,Nlat)/=0)stop "dmft_get_weiss_superc_rank2: Ntot%Nsites != 0"
    call assert_shape(Gloc, [Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_superc_rank2","Gloc")
    call assert_shape(Floc, [Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_superc_rank2","Floc")
    call assert_shape(Sigma,[Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_superc_rank2","Sigma")
    call assert_shape(Self, [Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_superc_rank2","Self")
    call assert_shape(Weiss,[Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_superc_rank2","Weiss")
    call assert_shape(Theta,[Ntot,Nlat*Nso,Lfreq],"dmft_get_weiss_superc_rank2","Theta")
    !
    Weiss = zero
    Theta = zero
    !
    allocate(calG0(2,Nlat,Nso,Nso))
    calG0 = zero
    !
    MPIloop: do i=1+mpi_rank,Lfreq,mpi_size
       do ilat=1,Nlat        
          call dmft_weiss_superc(&
               get_block(ilat,Gloc(:,:,i)) ,get_block(ilat,Floc(:,:,i)),&
               get_block(ilat,Sigma(:,:,i)),get_block(ilat,Self(:,:,i)),&
               calG0(1,ilat,:,:),calG0(2,ilat,:,:) )        
       enddo
       Weiss(:,:,i) =  from_rank3(calG0(1,:,:,:))
       Theta(:,:,i) =  from_rank3(calG0(2,:,:,:))
    enddo MPIloop
    call gf_reduce(Weiss)
    call gf_reduce(Theta)
    deallocate(calG0)
    !
  end subroutine dmft_get_weiss_superc_rank2



  !G/Sigma shape: [Nlat,Nso,Nso][:]
  !##################################################################
  subroutine dmft_get_weiss_superc_rank3(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:),intent(inout) :: Weiss
    complex(8),dimension(:,:,:,:),intent(inout) :: Theta
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nlat  = size(Gloc,1)
    Nso   = size(Gloc,2)
    Lfreq = size(Gloc,4)    
    Ntot  = Nlat*Nso
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nso,Nso]:"
    call assert_shape(Sigma,[Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_superc_rank3","Sigma")
    call assert_shape(Gloc, [Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_superc_rank3","Gloc")
    call assert_shape(Weiss,[Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_superc_rank3","Weiss")
    call assert_shape(Self, [Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_superc_rank3","Self")
    call assert_shape(Floc, [Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_superc_rank3","Floc")
    call assert_shape(Theta,[Nlat,Nso,Nso,Lfreq],"dmft_get_weiss_superc_rank3","Theta")
    !
    Weiss = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_superc(&
            from_rank3(Gloc(:,:,:,i)) , from_rank3(Floc(:,:,:,i)), &
            from_rank3(Sigma(:,:,:,i)), from_rank3(Self(:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Weiss(:,:,:,i) = to_rank3(Ttmp(1,:,:))
       Theta(:,:,:,i) = to_rank3(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Weiss)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_weiss_superc_rank3


  !G/Sigma shape: [Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_superc_rank4(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Theta
    !
    !MPI setup:
    call set_gf_mpi()
    !
    Nspin = size(Gloc,1)
    Norb  = size(Gloc,3)
    Lfreq = size(Gloc,5)    
    Nso   = Nspin*Norb
    !
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Sigma")
    call assert_shape(Gloc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Gloc")
    call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Weiss")
    call assert_shape(Self, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Self")
    call assert_shape(Floc, [Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Floc")
    call assert_shape(Theta,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank4","Theta")
    !
    Weiss = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_superc(&
            from_rank4(Gloc(:,:,:,:,i)) , from_rank4(Floc(:,:,:,:,i)), &
            from_rank4(Sigma(:,:,:,:,i)), from_rank4(Self(:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Weiss(:,:,:,:,i)  = to_rank4(Ttmp(1,:,:))
       Theta(:,:,:,:,i)  = to_rank4(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Weiss)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_weiss_superc_rank4


  !G/Sigma shape: [Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_superc_rank5(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Theta
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
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Sigma")
    call assert_shape(Gloc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Gloc")
    call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Weiss")
    call assert_shape(Self, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Self")
    call assert_shape(Floc, [Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Floc")
    call assert_shape(Theta,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank5","Theta")
    !
    Weiss = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_superc(&
            from_rank5(Gloc(:,:,:,:,:,i)) , from_rank5(Floc(:,:,:,:,:,i)), &
            from_rank5(Sigma(:,:,:,:,:,i)), from_rank5(Self(:,:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Weiss(:,:,:,:,:,i) = to_rank5(Ttmp(1,:,:))
       Theta(:,:,:,:,:,i) = to_rank5(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Weiss)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_weiss_superc_rank5

  !G/Sigma shape: [Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_superc_rank6(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Theta
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
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Sigma")
    call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Gloc")
    call assert_shape(Weiss,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Weiss")
    call assert_shape(Self, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Self")
    call assert_shape(Floc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Floc")
    call assert_shape(Theta,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank6","Theta")
    !
    Weiss = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_superc(&
            from_rank6(Gloc(:,:,:,:,:,:,i)) , from_rank6(Floc(:,:,:,:,:,:,i)), &
            from_rank6(Sigma(:,:,:,:,:,:,i)), from_rank6(Self(:,:,:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Weiss(:,:,:,:,:,:,i) = to_rank6(Ttmp(1,:,:))
       Theta(:,:,:,:,:,:,i) = to_rank6(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Weiss)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_weiss_superc_rank6

#if __GFORTRAN__ &&  __GNUC__ > 8
  !G/Sigma shape: [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb][:]
  !##################################################################
  subroutine dmft_get_weiss_superc_rank7(Gloc,Floc,Sigma,Self,Weiss,Theta)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Gloc
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Floc
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Sigma
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in)    :: Self
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Weiss
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(inout) :: Theta
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
    if(mpi_master)write(*,"(A)")"Get Weiss function [Nlat,Nlat,Nspin,Nspin,Norb,Norb]:"
    call assert_shape(Sigma,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank7","Sigma")
    call assert_shape(Gloc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank7","Gloc")
    call assert_shape(Weiss,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank7","Weiss")
    call assert_shape(Self, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank7","Self")
    call assert_shape(Floc, [Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank7","Floc")
    call assert_shape(Theta,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_get_weiss_superc_rank7","Theta")
    !
    Weiss = zero
    Theta = zero
    !
    if(allocated(Ttmp))deallocate(Ttmp)
    allocate(Ttmp(2,Ntot,Ntot))
    Ttmp = zero
    !
    MPIloop:do i=1+mpi_rank,Lfreq,mpi_size
       call dmft_weiss_superc(&
            from_rank7(Gloc(:,:,:,:,:,:,:,i)) , from_rank7(Floc(:,:,:,:,:,:,:,i)), &
            from_rank7(Sigma(:,:,:,:,:,:,:,i)), from_rank7(Self(:,:,:,:,:,:,:,i)), &
            Ttmp(1,:,:), Ttmp(2,:,:))
       Weiss(:,:,:,:,:,:,:,i) = to_rank7(Ttmp(1,:,:))
       Theta(:,:,:,:,:,:,:,i) = to_rank7(Ttmp(2,:,:))
    enddo MPIloop
    call gf_reduce(Weiss)
    call gf_reduce(Theta)
    deallocate(Ttmp)
    !
  end subroutine dmft_get_weiss_superc_rank7
#endif


end module SC_WEISS




