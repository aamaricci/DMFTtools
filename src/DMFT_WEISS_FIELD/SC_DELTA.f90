module SC_DELTA
  USE DMFT_CTRL_VARS
  USE SC_COMMON
  implicit none
  private




  interface dmft_delta
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_cluster
     module procedure :: dmft_get_delta_normal_ineq
     module procedure :: dmft_get_delta_normal_bethe
     module procedure :: dmft_get_delta_normal_bethe_ineq
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_ineq
  end interface dmft_delta


  public :: dmft_delta


contains




  subroutine dmft_get_delta_normal_main(Gloc,Smats,Delta,Hloc)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    !aux
    complex(8),dimension(:,:,:),allocatable       :: zeta_site ![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: Smats_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: invGloc_site![Nspin*Norb][Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable       :: calG0_site![Nspin*Norb][Nspin*Norb][Lmats]
    integer                                       :: Nspin,Norb,Nso,Lmats
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
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
    if(mpi_master)then
       !Retrieve parameters:
       call get_ctrl_var(beta,"BETA")
       call get_ctrl_var(xmu,"XMU")
       !
       !Testing part:
       Nspin = size(Gloc,1)
       Norb  = size(Gloc,3)
       Lmats = size(Gloc,5)
       Nso   = Nspin*Norb
       call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_main","Gloc")
       call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_main","Smats")
       call assert_shape(Delta,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_main","Delta")
       call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_main","Hloc")
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
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Hloc,Nspin,Norb) - &
               nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
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
    endif
#ifdef _MPI    
    if(check_MPI())call MPI_BCAST(Delta,size(Delta),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpi_ierr)
#endif
    !
  end subroutine dmft_get_delta_normal_main


  subroutine dmft_get_delta_normal_cluster(Gloc,Smats,Delta,Hloc)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc  ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Delta ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Hloc  ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    !aux
    complex(8),dimension(:,:,:),allocatable           :: zeta_site     ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: Smats_site    ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: invGloc_site  ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: calG0_site    ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
    integer                                           :: Nlat,Nspin,Norb,Nlso,Lmats
    integer                                           :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
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
    if(mpi_master)then
       !Retrieve parameters:
       call get_ctrl_var(beta,"BETA")
       call get_ctrl_var(xmu,"XMU")
       !
       !Testing part:
       Nlat  = size(Gloc,1)
       Nspin = size(Gloc,3)
       Norb  = size(Gloc,5)
       Lmats = size(Gloc,7)
       Nlso  = Nlat*Nspin*Norb
       call assert_shape(Gloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_cluster","Gloc")
       call assert_shape(Smats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_cluster","Smats")
       call assert_shape(Delta,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_cluster","Delta")
       call assert_shape(Hloc, [Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_cluster","Hloc")
       !
       if(allocated(wm))deallocate(wm)
       allocate(wm(Lmats))
       allocate(zeta_site(Nlso,Nlso,Lmats))
       allocate(Smats_site(Nlso,Nlso,Lmats))
       allocate(invGloc_site(Nlso,Nlso,Lmats))
       allocate(calG0_site(Nlso,Nlso,Lmats))
       !
       wm = pi/beta*(2*arange(1,Lmats)-1)
       !Dump the Gloc and the Smats into a [Nlat*Norb*Nspin]^2 matrix and create the zeta_site
       do i=1,Lmats
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nlso) - nnn2lso_reshape(Hloc,Nlat,Nspin,Norb) - &
               nnn2lso_reshape(Smats(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
          invGloc_site(:,:,i) = nnn2lso_reshape(Gloc(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
          Smats_site(:,:,i)   = nnn2lso_reshape(Smats(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
       enddo
       !
       !Invert the ilat-th site [Nlat*Norb*Nspin]**2 Gloc block matrix 
       do i=1,Lmats
          call inv(invGloc_site(:,:,i))
       enddo
       !
       ! [Delta]_ilat = [Zeta-Hloc-Sigma]_ilat - [Gloc]_ilat^-1
       calG0_site(:,:,1:Lmats) = zeta_site(:,:,1:Lmats) - invGloc_site(:,:,1:Lmats)
       !
       !Dump back the [Nlat*Norb*Nspin]**2 block of the ilat-th site into the 
       !output structure of [Nlat,Nlat,Nspin,Nspin,Norb,Norb] matrix
       do i=1,Lmats
          Delta(:,:,:,:,:,:,i) = lso2nnn_reshape(calG0_site(:,:,i),Nlat,Nspin,Norb)
       enddo
    endif
#ifdef _MPI    
    if(check_MPI())call MPI_BCAST(Delta,size(Delta),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpi_ierr)
#endif
    !
  end subroutine dmft_get_delta_normal_cluster

  subroutine dmft_get_delta_normal_ineq(Gloc,Smats,Delta,Hloc)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
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
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_ineq_mpi","Gloc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_ineq_mpi","Smats")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_ineq_mpi","Delta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_ineq_mpi","Hloc")
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
          zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nso) - &
               nn2so_reshape(Hloc(ilat,:,:,:,:),Nspin,Norb) - &
               nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
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
#ifdef _MPI    
    if(check_MPI())then
       call Mpi_AllReduce(Delta_tmp, Delta, size(Delta), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Delta=Delta_tmp
    endif
#else
    Delta=Delta_tmp
#endif
    !
  end subroutine dmft_get_delta_normal_ineq



  subroutine dmft_get_delta_normal_bethe(Gloc,Delta,Hloc,Wbands)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    real(8),dimension(:),intent(in)               :: Wbands ![Nspin*Norb]
    !aux
    integer                                       :: Nspin,Norb,Nso,Lmats
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
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
    if(mpi_master)then
       !Testing part:
       Nspin = size(Gloc,1)
       Norb  = size(Gloc,3)
       Lmats = size(Gloc,5)
       Nso   = Nspin*Norb
       call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_main","Gloc")
       call assert_shape(Delta,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_main","Delta")
       call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_main","Hloc")
       !
       if(allocated(wm))deallocate(wm)
       allocate(wm(Lmats))
       !
       wm = pi/beta*(2*arange(1,Lmats)-1)
       !
       !\calG0^{-1}_aa = d**2/4*Gmats
       Delta=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             Delta(ispin,ispin,iorb,iorb,:) = 0.25d0*Wbands(iorb+(ispin-1)*Norb)**2*Gloc(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
       deallocate(wm)
       !
    endif
#ifdef _MPI    
    if(check_MPI())call MPI_BCAST(Delta,size(Delta),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpi_ierr)
#endif
  end subroutine dmft_get_delta_normal_bethe




  subroutine dmft_get_delta_normal_bethe_ineq(Gloc,Delta,Hloc,Wbands)
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    real(8),dimension(:),intent(in)                 :: Wbands       ! [Nlat*Nspin*Norb]
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Delta_tmp    ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
    integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
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
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_bethe_ineq","Gloc")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_normal_bethe_ineq","Delta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_normal_ineq","Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Delta_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta_tmp = zero
    Delta     = zero
    MPIloop: do ilat=1+mpi_rank,Nlat,mpi_size
       !\calG0^{-1}_aa_i = d**2/4*Gmats_aa_i
       Delta=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             io       = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
             Delta_tmp(ilat,ispin,ispin,iorb,iorb,:) = 0.25d0*Wbands(io)**2*Gloc(ilat,ispin,ispin,iorb,iorb,:)
          enddo
       enddo
    end do MPIloop
    !
    !
#ifdef _MPI    
    if(check_MPI())then
       call Mpi_AllReduce(Delta_tmp, Delta, size(Delta), MPI_Double_Complex, MPI_Sum, MPI_COMM_WORLD, mpi_ierr)
    else
       Delta=Delta_tmp
    endif
#else
    Delta=Delta_tmp
#endif
    !
    !
    deallocate(wm)
    !
  end subroutine dmft_get_delta_normal_bethe_ineq




  subroutine dmft_get_delta_superc_main(Gloc,Floc,Smats,SAmats,Delta,Fdelta,Hloc)
    integer                                         :: MPI_COMM_WORLD
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc         ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Floc         ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats        !
    complex(8),dimension(:,:,:,:,:),intent(in)      :: SAmats        !
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta        !
    complex(8),dimension(:,:,:,:,:),intent(inout) :: fDelta        !
    complex(8),dimension(:,:,:,:),intent(in)        :: Hloc         ! [Nspin][Nspin][Norb][Norb]
    !aux
    complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    integer                                         :: Nspin,Norb,Nso,Nso2,Lmats
    integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
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
    if(mpi_master)then
       !Retrieve parameters:
       call get_ctrl_var(beta,"BETA")
       call get_ctrl_var(xmu,"XMU")
       !
       !Testing part:
       Nspin = size(Gloc,1)
       Norb  = size(Gloc,3)
       Lmats = size(Gloc,5)
       Nso   = Nspin*Norb
       Nso2  = 2*Nso
       call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_main_mpi","Gloc")
       call assert_shape(Floc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_main_mpi","Floc")
       call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_main_mpi","Smats")
       call assert_shape(SAmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_main_mpi","SAmats")
       call assert_shape(Delta,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_main_mpi","Delta")
       call assert_shape(fDelta,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_main_mpi","fDelta")
       call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_delta_superc_main_mpi","Hloc")
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
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - &
                        Hloc(ispin,jspin,iorb,jorb) - Smats(ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =  - SAmats(ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =  - SAmats(ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + &
                        Hloc(ispin,jspin,iorb,jorb) + conjg(Smats(ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Floc(ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Floc(ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = SAmats(ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = SAmats(ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(ispin,jspin,iorb,jorb,:))
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
                   Delta(ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   fDelta(ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    endif
#ifdef _MPI    
    if(check_MPI())then
       call MPI_BCAST(Delta,size(Delta),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpi_ierr)
       call MPI_BCAST(fDelta,size(fDelta),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpi_ierr)
    endif
#endif
    !
    !
  end subroutine dmft_get_delta_superc_main



  subroutine dmft_get_delta_superc_ineq(Gloc,Floc,Smats,SAmats,Delta,fDelta,Hloc)
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Floc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats        !
    complex(8),dimension(:,:,:,:,:,:),intent(in)      :: SAmats        ! 
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Delta        !
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: fDelta        ! 
    complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
    !aux
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Delta_tmp        !
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: fDelta_tmp        ! 
    complex(8),dimension(:,:,:),allocatable           :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
    complex(8),dimension(:,:,:),allocatable           :: Smats_site   !
    complex(8),dimension(:,:,:),allocatable           :: invGloc_site !
    complex(8),dimension(:,:,:),allocatable           :: calG0_site   !
    integer                                           :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
    integer                                           :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
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
    !Testing part:
    Nlat  = size(Gloc,1)
    Nspin = size(Gloc,2)
    Norb  = size(Gloc,4)
    Lmats = size(Gloc,6)
    Nso   = Nspin*Norb
    Nso2  = 2*Nso
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_ineq_mpi","Gloc")
    call assert_shape(Floc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_ineq_mpi","Floc")
    call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_ineq_mpi","Smats")
    call assert_shape(SAmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_ineq_mpi","SAmats")
    call assert_shape(Delta,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_ineq_mpi","Delta")
    call assert_shape(fDelta,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_delta_superc_ineq_mpi","fDelta")
    call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_delta_superc_ineq_mpi","Hloc")
    !
    if(allocated(wm))deallocate(wm)
    allocate(wm(Lmats))
    allocate(Delta_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(fDelta_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(zeta_site(Nso2,Nso2,Lmats))
    allocate(Smats_site(Nso2,Nso2,Lmats))
    allocate(invGloc_site(Nso2,Nso2,Lmats))
    allocate(calG0_site(Nso2,Nso2,Lmats))
    !
    wm = pi/beta*(2*arange(1,Lmats)-1)
    Delta_tmp   = zero
    fDelta_tmp   = zero
    Delta       = zero
    fDelta       = zero
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
                   zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io,jo+Nso,:)       =-SAmats(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo,:)       =-SAmats(ilat,ispin,jspin,iorb,jorb,:)
                   zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(ilat,ispin,jspin,iorb,jorb,:))
                   !
                   invGloc_site(io,jo,:)        = Gloc(ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io,jo+Nso,:)    = Floc(ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo,:)    = Floc(ilat,ispin,jspin,iorb,jorb,:)
                   invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(ilat,ispin,jspin,iorb,jorb,:))
                   !
                   Smats_site(io,jo,:)          = Smats(ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io,jo+Nso,:)      = SAmats(ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo,:)      = SAmats(ilat,ispin,jspin,iorb,jorb,:)
                   Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(ilat,ispin,jspin,iorb,jorb,:))
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
                   Delta_tmp(ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                   fDelta_tmp(ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
                enddo
             enddo
          enddo
       enddo
    end do MPIloop
    !
#ifdef _MPI    
    if(check_MPI())then
       call MPI_Allreduce(Delta_tmp, Delta, size(Delta), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
       call MPI_Allreduce(fDelta_tmp, fDelta, size(fDelta), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    else
       Delta=Delta_tmp
       fDelta=fDelta_tmp
    endif
#else
    Delta=Delta_tmp
    fDelta=fDelta_tmp
#endif
    !
  end subroutine dmft_get_delta_superc_ineq

end module SC_DELTA
