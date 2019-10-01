subroutine dmft_get_weiss_normal_main_mpi(MpiComm,Gloc,Smats,Weiss,Hloc)
  integer                                       :: MpiComm
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
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
  mpi_size  = get_size_MPI(MpiComm)
  mpi_rank =  get_rank_MPI(MpiComm)
  mpi_master= get_master_MPI(MpiComm)
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
     !
     call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_main","Gloc")
     call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_main","Smats")
     call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_main","Weiss")
     call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_main","Hloc")
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
  endif
  call Bcast_Mpi(MpiComm,Weiss)
  !
end subroutine dmft_get_weiss_normal_main_mpi



subroutine dmft_get_weiss_normal_cluster_mpi(MpiComm,Gloc,Smats,Weiss,Hloc)
  integer                                       :: MpiComm
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Gloc  ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Hloc  ! [Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  !aux
  complex(8),dimension(:,:,:),allocatable           :: zeta_site ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable           :: Smats_site![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable           :: invGloc_site![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable           :: calG0_site![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lmats]
  integer                                           :: Nlat,Nspin,Norb,Nlso,Lmats
  integer                                           :: i,ilat,jlat,iorb,jorb,ispin,jspin,io,jo
  !
  !MPI setup:
  mpi_size  = get_size_MPI(MpiComm)
  mpi_rank =  get_rank_MPI(MpiComm)
  mpi_master= get_master_MPI(MpiComm)
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
     Nlso   = Nspin*Norb
     call assert_shape(Gloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_cluster","Gloc")
     call assert_shape(Smats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_cluster","Smats")
     call assert_shape(Weiss,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_cluster","Weiss")
     call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_cluster","Hloc")
     !
     if(allocated(wm))deallocate(wm)
     allocate(wm(Lmats))
     allocate(zeta_site(Nlso,Nlso,Lmats))
     allocate(Smats_site(Nlso,Nlso,Lmats))
     allocate(invGloc_site(Nlso,Nlso,Lmats))
     allocate(calG0_site(Nlso,Nlso,Lmats))
     !
     wm = pi/beta*(2*arange(1,Lmats)-1)
     !Dump the Gloc and the Smats into a [Norb*Nspin]^2 matrix and create the zeta_site
     do i=1,Lmats
        zeta_site(:,:,i)    = (xi*wm(i)+xmu)*eye(Nlso) - nnn2lso_reshape(Hloc,Nlat,Nspin,Norb) -&
                              nnn2lso_reshape(Smats(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
        invGloc_site(:,:,i) = nnn2lso_reshape(Gloc(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
        Smats_site(:,:,i)   = nnn2lso_reshape(Smats(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
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
        Weiss(:,:,:,:,:,:,i) = lso2nnn_reshape(calG0_site(:,:,i),Nlat,Nspin,Norb)
     enddo
  endif
  call Bcast_Mpi(MpiComm,Weiss)
  !
end subroutine dmft_get_weiss_normal_cluster_mpi


subroutine dmft_get_weiss_normal_ineq_mpi(MpiComm,Gloc,Smats,Weiss,Hloc)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
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
  mpi_size  = get_size_MPI(MpiComm)
  mpi_rank =  get_rank_MPI(MpiComm)
  mpi_master= get_master_MPI(MpiComm)
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
  call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_ineq_mpi","Gloc")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_ineq_mpi","Smats")
  call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_ineq_mpi","Weiss")
  call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_ineq_mpi","Hloc")
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
end subroutine dmft_get_weiss_normal_ineq_mpi
