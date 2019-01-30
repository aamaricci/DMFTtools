subroutine dmft_get_weiss_normal_bethe_mpi(MpiComm,Gloc,Weiss,Hloc,Wbands)
  integer                                       :: MpiComm
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
  real(8),dimension(:),intent(in)               :: Wbands ![Nspin*Norb]
  !aux
  complex(8),dimension(size(Gloc,5))            :: invWeiss ![Lmats]
  integer                                       :: Nspin,Norb,Nso,Lmats
  integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
  !
  !MPI setup:
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
     call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_main","Gloc")
     call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_main","Weiss")
     call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_main","Hloc")
     !
     if(allocated(wm))deallocate(wm)
     allocate(wm(Lmats))
     !
     wm = pi/beta*(2*arange(1,Lmats)-1)
     !
     !\calG0^{-1}_aa = iw + mu - H_0 - d**2/4*Gmats
     Weiss=zero
     do ispin=1,Nspin
        do iorb=1,Norb
           invWeiss = (xi*wm(:)+xmu) - Hloc(ispin,ispin,iorb,iorb) - &
                0.25d0*Wbands(iorb+(ispin-1)*Norb)**2*Gloc(ispin,ispin,iorb,iorb,:)
           Weiss(ispin,ispin,iorb,iorb,:) = one/invWeiss
        enddo
     enddo
     !
     deallocate(wm)
     !
  endif
  call Bcast_Mpi(MpiComm,Weiss)
end subroutine dmft_get_weiss_normal_bethe_mpi




subroutine dmft_get_weiss_normal_bethe_ineq_mpi(MpiComm,Gloc,Weiss,Hloc,Wbands)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
  real(8),dimension(:),intent(in)                 :: Wbands       ! [Nlat*Nspin*Norb]
  !aux
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Weiss_tmp    ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(size(Gloc,6))              :: invWeiss    ![Lmats]
  integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
  integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
  !
  !MPI setup:
  mpi_size  = get_size_MPI(MpiComm)
  mpi_rank =  get_rank_MPI(MpiComm)
  mpi_master= get_master_MPI(MpiComm)
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
  call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_bethe_ineq","Gloc")
  call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_bethe_ineq","Weiss")
  call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_ineq","Hloc")
  !
  if(allocated(wm))deallocate(wm)
  allocate(wm(Lmats))
  allocate(Weiss_tmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  !
  wm = pi/beta*(2*arange(1,Lmats)-1)
  Weiss_tmp = zero
  Weiss     = zero
  MPIloop: do ilat=1+mpi_rank,Nlat,mpi_size
     !\calG0^{-1}_aa_i = iw + mu - H_0 - d**2/4*Gmats_aa_i
     Weiss=zero
     do ispin=1,Nspin
        do iorb=1,Norb
           io       = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
           invWeiss = (xi*wm(:)+xmu) - Hloc(ilat,ispin,ispin,iorb,iorb) - 0.25d0*Wbands(io)**2*Gloc(ilat,ispin,ispin,iorb,iorb,:)
           Weiss_tmp(ilat,ispin,ispin,iorb,iorb,:) = one/invWeiss
        enddo
     enddo
  end do MPIloop
  !
  !
  call Mpi_AllReduce(Weiss_tmp, Weiss, size(Weiss), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  !
  !
  deallocate(wm)
  !
end subroutine dmft_get_weiss_normal_bethe_ineq_mpi
