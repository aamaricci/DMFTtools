subroutine build_hk_model_kgrid_mpi(MpiComm,Hk,hk_model,Norb,kgrid,wdos)
  integer                                       :: MpiComm
  integer                                       :: Norb
  integer                                       :: Nktot
  real(8),dimension(:,:)                        :: kgrid ![Nktot][Ndim]
  real(8)                                       :: kvec(size(kgrid,2))
  integer                                       :: ik
  complex(8),dimension(Norb,Norb,size(kgrid,1)) :: hk,haux
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)                     :: kpoint
       integer                                  :: N
       complex(8),dimension(N,N)                :: hk_model
     end function hk_model
  end interface
  logical,optional                              :: wdos
  logical                                       :: wdos_
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nktot  = size(kgrid,1)
  !
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = hk_model(kgrid(ik,:),Norb)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
  if(wdos_)then
     allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(MpiComm,Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
     if(mpi_master)call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hk_model_kgrid_mpi



!##################################################################
!##################################################################
!##################################################################



subroutine build_hk_model_Nkvec_mpi(MpiComm,Hk,hk_model,Norb,Nkvec,wdos)
  integer                                        :: MpiComm
  integer,dimension(:),intent(in)                :: Nkvec
  integer                                        :: Norb
  real(8),dimension(product(Nkvec),size(Nkvec))  :: kgrid ![Nk][Ndim]
  integer                                        :: Nktot
  integer                                        :: ik
  complex(8),dimension(Norb,Norb,product(Nkvec)) :: hk,haux
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
!
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = hk_model(kgrid(ik,:),Norb)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
  if(wdos_)then
     allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hk_model_Nkvec_mpi





!##################################################################
!##################################################################
!##################################################################





subroutine build_hkr_model_kgrid_mpi(MpiComm,hk,hkr_model,Nlat,Norb,kgrid,pbc,wdos)
  integer                                                 :: MpiComm
  integer                                                 :: Nlat,Norb
  logical                                                 :: pbc
  real(8),dimension(:,:)                                  :: kgrid ![Nktot][Ndim]
  integer                                                 :: Nktot
  integer                                                 :: ik
  complex(8),dimension(Nlat*Norb,Nlat*Norb,size(kgrid,1)) :: hk,haux
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                               :: kpoint
       integer                                            :: Nlat,Norb
       logical                                            :: pbc
       complex(8),dimension(Nlat*Norb,Nlat*Norb)          :: hkr_model
     end function hkr_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nktot  = size(kgrid,1)
  !
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
  if(wdos_)then
     allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hkr_model_kgrid_mpi





!##################################################################
!##################################################################
!##################################################################






subroutine build_hkr_model_nkvec_mpi(MpiComm,hk,hkr_model,Nlat,Norb,Nkvec,pbc,wdos)
  integer                                                  :: MpiComm
  integer                                                  :: Nlat,Norb
  logical                                                  :: pbc
  integer,dimension(:),intent(in)                          :: Nkvec
  real(8),dimension(product(Nkvec),size(Nkvec))            :: kgrid ![Nk][Ndim]
  integer                                                  :: Nktot
  integer                                                  :: ik
  complex(8),dimension(Nlat*Norb,Nlat*Norb,product(Nkvec)) :: hk,haux
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                                :: kpoint
       integer                                             :: Nlat,Norb
       logical                                             :: pbc
       complex(8),dimension(Nlat*Norb,Nlat*Norb)           :: hkr_model
     end function hkr_model
  end interface
  logical,optional                                         :: wdos
  logical                                                  :: wdos_
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  !
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
  if(wdos_)then
     allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hkr_model_nkvec_mpi





!##################################################################
!##################################################################
!##################################################################




subroutine build_hk_path_mpi(MpiComm,hk,hk_model,Norb,kpath,Nkpath)
  integer                                        :: MpiComm
  integer                                                   :: Norb
  real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
  integer                                                   :: Nkpath
  integer                                                   :: Npts,Nktot
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  complex(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)  :: hk,haux
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  !
  Npts  =  size(kpath,1)          !# of k-points along the path
  Nktot = (Npts-1)*Nkpath
  !
  call TB_build_kgrid(kpath,Nkpath,kgrid)
  !
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
end subroutine build_hk_path_mpi








!##################################################################
!##################################################################
!##################################################################







subroutine build_hkR_path_mpi(MpiComm,hk,hkr_model,Nlat,Norb,kpath,Nkpath,pbc)
  integer                                        :: MpiComm
  integer                                                             :: Nlat,Norb
  logical                                                             :: pbc
  real(8),dimension(:,:)                                              :: kpath ![Npts][Ndim]
  integer                                                             :: Nkpath
  integer                                                             :: Npts,Nktot
  integer                                                             :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2))           :: kgrid
  complex(8),dimension(Nlat*Norb,Nlat*Norb,(size(kpath,1)-1)*Nkpath)  :: hk,haux
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                      :: kpoint
       integer                                   :: Nlat,Norb
       logical                                   :: pbc
       complex(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
     end function hkr_model
  end interface
  !
  Npts  =  size(kpath,1)          !# of k-points along the path
  Nktot = (Npts-1)*Nkpath
  !
  call TB_build_kgrid(kpath,Nkpath,kgrid)
  !
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
end subroutine build_hkR_path_mpi




!##################################################################
!##################################################################
!##################################################################


subroutine build_hk_w90_mpi(MpiComm,Hk,Nlso,Nkvec,Kpts_grid,wdos)
  integer                                                :: MpiComm
  integer                                                :: Nlso
  integer,dimension(:),intent(in)                        :: Nkvec
  complex(8),dimension(Nlso,Nlso,product(Nkvec))         :: Hk,haux
  !
  real(8),dimension(product(Nkvec),size(Nkvec)),optional :: Kpts_grid ![Nk][Ndim]
  logical,optional                                       :: wdos
  logical                                                :: wdos_
  !
  real(8),dimension(product(Nkvec),size(Nkvec))          :: kgrid ![Nk][Ndim]
  integer                                                :: Nk
  logical                                                :: IOfile
  integer                                                :: i,j,ik
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nk =product(Nkvec)
  !
  if(.not.TB_w90%status)stop "build_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
  !
  call TB_build_kgrid(Nkvec,Kgrid,.true.) !check bk_1,2,3 vectors have been set
  if(present(Kpts_grid))Kpts_grid=Kgrid
  !
  do ik=1+mpi_rank,Nk,mpi_size
     haux(:,:,ik) = w90_hk_model(Kgrid(ik,:),Nlso)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
  if(wdos_)then
     allocate(dos_Greal(TB_w90%Nlat,TB_w90%Nspin,TB_w90%Nspin,TB_w90%Norb,TB_w90%Norb,dos_Lreal))
     allocate(dos_wtk(Nk))
     dos_wtk=1d0/Nk
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(TB_w90%Nlat,TB_w90%Nspin,TB_w90%Nspin,TB_w90%Norb,TB_w90%Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal,trim(dos_file),iprint=1)
  endif
end subroutine build_hk_w90_mpi




subroutine build_hk_w90_path_mpi(MpiComm,Hk,Nlso,Kpath,Nkpath)
  integer                                                   :: MpiComm
  integer                                                   :: Nlso
  real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
  integer                                                   :: Nkpath
  integer                                                   :: Npts,Nktot
  complex(8),dimension(Nlso,Nlso,(size(kpath,1)-1)*Nkpath)  :: Hk,haux
  !
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  integer                                                   :: Nk
  logical                                                   :: IOfile
  integer                                                   :: i,j,ik
  !
  if(.not.TB_w90%status)stop "build_hk_w90_path: TB_w90 structure not allocated. Call setup_w90 first."
  !
  Npts  =  size(kpath,1)          !# of k-points along the path
  Nktot = (Npts-1)*Nkpath
  !
  call TB_build_kgrid(kpath,Nkpath,kgrid)
  !
  do ik=1+mpi_rank,Nktot,mpi_size
     haux(:,:,ik) = w90_hk_model(Kgrid(ik,:),Nlso)
  enddo
  Hk = zero
  call AllReduce_MPI(MpiComm,Haux,Hk)
  !
end subroutine build_hk_w90_path_mpi

















! subroutine build_hk_model_kgrid_d(Hk,hk_model,Norb,kgrid,wdos)
!   integer                                    :: Norb
!   integer                                    :: Nktot
!   real(8),dimension(:,:)                     :: kgrid ![Nktot][Ndim]
!   real(8)                                    :: kvec(size(kgrid,2))
!   integer                                    :: ik
!   real(8),dimension(Norb,Norb,size(kgrid,1)) :: hk
!   interface 
!      function hk_model(kpoint,N)
!        real(8),dimension(:)                  :: kpoint
!        integer                               :: N
!        real(8),dimension(N,N)                :: hk_model
!      end function hk_model
!   end interface
!   logical,optional                           :: wdos
!   logical                                    :: wdos_
!   wdos_=.false.;if(present(wdos))wdos_=wdos
!   !
!   Nktot  = size(kgrid,1)
!   !
!   do ik=1,Nktot
!      Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
!   enddo
!   !
!   if(wdos_)then
!      allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
!      allocate(dos_wtk(Nktot))
!      dos_wtk=1d0/Nktot
!      call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
!      call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
!   endif
! end subroutine build_hk_model_kgrid_d

! subroutine build_hk_model_nkvec_d(Hk,hk_model,Norb,Nkvec,wdos)
!   integer,dimension(:),intent(in)               :: Nkvec
!   integer                                       :: Norb
!   real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
!   integer                                       :: Nktot
!   integer                                       :: ik
!   real(8),dimension(Norb,Norb,product(Nkvec))   :: hk
!   interface 
!      function hk_model(kpoint,N)
!        real(8),dimension(:)    :: kpoint
!        integer                 :: N
!        real(8),dimension(N,N)  :: hk_model
!      end function hk_model
!   end interface
!   logical,optional                           :: wdos
!   logical                                    :: wdos_
!   wdos_=.false.;if(present(wdos))wdos_=wdos
!   !
!   call TB_build_kgrid(Nkvec,kgrid,.true.)
!   !
!   Nktot  = product(Nkvec)
!   do ik=1,Nktot
!      Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
!   enddo
!   !
!   if(wdos_)then
!      allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
!      allocate(dos_wtk(Nktot))
!      dos_wtk=1d0/Nktot
!      call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
!      call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
!   endif
! end subroutine build_hk_model_nkvec_d

! subroutine build_hkr_model_kgrid_d(hk,hkr_model,Nlat,Norb,kgrid,pbc,wdos)
!   integer                                              :: Nlat,Norb
!   logical                                              :: pbc
!   real(8),dimension(:,:)                               :: kgrid ![Nktot][Ndim]
!   integer                                              :: Nktot
!   integer                                              :: ik
!   real(8),dimension(Nlat*Norb,Nlat*Norb,size(kgrid,1)) :: hk
!   interface 
!      function hkr_model(kpoint,Nlat,Norb,pbc)
!        real(8),dimension(:)                            :: kpoint
!        integer                                         :: Nlat,Norb
!        logical                                         :: pbc
!        real(8),dimension(Nlat*Norb,Nlat*Norb)          :: hkr_model
!      end function hkr_model
!   end interface
!   logical,optional                           :: wdos
!   logical                                    :: wdos_
!   wdos_=.false.;if(present(wdos))wdos_=wdos
!   !
!   Nktot  = size(kgrid,1)
!   !
!   do ik=1,Nktot
!      Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
!   enddo
!   !
!   if(wdos_)then
!      allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
!      allocate(dos_wtk(Nktot))
!      dos_wtk=1d0/Nktot
!      call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
!      call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
!   endif
! end subroutine build_hkr_model_kgrid_d

! subroutine build_hkr_model_nkvec_d(hk,hkr_model,Nlat,Norb,Nkvec,pbc,wdos)
!   integer                                               :: Nlat,Norb
!   logical                                               :: pbc
!   integer,dimension(:),intent(in)                       :: Nkvec
!   real(8),dimension(product(Nkvec),size(Nkvec))         :: kgrid ![Nk][Ndim]
!   integer                                               :: Nktot
!   integer                                               :: ik
!   real(8),dimension(Nlat*Norb,Nlat*Norb,product(Nkvec)) :: hk
!   interface 
!      function hkr_model(kpoint,Nlat,Norb,pbc)
!        real(8),dimension(:)                             :: kpoint
!        integer                                          :: Nlat,Norb
!        logical                                          :: pbc
!        real(8),dimension(Nlat*Norb,Nlat*Norb)           :: hkr_model
!      end function hkr_model
!   end interface
!   logical,optional                           :: wdos
!   logical                                    :: wdos_
!   wdos_=.false.;if(present(wdos))wdos_=wdos
!   !
!   call TB_build_kgrid(Nkvec,kgrid,.true.)
!   !
!   Nktot  = product(Nkvec)
!   !
!   do ik=1,Nktot
!      Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
!   enddo
!   !
!   if(wdos_)then
!      allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
!      allocate(dos_wtk(Nktot))
!      dos_wtk=1d0/Nktot
!      call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
!      call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
!   endif
! end subroutine build_hkr_model_nkvec_d

! subroutine build_hk_path_d(hk,hk_model,Norb,kpath,Nkpath)
!   integer                                                   :: Norb
!   real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
!   integer                                                   :: Nkpath
!   integer                                                   :: Npts,Nktot
!   integer                                                   :: ik,i
!   real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
!   real(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)     :: hk
!   interface 
!      function hk_model(kpoint,N)
!        real(8),dimension(:)   :: kpoint
!        integer                :: N
!        real(8),dimension(N,N) :: hk_model
!      end function hk_model
!   end interface
!   !
!   Npts  =  size(kpath,1)          !# of k-points along the path
!   Nktot = (Npts-1)*Nkpath
!   !
!   call TB_build_kgrid(kpath,Nkpath,kgrid)
!   !
!   do ik=1,Nktot
!      Hk(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
!   enddo
!   !
! end subroutine build_hk_path_d

! subroutine build_hkR_path_d(hk,hkr_model,Nlat,Norb,kpath,Nkpath,pbc)
!   integer                                                          :: Nlat,Norb
!   logical                                                          :: pbc
!   real(8),dimension(:,:)                                           :: kpath ![Npts][Ndim]
!   integer                                                          :: Nkpath
!   integer                                                          :: Npts,Nktot
!   integer                                                          :: ik,i
!   real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2))        :: kgrid
!   real(8),dimension(Nlat*Norb,Nlat*Norb,(size(kpath,1)-1)*Nkpath) :: hk
!   interface 
!      function hkr_model(kpoint,Nlat,Norb,pbc)
!        real(8),dimension(:)                   :: kpoint
!        integer                                :: Nlat,Norb
!        logical                                :: pbc
!        real(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
!      end function hkr_model
!   end interface
!   !
!   Npts  =  size(kpath,1)          !# of k-points along the path
!   Nktot = (Npts-1)*Nkpath
!   !
!   call TB_build_kgrid(kpath,Nkpath,kgrid)
!   !
!   do ik=1,Nktot
!      Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
!   enddo
!   !
! end subroutine build_hkR_path_d
