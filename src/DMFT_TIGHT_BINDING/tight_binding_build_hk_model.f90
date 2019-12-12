subroutine build_hk_model_kgrid_d(Hk,hk_model,Norb,kgrid,wdos)
  integer                                    :: Norb
  integer                                    :: Nktot
  real(8),dimension(:,:)                     :: kgrid ![Nktot][Ndim]
  real(8)                                    :: kvec(size(kgrid,2))
  integer                                    :: ik
  real(8),dimension(Norb,Norb,size(kgrid,1)) :: hk
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)                  :: kpoint
       integer                               :: N
       real(8),dimension(N,N)                :: hk_model
     end function hk_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nktot  = size(kgrid,1)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hk_model_kgrid_d

subroutine build_hk_model_kgrid_c(Hk,hk_model,Norb,kgrid,wdos) 
  integer                                       :: Norb
  integer                                       :: Nktot
  real(8),dimension(:,:)                        :: kgrid ![Nktot][Ndim]
  real(8)                                       :: kvec(size(kgrid,2))
  integer                                       :: ik
  complex(8),dimension(Norb,Norb,size(kgrid,1)) :: hk
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)                     :: kpoint
       integer                                  :: N
       complex(8),dimension(N,N)                :: hk_model
     end function hk_model
  end interface
  logical,optional                              :: wdos
  logical                                       :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nktot  = size(kgrid,1)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hk_model_kgrid_c




!##################################################################
!##################################################################
!##################################################################




subroutine build_hk_model_nkvec_d(Hk,hk_model,Norb,Nkvec,wdos)
  integer,dimension(:),intent(in)               :: Nkvec
  integer                                       :: Norb
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
  integer                                       :: Nktot
  integer                                       :: ik
  real(8),dimension(Norb,Norb,product(Nkvec))   :: hk
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)    :: kpoint
       integer                 :: N
       real(8),dimension(N,N)  :: hk_model
     end function hk_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  do ik=1,Nktot
     Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hk_model_nkvec_d

subroutine build_hk_model_Nkvec_c(Hk,hk_model,Norb,Nkvec,wdos) 
  integer,dimension(:),intent(in)               :: Nkvec
  integer                                        :: Norb
  real(8),dimension(product(Nkvec),size(Nkvec))  :: kgrid ![Nk][Ndim]
  integer                                        :: Nktot
  integer                                        :: ik
  complex(8),dimension(Norb,Norb,product(Nkvec)) :: hk
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  do ik=1,Nktot
     Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hk_model_Nkvec_c





!##################################################################
!##################################################################
!##################################################################





subroutine build_hkr_model_kgrid_d(hk,hkr_model,Nlat,Norb,kgrid,pbc,wdos)
  integer                                              :: Nlat,Norb
  logical                                              :: pbc
  real(8),dimension(:,:)                               :: kgrid ![Nktot][Ndim]
  integer                                              :: Nktot
  integer                                              :: ik
  real(8),dimension(Nlat*Norb,Nlat*Norb,size(kgrid,1)) :: hk
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                            :: kpoint
       integer                                         :: Nlat,Norb
       logical                                         :: pbc
       real(8),dimension(Nlat*Norb,Nlat*Norb)          :: hkr_model
     end function hkr_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nktot  = size(kgrid,1)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hkr_model_kgrid_d

subroutine build_hkr_model_kgrid_c(hk,hkr_model,Nlat,Norb,kgrid,pbc,wdos)
  integer                                                 :: Nlat,Norb
  logical                                                 :: pbc
  real(8),dimension(:,:)                                  :: kgrid ![Nktot][Ndim]
  integer                                                 :: Nktot
  integer                                                 :: ik
  complex(8),dimension(Nlat*Norb,Nlat*Norb,size(kgrid,1)) :: hk
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
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nktot  = size(kgrid,1)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hkr_model_kgrid_c





!##################################################################
!##################################################################
!##################################################################






subroutine build_hkr_model_nkvec_d(hk,hkr_model,Nlat,Norb,Nkvec,pbc,wdos)
  integer                                               :: Nlat,Norb
  logical                                               :: pbc
  integer,dimension(:),intent(in)                       :: Nkvec
  real(8),dimension(product(Nkvec),size(Nkvec))         :: kgrid ![Nk][Ndim]
  integer                                               :: Nktot
  integer                                               :: ik
  real(8),dimension(Nlat*Norb,Nlat*Norb,product(Nkvec)) :: hk
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                             :: kpoint
       integer                                          :: Nlat,Norb
       logical                                          :: pbc
       real(8),dimension(Nlat*Norb,Nlat*Norb)           :: hkr_model
     end function hkr_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(one*Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hkr_model_nkvec_d

subroutine build_hkr_model_nkvec_c(hk,hkr_model,Nlat,Norb,Nkvec,pbc,wdos)
  integer                                                  :: Nlat,Norb
  logical                                                  :: pbc
  integer,dimension(:),intent(in)                          :: Nkvec
  real(8),dimension(product(Nkvec),size(Nkvec))            :: kgrid ![Nk][Ndim]
  integer                                                  :: Nktot
  integer                                                  :: ik
  complex(8),dimension(Nlat*Norb,Nlat*Norb,product(Nkvec)) :: hk
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                                :: kpoint
       integer                                             :: Nlat,Norb
       logical                                             :: pbc
       complex(8),dimension(Nlat*Norb,Nlat*Norb)           :: hkr_model
     end function hkr_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
     allocate(dos_wtk(Nktot))
     dos_wtk=1d0/Nktot
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
  endif
end subroutine build_hkr_model_nkvec_c





!##################################################################
!##################################################################
!##################################################################




subroutine build_hk_path_d(hk,hk_model,Norb,kpath,Nkpath)
  integer                                                   :: Norb
  real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
  integer                                                   :: Nkpath
  integer                                                   :: Npts,Nktot
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  real(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)     :: hk
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)   :: kpoint
       integer                :: N
       real(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  !
  Npts  =  size(kpath,1)          !# of k-points along the path
  Nktot = (Npts-1)*Nkpath
  !
  call TB_build_kgrid(kpath,Nkpath,kgrid)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
  enddo
  !
end subroutine build_hk_path_d

subroutine build_hk_path_c(hk,hk_model,Norb,kpath,Nkpath) 
  integer                                                   :: Norb
  real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
  integer                                                   :: Nkpath
  integer                                                   :: Npts,Nktot
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  complex(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)  :: hk
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
  do ik=1,Nktot
     Hk(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
  enddo
  !
end subroutine build_hk_path_c








!##################################################################
!##################################################################
!##################################################################







subroutine build_hkR_path_d(hk,hkr_model,Nlat,Norb,kpath,Nkpath,pbc)
  integer                                                          :: Nlat,Norb
  logical                                                          :: pbc
  real(8),dimension(:,:)                                           :: kpath ![Npts][Ndim]
  integer                                                          :: Nkpath
  integer                                                          :: Npts,Nktot
  integer                                                          :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2))        :: kgrid
  real(8),dimension(Nlat*Norb,Nlat*Norb,(size(kpath,1)-1)*Nkpath) :: hk
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                   :: kpoint
       integer                                :: Nlat,Norb
       logical                                :: pbc
       real(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
     end function hkr_model
  end interface
  !
  Npts  =  size(kpath,1)          !# of k-points along the path
  Nktot = (Npts-1)*Nkpath
  !
  call TB_build_kgrid(kpath,Nkpath,kgrid)
  !
  do ik=1,Nktot
     Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  !
end subroutine build_hkR_path_d

subroutine build_hkR_path_c(hk,hkr_model,Nlat,Norb,kpath,Nkpath,pbc)
  integer                                                             :: Nlat,Norb
  logical                                                             :: pbc
  real(8),dimension(:,:)                                              :: kpath ![Npts][Ndim]
  integer                                                             :: Nkpath
  integer                                                             :: Npts,Nktot
  integer                                                             :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2))           :: kgrid
  complex(8),dimension(Nlat*Norb,Nlat*Norb,(size(kpath,1)-1)*Nkpath) :: hk
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
  do ik=1,Nktot
     Hk(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
  enddo
  !
end subroutine build_hkR_path_c









!##################################################################
!##################################################################
!##################################################################







subroutine build_Hij_Nrvec(Hij,ts_model,Nso,Nrvec,Links,pbc,wdos)
  integer                                                     :: Nso
  integer,dimension(:),intent(in)                             :: Nrvec
  integer,dimension(:,:),intent(in)                           :: Links ![Nlink][dim]
  logical,optional                                            :: pbc
  logical                                                     :: pbc_
  integer,dimension(product(Nrvec),size(Nrvec))               :: RNgrid
  integer                                                     :: Nlat,Nlink
  integer                                                     :: ilat,jlat
  integer,dimension(size(Nrvec))                              :: Ri,Rj
  integer                                                     :: i,ilink
  complex(8),dimension(Nso,Nso,product(Nrvec),product(Nrvec)) :: Hij
  integer                                                     :: ir,ix,iy,iz,Nr(3)
  !
  interface
     function ts_model(link,Nso)
       integer                       :: link
       integer                       :: Nso
       complex(8),dimension(Nso,Nso) :: ts_model
     end function ts_model
  end interface
  logical,optional                           :: wdos
  logical                                    :: wdos_
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  pbc_ = .true. ; if(present(pbc))pbc_=pbc
  !
  if(size(Nrvec)/=size(Links,2))stop "TB_build_Hij ERROR: size(Nvec) != size(Links,2)"
  !
  Nlat  = product(Nrvec)
  Nlink = size(Links,1)
  !
  call TB_build_CoordGrid(Nrvec,RNgrid)
  !
  Hij = zero
  !
  do_lattice: do ilat = 1,Nlat
     Hij(:,:,ilat,ilat) = Hij(:,:,ilat,ilat) + ts_model(0,Nso)
     Ri = RNgrid(ilat,:)
     do_links: do ilink=1,Nlink
        Rj = Ri + Links(ilink,:)
        if(pbc_)then
           do i=1,size(Nrvec)
              if( Rj(i)==0 )Rj(i)=Nrvec(i)
              if( Rj(i)==Nrvec(i)+1)Rj(i)=1
           enddo
        endif
        jlat = TB_find_IndxCoord(Rj,Nrvec)
        if(jlat==0)cycle do_links
        !
        !>Build Hij
        Hij(:,:,ilat,jlat) = Hij(:,:,ilat,jlat) + ts_model(ilink,Nso)
     enddo do_links
  enddo do_lattice
  return
end subroutine build_Hij_Nrvec









!##################################################################
!##################################################################
!##################################################################










subroutine build_hk_w90(Hk,Nlso,Nkvec,Kpts_grid,wdos)
  integer                                                :: Nlso
  integer,dimension(:),intent(in)                        :: Nkvec
  complex(8),dimension(Nlso,Nlso,product(Nkvec))         :: Hk
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
  wdos_=.false.;if(present(wdos))wdos_=wdos
  !
  Nk =product(Nkvec)
  !
  if(.not.TB_w90%status)stop "build_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
  !
  call TB_build_kgrid(Nkvec,Kgrid,.true.) !check bk_1,2,3 vectors have been set
  if(present(Kpts_grid))Kpts_grid=Kgrid
  !
  do ik = 1,Nk
     Hk(:,:,ik) = w90_hk_model(Kgrid(ik,:),Nlso)
  enddo
  !
  if(wdos_)then
     allocate(dos_Greal(TB_w90%Nlat,TB_w90%Nspin,TB_w90%Nspin,TB_w90%Norb,TB_w90%Norb,dos_Lreal))
     allocate(dos_wtk(Nk))
     dos_wtk=1d0/Nk
     call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(TB_w90%Nlat,TB_w90%Nspin,TB_w90%Nspin,TB_w90%Norb,TB_w90%Norb,dos_Lreal))
     call dmft_print_gf_realaxis(dos_Greal,trim(dos_file),iprint=1)
  endif
end subroutine build_hk_w90




subroutine build_hk_w90_path(Hk,Nlso,Kpath,Nkpath)
  integer                                                   :: Nlso
  real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
  integer                                                   :: Nkpath
  integer                                                   :: Npts,Nktot
  complex(8),dimension(Nlso,Nlso,(size(kpath,1)-1)*Nkpath)  :: Hk
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
  do ik = 1,Nktot
     Hk(:,:,ik) = w90_hk_model(Kgrid(ik,:),Nlso)
  enddo
  !
end subroutine build_hk_w90_path





