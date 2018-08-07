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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
  wdos_=.true.;if(present(wdos))wdos_=wdos
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
        print*,ilat,jlat
        !>Build Hij
        Hij(:,:,ilat,jlat) = Hij(:,:,ilat,jlat) + ts_model(ilink,Nso)
     enddo do_links
     print*,""
  enddo do_lattice
  return
end subroutine build_Hij_Nrvec





! subroutine build_hk_model_1d_d(hk,hk_model,Nkvec,check_bk)
!   interface 
!      function hk_model(kpoint)
!        real(8),dimension(:) :: kpoint
!        real(8)              :: hk_model
!      end function hk_model
!   end interface
!   integer,dimension(:),intent(in)       :: Nkvec
!   logical,intent(in),optional           :: check_bk
!   logical                               :: check_bk_
!   real(8),dimension(1,1,product(Nkvec)) :: hk1
!   real(8),dimension(product(Nkvec))     :: hk
!   !
!   check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
!   !
!   call build_hk_model_Norb_d(hk1,hk_model1,1,Nkvec,check_bk_)
!   hk = hk1(1,1,:)
! contains
!   function hk_model1(kpoint,N)
!     real(8),dimension(:)    :: kpoint
!     integer                 :: N
!     real(8),dimension(1,1)  :: hk_model1
!     hk_model1 = hk_model(kpoint)
!   end function hk_model1
! end subroutine build_hk_model_1_d

! subroutine build_hk_model_1_c(hk,hk_model,Nkvec,check_bk)
!   interface 
!      function hk_model(kpoint)
!        real(8),dimension(:)      :: kpoint
!        complex(8)                :: hk_model
!      end function hk_model
!   end interface
!   integer,dimension(:),intent(in)          :: Nkvec
!   logical,intent(in),optional              :: check_bk
!   logical                                  :: check_bk_
!   complex(8),dimension(1,1,product(Nkvec)) :: hk1
!   complex(8),dimension(product(Nkvec))     :: hk
!   !
!   check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
!   !
!   call build_hk_model_Norb_c(hk1,hk_model1,1,Nkvec,check_bk_)
!   hk = hk1(1,1,:)
! contains
!   function hk_model1(kpoint,N)
!     real(8),dimension(:)      :: kpoint
!     integer                   :: N
!     complex(8),dimension(1,1) :: hk_model1
!     hk_model1 = hk_model(kpoint)
!   end function hk_model1
! end subroutine build_hk_model_1_c










! subroutine build_hkr_model_1_d(hk,hkr_model,Nlat,Norb,Nkvec,pbc,check_bk) 
!   interface 
!      function hkr_model(kpoint,Nlat,pbc)
!        real(8),dimension(:)         :: kpoint
!        integer                      :: Nlat
!        real(8),dimension(Nlat,Nlat) :: hkr_model
!        logical                      :: pbc
!      end function hkr_model
!   end interface
!   integer                                     :: Nlat,Norb
!   integer,dimension(:),intent(in)             :: Nkvec
!   logical                                     :: pbc
!   logical,intent(in),optional                 :: check_bk
!   logical                                     :: check_bk_
!   integer                                     :: Nktot,Nk(3)
!   integer                                     :: ik,ix,iy,iz
!   real(8)                                     :: kx,ky,kz,kvec(3)
!   real(8),dimension(Nlat,Nlat,product(Nkvec)) :: hk
!   !
!   check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
!   !
!   Nk=1
!   do ik=1,size(Nkvec)
!      Nk(ik)=Nkvec(ik)
!   enddo
!   Nktot  = product(Nk)
!   if(Nktot/=product(Nkvec))stop "TB_build_model ERROR: product(Nkvec) != Nktot"
!   !
!   call print_bk()
!   if(check_bk_.AND..not.set_bkvec)stop "TB_build_model ERROR: bk vectors not set"
!   !
!   ik=0
!   do iz=1,Nk(3)
!      kz = dble(iz-1)/Nk(3)
!      do iy=1,Nk(2)
!         ky = dble(iy-1)/Nk(2)
!         do ix=1,Nk(1)
!            kx = dble(ix-1)/Nk(1)
!            kvec = kx*bk_x + ky*bk_y + kz*bk_z
!            ik = ik+1
!            Hk(:,:,ik) = hkr_model(kvec,Nlat,pbc)
!         enddo
!      enddo
!   enddo
!   !
! end subroutine build_hkr_model_1_d


! subroutine build_hkr_model_1_c(hk,hkr_model,Nlat,Norb,Nkvec,pbc,check_bk) 
!   interface 
!      function hkr_model(kpoint,Nlat,pbc)
!        real(8),dimension(:)            :: kpoint
!        integer                         :: Nlat
!        complex(8),dimension(Nlat,Nlat) :: hkr_model
!        logical                         :: pbc
!      end function hkr_model
!   end interface
!   integer                                        :: Nlat,Norb
!   integer,dimension(:),intent(in)                :: Nkvec
!   logical                                        :: pbc
!   logical,intent(in),optional                    :: check_bk
!   logical                                        :: check_bk_
!   integer                                        :: Nktot,Nk(3)
!   integer                                        :: ik,ix,iy,iz
!   real(8)                                        :: kx,ky,kz,kvec(3)
!   complex(8),dimension(Nlat,Nlat,product(Nkvec)) :: hk
!   !
!   check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
!   !
!   Nk=1
!   do ik=1,size(Nkvec)
!      Nk(ik)=Nkvec(ik)
!   enddo
!   Nktot  = product(Nk)
!   if(Nktot/=product(Nkvec))stop "TB_build_model ERROR: product(Nkvec) != Nktot"
!   !
!   call print_bk()
!   if(check_bk_.AND..not.set_bkvec)stop "TB_build_model ERROR: bk vectors not set"
!   !
!   ik=0
!   do iz=1,Nk(3)
!      kz = dble(iz-1)/Nk(3)
!      do iy=1,Nk(2)
!         ky = dble(iy-1)/Nk(2)
!         do ix=1,Nk(1)
!            kx = dble(ix-1)/Nk(1)
!            kvec = kx*bk_x + ky*bk_y + kz*bk_z
!            ik = ik+1
!            Hk(:,:,ik) = hkr_model(kvec,Nlat,pbc)
!         enddo
!      enddo
!   enddo
!   !
! end subroutine build_hkr_model_1_c










! subroutine build_hkpath_model_1_d(hk,hk_model,kpath,Nkpath)
!   interface 
!      function hk_model(kpoint)
!        real(8),dimension(:) :: kpoint
!        real(8)              :: hk_model
!      end function hk_model
!   end interface
!   real(8),dimension(:,:)                                    :: kpath
!   integer                                                   :: Nkpath
!   integer                                                   :: Ndim,Npts,Nk
!   integer                                                   :: ik,i
!   real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
!   real(8),dimension(size(kpath,2))                          :: kvec
!   real(8),dimension((size(kpath,1)-1)*Nkpath)               :: hk
!   Npts=size(kpath,1)          !# of k-points along the path
!   Ndim=size(kpath,2)          !# of dimension of the k-vectors
!   Nk  = (Npts-1)*Nkpath
!   if(Ndim>3)stop "TB_build_model error: Ndim > 3"
!   do i=1,Ndim
!      kgrid(:,i) = kgrid_from_path(kpath,Npts,Nkpath,i)
!   enddo
!   do ik=1,Nk
!      Hk(ik) = hk_model(kgrid(ik,:))
!   enddo
! end subroutine build_hkpath_model_1_d

! subroutine build_hkpath_model_1_c(hk,hk_model,kpath,Nkpath) 
!   interface 
!      function hk_model(kpoint)
!        real(8),dimension(:) :: kpoint
!        complex(8)           :: hk_model
!      end function hk_model
!   end interface
!   real(8),dimension(:,:)                                    :: kpath
!   integer                                                   :: Nkpath
!   integer                                                   :: Ndim,Npts,Nk
!   integer                                                   :: ik,i
!   real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
!   real(8),dimension(size(kpath,2))                          :: kvec
!   complex(8),dimension((size(kpath,1)-1)*Nkpath)            :: hk
!   Npts=size(kpath,1)          !# of k-points along the path
!   Ndim=size(kpath,2)          !# of dimension of the k-vectors
!   Nk  = (Npts-1)*Nkpath
!   if(Ndim>3)stop "TB_build_model error: Ndim > 3"
!   do i=1,Ndim
!      kgrid(:,i) = kgrid_from_path(kpath,Npts,Nkpath,i)
!   enddo
!   do ik=1,Nk
!      Hk(ik) = hk_model(kgrid(ik,:))
!   enddo
! end subroutine build_hkpath_model_1_c





















!> this is an example taken from vasp!
!>>
!The second method generates k-meshes automatically, and requires only the input of subdivisions of the Brillouin zone in each
!direction and the origin ('shift') for the k-mesh. There are three possible input formats. The simplest one is only supported
!by VASP.4.5 and newer versions:
!
!   Automatic mesh
! 0              ! number of k-points = 0 ->automatic generation scheme 
! Auto           ! fully automatic
! 10             ! length (l)

!As before, the first line is treated as a comment. On the second line a number smaller or equal 0 must be specified. In the
!previous section, this value supplied the number of k-points, a zero in this line activates the automatic generation scheme.
!The fully automatic scheme, selected by the first character in the third line ('a'), generates $ \Gamma $ centered
!Monkhorst-Pack grids, where the numbers of subdivisions along each reciprocal lattice vector are given by

!   $\displaystyle N_1 = {\rm max}(1,l*\vert{\vec b}_1\vert+0.5)
!   $
!   $\displaystyle N_2 = {\rm max}(1,l*\vert{\vec b}_2\vert+0.5)
!   $
!   $\displaystyle N_3 = {\rm max}(1,l*\vert{\vec b}_3\vert+0.5).
!   $
!$ {\vec b}_i$ are the reciprocal lattice vectors, and  $ \vert{\vec b}_i\vert$ is their norm. VASP generates a equally
!spaced k-point grid with the coordinates:
!$\displaystyle {\vec k} = {\vec b}_1 \frac{n_1}{N_1} + {\vec b}_2 \frac{n_2}{N_2......rac{n_3}{N_3} ,
!\qquad n_1=0...,N_1-1 \quad n_2=0...,N_2-1 \quad n_3=0...,N_3-1 $
!Symmetry is used to map equivalent k-points to each other, which can reduce the total number of k-points significantly.
!Useful values for the length vary between 10 (large gap insulators) and 100 (d-metals).
!A slightly enhanced version, allows to supply the numbers for the subdivisions $ N_1$, $ N_2$ and $ N_3$ manually:

!Automatic mesh
! 0             ! number of k-points = 0 ->automatic generation scheme 
!Gamma          ! generate a Gamma centered grid
! 4  4  4       ! subdivisions N_1, N_2 and N_3 along recipr. l. vectors
! 0 . 0. 0.     ! optional shift of the mesh (s_1, s_2, s_3)
!
!In this case, the third linemight start with 'G' or 'g' --for generating meshes with their origin at the $ \Gamma $
!point (as above)-- or 'M' or 'm', which selects the original Monkhorst-Pack scheme.
!In the latter case k-point grids, with even ( $ {\rm mod}(N_i,2)=0)$ subdivisions are shifted off $ \Gamma $:
!
! $\displaystyle {\vec k} = {\vec b}_1 \frac{n_1+1/2}{N_1} + {\vec b}_2 \frac{n_2+1/2}{N_2} + {\vec b}_3 \frac{n_3+1/2}{N_3}$
!
!The fifth line is optional, and supplies an additional shift of the k-mesh (compared to the origin used in the Gamma
! centered or Monkhorst-Pack case). Usually the shift is zero, since the two most important cases are covered by the flags
! 'M' or 'm', 'G' or 'g'. The shift must be given in multiples of the length of the reciprocal lattice vectors, and the
! generated grids are then ('G' case):
!
!   $\displaystyle {\vec k} = {\vec b}_1 \frac{n_1+s_1}{N_1} + {\vec b}_2 \frac{n_2+s_2}{N_2} + {\vec b}_3 \frac{n_3+s_3}{N_3}.$
!
!   and ('M' case):
!
!   $\displaystyle {\vec k} = {\vec b}_1 \frac{n_1+s_1+1/2}{N_1} + {\vec b}_2 \frac{n_2+s_2+1/2}{N_2} + {\vec b}_3 \frac{n_3+s_3+1/2}{N_3}.$
!
!The selection 'M' without shift, is obviously equivalent to 'G' with a shift of 0.5 0.5 0.5, and vice versa.
!If the third line does not start with 'M', 'm', 'G' or 'g' an alternative input mode is selected. this mode is mainly
!for experts, and should not be used for casual VASP users. Here one can provide directly the generating basis vectors
!for the k-point mesh (in cartesian or reciprocal coordinates). The input file has the following format:
!
!Automatic generation
! 0 
!   Cartesian
! 0 .25 0.00 0.00
! 0 .00 0.25 0.00
! 0 .00 0.00 0.25
! 0 .00 0.00 0.00
!The entry in the third line switches between cartesian and reciprocal coordinates (as in the explicit input format described first
! - key characters 'C', 'c', 'K' or 'k' switch to cartesian coordinates). On the fifth, sixth and seventh line the generating basis
! vectors must be given and the eighth line supplies the shift (if one likes to shift the k-mesh off $ \Gamma $, default is to take
! the origin at $ \Gamma $, the shift is given in multiples of the generating basis vectors, only (0,0,0) and (1/2,1/2,1/2) and
! arbitrary combinations are usually usefull). This method can always be replaced by an appropriate Monkhorst-Pack setting. For
! instance for a fcc lattice the setting:
!   cart
! 0 .25 0 0
! 0 0 .25 0
! 0 0 0 .25
! 0 .5 0.5 0.5
!   is equivalent to
!   Monkhorst-pack
! 4 4 4 
! 0 0 0 
!   This input scheme is especially interesting to build meshes, which are based on the conventional cell (for instance sc
! for fcc and bcc), or the primitive cell if a large super cell is used. In the example above the k-point mesh is based on
! the sc-cell. (for the second input file the tetrahedron method can not be used because the shift breaks the symmetry
! (see below), whereas the first input file can be used together with the tetrahedron method).
!
!   Mind: The division scheme (or the generating basis of the k-mesh) must lead to a k-mesh which belongs to the same class
! of Bravais lattice as the reciprocal unit cell (Brillouin zone). 
