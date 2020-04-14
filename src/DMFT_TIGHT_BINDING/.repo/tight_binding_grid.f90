subroutine build_kgrid(Nkvec,kgrid,check_bk)
  integer,dimension(:)                          :: Nkvec ! Nk=product(Nkvec);Ndim=size(Nkvec)
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
  logical,intent(in),optional                   :: check_bk
  logical                                       :: check_bk_
  real(8),dimension(size(Nkvec))                :: kvec  ![Ndim]
  integer                                       :: ik,ix,iy,iz,Nk(3),ndim
  real(8)                                       :: kx,ky,kz
  !
  check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
  !
  ndim = size(Nkvec)          !dimension of the grid to be built
  !
  Nk=1
  do ik=1,size(Nkvec)
     Nk(ik)=Nkvec(ik)
  enddo
  if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
  !
  call print_bk()
  if(check_bk_.AND..not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
  !
  ik=0
  do iz=1,Nk(3)
     kz = dble(iz-1)/Nk(3)
     do iy=1,Nk(2)
        ky = dble(iy-1)/Nk(2)
        do ix=1,Nk(1)
           kx = dble(ix-1)/Nk(1)
           kvec = kx*bk_x(:ndim) + ky*bk_y(:ndim) + kz*bk_z(:ndim)
           ik=ik+1
           kgrid(ik,:)=kvec
        enddo
     enddo
  enddo
end subroutine build_kgrid

subroutine build_kgrid_generic(Nkvec,kgrid_x,kgrid_y,kgrid_z,check_bk)
  integer,dimension(:)                          :: Nkvec   ! Nk=product(Nkvec);Ndim=size(Nkvec)
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_x ![Nk][Ndim]
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_y ![Nk][Ndim]
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_z ![Nk][Ndim]
  logical,intent(in),optional                   :: check_bk
  logical                                       :: check_bk_
  real(8),dimension(size(Nkvec))                :: kvec    ![Ndim]
  integer                                       :: ik,ix,iy,iz,Nk(3),ndim
  real(8)                                       :: kx,ky,kz
  !
  check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
  !
  ndim = size(Nkvec)          !dimension of the grid to be built
  !
  Nk=1
  do ik=1,size(Nkvec)
     Nk(ik)=Nkvec(ik)
  enddo
  if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
  !
  call print_bk()
  if(check_bk_.AND..not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
  !
  ik=0
  do iz=1,Nk(3)
     kz = dble(iz-1)/Nk(3)
     do iy=1,Nk(2)
        ky = dble(iy-1)/Nk(2)
        do ix=1,Nk(1)
           kx = dble(ix-1)/Nk(1)
           ik=ik+1
           kgrid_x(ik,:)=kx*bk_x(:ndim)
           kgrid_y(ik,:)=ky*bk_y(:ndim)
           kgrid_z(ik,:)=kz*bk_z(:ndim)
        enddo
     enddo
  enddo
end subroutine build_kgrid_generic


subroutine build_rgrid(Nrvec,Rgrid,check_ei)
  integer,dimension(:)                          :: Nrvec ! Nr=product(Nrvec);Ndim=size(Nrvec)
  real(8),dimension(product(Nrvec),size(Nrvec)) :: Rgrid ![Nk][Ndim]
  logical,intent(in),optional                   :: check_ei
  logical                                       :: check_ei_
  real(8),dimension(size(Nrvec))                :: Rvec  ![Ndim]
  integer                                       :: ir,ix,iy,iz,Nr(3)
  !
  check_ei_=.false.;if(present(check_ei))check_ei_=check_ei
  !
  Nr=1
  do ir=1,size(Nrvec)
     Nr(ir)=Nrvec(ir)
  enddo
  if(product(Nr)/=product(Nrvec))stop "TB_build_Rgrid ERROR: product(Nrvec) != product(Nr)"
  !
  call print_ei()
  if(check_ei_.AND..not.set_eivec)stop "TB_build_Rgrid ERROR: Ei vectors not set"
  !
  ir=0
  do iz=1,Nr(3)
     do iy=1,Nr(2)
        do ix=1,Nr(1)
           Rvec = ix*ei_x + iy*ei_y + iz*ei_z
           ir=ir+1
           Rgrid(ir,:)=Rvec
        enddo
     enddo
  enddo
end subroutine build_rgrid



subroutine TB_write_grid(Grid,file_grid)
  real(8),dimension(:,:) :: Grid ![Nk/Nr,Ndim]
  character(len=*)       :: file_grid
  integer :: i,j,unit
  open(free_unit(unit),file=reg(file_grid),status="unknown",action="write",position="rewind")
  do i=1,size(Grid,1)
     write(unit,'(3F15.7)') (Grid(i,j),j=1,size(Grid,2))
  enddo
  close(unit)
end subroutine TB_write_grid




subroutine TB_build_CoordGrid(Nrvec,RNgrid)
  integer,dimension(:)                          :: Nrvec  ! Nr=product(Nrvec);Ndim=size(Nrvec)
  integer,dimension(product(Nrvec),size(Nrvec)) :: RNgrid ![Nr][Ndim]
  integer                                       :: ir,ix,iy,iz,Nr(3),Rvec(3)
  !
  Nr=1
  do ir=1,size(Nrvec)
     Nr(ir)=Nrvec(ir)
  enddo
  !
  ir=0
  do iz=1,Nr(3)
     do iy=1,Nr(2)
        do ix=1,Nr(1)
           ir=ir+1
           Rvec = [ix,iy,iz]
           RNgrid(ir,:) = Rvec(1:size(Nrvec))
        enddo
     enddo
  enddo
  ir=0
end subroutine TB_build_CoordGrid


function TB_find_IndxCoord(RNvec,Nrvec) result(ir)
  integer,dimension(:)           :: Nrvec
  integer,dimension(size(Nrvec)) :: RNvec
  integer                        :: Nr(3),Rvec(3)
  integer                        :: ir,ix,iy,iz
  logical :: bool
  !
  Nr=1
  Rvec=1
  do ir=1,size(Nrvec)
     Nr(ir)=Nrvec(ir)
     Rvec(ir)=RNvec(ir)
  enddo
  !
  ir=0
  do iz=1,Nr(3)
     do iy=1,Nr(2)
        do ix=1,Nr(1)
           ir=ir+1
           bool = Rvec(1)==ix.AND.Rvec(2)==iy.AND.Rvec(3)==iz
           if(bool)return
        enddo
     enddo
  enddo
  ir=0
end function TB_find_IndxCoord




subroutine kgrid_from_path_grid(kpath,Nkpath,kgrid)
  real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
  integer                                                   :: Nkpath
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid ![(Npts-1)*Nkpath][Ndim]
  real(8),dimension(size(kpath,2))                          :: kstart,kstop,kpoint,kdiff
  integer                                                   :: ipts,ik,ic,dim,Npts
  Npts = size(kpath,1)
  ic=0
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/dble(Nkpath)
     do ik=1,Nkpath
        ic=ic+1
        kpoint = kstart + (ik-1)*kdiff
        kgrid(ic,:)=kpoint
     enddo
  enddo
end subroutine kgrid_from_path_grid

subroutine kgrid_from_path_dim(kpath,Nkpath,dim,kxpath)
  real(8),dimension(:,:)                      :: kpath ![Npts][Ndim]
  integer                                     :: Nkpath
  integer                                     :: dim
  real(8),dimension((size(kpath,1)-1)*Nkpath) :: kxpath ![(Npts-1)*Nkpath][Ndim]
  real(8),dimension(size(kpath,2))            :: kstart,kstop,kpoint,kdiff
  integer                                     :: ipts,ik,ic,Npts
  Npts = size(kpath,1)
  if(dim>size(kpath,2))stop "TB_build_grid ERROR: dim > Ndim"
  ic=0
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/dble(Nkpath)
     do ik=1,Nkpath
        ic=ic+1
        kpoint = kstart + (ik-1)*kdiff
        kxpath(ic)=kpoint(dim)
     enddo
  enddo
end subroutine kgrid_from_path_dim
