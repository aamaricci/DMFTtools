module DMFT_TIGHT_BINDING
  USE SF_CONSTANTS, only: pi,pi2,xi,one,zero
  USE SF_IOTOOLS
  USE SF_LINALG, only: eigh
  USE SF_COLORS
  implicit none
  private


  interface TB_build_model
     module procedure build_hk_model_kgrid_d
     module procedure build_hk_model_kgrid_c
     module procedure build_hk_model_nkvec_d
     module procedure build_hk_model_nkvec_c
     !
     module procedure build_hkR_model_kgrid_d
     module procedure build_hkR_model_kgrid_c
     module procedure build_hkR_model_nkvec_d
     module procedure build_hkR_model_nkvec_c
     !
     module procedure build_hk_path_d
     module procedure build_hk_path_c
     !
     module procedure build_hkR_path_d
     module procedure build_hkR_path_c
  end interface TB_build_model



  interface TB_solve_model
     module procedure solve_Hk_along_BZpath
     module procedure solve_HkR_along_BZpath
  end interface TB_solve_model



  interface TB_write_hk
     module procedure write_hk_w90_func
     module procedure write_hk_w90_array
  end interface TB_write_hk

  interface TB_read_hk
     module procedure read_hk_w90_array
  end interface TB_read_hk


  interface TB_write_Hloc
     module procedure write_Hloc_1
     module procedure write_Hloc_2
  end interface TB_write_Hloc


  interface TB_read_Hloc
     module procedure read_Hloc_1
     module procedure read_Hloc_2
  end interface TB_read_Hloc


  interface TB_build_kgrid
     module procedure ::   build_kgrid
  !   module procedure :: f_build_kgrid
     module procedure ::   kgrid_from_path_grid
  !   module procedure :: f_kgrid_from_path_grid
     module procedure ::   kgrid_from_path_dim
  !   module procedure :: f_kgrid_from_path_dim
  end interface TB_build_kgrid





  !Some special points in the BZ:
  !we do everything in 3d.
  real(8),dimension(3),public,parameter :: kpoint_gamma=[0,0,0]*pi
  real(8),dimension(3),public,parameter :: kpoint_x1=[1,0,0]*pi
  real(8),dimension(3),public,parameter :: kpoint_x2=[0,1,0]*pi
  real(8),dimension(3),public,parameter :: kpoint_x3=[0,0,1]*pi
  real(8),dimension(3),public,parameter :: kpoint_m1=[1,1,0]*pi
  real(8),dimension(3),public,parameter :: kpoint_m2=[0,1,1]*pi
  real(8),dimension(3),public,parameter :: kpoint_m3=[1,0,1]*pi
  real(8),dimension(3),public,parameter :: kpoint_r=[1,1,1]*pi


  real(8),dimension(3),save             :: bk_x=[1d0,0d0,0d0]*pi2
  real(8),dimension(3),save             :: bk_y=[0d0,1d0,0d0]*pi2
  real(8),dimension(3),save             :: bk_z=[0d0,0d0,1d0]*pi2

  logical,save                          :: set_bkvec=.false.

  public :: TB_set_bk
  public :: TB_build_model
  public :: TB_solve_model
  public :: TB_write_hk
  public :: TB_read_hk
  public :: TB_write_hloc
  public :: TB_read_hloc
  public :: TB_build_kgrid




contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  
  !-------------------------------------------------------------------------------------------
  subroutine TB_set_bk(bkx,bky,bkz)
    real(8),dimension(:),intent(in)          :: bkx
    real(8),dimension(:),intent(in),optional :: bky
    real(8),dimension(:),intent(in),optional :: bkz
    bk_x = 0d0
    bk_x = bkx
    if(present(bky))then
       bk_y = 0d0
       bk_y = bky
    endif
    if(present(bkz))then
       bk_z = 0d0
       bk_z = bkz
    endif
    set_bkvec=.true.
  end subroutine TB_set_bk


  subroutine print_bk(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    unit=6
    if(present(pfile))open(free_unit(unit),file=reg(pfile))
    write(unit,"(A)")"Using Reciprocal Lattice vectors:"
    write(unit,"(A,3F9.4,A1)")"bk_x = [",(bk_x(i),i=1,3),"]"
    write(unit,"(A,3F9.4,A1)")"bk_y = [",(bk_y(i),i=1,3),"]"
    write(unit,"(A,3F9.4,A1)")"bk_z = [",(bk_z(i),i=1,3),"]"
    if(present(pfile))close(unit)
  end subroutine print_bk


  !-------------------------------------------------------------------------------------------
  !PURPOSE:  build the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) Hamiltonian matrix
  ! from the function user defined hk_model procedure.
  !-------------------------------------------------------------------------------------------
  include "tight_binding_build_hk_model.f90"



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  solve the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) along a given linear 
  ! path in the Brillouin Zone. A GNUPLOT script to plot the bands together with their
  ! character is generated.
  !-------------------------------------------------------------------------------------------
  include "tight_binding_solve_hk.f90"



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  read/write the Hamiltonian matrix H(k)
  !-------------------------------------------------------------------------------------------
  include "tight_binding_io_w90hk.f90"



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  read/write the local part of the Hamiltonian to a file
  !-------------------------------------------------------------------------------------------
  include "tight_binding_io_hloc.f90"



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  construct a grid of k-points with minimum information
  !-------------------------------------------------------------------------------------------
  subroutine build_kgrid(Nkvec,kgrid,check_bk)
    integer,dimension(:)                          :: Nkvec ! Nk=product(Nkvec);Ndim=size(Nkvec)
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
    logical,intent(in),optional                   :: check_bk
    logical                                       :: check_bk_
    real(8),dimension(size(Nkvec))                :: kvec  ![Ndim]
    integer                                       :: ik,ix,iy,iz,Nk(3)
    real(8)                                       :: kx,ky,kz
    !
    check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
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
             kvec = kx*bk_x + ky*bk_y + kz*bk_z
             ik=ik+1
             kgrid(ik,:)=kvec
          enddo
       enddo
    enddo
  end subroutine build_kgrid

  function f_build_kgrid(Nkvec,check_bk) result(kgrid)
    integer,dimension(:)                          :: Nkvec ! Nk=product(Nkvec);Ndim=size(Nkvec)
    logical,intent(in),optional                   :: check_bk
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
    logical                                       :: check_bk_
    real(8),dimension(size(Nkvec))                :: kvec  ![Ndim]
    integer                                       :: ik,ix,iy,iz,Nk(3)
    real(8)                                       :: kx,ky,kz
    !
    check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
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
             kvec = kx*bk_x + ky*bk_y + kz*bk_z
             ik=ik+1
             kgrid(ik,:)=kvec
          enddo
       enddo
    enddo
  end function f_build_kgrid





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

  function f_kgrid_from_path_grid(kpath,Nkpath) result(kgrid)
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
  end function f_kgrid_from_path_grid





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

  function f_kgrid_from_path_dim(kpath,Nkpath,dim) result(kxpath)
    real(8),dimension(:,:)                      :: kpath ![Npts][Ndim]
    integer                                     :: Nkpath
    real(8),dimension((size(kpath,1)-1)*Nkpath) :: kxpath ![(Npts-1)*Nkpath][Ndim]
    real(8),dimension(size(kpath,2))            :: kstart,kstop,kpoint,kdiff
    integer                                     :: ipts,ik,ic,dim,Npts
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
  end function f_kgrid_from_path_dim




END MODULE DMFT_TIGHT_BINDING


































! !-------------------------------------------------------------------------------------------
! !PURPOSE: reduce or expand the Hamiltonian in different format (to be expanded).
! !-------------------------------------------------------------------------------------------
! function shrink_Hkr(Hkr,Nlat,Nk,Norb) result(Hkr_)
!   complex(8),dimension(Nlat,Nlat,Nk,Norb,Norb) :: Hkr
!   integer                                      :: Norb,Nlat,Nk
!   integer                                      :: io,jo,ilat,iorb,jorb,jlat,ik
!   complex(8),dimension(Nlat*Norb,Nlat*Norb,Nk) :: Hkr_
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do iorb=1,Norb
!            do jorb=1,Norb
!               io=iorb + (ilat-1)*Norb
!               jo=jorb + (jlat-1)*Norb
!               do ik=1,Nk
!                  Hkr_(io,jo,ik)=Hkr(ilat,jlat,ik,iorb,jorb)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function shrink_Hkr
! !
! function expand_Hkr(Hkr_,Nlat,Nk,Norb) result(Hkr)
!   complex(8),dimension(Nlat,Nlat,Nk,Norb,Norb) :: Hkr
!   integer                                      :: Norb,Nlat,Nk
!   integer                                      :: io,jo,ilat,iorb,jorb,jlat,ik
!   complex(8),dimension(Nlat*Norb,Nlat*Norb,Nk) :: Hkr_
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do iorb=1,Norb
!            do jorb=1,Norb
!               io=iorb + (ilat-1)*Norb
!               jo=jorb + (jlat-1)*Norb
!               do ik=1,Nk
!                  Hkr(ilat,jlat,ik,iorb,jorb)=Hkr_(io,jo,ik)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
! end function expand_Hkr











! function kgrid(Nk,start,len)
!   integer               :: Nk,i
!   real(8),optional      :: start,len
!   real(8)               :: start_,len_
!   real(8),dimension(Nk) :: kgrid
!   start_=-pi ;if(present(start))start_=start
!   len_=2d0*pi;if(present(len))len_=len
!   do i=1,Nk
!      kgrid(i) = start_ + len_*(i-1)/Nk
!   enddo
! end function kgrid
!



!-------------------------------------------------------------------------------------------
!PURPOSE:  obtain the coordinates ix,iy,iz from the lattice index ik
!-------------------------------------------------------------------------------------------
! function indx2ix(ik,ndim) result(ix)
!   integer              :: ik
!   integer              :: ix
!   integer              :: nx_,ny_,nz_
!   integer,dimension(3) :: ndim
!   nx_=ndim(1)
!   ny_=ndim(2)
!   nz_=ndim(3)
!   ix=int(ik-1)/ny_/nz_+1
! end function indx2ix
! !
! function indx2iy(ik,ndim) result(iy)
!   integer              :: ik
!   integer              :: iy
!   integer              :: nx_,ny_,nz_
!   integer,dimension(3) :: ndim
!   nx_=ndim(1)
!   ny_=ndim(2)
!   nz_=ndim(3)
!   iy=mod(int(ik-1)/nz_,ny_)+1
! end function indx2iy
! !
! function indx2iz(ik,ndim) result(iz)
!   integer              :: ik
!   integer              :: iz
!   integer              :: nx_,ny_,nz_
!   integer,dimension(3) :: ndim
!   nx_=ndim(1)
!   ny_=ndim(2)
!   nz_=ndim(3)
!   iz=mod(ik-1,nz_)+1
! end function indx2iz
! !
! subroutine coord2indx(ik,ix,iy,iz,ndim)
!   integer              :: ix,iy,iz
!   integer,dimension(3) :: ndim
!   integer              :: nx,ny,nz
!   integer              :: ik
!   nx=ndim(1)
!   ny=ndim(2)
!   nz=ndim(3)
!   ik = (ix-1)*ny*nz+(iy-1)*nz+iz
! end subroutine coord2indx
! !
! subroutine indx2coord(ik,ix,iy,iz,ndim)
!   integer              :: ix,iy,iz
!   integer,dimension(3) :: ndim
!   integer              :: ik
!   ix=indx2ix(ik,ndim)
!   iy=indx2iy(ik,ndim)
!   iz=indx2iz(ik,ndim)
! end subroutine indx2coord

