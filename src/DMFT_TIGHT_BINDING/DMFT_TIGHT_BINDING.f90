module DMFT_TIGHT_BINDING
  USE SF_CONSTANTS, only: pi,pi2,xi,one,zero
  USE SF_IOTOOLS
  USE SF_LINALG, only: eigh,det
  USE SF_COLORS
  USE SF_TIMER, only:start_timer,stop_timer,eta
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
     !
     module procedure build_Hij_Nrvec
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
     module procedure ::   kgrid_from_path_grid
     module procedure ::   kgrid_from_path_dim
  end interface TB_build_kgrid


  interface TB_build_Rgrid
     module procedure ::   build_Rgrid
  end interface TB_build_Rgrid


  interface TB_print_bk
     module procedure :: print_bk
  end interface TB_print_bk

  interface TB_print_ei
     module procedure :: print_ei
  end interface TB_print_ei



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


  real(8),dimension(3),save             :: ei_x=[1d0,0d0,0d0]
  real(8),dimension(3),save             :: ei_y=[0d0,1d0,0d0]
  real(8),dimension(3),save             :: ei_z=[0d0,0d0,1d0]

  real(8),dimension(3),save             :: bk_x=[1d0,0d0,0d0]*pi2
  real(8),dimension(3),save             :: bk_y=[0d0,1d0,0d0]*pi2
  real(8),dimension(3),save             :: bk_z=[0d0,0d0,1d0]*pi2


  logical,save                          :: set_eivec=.false.
  logical,save                          :: set_bkvec=.false.



  public :: TB_set_ei
  public :: TB_set_bk
  !
  public :: TB_get_ei
  public :: TB_get_bk
  !
  public :: TB_reset_ei
  public :: TB_reset_bk
  !
  public :: TB_build_ei
  public :: TB_build_bk
  !
  public :: TB_print_ei
  public :: TB_print_bk
  !
  public :: TB_reciprocal_basis
  !
  public :: TB_build_kgrid
  public :: TB_build_Rgrid
  !
  public :: TB_build_CoordGrid
  public :: TB_find_IndxCoord
  !
  public :: TB_build_model
  public :: TB_solve_model
  !
  public :: TB_write_hk
  public :: TB_read_hk
  !
  public :: TB_write_hloc
  public :: TB_read_hloc






contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  
  !-------------------------------------------------------------------------------------------
  subroutine TB_reset_ei
    ei_x=[1d0,0d0,0d0]
    ei_y=[0d0,1d0,0d0]
    ei_z=[0d0,0d0,1d0]
    set_eivec=.false.
  end subroutine TB_reset_ei

  
  subroutine TB_reset_bk
    bk_x=[1d0,0d0,0d0]*pi2
    bk_y=[0d0,1d0,0d0]*pi2
    bk_z=[0d0,0d0,1d0]*pi2
    set_bkvec=.false.
  end subroutine TB_reset_bk


  subroutine TB_set_ei(eix,eiy,eiz)
    real(8),dimension(:),intent(in)                  :: eix
    real(8),dimension(size(eix)),intent(in),optional :: eiy
    real(8),dimension(size(eix)),intent(in),optional :: eiz
    ei_x = 0d0
    ei_x(1:size(eix)) = eix
    !
    if(present(eiy))then
       ei_y = 0d0
       ei_y(1:size(eiy)) = eiy
    endif
    !
    if(present(eiz))then
       ei_z = 0d0
       ei_z(1:size(eiz)) = eiz
    endif
    !
    set_eivec=.true.
  end subroutine TB_set_ei


  subroutine TB_set_bk(bkx,bky,bkz)
    real(8),dimension(:),intent(in)                  :: bkx
    real(8),dimension(size(bkx)),intent(in),optional :: bky
    real(8),dimension(size(bkx)),intent(in),optional :: bkz
    bk_x = 0d0
    bk_x(1:size(bkx)) = bkx
    !
    if(present(bky))then
       bk_y = 0d0
       bk_y(1:size(bky)) = bky
    endif
    !
    if(present(bkz))then
       bk_z = 0d0
       bk_z(1:size(bkz)) = bkz
    endif
    !
    set_bkvec=.true.
  end subroutine TB_set_bk





  subroutine TB_get_bk(bkx,bky,bkz)
    real(8),dimension(:),intent(inout)                  :: bkx
    real(8),dimension(size(bkx)),intent(inout),optional :: bky
    real(8),dimension(size(bkx)),intent(inout),optional :: bkz
    real(8),dimension(3)                                :: b1,b2,b3
    !
    call TB_reciprocal_basis(a1=ei_x,a2=ei_y,a3=ei_z,b1=b1,b2=b2,b3=b3)
    !
    bkx = b1(1:size(bkx))
    if(present(bky))bky = b2(1:size(bky))
    if(present(bkz))bkz = b3(1:size(bkz))
    !
  end subroutine TB_get_bk

  subroutine TB_get_ei(eix,eiy,eiz)
    real(8),dimension(:),intent(inout)                  :: eix
    real(8),dimension(size(eix)),intent(inout),optional :: eiy
    real(8),dimension(size(eix)),intent(inout),optional :: eiz
    real(8),dimension(3)                                :: a1,a2,a3
    !
    call TB_reciprocal_basis(a1=bk_x,a2=bk_y,a3=bk_z,b1=a1,b2=a2,b3=a3)
    !
    eix = a1(1:size(eix))
    if(present(eiy))eiy = a2(1:size(eiy))
    if(present(eiz))eiz = a3(1:size(eiz))
    !
  end subroutine TB_get_ei





  subroutine TB_build_bk(verbose)
    logical,optional :: verbose
    logical          :: verbose_
    if(.not.set_eivec)stop "TB_build_bk ERROR: Direct basis not set, set_eivec=F"
    verbose_=.false.;if(present(verbose))verbose_=verbose
    call TB_reciprocal_basis(&
         a1=ei_x,a2=ei_y,a3=ei_z,&
         b1=bk_x,b2=bk_y,b3=bk_z)
    set_bkvec=.true.
    if(verbose_)call print_bk
  end subroutine TB_build_bk

  subroutine TB_build_ei(verbose)
    logical,optional :: verbose
    logical          :: verbose_
    if(.not.set_bkvec)stop "TB_build_ei ERROR: Reciprocal basis not set, set_bkvec=F"
    verbose_=.false.;if(present(verbose))verbose_=verbose
    !note that we exchange the direct basis with the reciprocal
    call TB_reciprocal_basis(&
         b1=ei_x,b2=ei_y,b3=ei_z,&
         a1=bk_x,a2=bk_y,a3=bk_z)
    set_eivec=.true.
    if(verbose_)call print_ei
  end subroutine TB_build_ei





  subroutine print_ei(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    unit=6
    if(present(pfile))open(free_unit(unit),file=reg(pfile))
    write(unit,"(A)")"Using Direct Lattice vectors:"
    write(unit,"(A,3F8.4,A1)")"ei_x = [",(ei_x(i),i=1,3),"]"
    write(unit,"(A,3F8.4,A1)")"ei_y = [",(ei_y(i),i=1,3),"]"
    write(unit,"(A,3F8.4,A1)")"ei_z = [",(ei_z(i),i=1,3),"]"
    if(present(pfile))close(unit)
  end subroutine print_ei

  subroutine print_bk(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    unit=6
    if(present(pfile))open(free_unit(unit),file=reg(pfile))
    write(unit,"(A)")"Using Reciprocal Lattice vectors:"
    write(unit,"(A,3F8.4,A1)")"bk_x = [",(bk_x(i),i=1,3),"]"
    write(unit,"(A,3F8.4,A1)")"bk_y = [",(bk_y(i),i=1,3),"]"
    write(unit,"(A,3F8.4,A1)")"bk_z = [",(bk_z(i),i=1,3),"]"
    if(present(pfile))close(unit)
  end subroutine print_bk





  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Build the reciprocal space {b1[,b2,b3]} given the direct basis {a1[,a2,a3]}
  !  in dimensions up to 3
  !-------------------------------------------------------------------------------------------
  subroutine TB_reciprocal_basis(a1,a2,a3, b1,b2,b3)
    !This routine generates the reciprocal lattice vectors b1,b2,b3
    !given the space vectors a1,a2,a3.
    !the vectors are in units of the lattice constant (a=1) 
    real(8),dimension(:),intent(in)                 :: a1
    real(8),dimension(size(a1)),intent(in),optional :: a2
    real(8),dimension(size(a1)),intent(in),optional :: a3
    real(8),dimension(size(a1))                     :: b1
    real(8),dimension(size(a1)),optional            :: b2
    real(8),dimension(size(a1)),optional            :: b3
    !
    real(8),dimension(3)                            :: ar1,ar2,ar3
    real(8),dimension(3)                            :: bk1,bk2,bk3
    real(8)                                         :: den, s
    integer                                         :: iperm, i, j, k, l, ipol
    integer                                         :: N
    real(8),dimension(3,3)                          :: Mat
    !
    N = size(a1)
    !
    ar1=[1d0,0d0,0d0]
    ar2=[0d0,1d0,0d0]
    ar3=[0d0,0d0,1d0]
    !
    ar1(:N)=a1
    if(present(a2))ar2(:N)=a2
    if(present(a3))ar3(:N)=a3
    !
    den = det(transpose(reshape([ar1,ar2,ar3],shape=[3,3])))
    !
    !    here we compute the reciprocal vectors
    i = 1
    j = 2
    k = 3
    do ipol = 1, 3
       bk1(ipol) = (ar2(j)*ar3(k) - ar2(k)*ar3 (j) )/den*pi2
       bk2(ipol) = (ar3(j)*ar1(k) - ar3(k)*ar1 (j) )/den*pi2
       bk3(ipol) = (ar1(j)*ar2(k) - ar1(k)*ar2 (j) )/den*pi2
       l = i
       i = j
       j = k
       k = l
    enddo
    b1=bk1(:N)
    if(present(b2))b2=bk2(:N)
    if(present(b3))b3=bk3(:N)
    return
  end subroutine TB_reciprocal_basis




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

