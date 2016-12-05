module DMFT_TIGHT_BINDING
  USE SF_CONSTANTS, only: pi,pi2,xi,one,zero
  USE SF_IOTOOLS, only:free_unit,reg,txtfy,free_units
  USE SF_LINALG, only: eigh,eigh_jacobi
  USE SF_COLORS
  implicit none
  private


  interface TB_build_model
     module procedure build_hk_model_Norb_d
     module procedure build_hk_model_Norb_c
     module procedure build_hk_model_1_d
     module procedure build_hk_model_1_c
     module procedure build_hkpath_model_Norb_d
     module procedure build_hkpath_model_Norb_c
     module procedure build_hkpath_model_1_d 
     module procedure build_hkpath_model_1_c
     module procedure build_hkr_model_Norb_d
     module procedure build_hkr_model_1_d
     module procedure build_hkr_model_Norb_c
     module procedure build_hkr_model_1_c
  end interface TB_build_model
  public :: TB_build_model  



  interface TB_solve_path
     module procedure solve_Hk_along_BZpath
     module procedure solve_HkR_along_BZpath
  end interface TB_solve_path
  public :: TB_solve_path



  interface write_hk_w90
     module procedure write_hk_w90_1
     module procedure write_hk_w90_2
     module procedure write_hk_w90_3
     module procedure write_hk_w90_4
  end interface write_hk_w90
  public :: write_hk_w90

  interface read_hk_w90
     module procedure read_hk_w90_1
     module procedure read_hk_w90_2
  end interface read_hk_w90
  public :: read_hk_w90




  interface write_Hloc
     module procedure write_Hloc_1
     module procedure write_Hloc_2
  end interface write_Hloc
  public :: write_hloc


  interface read_Hloc
     module procedure read_Hloc_1
     module procedure read_Hloc_2
  end interface read_Hloc
  public :: read_hloc




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


contains



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  build the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) Hamiltonian matrix
  ! from the function user defined hk_model procedure.
  ! > multi-orbital/lattice sites Norb: real,complex
  ! > single orbital/lattice site: real,complex
  ! > multi-orbital/lattice sites Norb: real,complex
  ! > single orbital/lattice site: real,complex
  !-------------------------------------------------------------------------------------------
  include "tight_binding_build_hk_model.f90"
  include "tight_binding_build_hkR_model.f90"




  !-------------------------------------------------------------------------------------------
  !PURPOSE:  solve the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) along a given linear 
  ! path in the Brillouin Zone. A GNUPLOT script to plot the bands together with their
  ! character is generated.
  !-------------------------------------------------------------------------------------------
  include "tight_binding_solve_hk_path.f90"
  include "tight_binding_solve_hkR_path.f90"



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  construct a grid of k-points with minimum information
  !-------------------------------------------------------------------------------------------
  function kgrid(Nk,start,len)
    integer               :: Nk,i
    real(8),optional      :: start,len
    real(8)               :: start_,len_
    real(8),dimension(Nk) :: kgrid
    start_=-pi ;if(present(start))start_=start
    len_=2d0*pi;if(present(len))len_=len
    do i=1,Nk
       kgrid(i) = start_ + len_*(i-1)/Nk
    enddo
  end function kgrid
  !
  function kgrid_from_path(kpath,Npts,Nk,dim) result(kxpath)
    real(8),dimension(Npts,3)      :: kpath
    real(8),dimension((Npts-1)*Nk) :: kxpath
    real(8),dimension(3)           :: kstart,kstop,kpoint,kdiff
    integer                        :: ipts,ik,ic,dim,Nk,Npts
    if(dim>3)stop "error kigrid_from_path: dim > 3"
    ic=0
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nk)
       do ik=1,Nk
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          kxpath(ic)=kpoint(dim)
       enddo
    enddo
  end function kgrid_from_path








  !-------------------------------------------------------------------------------------------
  !PURPOSE:  obtain the coordinates ix,iy,iz from the lattice index ik
  !-------------------------------------------------------------------------------------------
  function indx2ix(ik,ndim) result(ix)
    integer              :: ik
    integer              :: ix
    integer              :: nx_,ny_,nz_
    integer,dimension(3) :: ndim
    nx_=ndim(1)
    ny_=ndim(2)
    nz_=ndim(3)
    ix=int(ik-1)/ny_/nz_+1
  end function indx2ix
  !
  function indx2iy(ik,ndim) result(iy)
    integer              :: ik
    integer              :: iy
    integer              :: nx_,ny_,nz_
    integer,dimension(3) :: ndim
    nx_=ndim(1)
    ny_=ndim(2)
    nz_=ndim(3)
    iy=mod(int(ik-1)/nz_,ny_)+1
  end function indx2iy
  !
  function indx2iz(ik,ndim) result(iz)
    integer              :: ik
    integer              :: iz
    integer              :: nx_,ny_,nz_
    integer,dimension(3) :: ndim
    nx_=ndim(1)
    ny_=ndim(2)
    nz_=ndim(3)
    iz=mod(ik-1,nz_)+1
  end function indx2iz
  !
  subroutine coord2indx(ik,ix,iy,iz,ndim)
    integer              :: ix,iy,iz
    integer,dimension(3) :: ndim
    integer              :: nx,ny,nz
    integer              :: ik
    nx=ndim(1)
    ny=ndim(2)
    nz=ndim(3)
    ik = (ix-1)*ny*nz+(iy-1)*nz+iz
  end subroutine coord2indx
  !
  subroutine indx2coord(ik,ix,iy,iz,ndim)
    integer              :: ix,iy,iz
    integer,dimension(3) :: ndim
    integer              :: ik
    ix=indx2ix(ik,ndim)
    iy=indx2iy(ik,ndim)
    iz=indx2iz(ik,ndim)
  end subroutine indx2coord








  !-------------------------------------------------------------------------------------------
  !PURPOSE:  write the Hamiltonian matrix H(k)
  !-------------------------------------------------------------------------------------------
  subroutine write_hk_w90_1(file,No,Nd,Np,Nineq,hk_model,kxgrid,kygrid,kzgrid)
    character(len=*)               :: file
    integer                        :: No,Nd,Np,Nineq
    integer                        :: Nktot,Nkx,Nky,Nkz,unit
    integer                        :: ik,ix,iy,iz,iorb,jorb
    real(8),dimension(:)           :: kxgrid,kygrid,kzgrid
    real(8)                        :: kx,ky,kz
    complex(8),dimension(No,No)    :: hk
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    Nkx = size(kxgrid)
    Nky = size(kygrid)
    Nkz = size(kzgrid)
    Nktot = Nkx*Nky*Nkz
    unit=free_unit()
    open(unit,file=reg(file))
    write(unit,'(1A1,1A12,1x,3(A2,1x))')"#",reg(txtfy(Nktot)),reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
    do ik=1,Nktot
       call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
       kx=kxgrid(ix)
       ky=kygrid(iy)
       kz=kzgrid(iz)
       Hk(:,:) = hk_model([kx,ky,kz],No)
       write(unit,"(3(F15.9,1x))")kx,ky,kz
       do iorb=1,No
          write(unit,"(20(2F15.9,1x))")(Hk(iorb,jorb),jorb=1,No)
       enddo
    enddo
    close(unit)
  end subroutine write_hk_w90_1
  !
  subroutine write_hk_w90_2(file,No,Nd,Np,Nineq,hk_model,kxgrid,kygrid,kzgrid)
    character(len=*)               :: file
    integer                        :: No,Nd,Np,Nineq
    integer                        :: Nktot,Nkx,Nky,Nkz,unit
    integer                        :: ik,ix,iy,iz,iorb,jorb
    real(8),dimension(:)           :: kxgrid,kygrid,kzgrid
    real(8)                        :: kx,ky,kz
    complex(8)                     :: hk
    interface 
       function hk_model(kpoint)
         real(8),dimension(:)      :: kpoint
         complex(8)                :: hk_model
       end function hk_model
    end interface
    Nkx = size(kxgrid)
    Nky = size(kygrid)
    Nkz = size(kzgrid)
    Nktot = Nkx*Nky*Nkz
    unit=free_unit()
    open(unit,file=reg(file))
    write(unit,'(1A1,1A12,1x,3(A2,1x))')"#",reg(txtfy(Nktot)),reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
    do ik=1,Nktot
       call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
       kx=kxgrid(ix)
       ky=kygrid(iy)
       kz=kzgrid(iz)
       Hk = hk_model([kx,ky,kz])
       write(unit,"(3(F15.9,1x))")kx,ky,kz
       write(unit,"(20(2F15.9,1x))")Hk
    enddo
  end subroutine write_hk_w90_2
  !
  subroutine write_hk_w90_3(file,No,Nd,Np,Nineq,hk,kxgrid,kygrid,kzgrid)
    character(len=*)                :: file
    integer                         :: No,Nd,Np,Nineq
    integer                         :: Nktot,Nkx,Nky,Nkz,unit
    integer                         :: ik,ix,iy,iz,i,j
    real(8),dimension(:)            :: kxgrid,kygrid,kzgrid
    real(8)                         :: kx,ky,kz
    complex(8),dimension(:,:,:)     :: Hk    
    unit=free_unit()
    Nkx = size(kxgrid)
    Nky = size(kygrid)
    Nkz = size(kzgrid)
    Nktot = Nkx*Nky*Nkz
    if(size(Hk,1)/=No)stop "write_hk_f90: error in dimension Hk,1"
    if(size(Hk,2)/=No)stop "write_hk_f90: error in dimension Hk,2"
    if(size(Hk,3)/=Nktot)stop "write_hk_f90: error in dimension Hk,3"
    open(unit,file=reg(file))
    write(unit,'(1A1,1A12,1x,3(A2,1x))')"#",reg(txtfy(Nktot)),reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
    do ik=1,Nktot
       call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
       kx=kxgrid(ix)
       ky=kygrid(iy)
       kz=kzgrid(iz)
       write(unit,"(3(F15.9,1x))")kx,ky,kz
       do i=1,No
          write(unit,"(20(2F15.9,1x))")(Hk(i,j,ik),j=1,No)
       enddo
    enddo
  end subroutine write_hk_w90_3
  !
  subroutine write_hk_w90_4(file,No,Nd,Np,Nineq,hk,kxgrid,kygrid,kzgrid)
    character(len=*)        :: file
    integer                 :: No,Nd,Np,Nineq
    integer                 :: Nktot,Nkx,Nky,Nkz,unit
    integer                 :: ik,ix,iy,iz,iorb,jorb
    real(8),dimension(:)    :: kxgrid,kygrid,kzgrid
    real(8)                 :: kx,ky,kz
    complex(8),dimension(:) :: Hk    
    unit=free_unit()
    Nkx = size(kxgrid)
    Nky = size(kygrid)
    Nkz = size(kzgrid)
    Nktot = Nkx*Nky*Nkz
    if(size(Hk,1)/=Nktot)stop "write_hk_f90: error in dimension Hk,1"
    open(unit,file=reg(file))
    write(unit,'(1A1,1A12,1x,3(A2,1x))')"#",reg(txtfy(Nktot)),reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
    do ik=1,Nktot
       call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
       kx=kxgrid(ix)
       ky=kygrid(iy)
       kz=kzgrid(iz)
       write(unit,"(3(F15.9,1x))")kx,ky,kz
       write(unit,"((2F15.9,1x))")Hk(ik)
    enddo
    close(unit)
  end subroutine write_hk_w90_4










  !-------------------------------------------------------------------------------------------
  !PURPOSE:  read the Hamiltonian matrix H(k)
  !-------------------------------------------------------------------------------------------
  subroutine read_hk_w90_1(file,No,Nd,Np,Nineq,hk,kxgrid,kygrid,kzgrid)
    character(len=*)            :: file
    integer                     :: No,Nd,Np,Nineq
    integer                     :: Nktot,Nk_,Nkx,Nky,Nkz,unit
    integer                     :: ik,ix,iy,iz,iorb,jorb
    real(8),dimension(:)        :: kxgrid,kygrid,kzgrid
    real(8)                     :: kx,ky,kz
    logical                     :: ioexist
    complex(8),dimension(:,:,:) :: Hk
    character(len=1)            :: achar
    inquire(file=reg(file),exist=ioexist)
    if(.not.ioexist)then
       write(*,*)"can not find file:"//reg(file)
       stop
    endif
    Nkx = size(kxgrid)
    Nky = size(kygrid)
    Nkz = size(kzgrid)
    Nktot = Nkx*Nky*Nkz
    unit=free_unit()
    open(unit,file=reg(file))
    read(unit,'(1A1,I12,1x,3(I2,1x))')achar,Nk_,Nd,Np,Nineq
    if(Nk_/=Nktot)stop "read_hk_f90: error in number of k-points: check kgrids and the header"
    if(size(Hk,3)/=Nktot)stop "read_hk_f90: error in size[Hk,3]"
    if(size(Hk,2)/=No)stop "read_hk_f90: error in size[Hk,2]"
    if(size(Hk,1)/=No)stop "read_hk_f90: error in size[Hk,1]"
    do ik=1,Nktot
       call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
       read(unit,"(3(F15.9,1x))")kx,ky,kz
       kxgrid(ix)=kx
       kygrid(iy)=ky
       kzgrid(iz)=kz
       do iorb=1,No
          read(unit,"(20(2F15.9,1x))")(Hk(iorb,jorb,ik),jorb=1,No)
       enddo
    enddo
    close(unit)
  end subroutine read_hk_w90_1
  !
  subroutine read_hk_w90_2(file,No,Nd,Np,Nineq,hk,kxgrid,kygrid,kzgrid)
    character(len=*)            :: file
    integer                     :: No,Nd,Np,Nineq
    integer                     :: Nktot,Nk_,Nkx,Nky,Nkz,unit
    integer                     :: ik,ix,iy,iz,iorb,jorb
    real(8),dimension(:)        :: kxgrid,kygrid,kzgrid
    real(8)                     :: kx,ky,kz
    logical                     :: ioexist
    complex(8),dimension(:)     :: Hk
    character(len=1)            :: achar
    inquire(file=reg(file),exist=ioexist)
    if(.not.ioexist)then
       write(*,*)"can not find file:"//reg(file)
       stop
    endif
    Nkx = size(kxgrid)
    Nky = size(kygrid)
    Nkz = size(kzgrid)
    Nktot = Nkx*Nky*Nkz
    unit=free_unit()
    open(unit,file=reg(file))
    read(unit,'(1A1,I12,1x,3(I2,1x))')achar,Nk_,Nd,Np,Nineq
    if(Nk_/=Nktot)stop "read_hk_f90: error in number of k-points: check kgrids and the header"
    if(size(Hk,1)/=Nktot)stop "read_hk_f90: error in dimension Hk,1"
    do ik=1,Nktot
       call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
       read(unit,"(3(F15.9,1x))")kx,ky,kz
       kxgrid(ix)=kx
       kygrid(iy)=ky
       kzgrid(iz)=kz
       read(unit,"(20(2F15.9,1x))")Hk(ik)
    enddo
    close(unit)
  end subroutine read_hk_w90_2











  !-------------------------------------------------------------------------------------------
  !PURPOSE:  write/read the local part of the Hamiltonian to a file
  !-------------------------------------------------------------------------------------------
  subroutine write_hloc_1(hloc,file) ![Nlso][Nlso]
    complex(8),dimension(:,:) :: Hloc
    character(len=*),optional :: file
    integer                   :: iorb,jorb,Ni,Nj,unit
    unit=6;
    if(present(file))then
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Ni=size(Hloc,1)
    Nj=size(Hloc,2)
    if(present(file))then
       do iorb=1,Ni
          write(unit,"(90F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       do iorb=1,Ni
          write(unit,"(90F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       close(unit)
    else
       do iorb=1,Ni
          write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
               ('(',dreal(Hloc(iorb,jorb)),',',dimag(Hloc(iorb,jorb)),')',jorb =1,Nj)
       enddo
       write(unit,*)""
    endif
  end subroutine write_hloc_1

  subroutine write_hloc_2(hloc,file) ![Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:,:,:) :: Hloc
    character(len=*),optional     :: file
    integer                       :: iorb,jorb,ispin,jspin
    integer                       :: Norb,Nspin,unit
    unit=6;
    if(present(file))then
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Nspin=size(Hloc,1)
    if(size(Hloc,2)/=Nspin)stop "write_hloc error: size[Hloc,2] != size[Hloc,1] == Nspin"
    Norb=size(Hloc,3)
    if(size(Hloc,4)/=Nspin)stop "write_hloc error: size[Hloc,4] != size[Hloc,3] == Norb"
    if(present(file))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       close(unit)
    else
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
                  (&
                  (&
                  '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
                  jorb =1,Norb),&
                  jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    endif
  end subroutine write_hloc_2


  subroutine read_hloc_1(hloc,file)
    complex(8),dimension(:,:)                    :: Hloc
    character(len=*)                             :: file
    integer                                      :: iorb,jorb,Ni,Nj,unit
    real(8),dimension(size(Hloc,1),size(Hloc,2)) :: reHloc,imHloc
    unit=free_unit()   
    open(unit,file=reg(file))
    Ni=size(Hloc,1)
    Nj=size(Hloc,2)
    do iorb=1,Ni
       read(unit,"(90F12.6)")(reHloc(iorb,jorb),jorb=1,Nj)
    enddo
    write(unit,*)""
    do iorb=1,Ni
       read(unit,"(90F12.6)")(imHloc(iorb,jorb),jorb=1,Nj)
    enddo
    close(unit)
    Hloc = dcmplx(reHloc,imHloc)
  end subroutine read_hloc_1

  subroutine read_hloc_2(hloc,file)
    complex(8),dimension(:,:,:,:)                                          :: Hloc
    character(len=*)                                                       :: file
    integer                                                                :: iorb,jorb,ispin,jspin,unit
    integer                                                                :: Nspin,Norb
    real(8),dimension(size(Hloc,1),size(Hloc,2),size(Hloc,3),size(Hloc,4)) :: reHloc,imHloc
    unit=free_unit()   
    open(unit,file=reg(file))
    Nspin=size(Hloc,1)
    if(size(Hloc,2)/=Nspin)stop "read_hloc error: size[Hloc,2] != size[Hloc,1] == Nspin"
    Norb=size(Hloc,3)
    if(size(Hloc,4)/=Nspin)stop "read_hloc error: size[Hloc,4] != size[Hloc,3] == Norb"
    do ispin=1,Nspin
       do iorb=1,Norb
          read(unit,"(1000F12.6)")((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    read(unit,*)
    do ispin=1,Nspin
       do iorb=1,Norb
          read(unit,"(1000F12.6)")((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    close(unit)
    Hloc = dcmplx(reHloc,imHloc)
  end subroutine read_hloc_2




  !-------------------------------------------------------------------------------------------
  !PURPOSE: reduce or expand the Hamiltonian in different format (to be expanded).
  !-------------------------------------------------------------------------------------------
  function shrink_Hkr(Hkr,Nlat,Nk,Norb) result(Hkr_)
    complex(8),dimension(Nlat,Nlat,Nk,Norb,Norb) :: Hkr
    integer                                      :: Norb,Nlat,Nk
    integer                                      :: io,jo,ilat,iorb,jorb,jlat,ik
    complex(8),dimension(Nlat*Norb,Nlat*Norb,Nk) :: Hkr_
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                io=iorb + (ilat-1)*Norb
                jo=jorb + (jlat-1)*Norb
                do ik=1,Nk
                   Hkr_(io,jo,ik)=Hkr(ilat,jlat,ik,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function shrink_Hkr
  !
  function expand_Hkr(Hkr_,Nlat,Nk,Norb) result(Hkr)
    complex(8),dimension(Nlat,Nlat,Nk,Norb,Norb) :: Hkr
    integer                                      :: Norb,Nlat,Nk
    integer                                      :: io,jo,ilat,iorb,jorb,jlat,ik
    complex(8),dimension(Nlat*Norb,Nlat*Norb,Nk) :: Hkr_
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                io=iorb + (ilat-1)*Norb
                jo=jorb + (jlat-1)*Norb
                do ik=1,Nk
                   Hkr(ilat,jlat,ik,iorb,jorb)=Hkr_(io,jo,ik)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function expand_Hkr



END MODULE DMFT_TIGHT_BINDING

