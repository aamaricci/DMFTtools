!> BUILD HK MODEL on a K-GRID (Norb, Norb=1)
!> TODO: find a better way to handle the kgrid, this version works only for cubic ones!
!
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


function build_hk_model_Norb_d(hk_model,Norb,kxgrid,kygrid,kzgrid) result(Hk)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)    :: kpoint
       real(8),dimension(N,N)  :: hk_model
     end function hk_model
  end interface
  integer                                                             :: Norb
  integer                                                             :: Nk,Nkx,Nky,Nkz
  integer                                                             :: ik,ix,iy,iz
  real(8),dimension(:)                                                :: kxgrid,kygrid,kzgrid
  real(8)                                                             :: kx,ky,kz
  real(8),dimension(Norb,Norb,size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(:,:,ik) = hk_model([kx,ky,kz],Norb)
  enddo
end function build_hk_model_Norb_d
!
function build_hk_model_Norb_c(hk_model,Norb,kxgrid,kygrid,kzgrid) result(Hk)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)       :: kpoint
       complex(8),dimension(N,N)  :: hk_model
     end function hk_model
  end interface
  integer                                                                :: Norb
  integer                                                                :: Nk,Nkx,Nky,Nkz
  integer                                                                :: ik,ix,iy,iz
  real(8),dimension(:)                                                   :: kxgrid,kygrid,kzgrid
  real(8)                                                                :: kx,ky,kz
  complex(8),dimension(Norb,Norb,size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(:,:,ik) = hk_model([kx,ky,kz],Norb)
  enddo
end function build_hk_model_Norb_c
!
!
function build_hk_model_1_d(hk_model,kxgrid,kygrid,kzgrid) result(Hk)
  interface 
     function hk_model(kpoint)
       real(8),dimension(:)   :: kpoint
       real(8)                :: hk_model
     end function hk_model
  end interface
  integer                                                      :: Nk,Nkx,Nky,Nkz
  integer                                                      :: ik,ix,iy,iz
  real(8),dimension(:)                                         :: kxgrid,kygrid,kzgrid
  real(8)                                                      :: kx,ky,kz
  real(8),dimension(size(kxgrid)*size(kygrid)*size(kzgrid))    :: hk
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(ik) = hk_model([kx,ky,kz])
  enddo
end function build_hk_model_1_d
!
function build_hk_model_1_c(hk_model,kxgrid,kygrid,kzgrid) result(Hk)
  interface 
     function hk_model(kpoint)
       real(8),dimension(:)      :: kpoint
       complex(8)                :: hk_model
     end function hk_model
  end interface
  integer                                                      :: Nk,Nkx,Nky,Nkz
  integer                                                      :: ik,ix,iy,iz
  real(8),dimension(:)                                         :: kxgrid,kygrid,kzgrid
  real(8)                                                      :: kx,ky,kz
  complex(8),dimension(size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(ik) = hk_model([kx,ky,kz])
  enddo
end function build_hk_model_1_c








!> BUILD HK MODEL on a K-PATH (Norb, Norb=1)
!> TODO: find a better way to handle the k-path
!
function build_hkpath_model_Norb_d(hk_model,Norb,kpath,Nkpath) result(Hk)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)   :: kpoint
       real(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  integer                                                   :: Norb
  real(8),dimension(:,:)                                    :: kpath
  integer                                                   :: Nkpath
  integer                                                   :: Ndim,Npts,Nk
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  real(8),dimension(size(kpath,2))                          :: kvec
  real(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)     :: hk
  Npts=size(kpath,1)          !# of k-points along the path
  Ndim=size(kpath,2)          !# of dimension of the k-vectors
  Nk  = (Npts-1)*Nkpath
  if(Ndim>3)stop "build_hkpath_model_Norb_d error: Ndim > 3"
  do i=1,Ndim
     kgrid(:,i) = kgrid_from_path(kpath,Npts,Nkpath,i)
  enddo
  do ik=1,Nk
     Hk(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
  enddo
end function build_hkpath_model_Norb_d
!
function build_hkpath_model_Norb_c(hk_model,Norb,kpath,Nkpath) result(Hk)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  integer                                                   :: Norb
  real(8),dimension(:,:)                                    :: kpath
  integer                                                   :: Nkpath
  integer                                                   :: Ndim,Npts,Nk
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  real(8),dimension(size(kpath,2))                          :: kvec
  complex(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)  :: hk
  Npts=size(kpath,1)          !# of k-points along the path
  Ndim=size(kpath,2)          !# of dimension of the k-vectors
  Nk  = (Npts-1)*Nkpath
  if(Ndim>3)stop "build_hkpath_model_Norb_c error: Ndim > 3"
  do i=1,Ndim
     kgrid(:,i) = kgrid_from_path(kpath,Npts,Nkpath,i)
  enddo
  do ik=1,Nk
     Hk(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
  enddo
end function build_hkpath_model_Norb_c
!
!
function build_hkpath_model_1_d(hk_model,kpath,Nkpath) result(Hk)
  interface 
     function hk_model(kpoint)
       real(8),dimension(:) :: kpoint
       real(8)              :: hk_model
     end function hk_model
  end interface
  real(8),dimension(:,:)                                    :: kpath
  integer                                                   :: Nkpath
  integer                                                   :: Ndim,Npts,Nk
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  real(8),dimension(size(kpath,2))                          :: kvec
  real(8),dimension((size(kpath,1)-1)*Nkpath)               :: hk
  Npts=size(kpath,1)          !# of k-points along the path
  Ndim=size(kpath,2)          !# of dimension of the k-vectors
  Nk  = (Npts-1)*Nkpath
  if(Ndim>3)stop "build_hk_model_Norb_d error: Ndim > 3"
  do i=1,Ndim
     kgrid(:,i) = kgrid_from_path(kpath,Npts,Nkpath,i)
  enddo
  do ik=1,Nk
     Hk(ik) = hk_model(kgrid(ik,:))
  enddo
end function build_hkpath_model_1_d
!
function build_hkpath_model_1_c(hk_model,kpath,Nkpath) result(Hk)
  interface 
     function hk_model(kpoint)
       real(8),dimension(:) :: kpoint
       complex(8)           :: hk_model
     end function hk_model
  end interface
  real(8),dimension(:,:)                                    :: kpath
  integer                                                   :: Nkpath
  integer                                                   :: Ndim,Npts,Nk
  integer                                                   :: ik,i
  real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
  real(8),dimension(size(kpath,2))                          :: kvec
  complex(8),dimension((size(kpath,1)-1)*Nkpath)            :: hk
  Npts=size(kpath,1)          !# of k-points along the path
  Ndim=size(kpath,2)          !# of dimension of the k-vectors
  Nk  = (Npts-1)*Nkpath
  if(Ndim>3)stop "build_hk_model_Norb_d error: Ndim > 3"
  do i=1,Ndim
     kgrid(:,i) = kgrid_from_path(kpath,Npts,Nkpath,i)
  enddo
  do ik=1,Nk
     Hk(ik) = hk_model(kgrid(ik,:))
  enddo
end function build_hkpath_model_1_c





















