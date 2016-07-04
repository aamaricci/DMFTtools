function build_hkr_model_Norb_d(hkr_model,Nlat,Norb,kxgrid,kygrid,kzgrid,pbc) result(Hk)
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                   :: kpoint
       integer                                :: Nlat,Norb
       logical                                :: pbc
       real(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
     end function hkr_model
  end interface
  integer                                                                       :: Nlat,Norb
  integer                                                                       :: Nk,Nkx,Nky,Nkz
  integer                                                                       :: ik,ix,iy,iz
  real(8),dimension(:)                                                          :: kxgrid,kygrid,kzgrid
  real(8)                                                                       :: kx,ky,kz
  real(8),dimension(Nlat*Norb,Nlat*Norb,size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  logical                                                                       :: pbc
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(:,:,ik) = hkr_model([kx,ky,kz],Nlat,Norb,pbc)
  enddo
end function build_hkr_model_Norb_d

function build_hkr_model_Norb_c(hkr_model,Nlat,Norb,kxgrid,kygrid,kzgrid,pbc) result(Hk)
  interface 
     function hkr_model(kpoint,Nlat,Norb,pbc)
       real(8),dimension(:)                      :: kpoint
       integer                                   :: Nlat,Norb
       logical                                   :: pbc
       complex(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
     end function hkr_model
  end interface
  integer                                                                          :: Nlat,Norb
  integer                                                                          :: Nk,Nkx,Nky,Nkz
  integer                                                                          :: ik,ix,iy,iz
  real(8),dimension(:)                                                             :: kxgrid,kygrid,kzgrid
  real(8)                                                                          :: kx,ky,kz
  complex(8),dimension(Nlat*Norb,Nlat*Norb,size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  logical                                                                          :: pbc
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(:,:,ik) = hkr_model([kx,ky,kz],Nlat,Norb,pbc)
  enddo
end function build_hkr_model_Norb_c
!
!
function build_hkr_model_1_d(hkr_model,Nlat,kxgrid,kygrid,kzgrid,pbc) result(Hk)
  interface 
     function hkr_model(kpoint,Nlat,pbc)
       real(8),dimension(:)         :: kpoint
       integer                      :: Nlat
       real(8),dimension(Nlat,Nlat) :: hkr_model
       logical                      :: pbc
     end function hkr_model
  end interface
  integer                                                             :: Nlat,Nk,Nkx,Nky,Nkz
  integer                                                             :: ik,ix,iy,iz
  real(8),dimension(:)                                                :: kxgrid,kygrid,kzgrid
  real(8)                                                             :: kx,ky,kz
  real(8),dimension(Nlat,Nlat,size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  logical                                                             :: pbc
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(:,:,ik) = hkr_model([kx,ky,kz],Nlat,pbc)
  enddo
end function build_hkr_model_1_d
!
function build_hkr_model_1_c(hkr_model,Nlat,kxgrid,kygrid,kzgrid,pbc) result(Hk)
  interface 
     function hkr_model(kpoint,Nlat,pbc)
       real(8),dimension(:)            :: kpoint
       integer                         :: Nlat
       complex(8),dimension(Nlat,Nlat) :: hkr_model
       logical                         :: pbc
     end function hkr_model
  end interface
  integer                                                                :: Nlat,Nk,Nkx,Nky,Nkz
  integer                                                                :: ik,ix,iy,iz
  real(8),dimension(:)                                                   :: kxgrid,kygrid,kzgrid
  real(8)                                                                :: kx,ky,kz
  complex(8),dimension(Nlat,Nlat,size(kxgrid)*size(kygrid)*size(kzgrid)) :: hk
  logical                                                                :: pbc
  Nkx = size(kxgrid)
  Nky = size(kygrid)
  Nkz = size(kzgrid)
  Nk  = Nkx*Nky*Nkz
  do ik=1,Nk
     call indx2coord(ik,ix,iy,iz,[Nkx,Nky,Nkz])
     kx=kxgrid(ix)
     ky=kygrid(iy)
     kz=kzgrid(iz)
     Hk(:,:,ik) = hkr_model([kx,ky,kz],Nlat,pbc)
  enddo
end function build_hkr_model_1_c
