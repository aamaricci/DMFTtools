subroutine dmft_get_gloc_matsubara_superc_main(Hk,Wtk,Gmats,Smats,hk_symm)
  complex(8),dimension(:,:,:,:),intent(in)        :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
  real(8),dimension(size(Hk,4)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
  logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nso,Nso,Lk],'dmft_get_gloc_matsubara_superc_main',"Hk")
  call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Smats")
  call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Gmats")
  !
  allocate(Gkmats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(zeta_mats(2,2,Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  do i=1,Lmats
     zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
     zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
     zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
     zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
  enddo
  !
  !invert (Z-Hk) for each k-point
  write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
  call start_timer
  Gmats=zero
  do ik=1,Lk
     call invert_gk_superc(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
     Gmats = Gmats + Gkmats*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_main


subroutine dmft_get_gloc_matsubara_superc_dos(Ebands,Dbands,Hloc,Gmats,Smats)
  real(8),dimension(:,:,:),intent(in)                           :: Ebands    ![2][Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in)   :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(2,size(Ebands,1)),intent(in)                :: Hloc      ![2][Nspin*Norb]
  complex(8),dimension(:,:,:,:,:,:),intent(in)                  :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout)               :: Gmats     !as Smats
  !allocatable arrays
  complex(8)                                                    :: gktmp(2),cdet
  complex(8)                                                    :: zeta_11,zeta_12,zeta_22 
  complex(8),dimension(:,:,:,:,:),allocatable                   :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                                       :: beta
  real(8)                                                       :: xmu,eps
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Ebands,3)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Ebands,[2,Nso,Lk],'dmft_get_gloc_matsubara_superc_dos',"Ebands")
  call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Smats")
  call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Gmats")
  !
  allocate(zeta_mats(2,2,Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  do i=1,Lmats
     zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
     zeta_mats(1,2,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
     zeta_mats(2,1,:,:,i) =                         -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
     zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
  enddo
  !
  !invert (Z-Hk) for each k-point
  write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
  call start_timer
  Gmats=zero
  do i=1,Lmats
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           zeta_11 = zeta_mats(1,1,io,io,i)
           zeta_12 = zeta_mats(1,2,io,io,i)
           zeta_12 = zeta_mats(2,2,io,io,i)
           do ik=1,Lk
              !
              cdet = (zeta_11-Hloc(1,io)-Ebands(1,io,ik))*(zeta_22-Hloc(2,io)-Ebands(2,io,ik)) - zeta_12**2
              gktmp(1)=-(zeta_22-Hloc(2,io)-Ebands(2,io,ik))/cdet
              gktmp(2)=  zeta_12/cdet
              Gmats(1,ispin,ispin,iorb,iorb,i) = Gmats(1,ispin,ispin,iorb,iorb,i) + gktmp(1)*Dbands(io,ik)
              Gmats(2,ispin,ispin,iorb,iorb,i) = Gmats(2,ispin,ispin,iorb,iorb,i) + gktmp(2)*Dbands(io,ik)
           enddo
        enddo
     enddo
     call eta(i,Lmats)
  enddo
  call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_dos



subroutine dmft_get_gloc_matsubara_superc_ineq(Hk,Wtk,Gmats,Smats,hk_symm)
  complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,4)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
  logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                           :: beta
  real(8)                                           :: xmu,eps
  real(8)                                           :: wini,wfin
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nlat  = size(Smats,2)
  Nspin = size(Smats,3)
  Norb  = size(Smats,5)
  Lmats = size(Smats,7)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_superc_ineq_main',"Hk")
  call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_ineq_main',"Smats")
  call assert_shape(Gmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_ineq_main',"Gmats")
  !
  write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
  !
  allocate(Gkmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))    
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  do ilat=1,Nlat
     !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
     !G22(iw) = -[G11[iw]]*
     !G21(iw) =   G12[w]
     do i=1,Lmats
        zeta_mats(1,1,ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) -        nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb)
        zeta_mats(1,2,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
        zeta_mats(2,1,ilat,:,:,i) =                         -        nn2so_reshape(Smats(2,ilat,:,:,:,:,i),Nspin,Norb)
        zeta_mats(2,2,ilat,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + conjg( nn2so_reshape(Smats(1,ilat,:,:,:,:,i),Nspin,Norb) )
     enddo
  enddo
  !
  call start_timer
  Gmats=zero
  do ik=1,Lk
     call invert_gk_superc_ineq(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
     Gmats = Gmats + Gkmats*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_ineq






subroutine dmft_get_gloc_matsubara_superc_gij(Hk,Wtk,Gmats,Fmats,Smats,hk_symm)
  complex(8),dimension(:,:,:,:)                       :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Fmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Smats           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  logical,optional                                  :: hk_symm(size(Hk,4))
  logical                                           :: hk_symm_(size(Hk,4))
  !
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkmats          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  integer                                           :: ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
  integer                                           :: Lk,Nlso,Lmats,Nlat,Nspin,Norb,Nso
  !
  real(8)                                           :: beta
  real(8)                                           :: xmu,eps
  real(8)                                           :: wini,wfin
  !
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nlat  = size(Smats,2)
  Nspin = size(Smats,3)
  Norb  = size(Smats,5)
  Lmats = size(Smats,7)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_superc_gij_main',"Hk")
  call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main',"Smats")
  call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main',"Gmats")
  call assert_shape(Fmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main',"Fmats")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  write(*,"(A)")"Get full Green's function (no print)"
  !
  allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
  allocate(Fkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Fkmats=zero
  allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats));zeta_mats=zero
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  zeta_mats=zero
  !SYMMETRIES in Matsubara-frequencies  [assuming a real order parameter]
  !G22(iw) = -[G11[iw]]*
  !G21(iw) =   G12[w]
  do ilat=1,Nlat
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
           zeta_mats(1,1,ilat,io,io,:) = xi*wm(:) + xmu
           zeta_mats(2,2,ilat,io,io,:) = xi*wm(:) - xmu
        enddo
     enddo
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 zeta_mats(1,1,ilat,io,jo,:) = zeta_mats(1,1,ilat,io,jo,:) - Smats(1,ilat,ispin,jspin,iorb,jorb,:)
                 zeta_mats(1,2,ilat,io,jo,:) = zeta_mats(1,2,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                 zeta_mats(2,1,ilat,io,jo,:) = zeta_mats(2,1,ilat,io,jo,:) - Smats(2,ilat,ispin,jspin,iorb,jorb,:)
                 zeta_mats(2,2,ilat,io,jo,:) = zeta_mats(2,2,ilat,io,jo,:) + conjg( Smats(1,ilat,ispin,jspin,iorb,jorb,:) )
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  call start_timer
  Gmats=zero
  do ik=1,Lk
     call invert_gk_superc_gij(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
     Gmats = Gmats + Gkmats*Wtk(ik)
     Fmats = Fmats + Fkmats*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_gij








