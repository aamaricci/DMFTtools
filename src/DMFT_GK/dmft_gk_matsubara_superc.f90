subroutine dmft_get_gk_matsubara_superc_main(Hk,Wtk,Gkmats,Smats)
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
  real(8),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_matsubara_superc_main',"Hk")
  call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Smats")
  call assert_shape(Gkmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Gkmats")
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
  Gkmats=zero
  call invert_gk_superc(zeta_mats,Hk,.false.,Gkmats)
end subroutine dmft_get_gk_matsubara_superc_main


subroutine dmft_get_gk_matsubara_superc_dos(Ebands,Dbands,Hloc,Gkmats,Smats)
  real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
  real(8),dimension(size(Ebands,2)),intent(in)    :: Dbands    ![Nspin*Norb]
  real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
  !allocatable arrays
  complex(8)                                      :: gktmp(2),cdet
  complex(8)                                      :: zeta_11,zeta_12,zeta_22 
  complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Ebands,[2,Nso],'dmft_get_gk_matsubara_superc_dos',"Ebands")
  call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Smats")
  call assert_shape(Gkmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_main',"Gkmats")
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
  Gkmats=zero
  do i=1,Lmats
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           zeta_11 = zeta_mats(1,1,io,io,i)
           zeta_12 = zeta_mats(1,2,io,io,i)
           zeta_12 = zeta_mats(2,2,io,io,i)
           !
           cdet = (zeta_11-Hloc(1,io)-Ebands(1,io))*(zeta_22-Hloc(2,io)-Ebands(2,io)) - zeta_12**2
           gktmp(1)=-(zeta_22-Hloc(2,io)-Ebands(2,io))/cdet
           gktmp(2)=  zeta_12/cdet
           Gkmats(1,ispin,ispin,iorb,iorb,i) = gktmp(1)*Dbands(io)
           Gkmats(2,ispin,ispin,iorb,iorb,i) = gktmp(2)*Dbands(io)
        enddo
     enddo
  enddo
end subroutine dmft_get_gk_matsubara_superc_dos



subroutine dmft_get_gk_matsubara_superc_ineq(Hk,Wtk,Gkmats,Smats)
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                           :: beta
  real(8)                                           :: xmu,eps
  real(8)                                           :: wini,wfin
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nlat  = size(Smats,2)
  Nspin = size(Smats,3)
  Norb  = size(Smats,5)
  Lmats = size(Smats,7)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_matsubara_superc_ineq_main',"Hk")
  call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_ineq_main',"Smats")
  call assert_shape(Gkmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_superc_ineq_main',"Gkmats")
  !  
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
  Gkmats=zero
  call invert_gk_superc_ineq(zeta_mats,Hk,.false.,Gkmats)
end subroutine dmft_get_gk_matsubara_superc_ineq








