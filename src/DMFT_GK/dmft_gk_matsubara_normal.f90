subroutine dmft_get_gk_matsubara_normal_main(Hk,Wtk,Gkmats,Smats)
  complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb]
  real(8),intent(in)                            :: Wtk       !
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
  !allocatable arrays
  complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                       :: beta
  real(8)                                       :: xmu,eps
  real(8)                                       :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Hk,[Nso,Nso],'dmft_get_gk_matsubara_normal_main',"Hk")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_main',"Smats")
  call assert_shape(Gkmats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_main',"Gkmats")
  !
  !Allocate and setup the Matsubara freq.
  allocate(zeta_mats(Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !
  do i=1,Lmats
     zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  Gkmats=zero
  call invert_gk_normal(zeta_mats,Hk,.false.,Gkmats) 
end subroutine dmft_get_gk_matsubara_normal_main





subroutine dmft_get_gk_matsubara_normal_dos(Ebands,Dbands,Hloc,Gkmats,Smats)
  real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb]
  real(8),dimension(size(Ebands)),intent(in)    :: Dbands    ![Nspin*Norb]
  real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
  !
  complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                       :: beta
  real(8)                                       :: xmu,eps
  real(8)                                       :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Nso   = Nspin*Norb    
  !Testing part:
  if(size(Ebands)/=Nso)stop "dmft_get_gk_matsubara_normal_dos_main ERROR: size(Ebands)!=Nso"
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_dos_main',"Smats")
  call assert_shape(Gkmats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_dos_main',"Gkmats")
  !
  !Allocate and setup the Matsubara freq.
  allocate(zeta_mats(Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !
  do i=1,Lmats
     zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  Gkmats=zero
  do i=1,Lmats
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           Gkmats(ispin,ispin,iorb,iorb,i) = Dbands(io)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io) )
        enddo
     enddo
  end do
end subroutine dmft_get_gk_matsubara_normal_dos






subroutine dmft_get_gk_matsubara_normal_ineq(Hk,Wtk,Gkmats,Smats,tridiag)
  complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
  real(8),intent(in)                              :: Wtk       !
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkmats     !as Smats
  logical,optional                                :: tridiag
  logical                                         :: tridiag_
  !allocatable arrays
  complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  real(8)                                         :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
  !
  Nlat  = size(Smats,1)
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_matsubara_normal_ineq_main',"Hk")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main',"Smats")
  call assert_shape(Gkmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gk_matsubara_normal_ineq_main',"Gkmats")
  !
  allocate(zeta_mats(Nlat,Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  do ilat=1,Nlat
     do i=1,Lmats
        zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - nn2so_reshape(Smats(ilat,:,:,:,:,i),Nspin,Norb)
     enddo
  enddo
  !
  !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
  Gkmats=zero
  if(.not.tridiag_)then
     call invert_gk_normal_ineq(zeta_mats,Hk,.false.,Gkmats)
  else
     call invert_gk_normal_tridiag(zeta_mats,Hk,.false.,Gkmats)
  endif
end subroutine dmft_get_gk_matsubara_normal_ineq

