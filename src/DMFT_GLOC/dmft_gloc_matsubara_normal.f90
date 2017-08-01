subroutine dmft_get_gloc_matsubara_normal_main(Hk,Wtk,Gmats,Smats,hk_symm)
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
  logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gkmats    !as Smats
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
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Hk,[Nso,Nso,Lk],'dmft_get_gloc_matsubara_normal_main',"Hk")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_main',"Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_main',"Gmats")
  !
  !Allocate and setup the Matsubara freq.
  allocate(Gkmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(zeta_mats(Nso,Nso,Lmats))
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !
  do i=1,Lmats
     zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  write(*,"(A)")"Get local Matsubara Green's function (no print)"
  call start_timer
  Gmats=zero
  do ik=1,Lk
     call invert_gk_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
     Gmats = Gmats + Gkmats*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_matsubara_normal_main





subroutine dmft_get_gloc_matsubara_normal_dos(Ebands,Dbands,Hloc,Gmats,Smats)
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                  :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout)               :: Gmats     !as Smats
  !
  complex(8)                                                  :: gktmp
  complex(8),dimension(:,:,:),allocatable                     :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                                     :: beta
  real(8)                                                     :: xmu,eps
  real(8)                                                     :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Lk    = size(Ebands,2)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_matsubara_normal_dos_main',"Ebands")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_dos_main',"Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_dos_main',"Gmats")
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
  write(*,"(A)")"Get local Matsubara Green's function (no print)"
  call start_timer
  Gmats=zero
  do i=1,Lmats
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           do ik=1,Lk
              gktmp = Dbands(io,ik)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io,ik) )
              Gmats(ispin,ispin,iorb,iorb,i) = Gmats(ispin,ispin,iorb,iorb,i) + gktmp
           enddo
        enddo
     enddo
     call eta(i,Lmats)
  end do
  call stop_timer
end subroutine dmft_get_gloc_matsubara_normal_dos






subroutine dmft_get_gloc_matsubara_normal_ineq(Hk,Wtk,Gmats,Smats,tridiag,hk_symm)
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  logical,optional                                :: tridiag
  logical                                         :: tridiag_
  logical,dimension(size(Hk,3)),optional          :: hk_symm
  logical,dimension((size(Hk,3)))                 :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
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
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nlat  = size(Smats,1)
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_normal_ineq_main',"Hk")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_ineq_main',"Smats")
  call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_ineq_main',"Gmats")
  !
  write(*,"(A)")"Get local Matsubara Green's function (no print)"
  if(.not.tridiag_)then
     write(*,"(A)")"Direct Inversion algorithm:"
  else
     write(*,"(A)")"Quantum Zipper algorithm:"
  endif
  !
  allocate(Gkmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
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
  call start_timer
  Gmats=zero
  if(.not.tridiag_)then
     do ik=1,Lk
        call invert_gk_normal_ineq(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
        Gmats = Gmats + Gkmats*Wtk(ik)
        call eta(ik,Lk)
     end do
  else
     do ik=1,Lk
        call invert_gk_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
        Gmats = Gmats + Gkmats*Wtk(ik)
        call eta(ik,Lk)
     end do
  endif
  call stop_timer
end subroutine dmft_get_gloc_matsubara_normal_ineq


subroutine dmft_get_gloc_matsubara_normal_gij(Hk,Wtk,Gmats,Smats,hk_symm)
  complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats     !      [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  logical,dimension(size(Hk,3)),optional            :: hk_symm
  logical,dimension((size(Hk,3)))                   :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:),allocatable         :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  real(8)                                         :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nlat  = size(Smats,1)
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_normal_gij_main',"Hk")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_gij_main',"Smats")
  call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_gij_main',"Gmats")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  write(*,"(A)")"Get full Green's function (no print)"
  !
  allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
  allocate(zeta_mats(Nlat,Nso,Nso,Lmats));zeta_mats=zero
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
  call start_timer
  Gmats=zero
  do ik=1,Lk
     call invert_gk_normal_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
     Gmats = Gmats + Gkmats*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_matsubara_normal_gij
