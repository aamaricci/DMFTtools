subroutine dmft_get_gloc_realaxis_normal_main(Hk,Wtk,Greal,Sreal,hk_symm)
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
  logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                       :: beta
  real(8)                                       :: xmu,eps
  real(8)                                       :: wini,wfin
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = size(Sreal,1)
  Norb  = size(Sreal,3)
  Lreal = size(Sreal,5)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Hk,[Nso,Nso,Lk],'dmft_get_gloc_realaxis_normal_main',"Hk")
  call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_main',"Sreal")
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_main',"Greal")
  !
  !Allocate and setup the Realaxis freq.
  allocate(Gkreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(zeta_real(Nso,Nso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do i=1,Lreal
     zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  write(*,"(A)")"Get local Realaxis Green's function (no print)"
  call start_timer
  Greal=zero
  do ik=1,Lk
     call invert_gk_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
     Greal = Greal + Gkreal*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_realaxis_normal_main


subroutine dmft_get_gloc_realaxis_normal_cluster(Hk,Wtk,Greal,Sreal,hk_symm)
  complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  logical,dimension(size(Hk,3)),optional            :: hk_symm   ![Lk]
  logical,dimension(size(Hk,3))                     :: hk_symm_  ![Lk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:),allocatable           :: zeta_real ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lreal]
  !
  real(8)                                           :: beta
  real(8)                                           :: xmu,eps
  real(8)                                           :: wini,wfin
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nlat  = size(Sreal,1)
  Nspin = size(Sreal,3)
  Norb  = size(Sreal,5)
  Lreal = size(Sreal,7)
  Lk    = size(Hk,3)
  Nlso  = Nlat*Nspin*Norb    
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_cluster',"Hk")
  call assert_shape(Sreal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_cluster',"Sreal")
  call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_cluster',"Greal")
  !
  !Allocate and setup the Realaxis freq.
  allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(zeta_real(Nlso,Nlso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do i=1,Lreal
     zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nlso) - nnn2lso_cluster_reshape(Sreal(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  write(*,"(A)")"Get local Realaxis Green's function (no print)"
  call start_timer
  Greal=zero
  do ik=1,Lk
     call invert_gk_normal_cluster(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
     Greal = Greal + Gkreal*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_realaxis_normal_cluster


subroutine dmft_get_gloc_realaxis_normal_dos(Ebands,Dbands,Hloc,Greal,Sreal)
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                  :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout)               :: Greal     !as Sreal
  !
  complex(8)                                                  :: gktmp
  complex(8),dimension(:,:,:),allocatable                     :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                                     :: xmu,eps
  real(8)                                                     :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  Nspin = size(Sreal,1)
  Norb  = size(Sreal,3)
  Lreal = size(Sreal,5)
  Lk    = size(Ebands,2)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_realaxis_normal_dos_main',"Ebands")
  call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos_main',"Sreal")
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos_main',"Greal")
  !
  !Allocate and setup the Realaxis freq.
  allocate(zeta_real(Nso,Nso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do i=1,Lreal
     zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
  enddo
  !
  !invert (Z-Hk) for each k-point
  write(*,"(A)")"Get local Realaxis Green's function (no print)"
  call start_timer
  Greal=zero
  do i=1,Lreal
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           do ik=1,Lk
              gktmp = Dbands(io,ik)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io,ik) )
              Greal(ispin,ispin,iorb,iorb,i) = Greal(ispin,ispin,iorb,iorb,i) + gktmp
           enddo
        enddo
     enddo
     call eta(i,Lreal)
  end do
  call stop_timer
end subroutine dmft_get_gloc_realaxis_normal_dos


subroutine dmft_get_gloc_realaxis_normal_ineq(Hk,Wtk,Greal,Sreal,tridiag,hk_symm)
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  logical,optional                                :: tridiag
  logical                                         :: tridiag_
  logical,dimension(size(Hk,3)),optional          :: hk_symm
  logical,dimension((size(Hk,3)))                 :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  real(8)                                         :: wini,wfin
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nlat  = size(Sreal,1)
  Nspin = size(Sreal,2)
  Norb  = size(Sreal,4)
  Lreal = size(Sreal,6)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_ineq_main',"Hk")
  call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_ineq_main',"Sreal")
  call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_ineq_main',"Greal")
  !
  write(*,"(A)")"Get local Realaxis Green's function (no print)"
  if(.not.tridiag_)then
     write(*,"(A)")"Direct Inversion algorithm:"
  else
     write(*,"(A)")"Quantum Zipper algorithm:"
  endif
  !
  allocate(Gkreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(zeta_real(Nlat,Nso,Nso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do ilat=1,Nlat
     do i=1,Lreal
        zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
     enddo
  enddo
  !
  !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
  call start_timer
  Greal=zero
  if(.not.tridiag_)then
     do ik=1,Lk
        call invert_gk_normal_ineq(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
        Greal = Greal + Gkreal*Wtk(ik)
        call eta(ik,Lk)
     end do
  else
     do ik=1,Lk
        call invert_gk_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
        Greal = Greal + Gkreal*Wtk(ik)
        call eta(ik,Lk)
     end do
  endif
  call stop_timer
end subroutine dmft_get_gloc_realaxis_normal_ineq


subroutine dmft_get_gloc_realaxis_normal_gij(Hk,Wtk,Greal,Sreal,hk_symm)
  complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  logical,dimension(size(Hk,3)),optional            :: hk_symm
  logical,dimension((size(Hk,3)))                   :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  real(8)                                         :: wini,wfin
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  Nlat  = size(Sreal,1)
  Nspin = size(Sreal,2)
  Norb  = size(Sreal,4)
  Lreal = size(Sreal,6)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_gij_main',"Hk")
  call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_gij_main',"Sreal")
  call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_gij_main',"Greal")
  !
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  write(*,"(A)")"Get full Green's function (no print)"
  !
  allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
  allocate(zeta_real(Nlat,Nso,Nso,Lreal));zeta_real=zero
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do ilat=1,Nlat
     do i=1,Lreal
        zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:,i),NSpin,Norb)
     enddo
  enddo
  !
  !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
  call start_timer
  Greal=zero
  do ik=1,Lk
     call invert_gk_normal_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
     Greal = Greal + Gkreal*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
end subroutine dmft_get_gloc_realaxis_normal_gij
