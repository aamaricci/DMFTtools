#define FROUTINE 'dmft_get_gloc_realaxis_normal_main'
subroutine dmft_get_gloc_realaxis_normal_main(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  integer,intent(in)                            :: iprint
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
  call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
  call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
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
  write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
  call start_timer
  Greal=zero
  do ik=1,Lk
     call invert_gk_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
     Greal = Greal + Gkreal*Wtk(ik)
     call eta(ik,Lk)
  end do
  call stop_timer
  call dmft_gloc_print_realaxis(wr,Greal,"Gloc",iprint)
end subroutine dmft_get_gloc_realaxis_normal_main
#undef FROUTINE


#define FROUTINE 'dmft_get_gloc_realaxis_normal_dos_main'
subroutine dmft_get_gloc_realaxis_normal_dos_main(Ebands,Dbands,Hloc,Greal,Sreal,iprint)
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                  :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout)               :: Greal     !as Sreal
  integer,intent(in)                                          :: iprint
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
  call assert_shape(Ebands,[Nso,Lk],FROUTINE,"Ebands")
  call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
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
  write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
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
  call dmft_gloc_print_realaxis(wr,Greal,"Gloc",iprint)
end subroutine dmft_get_gloc_realaxis_normal_dos_main
#undef FROUTINE



#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_main'
subroutine dmft_get_gloc_realaxis_normal_lattice_main(Hk,Wtk,Greal,Sreal,iprint,tridiag,hk_symm)
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  integer,intent(in)                              :: iprint    !
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
  call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
  call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Sreal")
  call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],FROUTINE,"Greal")
  !
  write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
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
        call invert_gk_normal_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
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
  call dmft_gloc_print_realaxis_lattice(wr,Greal,"LG",iprint)
end subroutine dmft_get_gloc_realaxis_normal_lattice_main
#undef FROUTINE







!##################################################################
!##################################################################
!##################################################################






#define FROUTINE 'dmft_get_gloc_realaxis_normal_1band'
subroutine dmft_get_gloc_realaxis_normal_1band(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
  complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
  real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
  complex(8),intent(in)                     :: Sreal(:)
  complex(8),intent(inout)                  :: Greal(size(Sreal))
  logical,optional                          :: hk_symm(size(Hk,1))
  logical                                   :: hk_symm_(size(Hk,1))
  integer                                   :: iprint
  !
  complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
  complex(8),dimension(1,1,1,1,size(Sreal)) :: Greal_
  complex(8),dimension(1,1,1,1,size(Sreal)) :: Sreal_
  !
  Hk_(1,1,:)        = Hk(:)
  Sreal_(1,1,1,1,:) = Sreal(:)
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  call dmft_get_gloc_realaxis_normal_main(Hk_,Wtk,Greal_,Sreal_,iprint,hk_symm_)
  Greal(:) = Greal_(1,1,1,1,:)
end subroutine dmft_get_gloc_realaxis_normal_1band
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_1band'
subroutine dmft_get_gloc_realaxis_normal_lattice_1band(Hk,Wtk,Greal,Sreal,iprint,tridiag,hk_symm)
  complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
  real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
  complex(8),dimension(:,:),intent(in)                            :: Sreal           ![Nlat][Lreal]
  complex(8),dimension(size(Sreal,1),size(Sreal,2)),intent(inout) :: Greal
  integer                                                         :: iprint
  logical,optional                                                :: tridiag
  logical                                                         :: tridiag_
  logical,optional                                                :: hk_symm(size(Hk,1))
  logical                                                         :: hk_symm_(size(Hk,1))
  !
  complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Greal_
  complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Sreal_
  !
  tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Sreal,[size(Hk,1),size(Sreal,2)],FROUTINE,"Sreal")
  Sreal_(:,1,1,1,1,:) = Sreal(:,:)
  call dmft_get_gloc_realaxis_normal_lattice_main(Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,hk_symm_)
  Greal(:,:) = Greal_(:,1,1,1,1,:)
end subroutine dmft_get_gloc_realaxis_normal_lattice_1band
#undef FROUTINE









!##################################################################
!##################################################################
!##################################################################





#define FROUTINE 'dmft_get_gloc_realaxis_normal_Nband'
subroutine dmft_get_gloc_realaxis_normal_Nband(Hk,Wtk,Greal,Sreal,iprint,hk_symm)
  complex(8),dimension(:,:,:),intent(in)                        :: Hk              ![Norb][Norb][Nk]
  real(8)                                                       :: Wtk(size(Hk,3)) ![Nk]
  complex(8),intent(in),dimension(:,:,:)                        :: Sreal(:,:,:)    ![Norb][Norb][Lreal]
  complex(8),intent(inout),dimension(:,:,:)                     :: Greal
  logical,optional                                              :: hk_symm(size(Hk,3))
  logical                                                       :: hk_symm_(size(Hk,3))
  integer                                                       :: iprint
  !
  complex(8),&
       dimension(1,1,size(Sreal,1),size(Sreal,2),size(Sreal,3)) :: Greal_
  complex(8),&
       dimension(1,1,size(Sreal,1),size(Sreal,2),size(Sreal,3)) :: Sreal_
  !
  integer                                                       :: Nspin,Norb,Nso,Lreal
  !
  Nspin = 1
  Norb  = size(Sreal,1)
  Lreal = size(Sreal,3)
  Nso   = Nspin*Norb
  call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Sreal,[Nso,Nso,Lreal],FROUTINE,"Sreal")
  call assert_shape(Greal,[Nso,Nso,Lreal],FROUTINE,"Greal")
  !
  Sreal_(1,1,:,:,:) = Sreal(:,:,:)
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  call dmft_get_gloc_realaxis_normal_main(Hk,Wtk,Greal_,Sreal_,iprint,hk_symm_)
  Greal(:,:,:) = Greal_(1,1,:,:,:)
end subroutine dmft_get_gloc_realaxis_normal_Nband
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_Nband'
subroutine dmft_get_gloc_realaxis_normal_lattice_Nband(Hk,Wtk,Greal,Sreal,iprint,tridiag,hk_symm)
  complex(8),dimension(:,:,:),intent(in)                                      :: Hk              ![Nlat*Norb][Nlat*Norb][Nk]
  real(8)                                                                     :: Wtk(size(Hk,3)) ![Nk]
  complex(8),intent(in)                                                       :: Sreal(:,:,:,:)  ![Nlat][Norb][Norb][Lreal]
  complex(8),intent(inout)                                                    :: Greal(:,:,:,:)
  logical,optional                                                            :: hk_symm(size(Hk,3))
  logical                                                                     :: hk_symm_(size(Hk,3))
  integer                                                                     :: iprint
  logical,optional                                                            :: tridiag
  logical                                                                     :: tridiag_
  !
  complex(8),&
       dimension(size(Sreal,1),1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Greal_ ![Nlat][1][1][Norb][Norb][Lreal]
  complex(8),&
       dimension(size(Sreal,1),1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Sreal_![Nlat][1][1][Norb][Norb][Lreal]
  integer                                                                     :: Nlat,Nspin,Norb,Nso,Nlso,Lreal
  !
  tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = 1    
  Nlat  = size(Sreal,1)
  Norb  = size(Sreal,2)
  Lreal = size(Sreal,4)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Sreal,[Nlat,Nso,Nso,Lreal],FROUTINE,"Sreal")
  call assert_shape(Greal,[Nlat,Nso,Nso,Lreal],FROUTINE,"Greal")
  !
  Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
  call dmft_get_gloc_realaxis_normal_lattice_main(Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,hk_symm_)
  Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
end subroutine dmft_get_gloc_realaxis_normal_lattice_Nband
#undef FROUTINE
