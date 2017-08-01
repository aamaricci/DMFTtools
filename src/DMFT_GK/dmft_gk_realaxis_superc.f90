subroutine dmft_get_gk_realaxis_superc_main(Hk,Wtk,Gkreal,Sreal)
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb]
  real(8),intent(in)                              :: Wtk       !
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
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
  Nspin = size(Sreal,2)
  Norb  = size(Sreal,4)
  Lreal = size(Sreal,6)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_realaxis_superc_main',"Hk")
  call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main',"Sreal")
  call assert_shape(Gkreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main',"Gkreal")
  !
  allocate(zeta_real(2,2,Nso,Nso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do i=1,Lreal
     zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
          nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
     zeta_real(1,2,:,:,i) =                                                  - &
          nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
     zeta_real(2,1,:,:,i) =                                                  - &
          nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
     zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
          conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
  enddo
  !
  !invert (Z-Hk) for each k-point
  Gkreal=zero
  call invert_gk_superc(zeta_real,Hk,.false.,Gkreal)
end subroutine dmft_get_gk_realaxis_superc_main


subroutine dmft_get_gk_realaxis_superc_dos(Ebands,Dbands,Hloc,Gkreal,Sreal)
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![2][Nspin*Norb]
  real(8),dimension(size(Ebands,2)),intent(in)   :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(2,size(Ebands,2)),intent(in)                :: Hloc      ![2][Nspin*Norb]
  complex(8),dimension(:,:,:,:,:,:),intent(in)                  :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout)               :: Gkreal     !as Sreal
  !allocatable arrays
  complex(8)                                                    :: gktmp(2),cdet
  complex(8)                                                    :: zeta_11,zeta_12,zeta_22 
  complex(8),dimension(:,:,:,:,:),allocatable                   :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                                       :: beta
  real(8)                                                       :: xmu,eps
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Sreal,2)
  Norb  = size(Sreal,4)
  Lreal = size(Sreal,6)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Ebands,[2,Nso],'dmft_get_gk_realaxis_superc_dos',"Ebands")
  call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main',"Sreal")
  call assert_shape(Gkreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main',"Gkreal")
  !
  allocate(zeta_real(2,2,Nso,Nso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do i=1,Lreal
     zeta_real(1,1,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
          nn2so_reshape(Sreal(1,:,:,:,:,i),Nspin,Norb)
     zeta_real(1,2,:,:,i) =                                                  - &
          nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
     zeta_real(2,1,:,:,i) =                                                  - &
          nn2so_reshape(Sreal(2,:,:,:,:,i),Nspin,Norb)
     zeta_real(2,2,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
          conjg( nn2so_reshape(Sreal(1,:,:,:,:,Lreal+1-i),Nspin,Norb) )
  enddo
  !
  !invert (Z-Hk) for each k-point
  Gkreal=zero
  do i=1,Lreal
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           zeta_11 = zeta_real(1,1,io,io,i)
           zeta_12 = zeta_real(1,2,io,io,i)
           zeta_12 = zeta_real(2,2,io,io,i)
           !
           cdet = (zeta_11-Hloc(1,io)-Ebands(1,io))*(zeta_22-Hloc(2,io)-Ebands(2,io)) - zeta_12**2
           gktmp(1)=-(zeta_22-Hloc(2,io)-Ebands(2,io))/cdet
           gktmp(2)=  zeta_12/cdet
           Gkreal(1,ispin,ispin,iorb,iorb,i) = Gkreal(1,ispin,ispin,iorb,iorb,i) + gktmp(1)*Dbands(io)
           Gkreal(2,ispin,ispin,iorb,iorb,i) = Gkreal(2,ispin,ispin,iorb,iorb,i) + gktmp(2)*Dbands(io)
        enddo
     enddo
  enddo
end subroutine dmft_get_gk_realaxis_superc_dos



subroutine dmft_get_gk_realaxis_superc_ineq(Hk,Wtk,Gkreal,Sreal)
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                           :: beta
  real(8)                                           :: xmu,eps
  real(8)                                           :: wini,wfin  
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  Nlat  = size(Sreal,2)
  Nspin = size(Sreal,3)
  Norb  = size(Sreal,5)
  Lreal = size(Sreal,7)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_realaxis_superc_ineq',"Hk")
  call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_ineq_main',"Sreal")
  call assert_shape(Gkreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_ineq_main',"Gkreal")
  !
  allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal))
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  do ilat=1,Nlat
     !SYMMETRIES in real-frequencies   [assuming a real order parameter]
     !G22(w)  = -[G11[-w]]*
     !G21(w)  =   G12[w]   
     do i=1,Lreal
        zeta_real(1,1,ilat,:,:,i) = (dcmplx(wr(i),eps)+ xmu)*eye(Nso)                - &
             nn2so_reshape(Sreal(1,ilat,:,:,:,:,i),Nspin,Norb)
        zeta_real(1,2,ilat,:,:,i) =                                                  - &
             nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
        zeta_real(2,1,ilat,:,:,i) =                                                  - &
             nn2so_reshape(Sreal(2,ilat,:,:,:,:,i),Nspin,Norb)
        zeta_real(2,2,ilat,:,:,i) = -conjg( dcmplx(wr(Lreal+1-i),eps)+xmu )*eye(Nso) + &
             conjg( nn2so_reshape(Sreal(1,ilat,:,:,:,:,Lreal+1-i),Nspin,Norb) )
     enddo
  enddo
  !
  Gkreal=zero
  call invert_gk_superc_ineq(zeta_real,Hk,.false.,Gkreal)
end subroutine dmft_get_gk_realaxis_superc_ineq

