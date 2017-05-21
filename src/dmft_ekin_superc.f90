!
!MAIN
subroutine dmft_kinetic_energy_superc_main(Hk,Wtk,Sigma,Self,Ekin,Eloc)
  complex(8),dimension(:,:,:)                     :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3))                   :: Wtk   ![Lk]
  complex(8),dimension(:,:,:)                     :: Sigma ![Nspin*Norb][Nspin*Norb][L]
  complex(8),dimension(:,:,:)                     :: Self ![Nspin*Norb][Nspin*Norb][L]
  real(8),dimension(size(Hk,1)),optional          :: Ekin,Eloc
  !
  integer                                         :: Lk,Nso,Liw
  integer                                         :: i,is,ik
  !
  integer                                         :: Norb,Nporb
  integer                                         :: Nspin  
  real(8)                                         :: beta
  real(8)                                         :: xmu
  !
  real(8),dimension(size(Hk,1),size(Hk,1))        :: Sigma_HF,Self_HF
  complex(8),dimension(size(Hk,1),size(Hk,1))     :: Ak,Bk,Ck,Hloc
  complex(8),dimension(size(Hk,1),size(Hk,1))     :: Gk,Tk
  complex(8),dimension(2*size(Hk,1),2*size(Hk,1)) :: GkNambu,HkNambu,HlocNambu,AkNambu
  real(8),dimension(2*size(Hk,1))                 :: NkNambu
  complex(8),dimension(2*size(Hk,1),2*size(Hk,1)) :: Evec
  real(8),dimension(2*size(Hk,1))                 :: Eval,Coef
  real(8)                                         :: spin_degeneracy
  !
  real(8),dimension(size(Hk,1))                   :: H0,Hl
  real(8),dimension(size(Hk,1))                   :: H0free,Hlfree
  real(8),dimension(size(Hk,1))                   :: Ekin_,Eloc_
  !
  !Retrieve parameters:
  call get_ctrl_var(Norb,"NORB")
  call get_ctrl_var(Nspin,"NSPIN")
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nso = size(Hk,1)
  Lk  = size(Hk,3)
  Liw = size(Sigma,3)
  !Testing:
  if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_main: Nso != Norb*Nspin [from Hk]"
  call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_superc_main","Hk")
  call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_main","Sigma")
  call assert_shape(Self,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_main","Self")
  !
  if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  wm = pi/beta*dble(2*arange(1,Liw)-1)
  !
  !Get HF part of the self-energy
  Sigma_HF = dreal(Sigma(:,:,Liw))
  Self_HF  = dreal(Self(:,:,Liw))
  !
  !Get the local Hamiltonian
  Hloc = sum(Hk(:,:,:),dim=3)/Lk
  where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
  if(size(Hloc,1)<16)then
     call print_hloc(Hloc)
  else
     call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
  endif
  !
  write(*,"(A)") "Kinetic energy computation"
  call start_timer()
  H0=0d0
  Hl=0d0
  !Get principal part: Tr[ Hk.(Gk-Tk) ]
  do ik=1,Lk
     Ak= Hk(:,:,ik)-Hloc
     do i=1,Liw
        Gknambu=zero
        Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)&
             - Sigma(:,:,i)        - Hk(:,:,ik)
        Gknambu(1:Nso,Nso+1:2*Nso)       = - Self(:,:,i)
        Gknambu(Nso+1:2*Nso,1:Nso)       = - Self(:,:,i)
        Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)&
             + conjg(Sigma(:,:,i)) + Hk(:,:,ik)
        call inv(Gknambu)
        Gk  = Gknambu(1:Nso,1:Nso) !Gk(iw)
        !
        Gknambu=zero
        Gknambu(1:Nso,1:Nso)             = (xi*wm(i) + xmu)*eye(Nso)&
             - Sigma_hf  - Hk(:,:,ik)
        Gknambu(1:Nso,Nso+1:2*Nso)       = - Self_hf
        Gknambu(Nso+1:2*Nso,1:Nso)       = - Self_hf
        Gknambu(Nso+1:2*Nso,Nso+1:2*Nso) = (xi*wm(i) - xmu)*eye(Nso)&
             + Sigma_hf  + Hk(:,:,ik)
        call inv(Gknambu)
        Tk = Gknambu(1:Nso,1:Nso) !G0k(iw)
        !
        Bk = matmul(Ak, Gk - Tk) !
        Ck = matmul(Hloc,Gk - Tk)
        do is=1,Nso
           H0(is) = H0(is) + Wtk(ik)*Bk(is,is)
           Hl(is) = Hl(is) + Wtk(ik)*Ck(is,is)
        enddo
     enddo
     call eta(ik,Lk)
  enddo
  call stop_timer
  spin_degeneracy=3d0-dble(Nspin) !2 if Nspin=1, 1 if Nspin=2
  H0=H0/beta*2.d0*spin_degeneracy
  Hl=Hl/beta*2.d0*spin_degeneracy
  !
  !
  H0free=0d0
  Hlfree=0d0
  HkNambu=zero
  HlocNambu=zero
  do ik=1,Lk
     HkNambu(1:Nso,1:Nso)               =  Hk(:,:,ik)-Hloc
     HkNambu(Nso+1:2*Nso,Nso+1:2*Nso)   = -Hk(:,:,ik)+Hloc
     HlocNambu(1:Nso,1:Nso)             =  Hloc
     HlocNambu(Nso+1:2*Nso,Nso+1:2*Nso) = -Hloc
     Evec(1:Nso,1:Nso)                 =  Hk(:,:,ik) + Sigma_hf
     Evec(1:Nso,Nso+1:2*Nso)           =             + Self_hf
     Evec(Nso+1:2*Nso,1:Nso)           =             + Self_hf
     Evec(Nso+1:2*Nso,Nso+1:2*Nso)     = -Hk(:,:,ik) - Sigma_hf
     call eigh(Evec,Eval)
     do is=1,2*Nso
        Coef = Evec(:,is)*conjg(Evec(:,is))
        NkNambu(is) = dot_product(Coef,fermi(Eval,beta))
     enddo
     AkNambu = matmul( diag(NkNambu) , HkNambu )
     do is=1,Nso
        H0free(is) = H0free(is) + Wtk(ik)*AkNambu(is,is) !Take only the 11 part
     enddo
     AkNambu = matmul( diag(NkNambu) , HlocNambu )
     do is=1,Nso
        Hlfree(is) = Hlfree(is) + Wtk(ik)*AkNambu(is,is)
     enddo
  enddo
  H0free=spin_degeneracy*H0free
  Hlfree=spin_degeneracy*Hlfree
  !
  Ekin_=H0+H0free
  Eloc_=Hl+Hlfree
  !
  call write_kinetic_info()
  call write_kinetic_value(Ekin_,Eloc_)
  if(present(Ekin))Ekin=Ekin_
  if(present(Eloc))Eloc=Eloc_
  !
  deallocate(wm)
end subroutine dmft_kinetic_energy_superc_main





!
!DOS
subroutine dmft_kinetic_energy_superc_dos(Ebands,Dbands,Hloc,Sigma,Self,Ekin,Eloc)
  real(8),dimension(:,:),intent(in)                           :: Ebands  ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands  ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc    ![Nspin*Norb]
  complex(8),dimension(:,:,:)                                 :: Sigma   ![Nspin*Norb][Nspin*Norb][L]
  complex(8),dimension(:,:,:)                                 :: Self    ![Nspin*Norb][Nspin*Norb][L]
  real(8),dimension(size(Ebands,1)),optional                  :: Ekin,Eloc
  !
  integer                                                     :: Lk,Nso,Liw
  integer                                                     :: i,is,ik
  !
  integer                                                     :: Norb,Nporb
  integer                                                     :: Nspin  
  real(8)                                                     :: beta
  real(8)                                                     :: xmu
  !
  real(8),dimension(size(Ebands,1),size(Ebands,1))            :: Sigma_HF,Self_HF
  complex(8)                                                  :: Ak,Bk,Ck
  complex(8)                                                  :: Gk,Tk
  !
  complex(8),dimension(2,2)                                   :: GkNambu
  complex(8),dimension(2,2)                                   :: HkNambu,HlocNambu,AkNambu
  real(8),dimension(2)                                        :: NkNambu
  complex(8),dimension(2,2)                                   :: Evec
  real(8),dimension(2)                                        :: Eval,Coef
  real(8)                                                     :: spin_degeneracy
  !
  real(8),dimension(size(Ebands,1))                           :: H0,Hl
  real(8),dimension(size(Ebands,1))                           :: H0free,Hlfree
  real(8),dimension(size(Ebands,1))                           :: Ekin_,Eloc_
  !
  real(8),dimension(:),allocatable                            :: wm
  !
  !Retrieve parameters:
  call get_ctrl_var(Norb,"NORB")
  call get_ctrl_var(Nspin,"NSPIN")
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nso = size(Ebands,1)
  Lk  = size(Ebands,2)
  Liw = size(Sigma,3)
  !Testing:
  if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_dos: Nso != Norb*Nspin [from Hk]"
  call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_dos","Sigma")
  call assert_shape(Self,[Nso,Nso,Liw],"dmft_kinetic_energy_superc_dos","Self")
  !
  if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  wm = pi/beta*dble(2*arange(1,Liw)-1)
  !
  !Get HF part of the self-energy
  Sigma_HF = dreal(Sigma(:,:,Liw))
  Self_HF  = dreal(Self(:,:,Liw))
  !
  !
  write(*,"(A)") "Kinetic energy computation"
  call start_timer()
  H0=0d0
  Hl=0d0
  !Get principal part: Tr[ Hk.(Gk-Tk) ]
  do ik=1,Lk
     do is=1,Nso
        Ak = Ebands(is,ik)
        do i=1,Liw
           Gknambu=zero
           Gknambu(1,1) = (xi*wm(i) + xmu)  - Sigma(is,is,i)        - Ebands(is,ik)
           Gknambu(1,2) =                   - Self(is,is,i)
           Gknambu(2,1) =                   - Self(is,is,i)
           Gknambu(2,2) = (xi*wm(i) - xmu)  + conjg(Sigma(is,is,i)) + Ebands(is,ik)
           call inv(Gknambu)
           Gk  = Gknambu(1,1) !Gk(iw)
           !
           Gknambu=zero
           Gknambu(1,1) = (xi*wm(i) + xmu)  - Sigma_hf(is,is)  - Ebands(is,ik)
           Gknambu(1,2) =                   - Self_hf(is,is)
           Gknambu(2,1) =                   - Self_hf(is,is)
           Gknambu(2,2) = (xi*wm(i) - xmu)  + Sigma_hf(is,is)  + Ebands(is,ik)
           call inv(Gknambu)
           Tk = Gknambu(1,1) !G0k(iw)
           !
           Bk = Ak*(Gk - Tk)
           Ck = Hloc(is)*(Gk - Tk)
           !
           H0(is) = H0(is) + Dbands(is,ik)*Bk
           Hl(is) = Hl(is) + Dbands(is,ik)*Ck
        enddo
     enddo
     call eta(ik,Lk)
  enddo
  call stop_timer
  spin_degeneracy=3d0-dble(Nspin) !2 if Nspin=1, 1 if Nspin=2
  H0=H0/beta*2.d0*spin_degeneracy
  Hl=Hl/beta*2.d0*spin_degeneracy
  !
  !
  H0free=0d0
  Hlfree=0d0
  do ik=1,Lk
     do is=1,Nso
        HkNambu=zero
        HlocNambu=zero
        HkNambu(1,1)   =  Ebands(is,ik)
        HkNambu(2,2)   = -Ebands(is,ik)
        HlocNambu(1,1) =  Hloc(is)
        HlocNambu(2,2) = -Hloc(is)
        Evec(1,1)      =  Ebands(is,ik) + Sigma_hf(is,is)
        Evec(1,2)      =                 + Self_hf(is,is)
        Evec(2,1)      =                 + Self_hf(is,is)
        Evec(2,2)      = -EBands(is,ik) - Sigma_hf(is,is)
        !
        call eigh(Evec,Eval)
        Coef = Evec(:,is)*conjg(Evec(:,is))
        NkNambu(is) = dot_product(Coef,fermi(Eval,beta))
        !
        AkNambu = matmul( diag(NkNambu) , HkNambu )
        H0free(is) = H0free(is) + Dbands(is,ik)*AkNambu(is,is) !Take only the 11 part
        !
        AkNambu = matmul( diag(NkNambu) , HlocNambu )
        Hlfree(is) = Hlfree(is) + Dbands(is,ik)*AkNambu(is,is)
        !
     enddo
  enddo
  H0free=spin_degeneracy*H0free
  Hlfree=spin_degeneracy*Hlfree
  !
  Ekin_=H0+H0free
  Eloc_=Hl+Hlfree
  !
  call write_kinetic_info()
  call write_kinetic_value(Ekin_,Eloc_)
  if(present(Ekin))Ekin=Ekin_
  if(present(Eloc))Eloc=Eloc_
  !
  deallocate(wm)
end subroutine dmft_kinetic_energy_superc_dos



!
!LATTICE (Nlat independent sites)
subroutine dmft_kinetic_energy_superc_lattice(Hk,Wtk,Sigma,Self,Ekin,Eloc)
  complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3))                                   :: Wtk       ! [Lk]
  complex(8),dimension(:,:,:,:)                                   :: Sigma  ! [Nlat][Nspin*Norb][Nspin*Norb][L]
  complex(8),dimension(:,:,:,:)                                   :: Self   ! [Nlat][Nspin*Norb][Nspin*Norb][L]
  real(8),dimension(size(Hk,1)),optional                          :: Ekin,Eloc
  !
  integer                                                         :: Lk,Nlso,Liw,Nso,Nlat
  integer                                                         :: ik
  integer                                                         :: i,iorb,ilat,ispin,io,is
  integer                                                         :: j,jorb,jlat,jspin,jo,js
  !
  integer                                                         :: Norb,Nporb
  integer                                                         :: Nspin  
  real(8)                                                         :: beta
  real(8)                                                         :: xmu
  !
  complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Sigma_HF ![Nlat][Nso][Nso]
  complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Self_HF  ![Nlat][Nso][Nso]
  complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Ak,Bk,Ck,Hloc,Hloc_tmp
  complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Gk,Tk
  complex(8),dimension(2*size(Hk,1),2*size(Hk,2))                 :: Gknambu,HkNambu,HlocNambu,AkNambu
  real(8),dimension(2*size(Hk,1))                                 :: NkNambu
  complex(8),dimension(2*size(Hk,1),2*size(Hk,2))                 :: Evec
  real(8),dimension(2*size(Hk,1))                                 :: Eval,Coef
  real(8)                                                         :: spin_degeneracy
  !
  real(8),dimension(size(Hk,1))                                   :: H0,Hl
  real(8),dimension(size(Hk,1))                                   :: H0free,Hlfree
  real(8),dimension(size(Hk,1))                                   :: Ekin_,Eloc_
  !
  !Retrieve parameters:
  call get_ctrl_var(Norb,"NORB")
  call get_ctrl_var(Nspin,"NSPIN")
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !Get generalized Lattice-Spin-Orbital index
  Nlso = size(Hk,1)
  Lk   = size(Hk,3)
  Nlat = size(Sigma,1)
  Nso  = size(Sigma,2)
  Liw  = size(Sigma,4)
  !Testing:
  if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_superc_lattice: Nso != Norb*Nspin [from Sigma]"
  call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_superc_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
  call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_lattice","Sigma")
  call assert_shape(Self,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_superc_lattice","Self")
  !
  !Allocate and setup the Matsubara freq.
  if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  wm = pi/beta*(2*arange(1,Liw)-1)
  !
  !Get HF part of the self-energy
  Sigma_HF = dreal(Sigma(:,:,:,Liw))![Nlat,Nso,Nso]
  Self_HF  = dreal(Self(:,:,:,Liw)) ![Nlat,Nso,Nso]
  !
  ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
  Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
  Hloc=zero
  do ilat=1,Nlat
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                 js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                 Hloc(is,js)=Hloc_tmp(is,js) 
              enddo
           enddo
        enddo
     enddo
  enddo
  where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
  if(size(Hloc,1)<16)then
     call print_hloc(Hloc)
  else
     call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
  endif
  !
  !Start the timer:
  write(*,"(A)") "Kinetic energy computation"
  call start_timer
  H0 = 0d0
  Hl = 0d0
  !Get principal part: Tr[ Hk.(Gk-Tk) ]
  do ik=1,Lk
     Ak    = Hk(:,:,ik) - Hloc(:,:)
     do i=1,Liw
        Gknambu=zero
        Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - &
             blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso) - Hk(:,:,ik)
        Gknambu(1:Nlso,Nlso+1:2*Nlso)        = - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
        Gknambu(Nlso+1:2*Nlso,1:Nlso)        = - blocks_to_matrix(Self(:,:,:,i),Nlat,Nso)
        Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + &
             conjg(blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso)) + Hk(:,:,ik)
        call inv(Gknambu(:,:))
        Gk = Gknambu(1:Nlso,1:Nlso)
        !
        Gknambu=zero
        Gknambu(1:Nlso,1:Nlso)               = (xi*wm(i) + xmu)*eye(Nlso) - &
             blocks_to_matrix(Sigma_HF,Nlat,Nso)  - Hk(:,:,ik)
        Gknambu(1:Nlso,Nlso+1:2*Nlso)        = - blocks_to_matrix(Self_HF,Nlat,Nso)
        Gknambu(Nlso+1:2*Nlso,1:Nlso)        = - blocks_to_matrix(Self_HF,Nlat,Nso)
        Gknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = (xi*wm(i) - xmu)*eye(Nlso) + &
             blocks_to_matrix(Sigma_HF,Nlat,Nso)  + Hk(:,:,ik)
        call inv(Gknambu(:,:))
        Tk = Gknambu(1:Nlso,1:Nlso)
        !
        Bk = matmul(Ak, Gk - Tk)
        Ck = matmul(Hloc, Gk - Tk)
        do is=1,Nlso
           H0(is) = H0(is) + Wtk(ik)*Bk(is,is)
           Hl(is) = Hl(is) + Wtk(ik)*Ck(is,is)
        enddo
     enddo
     call eta(ik,Lk)
  enddo
  call stop_timer
  spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
  H0 = H0/beta*2d0*spin_degeneracy!;print*,"Ekin_=",H0/Nlat
  Hl = Hl/beta*2d0*spin_degeneracy!;print*,"Eloc_=",Hl/Nlat
  !
  !
  !get tail subtracted contribution: Tr[ Hk.Tk ]
  H0free=0d0
  Hlfree=0d0
  HkNambu=zero
  HlocNambu=zero
  do ik=1,Lk
     HkNambu(1:Nlso,1:Nlso)                 =  Hk(:,:,ik)-Hloc
     HkNambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso)   = -Hk(:,:,ik)+Hloc
     HlocNambu(1:Nlso,1:Nlso)               =  Hloc
     HlocNambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = -Hloc
     Evec(1:Nlso,1:Nlso)                    =  Hk(:,:,ik) +  blocks_to_matrix(Sigma_HF,Nlat,Nso)
     Evec(1:Nlso,Nlso+1:2*Nlso)             =             +  blocks_to_matrix(Self_HF,Nlat,Nso)
     Evec(Nlso+1:2*Nlso,1:Nlso)             =             +  blocks_to_matrix(Self_HF,Nlat,Nso)
     Evec(Nlso+1:2*Nlso,Nlso+1:2*Nlso)      = -Hk(:,:,ik) -  blocks_to_matrix(Sigma_HF,Nlat,Nso)
     call eigh(Evec,Eval)
     NkNambu = fermi(Eval,beta)
     GkNambu = matmul(Evec,matmul(diag(NkNambu),conjg(transpose(Evec))))
     AkNambu = matmul(HkNambu  , GkNambu)
     do is=1,Nlso
        H0free(is) = H0free(is) + Wtk(ik)*AkNambu(is,is) !Take only the 11 part
     enddo
     ! H0free  = H0free + Wtk(ik)*trace_matrix( AkNambu(:Nlso,:Nlso) , Nlso) !take only the 11 part
     AkNambu = matmul(HlocNambu, GkNambu)
     do is=1,Nlso
        Hlfree(is) = Hlfree(is) + Wtk(ik)*AkNambu(is,is)
     enddo
     ! Hlfree  = Hlfree + Wtk(ik)*trace_matrix( AkNambu(:Nlso,:Nlso) , Nlso)
  enddo
  H0free=spin_degeneracy*H0free!;print*,"Efree=",H0free/Nlat
  Hlfree=spin_degeneracy*Hlfree!;print*,"Efree_loc=",Hlfree/Nlat
  !
  Ekin_=H0+H0free
  Eloc_=Hl+Hlfree
  !
  call write_kinetic_info()
  call write_kinetic_value(Ekin_,Eloc_,Nlat)
  if(present(Ekin))Ekin=Ekin_
  if(present(Eloc))Eloc=Eloc_
  !
  deallocate(wm)
end subroutine dmft_kinetic_energy_superc_lattice


