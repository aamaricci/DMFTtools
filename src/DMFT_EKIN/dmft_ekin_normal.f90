!
!MAIN
subroutine dmft_kinetic_energy_normal_main(Hk,Wtk,Sigma,Ekin,Eloc) 
  complex(8),dimension(:,:,:)                 :: Hk    ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3))               :: Wtk   ![Lk]
  complex(8),dimension(:,:,:)                 :: Sigma ![Nspin*Norb][Nspin*Norb][L]
  real(8),dimension(size(Hk,1)),optional      :: Ekin,Eloc
  !
  integer                                     :: Lk,Nso,Liw
  integer                                     :: i,ik,iso
  !
  integer                                     :: Norb,Nporb
  integer                                     :: Nspin  
  real(8)                                     :: beta
  real(8)                                     :: xmu
  !
  real(8),dimension(size(Hk,1),size(Hk,1))    :: Sigma_HF
  complex(8),dimension(size(Hk,1),size(Hk,1)) :: Ak,Bk,Ck,Dk,Hloc
  complex(8),dimension(size(Hk,1),size(Hk,1)) :: Gk,Tk
  !
  real(8),dimension(size(Hk,1))               :: Tail0,Tail1
  real(8),dimension(size(Hk,1))               :: Lail0,Lail1
  real(8)                                     :: spin_degeneracy
  !
  real(8),dimension(size(Hk,1))               :: H0,Hl
  !
  real(8),dimension(size(Hk,1))               :: Ekin_,Eloc_
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
  !if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_main: Nso != Norb*Nspin [from Hk]"
  call assert_shape(Hk,[Nso,Nso,Lk],"dmft_kinetic_energy_normal_main","Hk")
  call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_main","Sigma")
  !
  !Allocate and setup the Matsubara freq.
  if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  wm = pi/beta*dble(2*arange(1,Liw)-1)
  !
  !Get HF part of the self-energy
  Sigma_HF = dreal(Sigma(:,:,Liw))
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
  H0  = 0d0
  Hl  = 0d0
  !Get principal part: Tr[ Hk.(Gk-Tk) ]
  do ik=1,Lk
     Ak = Hk(:,:,ik) - Hloc(:,:)
     Bk =-Hk(:,:,ik) - Sigma_HF(:,:)
     do i=1,Liw
        Gk = (xi*wm(i)+xmu)*eye(Nso) - Sigma(:,:,i) - Hk(:,:,ik) 
        select case(Nso)
        case default
           call inv(Gk)
        case(1)
           Gk = 1d0/Gk
        end select
        Tk = eye(Nso)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
        Ck = matmul(Ak  ,Gk - Tk)
        Dk = matmul(Hloc,Gk - Tk)
        do iso=1,Nso
           H0(iso) = H0(iso) + Wtk(ik)*Ck(iso,iso)
           Hl(iso) = Hl(iso) + Wtk(ik)*Dk(iso,iso)
        enddo
     enddo
     call eta(ik,Lk)
  enddo
  call stop_timer()
  spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
  H0 = H0/beta*2*spin_degeneracy
  Hl = Hl/beta*2*spin_degeneracy
  !
  !get tail subtracted contribution: Tr[ Hk.Tk ]
  Tail0 = 0d0
  Tail1 = 0d0
  Lail0 = 0d0
  Lail1 = 0d0
  do ik=1,Lk
     Ak    = Hk(:,:,ik) - Hloc(:,:)
     Ck= matmul(Ak,(-Hk(:,:,ik)-Sigma_HF(:,:)))
     Dk= matmul(Hloc,(-Hk(:,:,ik)-Sigma_HF(:,:)))
     do iso=1,Nso
        Tail0(iso) = Tail0(iso) + 0.5d0*Wtk(ik)*Ak(iso,iso)
        Tail1(iso) = Tail1(iso) + 0.25d0*Wtk(ik)*Ck(iso,iso)
        Lail0(iso) = Lail0(iso) + 0.5d0*Wtk(ik)*Hloc(iso,iso)
        Lail1(iso) = Lail1(iso) + 0.25d0*Wtk(ik)*Dk(iso,iso)
     enddo
  enddo
  Tail0=Tail0*spin_degeneracy
  Tail1=Tail1*beta*spin_degeneracy
  Lail0=Lail0*spin_degeneracy
  Lail1=Lail1*beta*spin_degeneracy
  !
  Ekin_=H0+Tail0+Tail1
  Eloc_=Hl+Lail0+Lail1
  !
  call write_kinetic_info()
  call write_kinetic_value(Ekin_,Eloc_)
  if(present(Ekin))Ekin=Ekin_
  if(present(Eloc))Eloc=Eloc_
  !
  deallocate(wm)
end subroutine dmft_kinetic_energy_normal_main




!
!DOS
subroutine dmft_kinetic_energy_normal_dos(Ebands,Dbands,Hloc,Sigma,Ekin,Eloc) 
  real(8),dimension(:,:),intent(in)                           :: Ebands  ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands  ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc    ![Nspin*Norb]
  complex(8),dimension(:,:,:)                                 :: Sigma   ![Nspin*Norb][Nspin*Norb][L]
  real(8),dimension(size(Ebands,1)),optional                  :: Ekin,Eloc
  !
  integer                                                     :: Lk,Nso,Liw
  integer                                                     :: i,ik,iso
  !
  integer                                                     :: Norb,Nporb
  integer                                                     :: Nspin  
  real(8)                                                     :: beta
  real(8)                                                     :: xmu
  !
  real(8),dimension(size(Ebands,1),size(Ebands,1))            :: Sigma_HF
  !
  complex(8)                                                  :: Ak,Bk,Ck,Dk
  complex(8)                                                  :: Gk,Tk
  !
  real(8),dimension(size(Ebands,1))                           :: Tail0,Tail1
  real(8),dimension(size(Ebands,1))                           :: Lail0,Lail1
  real(8)                                                     :: spin_degeneracy
  !
  real(8),dimension(size(Ebands,1))                           :: H0,Hl
  !
  real(8),dimension(size(Ebands,1))                           :: Ekin_,Eloc_
  !
  real(8),dimension(:),allocatable                            :: wm
  !
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
  if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_dos: Nso != Norb*Nspin [from Hk]"
  call assert_shape(Sigma,[Nso,Nso,Liw],"dmft_kinetic_energy_normal_dos","Sigma")
  !
  !Allocate and setup the Matsubara freq.
  if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  wm = pi/beta*dble(2*arange(1,Liw)-1)
  !
  !Get HF part of the self-energy
  Sigma_HF = dreal(Sigma(:,:,Liw))
  !
  !
  write(*,"(A)") "Kinetic energy computation"
  call start_timer()
  H0  = 0d0
  Hl  = 0d0
  !Get principal part: Tr[ Hk.(Gk-Tk) ]
  do ik=1,Lk
     do iso=1,Nso
        Ak = Ebands(iso,ik)
        Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso) 
        do i=1,Liw
           Gk = (xi*wm(i)+xmu) - Sigma(iso,iso,i) - Ebands(iso,ik) - Hloc(iso)
           Gk = 1d0/Gk
           Tk = 1d0/(xi*wm(i)) - Bk/(xi*wm(i))**2
           Ck = Ak*(Gk - Tk)
           Dk = Hloc(iso)*(Gk - Tk)
           H0(iso) = H0(iso) + Dbands(iso,ik)*Ck
           Hl(iso) = Hl(iso) + Dbands(iso,ik)*Dk
        enddo
     enddo
     call eta(ik,Lk)
  enddo
  call stop_timer()
  spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
  H0 = H0/beta*2*spin_degeneracy
  Hl = Hl/beta*2*spin_degeneracy
  !
  !get tail subtracted contribution: Tr[ Hk.Tk ]
  Tail0 = 0d0
  Tail1 = 0d0
  Lail0 = 0d0
  Lail1 = 0d0
  do ik=1,Lk
     do iso=1,Nso
        Ak = Ebands(iso,ik)
        Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
        Ck= Ak*Bk
        Dk= Hloc(iso)*Bk
        Tail0(iso) = Tail0(iso) + 0.5d0*Dbands(iso,ik)*Ak
        Tail1(iso) = Tail1(iso) + 0.25d0*Dbands(iso,ik)*Ck
        Lail0(iso) = Lail0(iso) + 0.5d0*Dbands(iso,ik)*Hloc(iso)
        Lail1(iso) = Lail1(iso) + 0.25d0*Dbands(iso,ik)*Dk
     enddo
  enddo
  Tail0=Tail0*spin_degeneracy
  Tail1=Tail1*beta*spin_degeneracy
  Lail0=Lail0*spin_degeneracy
  Lail1=Lail1*beta*spin_degeneracy
  !
  Ekin_=H0+Tail0+Tail1
  Eloc_=Hl+Lail0+Lail1
  !
  call write_kinetic_info()
  call write_kinetic_value(Ekin_,Eloc_)
  if(present(Ekin))Ekin=Ekin_
  if(present(Eloc))Eloc=Eloc_
  !
  deallocate(wm)
end subroutine dmft_kinetic_energy_normal_dos


!
!LATTICE (Nlat independent sites)
subroutine dmft_kinetic_energy_normal_lattice(Hk,Wtk,Sigma,Ekin,Eloc)
  complex(8),dimension(:,:,:)                                     :: Hk        ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3))                                   :: Wtk       ! [Lk]
  complex(8),dimension(:,:,:,:)                                   :: Sigma     ! [Nlat][Nspin*Norb][Nspin*Norb][L]
  real(8),dimension(size(Hk,1)),optional                          :: Ekin,Eloc
  !aux
  integer                                                         :: Lk,Nlso,Liw,Nlat,Nso
  integer                                                         :: ik,iso
  integer                                                         :: i,iorb,ilat,ispin,io,is
  integer                                                         :: j,jorb,jlat,jspin,jo,js
  !
  integer                                                         :: Norb,Nporb
  integer                                                         :: Nspin  
  real(8)                                                         :: beta
  real(8)                                                         :: xmu
  !
  complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Sigma_HF
  complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Ak,Bk,Ck,Dk,Hloc,Hloc_tmp
  complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Tk
  complex(8),dimension(size(Hk,1),size(Hk,1))                     :: Gk
  !
  real(8),dimension(size(Hk,1))                                   :: Tail0,Tail1
  real(8),dimension(size(Hk,1))                                   :: Lail0,Lail1
  real(8)                                                         :: spin_degeneracy
  !
  real(8),dimension(size(Hk,1))                                   :: H0,Hl
  !
  real(8),dimension(size(Hk,1))                                   :: Ekin_,Eloc_
  !
  !
  !Retrieve parameters:
  call get_ctrl_var(Norb,"NORB")
  call get_ctrl_var(Nspin,"NSPIN")
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nlso = size(Hk,1)
  Lk   = size(Hk,3)
  Nlat = size(Sigma,1)
  Nso  = size(Sigma,2)
  Liw  = size(Sigma,4)
  !Testing:
  if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_normal_lattice: Nso != Norb*Nspin [from Sigma]"
  call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_normal_lattice","Hk") !implcitly test that Nlat*Nso=Nlso
  call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_normal_lattice","Sigma")
  !
  !Allocate and setup the Matsubara freq.
  if(allocated(wm))deallocate(wm);allocate(wm(Liw))
  wm = pi/beta*(2*arange(1,Liw)-1)
  !
  !Get HF part of the self-energy
  Sigma_HF(:,:,:) = dreal(Sigma(:,:,:,Liw))
  !
  ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
  Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
  Hloc=0d0
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
  H0           = 0d0
  Hl           = 0d0
  !Get principal part: Tr[ Hk.(Gk-Tk) ]
  do ik=1,Lk
     Ak    =  Hk(:,:,ik) - Hloc(:,:)
     Bk    = -Hk(:,:,ik) - blocks_to_matrix(Sigma_HF(:,:,:),Nlat,Nso) !Sigma_HF [Nlat,Nso,Nso]--> [Nlso,Nslo]
     do i=1,Liw
        Gk = (xi*wm(i)+xmu)*eye(Nlso) - blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso) - Hk(:,:,ik) !Sigma [Nlat,Nso,Nso,*]--> [Nlso,Nslo]
        call inv(Gk(:,:))
        Tk = eye(Nlso)/(xi*wm(i)) - Bk/(xi*wm(i))**2
        Ck = matmul(Ak,Gk(:,:) - Tk)
        Dk = matmul(Hloc,Gk(:,:) - Tk)
        do is=1,Nlso
           H0(is) = H0(is) + Wtk(ik)*Ck(is,is)
           Hl(is) = Hl(is) + Wtk(ik)*Dk(is,is)
        enddo
     enddo
     call eta(ik,Lk)
  enddo
  call stop_timer
  spin_degeneracy=3d0-Nspin !2 if Nspin=1, 1 if Nspin=2
  H0 = H0/beta*2d0*spin_degeneracy
  Hl = Hl/beta*2d0*spin_degeneracy
  !
  !get tail subtracted contribution: Tr[ Hk.Tk ]
  Tail0=0d0
  Tail1=0d0
  Lail0=0d0
  Lail1=0d0
  do ik=1,Lk
     Ak    =  Hk(:,:,ik) - Hloc(:,:)
     Bk    = -Hk(:,:,ik) - blocks_to_matrix(Sigma_HF(:,:,:),Nlat,Nso) !Sigma_HF [Nlat,Nso,Nso]--> [Nlso,Nslo]
     Ck= matmul(Ak,Bk)
     Dk= matmul(Hloc,Bk)
     do is=1,Nlso
        Tail0(is) = Tail0(is) + 0.5d0*Wtk(ik)*Ak(is,is)
        Tail1(is) = Tail1(is) + 0.25d0*Wtk(ik)*Ck(is,is)
        Lail0(is) = Lail0(is) + 0.5d0*Wtk(ik)*Hloc(is,is)
        Lail1(is) = Lail1(is) + 0.25d0*Wtk(ik)*Dk(is,is)
     enddo
  enddo
  Tail0=Tail0*spin_degeneracy
  Tail1=Tail1*beta*spin_degeneracy
  Lail0=Lail0*spin_degeneracy
  Lail1=Lail1*beta*spin_degeneracy
  !
  Ekin_=H0+Tail0+Tail1
  Eloc_=Hl+Lail0+Lail1
  !
  call write_kinetic_info()
  call write_kinetic_value(Ekin_,Eloc_,Nlat,Nso)
  if(present(Ekin))Ekin=Ekin_
  if(present(Eloc))Eloc=Eloc_
  !
  deallocate(wm)
end subroutine dmft_kinetic_energy_normal_lattice
