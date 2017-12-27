subroutine dmft_get_gloc_matsubara_superc_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,mpi_split,hk_symm)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
  real(8),dimension(size(Hk,4)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  character(len=*),optional                       :: mpi_split  
  character(len=1)                                :: mpi_split_ 
  logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
  logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
  complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
  !
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nso,Nso,Lk],'dmft_get_gloc_matsubara_superc_main_mpi',"Hk")
  call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main_mpi',"Smats")
  call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main_mpi',"Gmats")
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
  if(mpi_master)write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
  if(mpi_master)call start_timer
  Gmats=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_superc_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_superc_mpi(MpiComm,zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
        Gmats = Gmats + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lmats))
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_superc(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
        Gtmp = Gtmp + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_main_mpi


subroutine dmft_get_gloc_matsubara_superc_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Gmats,Smats)
  integer                                                               :: MpiComm
  real(8),dimension(:,:),intent(in)                                     :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in)           :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                          :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:,:),intent(in)                          :: Smats     ![2][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout)                       :: Gmats     !as Smats
  ! arrays
  complex(8),dimension(2,2,size(Ebands,1),size(Ebands,1),size(Smats,6)) :: zeta_mats ![2][2][Nspin*Norb][Nspin*Norb][Lmats]
  complex(8),dimension(2*size(Ebands,1),2*size(Ebands,1))               :: Gmatrix
  !
  complex(8),dimension(:,:,:,:,:,:),allocatable                         :: Gtmp    !as Smats
  !
  real(8)                                                               :: beta
  real(8)                                                               :: xmu,eps
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Ebands,2)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_matsubara_superc_dos',"Ebands")
  call assert_shape(Smats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Smats")
  call assert_shape(Gmats,[2,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_main',"Gmats")
  !
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  zeta_mats=zero
  do i=1,Lmats
     zeta_mats(1,1,:,:,i) = (xi*wm(i)+xmu)*eye(Nso) - diag(Hloc) -        nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb)
     zeta_mats(1,2,:,:,i) =                                      -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
     zeta_mats(2,1,:,:,i) =                                      -        nn2so_reshape(Smats(2,:,:,:,:,i),Nspin,Norb)
     zeta_mats(2,2,:,:,i) = (xi*wm(i)-xmu)*eye(Nso) + diag(Hloc) + conjg( nn2so_reshape(Smats(1,:,:,:,:,i),Nspin,Norb) )
  enddo
  !
  !invert (Z-Hk) for each k-point
  if(mpi_master)write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
  if(mpi_master)call start_timer
  allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lmats))
  Gtmp=zero
  Gmats=zero
  do ik = 1+mpi_rank, Lk, mpi_size
     !
     do i=1,Lmats
        Gmatrix  = zero
        Gmatrix(1:Nso,1:Nso)             = zeta_mats(1,1,:,:,i) - diag(Ebands(:,ik))
        Gmatrix(1:Nso,Nso+1:2*Nso)       = zeta_mats(1,2,:,:,i)
        Gmatrix(Nso+1:2*Nso,1:Nso)       = zeta_mats(2,1,:,:,i)
        Gmatrix(Nso+1:2*Nso,Nso+1:2*Nso) = zeta_mats(2,2,:,:,i) + diag(Ebands(:,ik))
        !
        call inv(Gmatrix)  ! PAY ATTENTION HERE: it is not guaranteed that Gloc is a symmetric matrix
        !
        do ispin=1,Nspin
           do iorb=1,Norb
              io = iorb + (ispin-1)*Norb
              Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + Gmatrix(io,io)*Dbands(io,ik)
              Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + Gmatrix(io,Nso+io)*Dbands(io,ik)
           enddo
        enddo
        !
     enddo
     !
     call eta(ik,Lk)
  enddo
  call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  if(mpi_master)call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_dos_mpi


subroutine dmft_get_gloc_matsubara_superc_ineq_mpi(MpiComm,Hk,Wtk,Gmats,Smats,mpi_split,hk_symm)
  integer                                           :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,4)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Smats     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  character(len=*),optional                         :: mpi_split  
  character(len=1)                                  :: mpi_split_ 
  logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
  logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  !
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
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
  call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_superc_ineq_main_mpi',"Hk")
  call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_ineq_main_mpi',"Smats")
  call assert_shape(Gmats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_ineq_main_mpi',"Gmats")
  !
  if(mpi_master)write(*,"(A)")"Get local Matsubara Superc Green's function (no print)"
  if(mpi_master)call start_timer
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
  if(mpi_master)call start_timer
  Gmats=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_superc_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_superc_ineq_mpi(MpiComm,zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
        Gmats = Gmats + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_superc_ineq(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats)
        Gtmp = Gtmp + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_ineq_mpi



subroutine dmft_get_gloc_matsubara_superc_gij_mpi(MpiComm,Hk,Wtk,Gmats,Fmats,Smats,mpi_split,hk_symm)
  integer                                           :: MpiComm
  complex(8),dimension(:,:,:,:)                       :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Gmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Fmats           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Smats           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  character(len=*),optional                         :: mpi_split  
  character(len=1)                                  :: mpi_split_ 
  logical,optional                                  :: hk_symm(size(Hk,4))
  logical                                           :: hk_symm_(size(Hk,4))
  !
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_mats       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkmats,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  integer                                           :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
  !
  !
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
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
  call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Hk")
  call assert_shape(Smats,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Smats")
  call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Gmats")
  call assert_shape(Fmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_superc_gij_main_mpi',"Fmats")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  if(mpi_master)write(*,"(A)")"Get full Green's function (no print)"
  !
  allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
  allocate(Fkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Fkmats=zero
  allocate(zeta_mats(2,2,Nlat,Nso,Nso,Lmats));zeta_mats=zero
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
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
  if(mpi_master)call start_timer
  Gmats=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_superc_gij_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_superc_gij(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
        Gmats = Gmats + Gkmats*Wtk(ik)
        Fmats = Fmats + Fkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
     allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Ftmp=zero
     do ik=1,Lk
        call invert_gk_superc_gij(zeta_mats,Hk(:,:,:,ik),hk_symm_(ik),Gkmats,Fkmats)
        Gtmp = Gtmp + Gkmats*Wtk(ik)
        Ftmp = Ftmp + Fkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp, Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     call Mpi_AllReduce(Ftmp, Fmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp,Ftmp)
     !
  end select
  if(mpi_master)call stop_timer
end subroutine dmft_get_gloc_matsubara_superc_gij_mpi










