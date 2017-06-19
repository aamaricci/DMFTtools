
subroutine dmft_get_gloc_realaxis_superc_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb][Nk]
  real(8),dimension(size(Hk,4)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  integer,intent(in)                              :: iprint
  character(len=*),optional                       :: mpi_split  
  character(len=1)                                :: mpi_split_ 
  logical,dimension(size(Hk,4)),optional          :: hk_symm   ![Nk]
  logical,dimension(size(Hk,4))                   :: hk_symm_  ![Nk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
  complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
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
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = size(Sreal,2)
  Norb  = size(Sreal,4)
  Lreal = size(Sreal,6)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nso,Nso,Lk],'dmft_get_gloc_realaxis_superc_main_mpi',"Hk")
  call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_mpi',"Sreal")
  call assert_shape(Greal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_mpi',"Greal")
  !
  allocate(Gkreal(2,Nspin,Nspin,Norb,Norb,Lreal))
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
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)call start_timer
  Greal=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_realaxis_superc_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_superc_mpi(MpiComm,zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
        Greal = Greal + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lreal))
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_superc(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
        Gtmp = Gtmp + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_realaxis(wr,Greal(1,:,:,:,:,:),"Gloc",iprint)
  if(mpi_master)call dmft_gloc_print_realaxis(wr,Greal(2,:,:,:,:,:),"Floc",iprint)
end subroutine dmft_get_gloc_realaxis_superc_main_mpi


subroutine dmft_get_gloc_realaxis_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                           :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,4)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  integer,intent(in)                                :: iprint
  character(len=*),optional                         :: mpi_split  
  character(len=1)                                  :: mpi_split_ 
  logical,dimension(size(Hk,4)),optional            :: hk_symm   ![Nk]
  logical,dimension(size(Hk,4))                     :: hk_symm_  ![Nk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
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
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nlat  = size(Sreal,2)
  Nspin = size(Sreal,3)
  Norb  = size(Sreal,5)
  Lreal = size(Sreal,7)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_superc_lattice_main_mpi',"Hk")
  call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_lattice_main_mpi',"Sreal")
  call assert_shape(Greal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_lattice_main_mpi',"Greal")
  !
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)call start_timer
  !
  allocate(Gkreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
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
  if(mpi_master)call start_timer
  Greal=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_realaxis_superc_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_superc_lattice_mpi(MpiComm,zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
        Greal = Greal + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_superc_lattice(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal)
        Gtmp = Gtmp + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_realaxis_lattice(wr,Greal(1,:,:,:,:,:,:),"LG",iprint)
  if(mpi_master)call dmft_gloc_print_realaxis_lattice(wr,Greal(2,:,:,:,:,:,:),"LF",iprint)
end subroutine dmft_get_gloc_realaxis_superc_lattice_main_mpi

subroutine dmft_get_gloc_realaxis_superc_gij_main_mpi(MpiComm,Hk,Wtk,Greal,Freal,Sreal,iprint,mpi_split,hk_symm)
  integer                                           :: MpiComm
  complex(8),dimension(:,:,:,:)                       :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  real(8)                                           :: Wtk(size(Hk,4)) ![Nk]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Greal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Freal           ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),intent(inout),dimension(:,:,:,:,:,:,:) :: Sreal           ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  integer                                           :: iprint
  character(len=*),optional                         :: mpi_split  
  character(len=1)                                  :: mpi_split_ 
  logical,optional                                  :: hk_symm(size(Hk,4))
  logical                                           :: hk_symm_(size(Hk,4))
  !
  complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real       ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal,Gtmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Fkreal,Ftmp          ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  integer                                           :: ik,Lk,Nlso,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,js
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
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  Nlat  = size(Sreal,2)
  Nspin = size(Sreal,3)
  Norb  = size(Sreal,5)
  Lreal = size(Sreal,7)
  Lk    = size(Hk,4)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[2,Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Hk")
  call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Sreal")
  call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Greal")
  call assert_shape(Freal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_main_gij_mpi',"Freal")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  if(mpi_master)write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
  !
  allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
  allocate(Fkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Fkreal=zero
  allocate(zeta_real(2,2,Nlat,Nso,Nso,Lreal));zeta_real=zero
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  !SYMMETRIES in real-frequencies   [assuming a real order parameter]
  !G22(w)  = -[G11[-w]]*
  !G21(w)  =   G12[w]
  do ilat=1,Nlat
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           js = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
           zeta_real(1,1,ilat,io,io,:) = dcmplx(wr(:),eps) + xmu
           zeta_real(2,2,ilat,io,io,:) = -conjg( dcmplx(wr(Lreal:1:-1),eps) + xmu )
        enddo
     enddo
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 zeta_real(1,1,ilat,io,jo,:) = zeta_real(1,1,ilat,io,jo,:) - &
                      Sreal(1,ilat,ispin,jspin,iorb,jorb,:)
                 zeta_real(1,2,ilat,io,jo,:) = zeta_real(1,2,ilat,io,jo,:) - &
                      Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                 zeta_real(2,1,ilat,io,jo,:) = zeta_real(2,1,ilat,io,jo,:) - &
                      Sreal(2,ilat,ispin,jspin,iorb,jorb,:)
                 zeta_real(2,2,ilat,io,jo,:) = zeta_real(2,2,ilat,io,jo,:) + &
                      conjg( Sreal(1,ilat,ispin,jspin,iorb,jorb,Lreal:1:-1) )
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  if(mpi_master)call start_timer
  Greal=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_normal_realaxis_gij_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_superc_gij_mpi(MpiComm,zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
        Greal = Greal + Gkreal*Wtk(ik)
        Freal = Freal + Fkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
     allocate(Ftmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Ftmp=zero
     do ik=1,Lk
        call invert_gk_superc_gij(zeta_real,Hk(:,:,:,ik),hk_symm_(ik),Gkreal,Fkreal)
        Gtmp = Gtmp + Gkreal*Wtk(ik)
        Ftmp = Ftmp + Fkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp, Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     call Mpi_AllReduce(Ftmp, Freal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp,Ftmp)
     !
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_realaxis_gij(wm,Greal,"Gij",iprint)
  if(mpi_master)call dmft_gloc_print_realaxis_gij(wm,Freal,"Fij",iprint)
end subroutine dmft_get_gloc_realaxis_superc_gij_main_mpi










subroutine dmft_get_gloc_realaxis_superc_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                       :: MpiComm
  complex(8),dimension(:,:),intent(in)            :: Hk              ![2][Nk]
  real(8),intent(in)                            :: Wtk(size(Hk,2))   ![Nk]
  complex(8),intent(in)                         :: Sreal(:,:)
  complex(8),intent(inout)                      :: Greal(2,size(Sreal,2))
  logical,optional                              :: hk_symm(size(Hk,2))
  logical                                       :: hk_symm_(size(Hk,2))
  integer                                       :: iprint
  character(len=*),optional                     :: mpi_split  
  character(len=1)                              :: mpi_split_ 
  !
  complex(8),dimension(2,1,1,size(Hk,2))            :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
  complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Greal_
  complex(8),dimension(2,1,1,1,1,size(Sreal,2)) :: Sreal_
  !
  call assert_shape(Sreal,[2,size(Sreal,2)],'dmft_get_gloc_realaxis_superc_1band_mpi',"Sreal")
  !
  Hk_(:,1,1,:)          = Hk
  Sreal_(:,1,1,1,1,:) = Sreal(:,:)
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  call dmft_get_gloc_realaxis_superc_main_mpi(MpiComm,Hk_,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
  Greal(:,:) = Greal_(:,1,1,1,1,:)
end subroutine dmft_get_gloc_realaxis_superc_1band_mpi


subroutine dmft_get_gloc_realaxis_superc_lattice_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                                  :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)                   :: Hk              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  real(8),intent(in)                                       :: Wtk(size(Hk,4)) ![Nk]
  complex(8),intent(in)                                    :: Sreal(:,:,:)
  complex(8),intent(inout)                                 :: Greal(2,size(Hk,2),size(Sreal,3))
  logical,optional                                         :: hk_symm(size(Hk,4))
  logical                                                  :: hk_symm_(size(Hk,4))
  integer                                                  :: iprint
  character(len=*),optional                                :: mpi_split  
  character(len=1)                                         :: mpi_split_ 
  !
  complex(8),dimension(2,size(Hk,2),1,1,1,1,size(Sreal,3)) :: Greal_
  complex(8),dimension(2,size(Hk,2),1,1,1,1,size(Sreal,3)) :: Sreal_
  !
  call assert_shape(Hk,[2,size(Hk,2),size(Hk,2),size(Hk,4)],'dmft_get_gloc_realaxis_superc_lattice_1band_mpi',"Hk")
  call assert_shape(Sreal,[2,size(Hk,2),size(Sreal,3)],'dmft_get_gloc_realaxis_superc_lattice_1band_mpi',"Sreal")
  !
  Sreal_(:,:,1,1,1,1,:) = Sreal(:,:,:)
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  call dmft_get_gloc_realaxis_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
  Greal(:,:,:) = Greal_(:,:,1,1,1,1,:)
end subroutine dmft_get_gloc_realaxis_superc_lattice_1band_mpi









subroutine dmft_get_gloc_realaxis_superc_Nband_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                                         :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)                          :: Hk              ![2][Norb][Norb][Nk]
  real(8),intent(in)                                              :: Wtk(size(Hk,4)) ![Nk]
  complex(8),intent(in)                                           :: Sreal(:,:,:,:)  ![2][Norb][Norb][Lreal]
  complex(8),intent(inout)                                        :: Greal(:,:,:,:)  !as Sreal
  logical,optional                                                :: hk_symm(size(Hk,4))
  logical                                                         :: hk_symm_(size(Hk,4))
  integer                                                         :: iprint
  character(len=*),optional                                       :: mpi_split  
  character(len=1)                                                :: mpi_split_ 
  !
  complex(8),&
       dimension(2,1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Greal_
  complex(8),&
       dimension(2,1,1,size(Sreal,2),size(Sreal,3),size(Sreal,4)) :: Sreal_
  integer                                                         :: Nspin,Norb,Nso,Lreal
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = 1
  Norb  = size(Sreal,2)
  Lreal = size(Sreal,4)
  Nso   = Nspin*Norb
  call assert_shape(Hk,[2,Nso,Nso,size(Hk,4)],'dmft_get_gloc_realaxis_superc_Nband_mpi',"Hk")
  call assert_shape(Sreal,[2,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_Nband_mpi',"Sreal")
  call assert_shape(Greal,[2,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_Nband_mpi',"Greal")
  !
  Sreal_(:,1,1,:,:,:) = Sreal(:,:,:,:)
  call dmft_get_gloc_realaxis_superc_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
  Greal(:,:,:,:) = Greal_(:,1,1,:,:,:)
end subroutine dmft_get_gloc_realaxis_superc_Nband_mpi



subroutine dmft_get_gloc_realaxis_superc_lattice_Nband_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                                                       :: MpiComm
  complex(8),dimension(:,:,:,:),intent(in)                                        :: Hk               ![2][Nlat*Norb][Nlat*Norb][Nk]
  real(8),intent(in)                                                            :: Wtk(size(Hk,4))  ![Nk]
  complex(8),intent(in)                                                         :: Sreal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lreal]
  complex(8),intent(inout)                                                      :: Greal(:,:,:,:,:) ![2][Nlat][Norb][Norb][Lreal]
  logical,optional                                                              :: hk_symm(size(Hk,4))
  logical                                                                       :: hk_symm_(size(Hk,4))
  integer                                                                       :: iprint
  character(len=*),optional                                                     :: mpi_split  
  character(len=1)                                                              :: mpi_split_ 
  !
  complex(8),&
       dimension(2,size(Sreal,2),1,1,size(Sreal,3),size(Sreal,4),size(Sreal,5)) :: Greal_
  complex(8),&
       dimension(2,size(Sreal,2),1,1,size(Sreal,3),size(Sreal,4),size(Sreal,5)) :: Sreal_
  integer                                                                       :: Nspin,Norb,Nso,Nlso,Lreal
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = 1
  Nlat  = size(Sreal,2)
  Norb  = size(Sreal,3)
  Lreal = size(Sreal,5)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  call assert_shape(Hk,[2,Nlso,Nlso,size(Hk,4)],'dmft_get_gloc_realaxis_superc_lattice_Nband_mpi',"Hk")
  call assert_shape(Sreal,[2,Nlat,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_lattice_Nband_mpi',"Sreal")
  call assert_shape(Greal,[2,Nlat,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_superc_lattice_Nband_mpi',"Greal")
  !
  Sreal_(:,:,1,1,:,:,:) = Sreal(:,:,:,:,:)
  call dmft_get_gloc_realaxis_superc_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
  Greal(:,:,:,:,:) = Greal_(:,:,1,1,:,:,:)
end subroutine dmft_get_gloc_realaxis_superc_lattice_Nband_mpi


