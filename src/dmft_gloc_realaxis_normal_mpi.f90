#define FROUTINE 'dmft_get_gloc_realaxis_normal_main_mpi'
subroutine dmft_get_gloc_realaxis_normal_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                       :: MpiComm
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  integer,intent(in)                            :: iprint
  character(len=*),optional                     :: mpi_split  
  character(len=1)                              :: mpi_split_ 
  logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
  logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gkreal    !as Sreal  
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
  complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal] 
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
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
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)call start_timer
  Greal=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_realaxis_normal_main_mpi: ! mpi_split_ in ['w','k'] "
  case ('w')
     do ik=1,Lk
        call invert_gk_normal_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
        Greal = Greal + Gkreal*Wtk(ik)
        call eta(ik,Lk)
     end do
  case ('k')
     allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_normal(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
        Gtmp = Gtmp + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_realaxis(wr,Greal,"Gloc",iprint)
end subroutine dmft_get_gloc_realaxis_normal_main_mpi
#undef FROUTINE


#define FROUTINE 'dmft_get_gloc_realaxis_normal_dos_main_mpi'
subroutine dmft_get_gloc_realaxis_normal_dos_main_mpi(MpiComm,Ebands,Dbands,Hloc,Greal,Sreal,iprint,mpi_split)
  integer                                                     :: MpiComm
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                  :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout)               :: Greal     !as Sreal
  integer,intent(in)                                          :: iprint
  character(len=*),optional                                   :: mpi_split  
  character(len=1)                                            :: mpi_split_ 
  !
  complex(8)                                                  :: gktmp
  complex(8),dimension(:,:,:),allocatable                     :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),allocatable                 :: Gtmp      !as Sreal
  !
  real(8)                                                     :: beta
  real(8)                                                     :: xmu,eps
  real(8)                                                     :: wini,wfin
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !Retrieve parameters:
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
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
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)call start_timer
  Greal=zero
  allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
  !
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_realaxis_normal_dos_main_mpi: ! mpi_split_ in ['w','k'] "
  case ('w')
     do i = 1+mpi_rank, Lreal, mpi_size
        do ispin=1,Nspin
           do iorb=1,Norb
              io = iorb + (ispin-1)*Norb
              do ik=1,Lk
                 gktmp = Dbands(io,ik)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io,ik) )
                 Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
              enddo
           enddo
        enddo
        if(mpi_master)call eta(i,Lreal)
     end do
     call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  case ('k')
     do ik = 1+mpi_rank, Lk, mpi_size
        do ispin=1,Nspin
           do iorb=1,Norb
              io = iorb + (ispin-1)*Norb
              do i=1,Lreal
                 gktmp = Dbands(io,ik)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io,ik) )
                 Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
              enddo
           enddo
        enddo
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_realaxis(wm,Greal,"Gloc",iprint)
end subroutine dmft_get_gloc_realaxis_normal_dos_main_mpi
#undef FROUTINE



#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_main_mpi'
subroutine dmft_get_gloc_realaxis_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,tridiag,mpi_split,hk_symm)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
  integer,intent(in)                              :: iprint    !
  logical,optional                                :: tridiag
  logical                                         :: tridiag_
  character(len=*),optional                       :: mpi_split  
  character(len=1)                                :: mpi_split_ 
  logical,dimension(size(Hk,3)),optional          :: hk_symm
  logical,dimension((size(Hk,3)))                 :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
  complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
  !
  real(8)                                         :: beta
  real(8)                                         :: xmu,eps
  real(8)                                         :: wini,wfin  
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  call get_ctrl_var(eps,"EPS")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
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
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)then
     if(.not.tridiag_)then
        write(*,"(A)")"Direct Inversion algorithm:"
     else
        write(*,"(A)")"Quantum Zipper algorithm:"
     endif
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
  if(mpi_master)call start_timer
  Greal=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_realaxis_normal_lattice_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     if(.not.tridiag_)then
        do ik=1,Lk
           call invert_gk_normal_lattice_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
           Greal = Greal + Gkreal*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     else
        do ik=1,Lk
           call invert_gk_normal_tridiag_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
           Greal = Greal + Gkreal*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     endif
     !
  case ('k')
     allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
     if(.not.tridiag_)then
        do ik=1+mpi_rank,Lk,mpi_size
           call invert_gk_normal_lattice(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
           Gtmp = Gtmp + Gkreal*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     else
        do ik=1+mpi_rank,Lk,mpi_size
           call invert_gk_normal_tridiag(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
           Gtmp = Gtmp + Gkreal*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     endif
     call Mpi_AllReduce(Gtmp,Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_realaxis_lattice(wr,Greal,"LG",iprint)
end subroutine dmft_get_gloc_realaxis_normal_lattice_main_mpi
#undef FROUTINE







!##################################################################
!##################################################################
!##################################################################







#define FROUTINE 'dmft_get_gloc_realaxis_normal_1band_mpi'
subroutine dmft_get_gloc_realaxis_normal_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,mpi_split,hk_symm)
  integer                                   :: MpiComm
  complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
  real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
  complex(8),intent(in)                     :: Sreal(:)
  complex(8),intent(inout)                  :: Greal(size(Sreal))
  logical,optional                          :: hk_symm(size(Hk,1))
  logical                                   :: hk_symm_(size(Hk,1))
  integer                                   :: iprint
  character(len=*),optional                 :: mpi_split  
  character(len=1)                          :: mpi_split_ 
  !
  complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
  complex(8),dimension(1,1,1,1,size(Sreal)) :: Greal_
  complex(8),dimension(1,1,1,1,size(Sreal)) :: Sreal_
  !
  Hk_(1,1,:)        = Hk(:)
  Sreal_(1,1,1,1,:) = Sreal(:)
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  call dmft_get_gloc_realaxis_normal_main_mpi(MpiComm,Hk_,Wtk,Greal_,Sreal_,iprint,mpi_split_,hk_symm_)
  Greal(:) = Greal_(1,1,1,1,:)
end subroutine dmft_get_gloc_realaxis_normal_1band_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_realaxis_normal_lattice_1band_mpi'
subroutine dmft_get_gloc_realaxis_normal_lattice_1band_mpi(MpiComm,Hk,Wtk,Greal,Sreal,iprint,tridiag,mpi_split,hk_symm)
  integer                                                         :: MpiComm
  complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
  real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
  complex(8),dimension(:,:),intent(in)                            :: Sreal           ![Nlat][Lreal]
  complex(8),dimension(size(Sreal,1),size(Sreal,2)),intent(inout) :: Greal
  integer                                                         :: iprint
  logical,optional                                                :: tridiag
  logical                                                         :: tridiag_
  character(len=*),optional                                       :: mpi_split  
  character(len=1)                                                :: mpi_split_ 
  logical,optional                                                :: hk_symm(size(Hk,1))
  logical                                                         :: hk_symm_(size(Hk,1))
  !
  complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Greal_
  complex(8),dimension(size(Sreal,1),1,1,1,1,size(Sreal,2))       :: Sreal_
  !
  tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Sreal,[size(Hk,1),size(Sreal,2)],FROUTINE,"Sreal")
  Sreal_(:,1,1,1,1,:) = Sreal(:,:)
  call dmft_get_gloc_realaxis_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Greal_,Sreal_,iprint,tridiag_,mpi_split_,hk_symm_)
  Greal(:,:) = Greal_(:,1,1,1,1,:)
end subroutine dmft_get_gloc_realaxis_normal_lattice_1band_mpi
#undef FROUTINE













!##################################################################
!##################################################################
!##################################################################


