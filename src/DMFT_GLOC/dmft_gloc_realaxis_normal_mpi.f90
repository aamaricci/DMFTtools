subroutine dmft_get_gloc_realaxis_normal_main_mpi(MpiComm,Hk,Wtk,Greal,Sreal,mpi_split,hk_symm)
  integer                                       :: MpiComm
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Greal     !as Sreal
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
  call assert_shape(Hk,[Nso,Nso,Lk],'dmft_get_gloc_realaxis_normal_main_mpi',"Hk")
  call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_main_mpi',"Sreal")
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_main_mpi',"Greal")
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
  if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
  if(mpi_master)call start_timer
  Greal=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_realaxis_normal_main_mpi: ! mpi_split_ in ['w','k'] "
  case ('w')
     do ik=1,Lk
        call invert_gk_normal_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)      
        Greal = Greal + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
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
end subroutine dmft_get_gloc_realaxis_normal_main_mpi


subroutine dmft_get_gloc_realaxis_normal_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Greal,Sreal,mpi_split)
  integer                                                     :: MpiComm
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                  :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:),intent(inout)               :: Greal     !as Sreal
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
  call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_realaxis_normal_dos_main_mpi',"Ebands")
  call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos_main_mpi',"Sreal")
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos_main_mpi',"Greal")
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
  if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
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
end subroutine dmft_get_gloc_realaxis_normal_dos_mpi


subroutine dmft_get_gloc_realaxis_normal_ineq_mpi(MpiComm,Hk,Wtk,Greal,Sreal,tridiag,mpi_split,hk_symm)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     !as Sreal
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
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_ineq_main_mpi',"Hk")
  call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_ineq_main_mpi',"Sreal")
  call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_ineq_main_mpi',"Greal")
  !
  if(mpi_master)write(*,"(A)")"Get local Realaxis Green's function (no print)"
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
     stop "dmft_get_gloc_realaxis_normal_ineq_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     if(.not.tridiag_)then
        do ik=1,Lk
           call invert_gk_normal_ineq_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
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
           call invert_gk_normal_ineq(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
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
end subroutine dmft_get_gloc_realaxis_normal_ineq_mpi


subroutine dmft_get_gloc_realaxis_normal_gij_mpi(MpiComm,Hk,Wtk,Greal,Sreal,mpi_split,hk_symm)
  integer                                           :: MpiComm
  complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  character(len=*),optional                         :: mpi_split  
  character(len=1)                                  :: mpi_split_ 
  logical,dimension(size(Hk,3)),optional            :: hk_symm
  logical,dimension((size(Hk,3)))                   :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkreal    !as Sreal
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
  complex(8),dimension(:,:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
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
  Nlat  = size(Sreal,1)
  Nspin = size(Sreal,2)
  Norb  = size(Sreal,4)
  Lreal = size(Sreal,6)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_realaxis_normal_gij_main_mpi',"Hk")
  call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_gij_main_mpi',"Sreal")
  call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_gij_main_mpi',"Greal")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  if(mpi_master)write(*,"(A)")"Get full Green's function (no print)"
  !
  allocate(Gkreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
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
     stop "dmft_get_gloc_realaxis_normal_gij_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_normal_gij_mpi(MpiComm,zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
        Greal = Greal + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero     
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_normal_gij(zeta_real,Hk(:,:,ik),hk_symm_(ik),Gkreal)
        Gtmp = Gtmp + Gkreal*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp, Greal, size(Greal), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
end subroutine dmft_get_gloc_realaxis_normal_gij_mpi









