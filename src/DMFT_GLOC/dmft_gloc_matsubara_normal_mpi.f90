subroutine dmft_get_gloc_matsubara_normal_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
  integer                                       :: MpiComm
  complex(8),dimension(:,:,:),intent(in)        :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
  real(8),dimension(size(Hk,3)),intent(in)      :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  integer,intent(in)                            :: iprint
  character(len=*),optional                     :: mpi_split  
  character(len=1)                              :: mpi_split_ 
  logical,dimension(size(Hk,3)),optional        :: hk_symm   ![Lk]
  logical,dimension(size(Hk,3))                 :: hk_symm_  ![Lk]
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gkmats    !as Smats  
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Smats
  complex(8),dimension(:,:,:),allocatable       :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_  =.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Hk,[Nso,Nso,Lk],"dmft_get_gloc_matsubara_normal_main_mpi","Hk")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_main_mpi","Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_main_mpi","Gmats")
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
  if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)call start_timer
  Gmats=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_normal_main_mpi: ! mpi_split_ in ['w','k'] "
  case ('w')
     do ik=1,Lk
        call invert_gk_normal_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)      
        Gmats = Gmats + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
  case ('k')
     allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_normal(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
        Gtmp = Gtmp + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_matsubara(wm,Gmats,"Gloc",iprint)
end subroutine dmft_get_gloc_matsubara_normal_main_mpi


subroutine dmft_get_gloc_matsubara_normal_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Gmats,Smats,iprint,mpi_split)
  integer                                                     :: MpiComm
  real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
  real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
  complex(8),dimension(:,:,:,:,:),intent(in)                  :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout)               :: Gmats     !as Smats
  integer,intent(in)                                          :: iprint
  character(len=*),optional                                   :: mpi_split  
  character(len=1)                                            :: mpi_split_ 
  !
  complex(8)                                                  :: gktmp
  complex(8),dimension(:,:,:),allocatable                     :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),allocatable                 :: Gtmp      !as Smats
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
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  !
  Nspin = size(Smats,1)
  Norb  = size(Smats,3)
  Lmats = size(Smats,5)
  Lk    = size(Ebands,2)
  Nso   = Nspin*Norb    
  !Testing part:
  call assert_shape(Ebands,[Nso,Lk],"dmft_get_gloc_matsubara_normal_dos_main_mpi","Ebands")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_dos_main_mpi","Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_gloc_matsubara_normal_dos_main_mpi","Gmats")
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
  if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)call start_timer
  Gmats=zero
  allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
  !
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_normal_dos_main_mpi: ! mpi_split_ in ['w','k'] "
  case ('w')
     do i = 1+mpi_rank, Lmats, mpi_size
        do ispin=1,Nspin
           do iorb=1,Norb
              io = iorb + (ispin-1)*Norb
              do ik=1,Lk
                 gktmp = Dbands(io,ik)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io,ik) )
                 Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
              enddo
           enddo
        enddo
        if(mpi_master)call eta(i,Lmats)
     end do
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  case ('k')
     do ik = 1+mpi_rank, Lk, mpi_size
        do ispin=1,Nspin
           do iorb=1,Norb
              io = iorb + (ispin-1)*Norb
              do i=1,Lmats
                 gktmp = Dbands(io,ik)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io,ik) )
                 Gtmp(ispin,ispin,iorb,iorb,i) = Gtmp(ispin,ispin,iorb,iorb,i) + gktmp
              enddo
           enddo
        enddo
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_matsubara(wm,Gmats,"Gloc",iprint)
end subroutine dmft_get_gloc_matsubara_normal_dos_mpi


subroutine dmft_get_gloc_matsubara_normal_ineq_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
  integer                                         :: MpiComm
  complex(8),dimension(:,:,:),intent(in)          :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)        :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Smats     ![Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gmats     !as Smats
  integer,intent(in)                              :: iprint    !
  logical,optional                                :: tridiag
  logical                                         :: tridiag_
  character(len=*),optional                       :: mpi_split  
  character(len=1)                                :: mpi_split_ 
  logical,dimension(size(Hk,3)),optional          :: hk_symm
  logical,dimension((size(Hk,3)))                 :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Smats
  complex(8),dimension(:,:,:,:),allocatable       :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
  !
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
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
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Hk")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Smats")
  call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_ineq_main_mpi',"Gmats")
  !
  if(mpi_master)write(*,"(A)")"Get local Matsubara Green's function (print mode:"//reg(txtfy(iprint))//")"
  if(mpi_master)then
     if(.not.tridiag_)then
        write(*,"(A)")"Direct Inversion algorithm:"
     else
        write(*,"(A)")"Quantum Zipper algorithm:"
     endif
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
  if(mpi_master)call start_timer
  Gmats=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_normal_main_ineq_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     if(.not.tridiag_)then
        do ik=1,Lk
           call invert_gk_normal_ineq_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
           Gmats = Gmats + Gkmats*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     else
        do ik=1,Lk
           call invert_gk_normal_tridiag_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
           Gmats = Gmats + Gkmats*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     endif
     !
  case ('k')
     allocate(Gtmp(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero
     if(.not.tridiag_)then
        do ik=1+mpi_rank,Lk,mpi_size
           call invert_gk_normal_ineq(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
           Gtmp = Gtmp + Gkmats*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     else
        do ik=1+mpi_rank,Lk,mpi_size
           call invert_gk_normal_tridiag(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
           Gtmp = Gtmp + Gkmats*Wtk(ik)
           if(mpi_master)call eta(ik,Lk)
        end do
     endif
     call Mpi_AllReduce(Gtmp,Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_matsubara_ineq(wm,Gmats,"LG",iprint)
end subroutine dmft_get_gloc_matsubara_normal_ineq_mpi








subroutine dmft_get_gloc_matsubara_normal_gij_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
  integer                                           :: MpiComm
  complex(8),dimension(:,:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
  real(8),dimension(size(Hk,3)),intent(in)          :: Wtk       ![Nk]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats     !      [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  integer,intent(in)                                :: iprint
  character(len=*),optional                         :: mpi_split  
  character(len=1)                                  :: mpi_split_ 
  logical,dimension(size(Hk,3)),optional            :: hk_symm
  logical,dimension((size(Hk,3)))                   :: hk_symm_
  !allocatable arrays
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gkmats    !as Smats
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp
  complex(8),dimension(:,:,:,:),allocatable         :: zeta_mats ![Nlat][Nspin*Norb][Nspin*Norb][Lmats]
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
  Nlat  = size(Smats,1)
  Nspin = size(Smats,2)
  Norb  = size(Smats,4)
  Lmats = size(Smats,6)
  Lk    = size(Hk,3)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  !Testing part:
  call assert_shape(Hk,[Nlso,Nlso,Lk],'dmft_get_gloc_matsubara_normal_gij_main_mpi',"Hk")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_gij_main_mpi',"Smats")
  call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_gij_main_mpi',"Gmats")
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  if(mpi_master)write(*,"(A)")"Get full Green's function (print mode:"//reg(txtfy(iprint))//")"
  !
  allocate(Gkmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gkmats=zero
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
  if(mpi_master)call start_timer
  Gmats=zero
  select case(mpi_split_)
  case default
     stop "dmft_get_gloc_matsubara_normal_gij_main_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     do ik=1,Lk
        call invert_gk_normal_gij_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
        Gmats = Gmats + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     !
  case ('k')
     allocate(Gtmp(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gtmp=zero     
     do ik=1+mpi_rank,Lk,mpi_size
        call invert_gk_normal_gij(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
        Gtmp = Gtmp + Gkmats*Wtk(ik)
        if(mpi_master)call eta(ik,Lk)
     end do
     call Mpi_AllReduce(Gtmp, Gmats, size(Gmats), MPI_Double_Complex, MPI_Sum, MpiComm, MPI_ierr)
     deallocate(Gtmp)
     !
  end select
  if(mpi_master)call stop_timer
  if(mpi_master)call dmft_gloc_print_matsubara_gij(wm,Gmats,"Gij",iprint)
end subroutine dmft_get_gloc_matsubara_normal_gij_mpi




