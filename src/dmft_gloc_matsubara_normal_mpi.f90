#define FROUTINE 'dmft_get_gloc_matsubara_normal_main_mpi'
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
  call assert_shape(Hk,[Nso,Nso,Lk],FROUTINE,"Hk")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
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
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
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
#undef FROUTINE


#define FROUTINE 'dmft_get_gloc_matsubara_normal_dos_main_mpi'
subroutine dmft_get_gloc_matsubara_normal_dos_main_mpi(MpiComm,Ebands,Dbands,Hloc,Gmats,Smats,iprint,mpi_split)
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
  call assert_shape(Ebands,[Nso,Lk],FROUTINE,"Ebands")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
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
  if(mpi_master)write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
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
end subroutine dmft_get_gloc_matsubara_normal_dos_main_mpi
#undef FROUTINE


#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_main_mpi'
subroutine dmft_get_gloc_matsubara_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
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
  call assert_shape(Hk,[Nlso,Nlso,Lk],FROUTINE,"Hk")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Smats")
  call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],FROUTINE,"Gmats")
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
     stop "dmft_get_gloc_matsubara_normal_main_lattice_mpi: ! mpi_split_ in ['w','k'] "
     !
  case ('w')
     if(.not.tridiag_)then
        do ik=1,Lk
           call invert_gk_normal_lattice_mpi(MpiComm,zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
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
           call invert_gk_normal_lattice(zeta_mats,Hk(:,:,ik),hk_symm_(ik),Gkmats)
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
  if(mpi_master)call dmft_gloc_print_matsubara_lattice(wm,Gmats,"LG",iprint)
end subroutine dmft_get_gloc_matsubara_normal_lattice_main_mpi
#undef FROUTINE












!##################################################################
!##################################################################
!##################################################################











#define FROUTINE 'dmft_get_gloc_matsubara_normal_1band_mpi'
subroutine dmft_get_gloc_matsubara_normal_1band_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
  integer                                   :: MpiComm
  complex(8),dimension(:),intent(in)        :: Hk              ![Nk]
  real(8),intent(in)                        :: Wtk(size(Hk))   ![Nk]
  complex(8),intent(in)                     :: Smats(:)
  complex(8),intent(inout)                  :: Gmats(size(Smats))
  logical,optional                          :: hk_symm(size(Hk,1))
  logical                                   :: hk_symm_(size(Hk,1))
  integer                                   :: iprint
  character(len=*),optional                 :: mpi_split  
  character(len=1)                          :: mpi_split_ 
  !
  complex(8),dimension(1,1,size(Hk))        :: Hk_             ![Norb*Nspin][Norb*Nspin][Nk]
  complex(8),dimension(1,1,1,1,size(Smats)) :: Gmats_
  complex(8),dimension(1,1,1,1,size(Smats)) :: Smats_
  !
  Hk_(1,1,:)        = Hk(:)
  Smats_(1,1,1,1,:) = Smats(:)
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  call dmft_get_gloc_matsubara_normal_main_mpi(MpiComm,Hk_,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
  Gmats(:) = Gmats_(1,1,1,1,:)
end subroutine dmft_get_gloc_matsubara_normal_1band_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_1band_mpi'
subroutine dmft_get_gloc_matsubara_normal_lattice_1band_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
  integer                                                         :: MpiComm
  complex(8),dimension(:,:,:),intent(in)                          :: Hk              ![Nlat][Nlat][Nk]
  real(8),intent(in)                                              :: Wtk(size(Hk,3)) ![Nk]
  complex(8),dimension(:,:),intent(in)                            :: Smats           ![Nlat][Lmats]
  complex(8),dimension(size(Smats,1),size(Smats,2)),intent(inout) :: Gmats
  integer                                                         :: iprint
  logical,optional                                                :: tridiag
  logical                                                         :: tridiag_
  character(len=*),optional                                       :: mpi_split  
  character(len=1)                                                :: mpi_split_ 
  logical,optional                                                :: hk_symm(size(Hk,1))
  logical                                                         :: hk_symm_(size(Hk,1))
  !
  complex(8),dimension(size(Smats,1),1,1,1,1,size(Smats,2))       :: Gmats_
  complex(8),dimension(size(Smats,1),1,1,1,1,size(Smats,2))       :: Smats_
  !
  tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  call assert_shape(Hk,[size(Hk,1),size(Hk,1),size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Smats,[size(Hk,1),size(Smats,2)],FROUTINE,"Smats")
  Smats_(:,1,1,1,1,:) = Smats(:,:)
  call dmft_get_gloc_matsubara_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,tridiag_,mpi_split_,hk_symm_)
  Gmats(:,:) = Gmats_(:,1,1,1,1,:)
end subroutine dmft_get_gloc_matsubara_normal_lattice_1band_mpi
#undef FROUTINE





!##################################################################
!##################################################################
!##################################################################







#define FROUTINE 'dmft_get_gloc_matsubara_normal_Nband_mpi'
subroutine dmft_get_gloc_matsubara_normal_Nband_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,mpi_split,hk_symm)
  integer                                                             :: MpiComm
  complex(8),dimension(:,:,:),intent(in)                              :: Hk              ![Norb][Norb][Nk]
  real(8)                                                             :: Wtk(size(Hk,3)) ![Nk]
  complex(8),intent(in),dimension(:,:,:)                              :: Smats(:,:,:)    ![Norb][Norb][Lmats]
  complex(8),intent(inout),dimension(:,:,:)                           :: Gmats
  logical,optional                                                    :: hk_symm(size(Hk,3))
  logical                                                             :: hk_symm_(size(Hk,3))
  integer                                                             :: iprint
  character(len=*),optional                                           :: mpi_split  
  character(len=1)                                                    :: mpi_split_ 
  !
  complex(8),dimension(1,1,size(Smats,1),size(Smats,2),size(Smats,3)) :: Gmats_
  complex(8),dimension(1,1,size(Smats,1),size(Smats,2),size(Smats,3)) :: Smats_
  integer                                                             :: Nspin,Norb,Nso,Lmats
  !
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = 1
  Norb  = size(Smats,1)
  Lmats = size(Smats,3)
  Nso   = Nspin*Norb
  call assert_shape(Hk,[Nso,Nso,size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Smats,[Nso,Nso,Lmats],FROUTINE,"Smats")
  call assert_shape(Gmats,[Nso,Nso,Lmats],FROUTINE,"Gmats")
  !
  Smats_(1,1,:,:,:) = Smats(:,:,:)
  call dmft_get_gloc_matsubara_normal_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,mpi_split_,hk_symm_)
  Gmats(:,:,:) = Gmats_(1,1,:,:,:)
end subroutine dmft_get_gloc_matsubara_normal_Nband_mpi
#undef FROUTINE

#define FROUTINE 'dmft_get_gloc_matsubara_normal_lattice_Nband_mpi'
subroutine dmft_get_gloc_matsubara_normal_lattice_Nband_mpi(MpiComm,Hk,Wtk,Gmats,Smats,iprint,tridiag,mpi_split,hk_symm)
  integer                                                                     :: MpiComm
  complex(8),dimension(:,:,:),intent(in)                                      :: Hk              ![Nlat*Norb][Nlat*Norb][Nk]
  real(8)                                                                     :: Wtk(size(Hk,3)) ![Nk]
  complex(8),intent(in)                                                       :: Smats(:,:,:,:)  ![Nlat][Norb][Norb][Lmats]
  complex(8),intent(inout)                                                    :: Gmats(:,:,:,:)
  logical,optional                                                            :: hk_symm(size(Hk,3))
  logical                                                                     :: hk_symm_(size(Hk,3))
  integer                                                                     :: iprint
  logical,optional                                                            :: tridiag
  logical                                                                     :: tridiag_
  character(len=*),optional                                                   :: mpi_split  
  character(len=1)                                                            :: mpi_split_ 
  !
  complex(8),&
       dimension(size(Smats,1),1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Gmats_ ![Nlat][1][1][Norb][Norb][Lmats]
  complex(8),&
       dimension(size(Smats,1),1,1,size(Smats,2),size(Smats,3),size(Smats,4)) :: Smats_![Nlat][1][1][Norb][Norb][Lmats]
  !
  integer                                                                     :: Nspin,Norb,Nso,Nlso,Lmats
  !
  !
  tridiag_=.false.;if(present(tridiag)) tridiag_=tridiag
  mpi_split_='w'    ;if(present(mpi_split)) mpi_split_=mpi_split
  hk_symm_=.false.;if(present(hk_symm)) hk_symm_=hk_symm
  !
  Nspin = 1    
  Nlat  = size(Smats,1)
  Norb  = size(Smats,2)
  Lmats = size(Smats,4)
  Nso   = Nspin*Norb
  Nlso  = Nlat*Nspin*Norb
  call assert_shape(Hk,[Nlso,Nlso,size(Hk,3)],FROUTINE,"Hk")
  call assert_shape(Smats,[Nlat,Nso,Nso,Lmats],FROUTINE,"Smats")
  call assert_shape(Gmats,[Nlat,Nso,Nso,Lmats],FROUTINE,"Gmats")
  !
  Smats_(:,1,1,:,:,:) = Smats(:,:,:,:)
  call dmft_get_gloc_matsubara_normal_lattice_main_mpi(MpiComm,Hk,Wtk,Gmats_,Smats_,iprint,tridiag_,mpi_split_,hk_symm_)
  Gmats(:,:,:,:) = Gmats_(:,1,1,:,:,:)
end subroutine dmft_get_gloc_matsubara_normal_lattice_Nband_mpi
#undef FROUTINE



