module GK_REALAXIS
  USE DMFT_CTRL_VARS
  USE GK_COMMON
  implicit none
  private
  interface dmft_gk_realaxis
     module procedure :: dmft_get_gk_realaxis_normal_main
     module procedure :: dmft_get_gk_realaxis_normal_dos
     module procedure :: dmft_get_gk_realaxis_normal_ineq
     module procedure :: dmft_get_gk_realaxis_superc_main
     module procedure :: dmft_get_gk_realaxis_superc_dos
     module procedure :: dmft_get_gk_realaxis_superc_ineq
  end interface dmft_gk_realaxis

  public :: dmft_gk_realaxis

contains



  subroutine dmft_get_gk_realaxis_normal_main(Hk,Gkreal,Sreal)
    complex(8),dimension(:,:),intent(in)          :: Hk        ![Nspin*Norb][Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal] 
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Hk,[Nso,Nso],'dmft_get_gk_realaxis_normal_main_mpi',"Hk")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_main_mpi',"Gkreal")
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
    Gkreal=zero
    call invert_gk_normal(zeta_real,Hk,.false.,Gkreal)      
  end subroutine dmft_get_gk_realaxis_normal_main


  subroutine dmft_get_gk_realaxis_normal_dos(Ebands,Dbands,Hloc,Gkreal,Sreal)
    real(8),dimension(:),intent(in)               :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands)),intent(in)    :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !
    complex(8),dimension(:,:,:),allocatable       :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),allocatable   :: Gtmp      !as Sreal
    !
    real(8)                                       :: beta
    real(8)                                       :: xmu,eps
    real(8)                                       :: wini,wfin
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !Retrieve parameters:
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Nso   = Nspin*Norb    
    !Testing part:
    if(size(Ebands)/=Nso)stop "dmft_get_gk_realaxis_normal_dos_main_mpi ERROR: size(Ebands)!=Nso"
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_dos_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_dos_main_mpi',"Gkreal")
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
    Gkreal=zero
    allocate(Gtmp(Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
    !
    do i = 1+mpi_rank, Lreal, mpi_size
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             Gtmp(ispin,ispin,iorb,iorb,i) = Dbands(io)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io) )
          enddo
       enddo
    end do
#ifdef _MPI    
    if(check_MPI())then
       Gkreal=zero
       call AllReduce_MPI(MPI_COMM_WORLD,Gtmp,Gkreal)
    else
       Gkreal=Gtmp
    endif
#else
    Gkreal=Gtmp
#endif
  end subroutine dmft_get_gk_realaxis_normal_dos



  subroutine dmft_get_gk_realaxis_normal_ineq(Hk,Gkreal,Sreal,tridiag)
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    logical,optional                                :: tridiag
    logical                                         :: tridiag_
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:),allocatable       :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    real(8)                                         :: wini,wfin  
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    tridiag_=.false.;if(present(tridiag))tridiag_=tridiag
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Lreal = size(Sreal,6)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hk,[Nlso,Nlso],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Sreal")
    call assert_shape(Gkreal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_normal_ineq_main_mpi',"Gkreal")
    !
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
    Gkreal=zero
    if(.not.tridiag_)then
       call invert_gk_normal_ineq(zeta_real,Hk,.false.,Gkreal)
    else
       call invert_gk_normal_tridiag(zeta_real,Hk,.false.,Gkreal)
    endif
  end subroutine dmft_get_gk_realaxis_normal_ineq






















  subroutine dmft_get_gk_realaxis_superc_main(Hk,Gkreal,Sreal)
    complex(8),dimension(:,:,:),intent(in)          :: Hk        ![2][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
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
    call assert_shape(Hk,[2,Nso,Nso],'dmft_get_gk_realaxis_superc_main_mpi',"Hk")
    call assert_shape(Sreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main_mpi',"Sreal")
    call assert_shape(Gkreal,[2,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_main_mpi',"Gkreal")
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
    real(8),dimension(:,:),intent(in)               :: Ebands    ![2][Nspin*Norb]
    real(8),dimension(size(Ebands,2)),intent(in)    :: Dbands    ![Nspin*Norb]
    real(8),dimension(2,size(Ebands,2)),intent(in)  :: Hloc      ![2][Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !allocatable arrays
    complex(8)                                      :: gktmp(2),cdet
    complex(8)                                      :: zeta_11,zeta_12,zeta_22 
    complex(8),dimension(:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nspin*Norb][Nspin*Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gtmp     ![2][Nspin][Nspin][Norb][Norb][Lreal]
    !
    real(8)                                         :: beta
    real(8)                                         :: xmu,eps
    !
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
    allocate(Gtmp(2,Nspin,Nspin,Norb,Norb,Lreal));Gtmp=zero
    !
    do i = 1+mpi_rank, Lmats, mpi_size
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
             Gtmp(1,ispin,ispin,iorb,iorb,i) = Gtmp(1,ispin,ispin,iorb,iorb,i) + gktmp(1)*Dbands(io)
             Gtmp(2,ispin,ispin,iorb,iorb,i) = Gtmp(2,ispin,ispin,iorb,iorb,i) + gktmp(2)*Dbands(io)
          enddo
       enddo
    enddo
#ifdef _MPI    
    if(check_MPI())then
       Gkreal=zero
       call AllReduce_MPI(MPI_COMM_WORLD,Gtmp,Gkreal)
    else
       Gkreal=Gtmp
    endif
#else
    Gkreal=Gtmp
#endif
  end subroutine dmft_get_gk_realaxis_superc_dos


  subroutine dmft_get_gk_realaxis_superc_ineq(Hk,Gkreal,Sreal)
    complex(8),dimension(:,:,:),intent(in)            :: Hk        ![2][Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: Sreal     ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gkreal     !as Sreal
    !allocatable arrays
    complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: Gtmp    !as Sreal
    complex(8),dimension(:,:,:,:,:,:),allocatable     :: zeta_real ![2][2][Nlat][Nspin*Norb][Nspin*Norb][Lreal]
    !
    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
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
    call assert_shape(Hk,[2,Nlso,Nlso],'dmft_get_gk_realaxis_superc_ineq_main_mpi',"Hk")
    call assert_shape(Sreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_ineq_main_mpi',"Sreal")
    call assert_shape(Gkreal,[2,Nlat,Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gk_realaxis_superc_ineq_main_mpi',"Gkreal")
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















end module GK_REALAXIS
