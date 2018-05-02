MODULE DMFT_EKIN
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS,   only:free_unit,str,reg
  USE SF_ARRAYS,    only:arange
  USE SF_TIMER
  USE SF_LINALG,    only:inv,eye,eigh,diag
  USE SF_MISC,      only:assert_shape
  USE SF_SPECIAL,   only:fermi
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private


  interface dmft_kinetic_energy
     module procedure :: dmft_kinetic_energy_normal_main
     module procedure :: dmft_kinetic_energy_normal_dos
     module procedure :: dmft_kinetic_energy_normal_lattice
     module procedure :: dmft_kinetic_energy_superc_main
     module procedure :: dmft_kinetic_energy_superc_dos
     module procedure :: dmft_kinetic_energy_superc_lattice
#ifdef _MPI
     module procedure :: dmft_kinetic_energy_normal_main_mpi
     module procedure :: dmft_kinetic_energy_normal_dos_mpi
     module procedure :: dmft_kinetic_energy_normal_lattice_mpi
     module procedure :: dmft_kinetic_energy_superc_main_mpi
     module procedure :: dmft_kinetic_energy_superc_dos_mpi
     module procedure :: dmft_kinetic_energy_superc_lattice_mpi
#endif
     !ADDITIONAL INTERFACES:
     !1-Band ([Nlat])[L] for Sigma & Self
     module procedure :: dmft_kinetic_energy_normal_1band
     module procedure :: dmft_kinetic_energy_normal_1band_lattice
     module procedure :: dmft_kinetic_energy_superc_1band
     module procedure :: dmft_kinetic_energy_superc_1band_lattice
#ifdef _MPI
     module procedure :: dmft_kinetic_energy_normal_1band_mpi
     module procedure :: dmft_kinetic_energy_normal_1band_lattice_mpi
     module procedure :: dmft_kinetic_energy_superc_1band_mpi
     module procedure :: dmft_kinetic_energy_superc_1band_lattice_mpi
#endif
     !NN-Shape ([Nlat])[Nspin][Nspin][Norb][Norb][L] for Sigma & Self
     module procedure :: dmft_kinetic_energy_normal_NN
     module procedure :: dmft_kinetic_energy_normal_NN_dos
     module procedure :: dmft_kinetic_energy_normal_NN_lattice
     module procedure :: dmft_kinetic_energy_superc_NN
     module procedure :: dmft_kinetic_energy_superc_NN_dos
     module procedure :: dmft_kinetic_energy_superc_NN_lattice
#ifdef _MPI
     module procedure :: dmft_kinetic_energy_normal_NN_mpi
     module procedure :: dmft_kinetic_energy_normal_NN_dos_mpi
     module procedure :: dmft_kinetic_energy_normal_NN_lattice_mpi
     module procedure :: dmft_kinetic_energy_superc_NN_mpi
     module procedure :: dmft_kinetic_energy_superc_NN_dos_mpi
     module procedure :: dmft_kinetic_energy_superc_NN_lattice_mpi
#endif
  end interface dmft_kinetic_energy


  interface print_Hloc
     module procedure print_Hloc_2
     module procedure print_Hloc_4
  end interface print_Hloc



  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block


  real(8),dimension(:),allocatable        :: wm



  !PUBLIC in DMFT
  public :: dmft_kinetic_energy

contains 




  !----------------------------------------------------------------------------------------
  !PURPOSE: Evaluate the Kinetic energy of the lattice model. (NORMAL / SUPERC)
  ! input: 
  ! 1. the Hamiltonian matrix H(k) [& -conjg(H(-k))]
  ! 2. DMFT normal&anomalous self-energy Sigma [& Self].
  !
  ! The main routines accept self-energy as:
  ! - Sigma & Self: [Nspin*Norb][Nspin*Norb][L]
  ! - Sigma & Self: [Nlat][Nspin*Norb][Nspin*Norb][L]
  !----------------------------------------------------------------------------------------
  include "dmft_ekin_normal.f90"
#ifdef _MPI
  include "dmft_ekin_normal_mpi.f90"
#endif

  include "dmft_ekin_superc.f90"
#ifdef _MPI
  include "dmft_ekin_superc_mpi.f90"
#endif














  !####################################################################
  !                    ADDITIONAL INTERFACES 
  !####################################################################
  !********************************************************************
  !                           1 BAND
  ! - Sigma & Self shape:       [L] (no Nspin[=1], no Norb[=1])
  ! - Sigma & Self shape: [Nlat][L] (no Nspin[=1], no Norb[=1])
  !********************************************************************
  subroutine dmft_kinetic_energy_normal_1band(Hk,Wtk,Sigma,Ekin,Eloc)
    complex(8),dimension(:)                   :: Hk
    real(8),dimension(size(Hk))               :: Wtk
    complex(8),dimension(:)                   :: Sigma
    real(8),optional                          :: Ekin,Eloc
    complex(8),dimension(1,1,size(Sigma))     :: Sigma_
    complex(8),dimension(1,1,size(Hk))        :: Hk_
    real(8),dimension(1)                      :: Ekin_,Eloc_
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call dmft_kinetic_energy_normal_main(Hk_,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_(1)
    if(present(Eloc))Eloc=Eloc_(1)
  end subroutine dmft_kinetic_energy_normal_1band

  subroutine dmft_kinetic_energy_normal_1band_lattice(Hk,Wtk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)                           :: Hk     ! [Nlat*1][Nlat*1][Nk] !Nspin*Norb=1
    real(8),dimension(size(Hk,3))                         :: wtk    ! [Nk]
    complex(8),dimension(:,:)                             :: Sigma  ! [Nlat][L]
    real(8),dimension(size(Hk,1)),optional                :: Ekin,Eloc
    complex(8),dimension(size(Sigma,1),1,1,size(Sigma,2)) :: Sigma_ ! [Nlat][1][1][L]
    integer                                               :: ilat
    integer                                               :: Nlat,Nlso,Nso,Lk
    real(8),dimension(size(Hk,1))                         :: Ekin_,Eloc_
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_normal_1band_lattice","Hk") 
    Sigma_(:,1,1,:) = Sigma(:,:)
    call dmft_kinetic_energy_normal_lattice(Hk,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_1band_lattice


#ifdef _MPI
  subroutine dmft_kinetic_energy_normal_1band_mpi(MpiComm,Hk,Wtk,Sigma,Ekin,Eloc)
    integer                                   :: MpiComm
    complex(8),dimension(:)                   :: Hk
    real(8),dimension(size(Hk))               :: Wtk
    complex(8),dimension(:)                   :: Sigma
    real(8),optional                          :: Ekin,Eloc
    complex(8),dimension(1,1,size(Sigma))     :: Sigma_
    complex(8),dimension(1,1,size(Hk))        :: Hk_
    real(8),dimension(1)                      :: Ekin_,Eloc_
    Sigma_(1,1,:) = Sigma
    Hk_(1,1,:)    = Hk
    call dmft_kinetic_energy_normal_main_mpi(MpiComm,Hk_,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_(1)
    if(present(Eloc))Eloc=Eloc_(1)
  end subroutine dmft_kinetic_energy_normal_1band_mpi

  subroutine dmft_kinetic_energy_normal_1band_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Ekin,Eloc)
    integer                                               :: MpiComm
    complex(8),dimension(:,:,:)                           :: Hk     ! [Nlat*1][Nlat*1][Nk]
    real(8),dimension(size(Hk,3))                         :: wtk    ! [Nk]
    complex(8),dimension(:,:)                             :: Sigma  ! [Nlat][L]
    real(8),dimension(size(Hk,1)),optional                :: Ekin,Eloc
    complex(8),dimension(size(Sigma,1),1,1,size(Sigma,2)) :: Sigma_ ! [Nlat][1][1][L]
    integer                                               :: ilat
    integer                                               :: Nlat,Nlso,Nso,Lk
    real(8),dimension(size(Hk,1))                         :: Ekin_,Eloc_
    Nlat = size(Sigma,1)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_normal_1band_lattice_mpi","Hk") !implictly test Nlat*1*=Nlso
    Sigma_(:,1,1,:) = Sigma(:,:)
    call dmft_kinetic_energy_normal_lattice_mpi(MpiComm,Hk,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_1band_lattice_mpi
#endif








  subroutine dmft_kinetic_energy_superc_1band(Hk,Wtk,Sigma,Self,Ekin,Eloc)
    real(8),dimension(:)                  :: Hk    ![Lk]
    real(8),dimension(size(Hk))           :: Wtk   ![Lk]
    complex(8),dimension(:)               :: Sigma ![L]
    complex(8),dimension(size(Sigma))     :: Self  ![L]
    real(8),optional                      :: Ekin,Eloc
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Sigma)) :: Self_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    real(8),dimension(1)                  :: Ekin_,Eloc_
    Sigma_(1,1,:)  = Sigma(:)
    Self_(1,1,:)   = Self(:)
    Hk_(1,1,:)     = Hk
    call dmft_kinetic_energy_superc_main(Hk_,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_(1)
    if(present(Eloc))Eloc=Eloc_(1)
  end subroutine dmft_kinetic_energy_superc_1band


  subroutine dmft_kinetic_energy_superc_1band_lattice(Hk,Wtk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)                           :: Hk     ! [Nlat*1][Nlat*1][Nk]
    real(8),dimension(size(Hk,3))                         :: wtk    ! [Nk]
    complex(8),dimension(:,:)                             :: Sigma  ! [Nlat][L]
    complex(8),dimension(size(Sigma,1),size(Sigma,2))     :: Self   ! [Nlat][L]
    real(8),dimension(size(Hk,1)),optional                :: Ekin,Eloc
    complex(8),dimension(size(Sigma,1),1,1,size(Sigma,2)) :: Sigma_ ! [Nlat][1][1][L]
    complex(8),dimension(size(Sigma,1),1,1,size(Sigma,2)) :: Self_ ! [Nlat][1][1][L]
    integer                                               :: i,iorb,ilat,ispin,io,is
    integer                                               :: j,jorb,jlat,jspin,jo,js
    integer                                               :: Nlat,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1))                         :: Ekin_,Eloc_
    Nlat = size(Sigma,1)
    Liw  = size(Sigma,2)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_superc_1band_lattice","Hk") !implictly test Nlat*1*1=Nlso
    Sigma_(:,1,1,:)  = Sigma
    Self_(:,1,1,:) = Self
    call dmft_kinetic_energy_superc_lattice(Hk,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_1band_lattice


#ifdef _MPI
  subroutine dmft_kinetic_energy_superc_1band_mpi(MpiComm,Hk,Wtk,Sigma,Self,Ekin,Eloc)
    integer                               :: MpiComm
    real(8),dimension(:)                  :: Hk
    real(8),dimension(size(Hk))           :: Wtk
    complex(8),dimension(:)               :: Sigma
    complex(8),dimension(size(Sigma))     :: Self
    real(8),optional                      :: Ekin,Eloc
    complex(8),dimension(1,1,size(Sigma)) :: Sigma_
    complex(8),dimension(1,1,size(Sigma)) :: Self_
    complex(8),dimension(1,1,size(Hk))    :: Hk_
    real(8),dimension(1)                  :: Ekin_,Eloc_
    Sigma_(1,1,:)  = Sigma(:)
    Self_(1,1,:)   = Self(:)
    Hk_(1,1,:)     = Hk
    call dmft_kinetic_energy_superc_main_mpi(MpiComm,Hk_,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_(1)
    if(present(Eloc))Eloc=Eloc_(1)
  end subroutine dmft_kinetic_energy_superc_1band_mpi

  subroutine dmft_kinetic_energy_superc_1band_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Self,Ekin,Eloc)
    integer                                               :: MpiComm
    complex(8),dimension(:,:,:)                           :: Hk     ! [Nlat*1][Nlat*1][Nk]
    real(8),dimension(size(Hk,3))                         :: wtk    ! [Nk]
    complex(8),dimension(:,:)                             :: Sigma  ! [Nlat][L]
    complex(8),dimension(size(Sigma,1),size(Sigma,2))     :: Self ! [Nlat][L]  
    real(8),dimension(size(Hk,1)),optional                :: Ekin,Eloc
    complex(8),dimension(size(Sigma,1),1,1,size(Sigma,2)) :: Sigma_ ! [Nlat][1][1][L]
    complex(8),dimension(size(Sigma,1),1,1,size(Sigma,2)) :: Self_ ! [Nlat][1][1][L]
    integer                                               :: i,iorb,ilat,ispin,io,is
    integer                                               :: j,jorb,jlat,jspin,jo,js
    integer                                               :: Nlat,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1))                         :: Ekin_,Eloc_
    Nlat = size(Sigma,1)
    Liw  = size(Sigma,2)
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    call assert_shape(Hk,[Nlat,Nlso,Lk],"dmft_kinetic_energy_superc_1band_lattice_mpi","Hk") !implictly test Nlat*1*1=Nlso
    Sigma_(:,1,1,:)  = Sigma
    Self_(:,1,1,:) = Self
    call dmft_kinetic_energy_superc_lattice_mpi(MpiComm,Hk,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_1band_lattice_mpi
#endif





  !********************************************************************
  !                             NN FORM
  ! - Sigma & Self shape:       [Nspin][Nspin][Norb][Norb][L]
  ! - Sigma & Self shape: [Nlat][Nspin][Nspin][Norb][Norb][L]
  !********************************************************************
  subroutine dmft_kinetic_energy_normal_NN(Hk,Wtk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Hk,1)),optional  :: Ekin,Eloc
    real(8),dimension(size(Hk,1))           :: Ekin_,Eloc_
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN","Sigma")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_normal_NN","Hk")
    allocate(Sigma_(Nso,Nso,Lmats))
    Sigma_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_main(Hk,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN

  subroutine dmft_kinetic_energy_normal_NN_dos(Ebands,Dbands,Hloc,Sigma,Ekin,Eloc)
    real(8),dimension(:,:)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1))                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:)                  :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                          :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Ebands,1)),optional       :: Ekin,Eloc
    real(8),dimension(size(Ebands,1))                :: Ekin_,Eloc_
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN","Sigma")
    allocate(Sigma_(Nso,Nso,Lmats))
    Sigma_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_dos(Ebands,Dbands,Hloc,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_dos

  subroutine dmft_kinetic_energy_normal_NN_lattice(Hk,Wtk,Sigma,Ekin,Eloc)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1)),optional    :: Ekin,Eloc
    real(8),dimension(size(Hk,1))             :: Ekin_,Eloc_
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"dmft_kinetic_energy_normal_NN_lattice","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"dmft_kinetic_energy_normal_NN_lattice","Sigma")
    allocate(Sigma_(Nlat,Nso,Nso,Liw))
    Sigma_ = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Sigma_(:,is,js,:) = Sigma(:,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_lattice(Hk,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_lattice


#ifdef _MPI
  subroutine dmft_kinetic_energy_normal_NN_mpi(MpiComm,Hk,Wtk,Sigma,Ekin,Eloc)
    integer                                 :: MpiComm
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Hk,1)),optional  :: Ekin,Eloc
    real(8),dimension(size(Hk,1))           :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN_mpi","Sigma")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_normal_NN_mpi","Hk")
    allocate(Sigma_(Nso,Nso,Lmats))
    Sigma_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_main_mpi(MpiComm,Hk,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_mpi

  subroutine dmft_kinetic_energy_normal_NN_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Sigma,Ekin,Eloc)
    integer                                          :: MpiComm
    real(8),dimension(:,:)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1))                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:)                  :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    integer                                          :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Ebands,1)),optional       :: Ekin,Eloc
    real(8),dimension(size(Ebands,1))                :: Ekin_,Eloc_
    !
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_normal_NN","Sigma")
    allocate(Sigma_(Nso,Nso,Lmats))
    Sigma_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) = Sigma(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_dos_mpi

  subroutine dmft_kinetic_energy_normal_NN_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Ekin,Eloc)
    integer                                   :: MpiComm
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1)),optional    :: Ekin,Eloc
    real(8),dimension(size(Hk,1))             :: Ekin_,Eloc_
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlso,Nlso,Lk],"dmft_kinetic_energy_normal_NN_lattice_mpi","Hk")
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"dmft_kinetic_energy_normal_NN_lattice_mpi","Sigma")
    allocate(Sigma_(Nlat,Nso,Nso,Liw))
    Sigma_ = zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Sigma_(:,is,js,:) = Sigma(:,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_normal_lattice_mpi(MpiComm,Hk,Wtk,Sigma_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_normal_NN_lattice_mpi
#endif








  subroutine dmft_kinetic_energy_superc_NN(Hk,Wtk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)         :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Hk,1)),optional  :: Ekin,Eloc
    real(8),dimension(size(Hk,1))           :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Self")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_superc_NN","Hk")
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) =  Sigma(ispin,jspin,iorb,jorb,:)
                Self_(is,js,:)= Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_main(Hk,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN

  subroutine dmft_kinetic_energy_superc_NN_dos(Ebands,Dbands,Hloc,Sigma,Self,Ekin,Eloc)
    real(8),dimension(:,:)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1))                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:)                  :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)                  :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                          :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Ebands,1)),optional       :: Ekin,Eloc
    real(8),dimension(size(Ebands,1))                :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN","Self")
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) =  Sigma(ispin,jspin,iorb,jorb,:)
                Self_(is,js,:)= Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_dos(Ebands,Dbands,Hloc,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_dos

  subroutine dmft_kinetic_energy_superc_NN_lattice(Hk,Wtk,Sigma,Self,Ekin,Eloc)
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)         :: Self   ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Self_  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                   :: i,iorb,ilat,ispin,is
    integer                                   :: j,jorb,jlat,jspin,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1)),optional    :: Ekin,Eloc
    real(8),dimension(size(Hk,1))             :: Ekin_,Eloc_
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    !Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"kinetic_energy_lattice_superc_2","Hk") !implictly check Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Sigma")
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Self")
    allocate(Sigma_(Nlat,Nso,Nso,Liw))
    allocate(Self_(Nlat,Nso,Nso,Liw))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(:,is,js,:) =  Sigma(:,ispin,jspin,iorb,jorb,:)
                Self_(:,is,js,:)= Self(:,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_lattice(Hk,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_lattice

#ifdef _MPI
  subroutine dmft_kinetic_energy_superc_NN_mpi(MpiComm,Hk,Wtk,Sigma,Self,Ekin,Eloc)
    integer                                 :: MpiComm
    complex(8),dimension(:,:,:)             :: Hk     ![Nspin*Norb][Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))           :: Wtk    ![Nk]
    complex(8),dimension(:,:,:,:,:)         :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)         :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                 :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Hk,1)),optional  :: Ekin,Eloc
    real(8),dimension(size(Hk,1))           :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN_mpi","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN_mpi","Self")
    call assert_shape(Hk,[Nso,Nso,size(Hk,3)],"dmft_kinetic_energy_superc_NN_mpi","Hk")
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) =  Sigma(ispin,jspin,iorb,jorb,:)
                Self_(is,js,:)= Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_main_mpi(MpiComm,Hk,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_mpi

  subroutine dmft_kinetic_energy_superc_NN_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Sigma,Self,Ekin,Eloc)
    integer                                          :: MpiComm
    real(8),dimension(:,:)                           :: Ebands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)) :: Dbands  ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1))                :: Hloc    ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:)                  :: Sigma  ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:)                  :: Self   ![Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Sigma_ ![Nspin*Norb][Nspin*Norb][L]
    complex(8),dimension(:,:,:),allocatable          :: Self_  ![Nspin*Norb][Nspin*Norb][L]
    integer                                          :: Nspin,Norb,Lmats,Nso,i,iorb,jorb,ispin,jspin,is,js
    real(8),dimension(size(Ebands,1)),optional       :: Ekin,Eloc
    real(8),dimension(size(Ebands,1))                :: Ekin_,Eloc_
    Nspin = size(Sigma,1)
    Norb  = size(Sigma,3)
    Lmats = size(Sigma,5)
    Nso   = Nspin*Norb
    call assert_shape(Sigma,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN_mpi","Sigma")
    call assert_shape(Self,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_kinetic_energy_superc_NN_mpi","Self")
    allocate(Sigma_(Nso,Nso,Lmats))
    allocate(Self_(Nso,Nso,Lmats))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(is,js,:) =  Sigma(ispin,jspin,iorb,jorb,:)
                Self_(is,js,:)= Self(ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_dos_mpi(MpiComm,Ebands,Dbands,Hloc,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_dos_mpi

  subroutine dmft_kinetic_energy_superc_NN_lattice_mpi(MpiComm,Hk,Wtk,Sigma,Self,Ekin,Eloc)
    integer                                   :: MpiComm
    complex(8),dimension(:,:,:)               :: Hk     ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk]
    real(8),dimension(size(Hk,3))             :: wtk    ! [Nk]
    complex(8),dimension(:,:,:,:,:,:)         :: Sigma  ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:,:,:)         :: Self   ! [Nlat][Nspin][Nspin][Norb][Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Sigma_ ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    complex(8),dimension(:,:,:,:),allocatable :: Self_  ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                                   :: i,iorb,ilat,ispin,is
    integer                                   :: j,jorb,jlat,jspin,js
    integer                                   :: Nlat,Nspin,Norb,Nlso,Nso,Lk,Liw
    real(8),dimension(size(Hk,1)),optional    :: Ekin,Eloc
    real(8),dimension(size(Hk,1))             :: Ekin_,Eloc_
    !Get generalized Lattice-Spin-Orbital index
    Nlat = size(Sigma,1)
    Nspin= size(Sigma,2)
    Norb = size(Sigma,4)
    Nso  = Nspin*Norb
    Nlso = size(Hk,1)
    Liw  = size(Sigma,6)
    Lk   = size(Hk,3)
    !Nlso = Nlat*Nspin*Norb
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"kinetic_energy_lattice_superc_2","Hk") !implictly check Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Sigma")
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Liw],"kinetic_energy_lattice_superc_2","Self")
    allocate(Sigma_(Nlat,Nso,Nso,Liw))
    allocate(Self_(Nlat,Nso,Nso,Liw))
    Sigma_=zero
    Self_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Sigma_(:,is,js,:) =  Sigma(:,ispin,jspin,iorb,jorb,:)
                Self_(:,is,js,:)= Self(:,ispin,jspin,iorb,jorb,:)
             enddo
          enddo
       enddo
    enddo
    call dmft_kinetic_energy_superc_lattice_mpi(MpiComm,Hk,Wtk,Sigma_,Self_,Ekin_,Eloc_)
    if(present(Ekin))Ekin=Ekin_
    if(present(Eloc))Eloc=Eloc_
  end subroutine dmft_kinetic_energy_superc_NN_lattice_mpi
#endif








  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8)                    :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !+-----------------------------------------------------------------------------+!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix,Nlat,Nso) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !+-----------------------------------------------------------------------------+!
  function select_block_Nlso(ip,Matrix,Nlat,Nso) result(Vblock)
    integer                                 :: ip
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nso,Nso)           :: Vblock
    integer                                 :: i,j
    Vblock=zero
    i = 1+(ip-1)*Nso
    j =       ip*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ip,Matrix,Nlat,Nspin,Norb) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block_nnn



  !+------------------------------------------------------------------+
  !PURPOSE: Print the local part of the Non-interacting Hamiltonian
  !+------------------------------------------------------------------+
  subroutine print_Hloc_2(hloc,file)
    complex(8),dimension(:,:) :: hloc
    character(len=*),optional :: file
    integer                   :: Ni,Nj,iorb,jorb,unit
    unit=6;
    if(present(file))then
       write(unit,"(A)")"print_Hloc on file :"//reg(file)
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Ni = size(hloc,1)
    Nj = size(hloc,2)
    if(present(file))then
       do iorb=1,Ni
          write(unit,"(10000F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       do iorb=1,Ni
          write(unit,"(10000F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       close(unit)
    else
       do iorb=1,Ni
          write(unit,"(10000(A1,F7.3,A1,F7.3,A1,2x))")&
               ('(',dreal(Hloc(iorb,jorb)),',',dimag(Hloc(iorb,jorb)),')',jorb =1,Nj)
       enddo
       write(unit,*)""
    endif
  end subroutine print_Hloc_2
  !
  subroutine print_Hloc_4(hloc,file)
    complex(8),dimension(:,:,:,:) :: hloc ![Nspin][Nspin][Norb][Norb]
    character(len=*),optional     :: file
    integer                       :: Nspin,Norb
    integer                       :: iorb,jorb,ispin,jspin,unit
    unit=6;
    if(present(file))then
       write(unit,"(A)")"print_Hloc on file :"//reg(file)
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Nspin = size(hloc,1)
    Norb  = size(hloc,3)
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"print_Hloc","Hloc")
    if(present(file))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(10000F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(10000F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    else
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(10000(A1,F7.3,A1,F7.3,A1,2x))")&
                  (&
                  (&
                  '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
                  jorb =1,Norb),&
                  jspin=1,Nspin)
          enddo
       enddo
    endif
  end subroutine print_Hloc_4





  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_info()
    !to be removed
  end subroutine write_kinetic_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : Write energies to file
  !+-------------------------------------------------------------------+
  subroutine write_kinetic_value(Ekin,Eloc,Nlat,Nso)
    real(8),dimension(:)               :: Ekin
    real(8),dimension(size(Ekin))      :: Eloc
    real(8),dimension(:,:),allocatable :: Ekin_,Eloc_
    integer,optional                   :: Nlat,Nso
    integer                            :: Nlso
    integer                            :: i,iso,ilat,unit
    

    if(.not.present(Nlat))then
       !
       Nso = size(Ekin)
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nso),&
            (str(2+Nso+i)//"<Eloc"//str(i)//">",i=1,Nso)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.dat")
       write(unit,"(90F15.9)")sum(Ekin),sum(Eloc),(Ekin(i),i=1,Nso),(Eloc(i),i=1,Nso)
       close(unit)
       !
    else
       !
       if(.not.present(Nso))stop "ERROR write_kinetic_value: Nlat present but Nso not present."
       !
       Nlso = size(Ekin)
       if(Nlso /= Nlat*Nso)stop "Error write_kinetic_value: Nlso != Nlat*Nso" 
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nso),&
            (str(2+Nso+i)//"<Eloc"//str(i)//">",i=1,Nso)
       close(unit)
       !
       !
       allocate(Ekin_(Nlat,Nso))
       allocate(Eloc_(Nlat,Nso))
       do ilat=1,Nlat
          do iso=1,Nso
             i = iso + (ilat-1)*Nso
             Ekin_(ilat,iso) = Ekin(i)
             Eloc_(ilat,iso) = Eloc(i)
          enddo
       enddo
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.dat")       
       write(unit,"(90F15.9)")sum(Ekin_)/Nlat,sum(Eloc_)/Nlat
       do ilat=1,Nlat
          write(unit,"(100000F15.9)")sum(Ekin_(ilat,:)),sum(Eloc_(ilat,:)),&
               (Ekin_(ilat,i),i=1,Nso),&
               (Eloc_(ilat,i),i=1,Nso)
       enddo
       close(unit)
    endif
  end subroutine write_kinetic_value









#ifdef _MPI
  function MPI_Get_size(comm) result(size)
    integer :: comm
    integer :: size,ierr
    call MPI_Comm_size(comm,size,ierr)
  end function MPI_Get_size



  function MPI_Get_rank(comm) result(rank)
    integer :: comm
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
  end function MPI_Get_rank



  function MPI_Get_master(comm) result(master)
    integer :: comm
    logical :: master
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    master=.false.
    if(rank==0)master=.true.
  end function MPI_Get_master
#endif



end MODULE DMFT_EKIN


