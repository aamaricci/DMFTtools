module DMFT_WEISS_FIELD
  USE SF_TIMER
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only:reg,txtfy
  USE SF_ARRAYS,    only:arange
  USE SF_LINALG,    only:eye,inv
  USE SF_MISC,      only:assert_shape
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private



  interface dmft_self_consistency
     module procedure :: dmft_sc_normal_main
     module procedure :: dmft_sc_normal_ineq
     module procedure :: dmft_sc_normal_bethe
     module procedure :: dmft_sc_normal_bethe_ineq
     module procedure :: dmft_sc_superc_main
     module procedure :: dmft_sc_superc_ineq
#ifdef _MPI
     module procedure :: dmft_sc_normal_main_mpi
     module procedure :: dmft_sc_normal_ineq_mpi
     module procedure :: dmft_sc_normal_bethe_mpi
     module procedure :: dmft_sc_normal_bethe_ineq_mpi
     module procedure :: dmft_sc_superc_main_mpi
     module procedure :: dmft_sc_superc_ineq_mpi
#endif
  end interface dmft_self_consistency


  interface dmft_weiss
     module procedure :: dmft_get_weiss_normal_main
     module procedure :: dmft_get_weiss_normal_ineq
     module procedure :: dmft_get_weiss_normal_bethe
     module procedure :: dmft_get_weiss_normal_bethe_ineq
     module procedure :: dmft_get_weiss_superc_main
     module procedure :: dmft_get_weiss_superc_ineq
#ifdef _MPI
     module procedure :: dmft_get_weiss_normal_main_mpi
     module procedure :: dmft_get_weiss_normal_ineq_mpi
     module procedure :: dmft_get_weiss_normal_bethe_mpi
     module procedure :: dmft_get_weiss_normal_bethe_ineq_mpi
     module procedure :: dmft_get_weiss_superc_main_mpi
     module procedure :: dmft_get_weiss_superc_ineq_mpi
#endif
  end interface dmft_weiss


  interface dmft_delta
     module procedure :: dmft_get_delta_normal_main
     module procedure :: dmft_get_delta_normal_ineq
     module procedure :: dmft_get_delta_normal_bethe
     module procedure :: dmft_get_delta_normal_bethe_ineq
     module procedure :: dmft_get_delta_superc_main
     module procedure :: dmft_get_delta_superc_ineq
#ifdef _MPI
     module procedure :: dmft_get_delta_normal_main_mpi
     module procedure :: dmft_get_delta_normal_ineq_mpi
     module procedure :: dmft_get_delta_normal_bethe_mpi
     module procedure :: dmft_get_delta_normal_bethe_ineq_mpi
     module procedure :: dmft_get_delta_superc_main_mpi
     module procedure :: dmft_get_delta_superc_ineq_mpi
#endif
  end interface dmft_delta


  public :: dmft_self_consistency
  public :: dmft_weiss
  public :: dmft_delta



  !##################################################################
  !##################################################################
  !##################################################################



  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block


  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape

  real(8),dimension(:),allocatable :: wm !Matsubara frequencies
  character(len=128)               :: suffix
  character(len=128)               :: weiss_suffix=".dat"
  integer                          :: Nlat,Nspin,Norb
  integer                          :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
  !
  integer                          :: mpi_ierr
  integer                          :: mpi_rank
  integer                          :: mpi_size
  logical                          :: mpi_master
  !
  real(8)                          :: beta
  real(8)                          :: xmu


contains


  !--------------------------------------------------------------------!
  !PURPOSE: Perform the DMFT local self-consistency using Weiss or Delta
  ! equations and given G_loc/F_loc and Sigma/Self
  ! INPUT:
  ! 1. GLOC/FLOC  : ([Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 2. Sigma/SELF : ([Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 4. Hloc       : local part of the non-interacting Hamiltonian
  !--------------------------------------------------------------------!
  include "dmft_self_consistency.f90" 
#ifdef _MPI
  include "dmft_self_consistency_mpi.f90"
#endif



  !--------------------------------------------------------------------!
  !PURPOSE: Get the local Weiss Field calG0 or using self-consistency 
  ! equations and given G_loc/F_loc and Sigma/Self
  ! INPUT:
  ! 1. GLOC/FLOC  : ([Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 2. Sigma/SELF : ([Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 4. Hloc       : local part of the non-interacting Hamiltonian
  !--------------------------------------------------------------------!
  !                    WEISS FIELD - NORMAL
  !--------------------------------------------------------------------!
  include "dmft_weiss_normal.f90"
  include "dmft_weiss_bethe.f90"
#ifdef _MPI
  include "dmft_weiss_normal_mpi.f90"
  include "dmft_weiss_bethe_mpi.f90"
#endif


  !--------------------------------------------------------------------!
  !                    WEISS FIELD - SUPERC
  !--------------------------------------------------------------------!
  include "dmft_weiss_superc.f90"
#ifdef _MPI
  include "dmft_weiss_superc_mpi.f90"
#endif


  !##################################################################
  !##################################################################
  !##################################################################







  !--------------------------------------------------------------------!
  !PURPOSE: Get the local Hybridization \Delta functino using 
  ! equations and given G_loc/F_loc and Sigma/Self
  ! INPUT:
  ! 1. GLOC/FLOC  : ([Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 2. Sigma/SELF : ([Nlat])[Nspin][Nspin][Norb][Norb][Lmats]
  ! 4. Hloc       : local part of the non-interacting Hamiltonian
  !--------------------------------------------------------------------!
  !                    WEISS FIELD - NORMAL
  !--------------------------------------------------------------------!
  include "dmft_delta_normal.f90"
  include "dmft_delta_bethe.f90"
#ifdef _MPI
  include "dmft_delta_normal_mpi.f90"
  include "dmft_delta_bethe_mpi.f90"
#endif


  !--------------------------------------------------------------------!
  !                    WEISS FIELD - SUPERC
  !--------------------------------------------------------------------!
  include "dmft_delta_superc.f90"
#ifdef _MPI
  include "dmft_delta_superc_mpi.f90"
#endif








  !####################################################################
  !                    computational ROUTINES
  !####################################################################
  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
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






  !--------------------------------------------------------------------!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !--------------------------------------------------------------------!
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








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso


end module DMFT_WEISS_FIELD
