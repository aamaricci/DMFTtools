module DMFT_CONVERGENCE
  USE SF_FONTS
  USE SF_IOTOOLS, only: reg
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private

  logical :: mpi_master
  integer :: mpi_rank,mpi_size,irank


  !CONVERGENCE
  public :: check_convergence,check_convergence_global,check_convergence_local


  interface check_convergence
     module procedure i0_check_convergence_relative
     module procedure i1_check_convergence_relative
     module procedure i2_check_convergence_relative
     module procedure i3_check_convergence_relative
     module procedure d0_check_convergence_relative_
     module procedure d0_check_convergence_relative
     module procedure d1_check_convergence_relative
     module procedure d2_check_convergence_relative
     module procedure d3_check_convergence_relative
     module procedure z0_check_convergence_relative
     module procedure z1_check_convergence_relative
     module procedure z2_check_convergence_relative
     module procedure z3_check_convergence_relative
  end interface check_convergence


  interface check_convergence_local
     module procedure i0_check_convergence_local
     module procedure i1_check_convergence_local
     module procedure i2_check_convergence_local
     module procedure d0_check_convergence_local
     module procedure d1_check_convergence_local
     module procedure d2_check_convergence_local
     module procedure z0_check_convergence_local
     module procedure z1_check_convergence_local
     module procedure z2_check_convergence_local
  end interface check_convergence_local

  interface check_convergence_global
     module procedure i0_check_convergence_global
     module procedure i1_check_convergence_global
     module procedure i2_check_convergence_global
     module procedure d0_check_convergence_global
     module procedure d1_check_convergence_global
     module procedure d2_check_convergence_global
     module procedure z0_check_convergence_global
     module procedure z1_check_convergence_global
     module procedure z2_check_convergence_global
  end interface check_convergence_global



contains





  !err = sum(NEW-OLD)/sum(NEW)
  include "error_convergence_relative.f90"

  !err = abs(NEW-OLD)
  include "error_convergence_local.f90"

  !err = sum((NEW-OLD)/NEW)
  include "error_convergence_global.f90"



  subroutine setup_mpi()
    mpi_master=.true.    
#ifdef _MPI    
    if(check_MPI())then
       mpi_master  = get_master_MPI()
    else
       mpi_master=.true.
    endif
#else
    mpi_master=.true.
#endif
  end subroutine setup_mpi





end module DMFT_CONVERGENCE
