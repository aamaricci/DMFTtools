module DMFT_GFIO
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,str,splot,sread,file_gunzip,file_gzip
  USE SF_MISC,      only: assert_shape
  USE SF_ARRAYS,    only: arange,linspace
  USE DMFT_CTRL_VARS

  implicit none
  private


  interface dmft_print_gf_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
  end interface dmft_print_gf_matsubara
  !
  interface dmft_read_gf_matsubara
     module procedure :: dmft_gf_read_matsubara_main
     module procedure :: dmft_gf_read_matsubara_ineq
     module procedure :: dmft_gij_read_matsubara
  end interface dmft_read_gf_matsubara


  interface dmft_print_gf_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
  end interface dmft_print_gf_realaxis
  !
  interface dmft_read_gf_realaxis
     module procedure :: dmft_gf_read_realaxis_main
     module procedure :: dmft_gf_read_realaxis_ineq
     module procedure :: dmft_gij_read_realaxis
  end interface dmft_read_gf_realaxis



  public :: set_gf_suffix
  !
  public :: dmft_print_gf_matsubara
  public :: dmft_print_gf_realaxis

  public :: dmft_read_gf_matsubara
  public :: dmft_read_gf_realaxis

  real(8),dimension(:),allocatable          :: wm !Matsubara frequencies
  real(8),dimension(:),allocatable          :: wr !Real frequencies

  character(len=128)                        :: suffix
  character(len=128)                        :: gf_suffix='.dat'
  integer                                   :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats
  integer                                   :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
  !
  real(8)                                   :: beta
  real(8)                                   :: wini,wfin 



contains



  subroutine set_gf_suffix(string)
    character(len=*) :: string
    gf_suffix=reg(string)
  end subroutine set_gf_suffix


  !> PRINT GF:
  include "print_gf.f90"


  !> READ GF:
  include "read_gf.f90"

end module DMFT_GFIO
