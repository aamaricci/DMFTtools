module DMFT_GFIO
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,str,splot,sread,file_gunzip,file_gzip,file_targz,file_untargz
  USE SF_MISC,      only: assert_shape
  USE SF_ARRAYS,    only: arange,linspace
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  !
  USE DMFT_CTRL_VARS

  implicit none
  private


  interface write_function_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif

  end interface write_function_matsubara

  interface dmft_write_function_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface dmft_write_function_matsubara

  interface write_gf_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface write_gf_matsubara

  interface dmft_write_gf_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface dmft_write_gf_matsubara

  interface print_function_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface print_function_matsubara

  interface dmft_print_function_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface dmft_print_function_matsubara

  interface print_gf_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface print_gf_matsubara

  interface dmft_print_gf_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
     module procedure :: dmft_gij_print_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_matsubara_cluster_ineq
#endif
  end interface dmft_print_gf_matsubara





  interface write_function_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface write_function_realaxis

  interface dmft_write_function_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface dmft_write_function_realaxis

  interface write_gf_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface write_gf_realaxis

  interface dmft_write_gf_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface dmft_write_gf_realaxis

  interface print_function_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface print_function_realaxis

  interface dmft_print_function_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface dmft_print_function_realaxis

  interface print_gf_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface print_gf_realaxis

  interface dmft_print_gf_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
     module procedure :: dmft_gij_print_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_realaxis_cluster_ineq
#endif
  end interface dmft_print_gf_realaxis





  interface read_function_matsubara
     module procedure :: dmft_gf_read_matsubara_main
     module procedure :: dmft_gf_read_matsubara_ineq
     module procedure :: dmft_gij_read_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_matsubara_cluster_ineq
#endif
  end interface read_function_matsubara

  interface dmft_read_function_matsubara
     module procedure :: dmft_gf_read_matsubara_main
     module procedure :: dmft_gf_read_matsubara_ineq
     module procedure :: dmft_gij_read_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_matsubara_cluster_ineq
#endif
  end interface dmft_read_function_matsubara

  interface read_gf_matsubara
     module procedure :: dmft_gf_read_matsubara_main
     module procedure :: dmft_gf_read_matsubara_ineq
     module procedure :: dmft_gij_read_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_matsubara_cluster_ineq
#endif
  end interface read_gf_matsubara

  interface dmft_read_gf_matsubara
     module procedure :: dmft_gf_read_matsubara_main
     module procedure :: dmft_gf_read_matsubara_ineq
     module procedure :: dmft_gij_read_matsubara
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_matsubara_cluster_ineq
#endif
  end interface dmft_read_gf_matsubara




  interface read_function_realaxis
     module procedure :: dmft_gf_read_realaxis_main
     module procedure :: dmft_gf_read_realaxis_ineq
     module procedure :: dmft_gij_read_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_realaxis_cluster_ineq
#endif     
  end interface read_function_realaxis

  interface dmft_read_function_realaxis
     module procedure :: dmft_gf_read_realaxis_main
     module procedure :: dmft_gf_read_realaxis_ineq
     module procedure :: dmft_gij_read_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_realaxis_cluster_ineq
#endif     
  end interface dmft_read_function_realaxis

  interface read_gf_realaxis
     module procedure :: dmft_gf_read_realaxis_main
     module procedure :: dmft_gf_read_realaxis_ineq
     module procedure :: dmft_gij_read_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_realaxis_cluster_ineq
#endif     
  end interface read_gf_realaxis

  interface dmft_read_gf_realaxis
     module procedure :: dmft_gf_read_realaxis_main
     module procedure :: dmft_gf_read_realaxis_ineq
     module procedure :: dmft_gij_read_realaxis
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_read_realaxis_cluster_ineq
#endif     
  end interface dmft_read_gf_realaxis


  public :: set_gf_suffix

  public :: write_function_matsubara
  public :: write_gf_matsubara
  public :: dmft_write_function_matsubara
  public :: dmft_write_gf_matsubara

  public :: print_function_matsubara
  public :: print_gf_matsubara
  public :: dmft_print_function_matsubara
  public :: dmft_print_gf_matsubara


  public :: write_function_realaxis
  public :: write_gf_realaxis
  public :: dmft_write_function_realaxis
  public :: dmft_write_gf_realaxis

  public :: print_function_realaxis
  public :: print_gf_realaxis
  public :: dmft_print_function_realaxis
  public :: dmft_print_gf_realaxis

  public :: read_function_matsubara
  public :: read_gf_matsubara
  public :: dmft_read_function_matsubara
  public :: dmft_read_gf_matsubara

  public :: read_function_realaxis
  public :: read_gf_realaxis
  public :: dmft_read_function_realaxis
  public :: dmft_read_gf_realaxis



  real(8),dimension(:),allocatable :: wm !Matsubara frequencies
  real(8),dimension(:),allocatable :: wr !Real frequencies

  character(len=128)               :: suffix
  character(len=128)               :: gf_suffix='.dat'
  integer                          :: Lk,Nineq,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats
  integer                          :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js,iineq
  !
  real(8)                          :: beta
  real(8)                          :: wini,wfin 


  integer                          :: mpi_ierr
  integer                          :: mpi_rank
  integer                          :: mpi_size
  logical                          :: mpi_master


contains



  subroutine set_gf_suffix(string)
    character(len=*) :: string
    gf_suffix=reg(string)
  end subroutine set_gf_suffix


  !> PRINT GF:
  !  
  !> MATSUBARA: single site
  subroutine dmft_gf_print_matsubara_main(Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    Nspin = size(Gmats,1)
    Norb  = size(Gmats,3)
    Lmats = size(Gmats,5)
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_print_matsubara",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"//str(gf_suffix)
                call splot(reg(suffix),wm,Gmats(ispin,ispin,iorb,iorb,:))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_iw"//reg(gf_suffix)
                   call splot(reg(suffix),wm,Gmats(ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
          !
       case default                  !print all off-diagonals
          write(*,"(A,1x,A)") reg(fname),"matsubara: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"//reg(gf_suffix)
                      call splot(reg(suffix),wm,Gmats(ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_print_matsubara_main


  !> MATSUBARA: ineq sites
  subroutine dmft_gf_print_matsubara_ineq(Gmats,fname,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Gmats
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    logical,optional                             :: itar
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,2)
    Norb  = size(Gmats,4)
    Lmats = size(Gmats,6)
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_print_matsubara_ineq",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin-orbital diagonal elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"//reg(gf_suffix)
                call splot(reg(suffix),wm,Gmats(:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_iw"//reg(gf_suffix)
                   call splot(reg(suffix),wm,Gmats(:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
          !
       case(3)                  !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. No Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"//reg(gf_suffix)
                      call splot(reg(suffix),wm,Gmats(:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case(4)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin-orbital diagonal elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nlat
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call splot(reg(suffix),wm,Gmats(ilat,ispin,ispin,iorb,iorb,:))
                enddo
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"
                if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                     pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
             enddo
          enddo
          !
       case(5)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   do ilat=1,Nlat
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                      call splot(reg(suffix),wm,Gmats(ilat,ispin,ispin,iorb,jorb,:))
                   enddo
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_iw"
                   if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                        pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      do ilat=1,Nlat
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                         call splot(reg(suffix),wm,Gmats(ilat,ispin,jspin,iorb,jorb,:))
                      enddo
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"
                      if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                           pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_print_matsubara_ineq



  !> MATSUBARA: cluster + ineq sites
#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_gf_print_matsubara_cluster_ineq(Gmats,fname,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in) :: Gmats
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    logical,optional                             :: itar
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    !
    Nineq = size(Gmats,1)
    Nlat  = size(Gmats,2)
    Nspin = size(Gmats,4)
    Norb  = size(Gmats,6)
    Lmats = size(Gmats,8)
    call assert_shape(Gmats,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_print_matsubara_ineq",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         do iineq=1,Nineq
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//&
                                 "_iw_"//reg(index)//str(iineq,pad)//reg(gf_suffix)
                            call splot(reg(suffix),wm,Gmats(iineq,ilat,jlat,ispin,ispin,iorb,jorb,:))
                         enddo
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_iw"
                         if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            do iineq=1,Nineq
                               suffix=reg(fname)//&
                                    "_i"//str(ilat)//str(jlat)//&
                                    "_l"//str(iorb)//str(jorb)//&
                                    "_s"//str(ispin)//str(jspin)//&
                                    "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                               call splot(reg(suffix),wm,Gmats(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:))
                            enddo
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_iw"
                            if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                                 pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_print_matsubara_cluster_ineq
#endif

  !> MATSUBARA: full Gij
  subroutine dmft_gij_print_matsubara(Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Gmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats/Lreal]
    real(8),dimension(size(Gmats,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,3)
    Norb  = size(Gmats,5)
    Lmats = size(Gmats,7)
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gij_print_matsubara",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(0)
          write(*,"(A,1x,A)")reg(fname),"matsubara: not written on file."
       case(1)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"matsubara: spin-orbital diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"//reg(gf_suffix)
                call splot(reg(suffix),wm,Gmats(:,:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       case(2)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin-orbital diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      suffix=reg(fname)//&
                           "_i"//str(ilat)//str(jlat)//&
                           "_l"//str(iorb)//str(iorb)//&
                           "_s"//str(ispin)//&
                           "_iw"//str(gf_suffix)
                      call splot(reg(suffix),wm,Gmats(ilat,jlat,ispin,ispin,iorb,iorb,:))
                   enddo
                enddo
             enddo
          enddo
       case(3)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"matsubara: spin diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_iw"//reg(gf_suffix)
                   call splot(reg(suffix),wm,Gmats(:,:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       case(4)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_iw"//str(gf_suffix)
                         call splot(reg(suffix),wm,Gmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       case(5)              !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"//reg(gf_suffix)
                      call splot(reg(suffix),wm,Gmats(:,:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
       case default                  !print all off-diagonals, many files
          write(*,"(A,1x,A)") reg(fname),"matsubara: write all elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_iw"//reg(gf_suffix)
                            call splot(reg(suffix),wm,Gmats(ilat,jlat,ispin,jspin,iorb,jorb,:))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine dmft_gij_print_matsubara





  !##################################################################
  !##################################################################
  !##################################################################







  !> REALAXIS: single site
  subroutine dmft_gf_print_realaxis_main(Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Greal
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !
    Nspin = size(Greal,1)
    Norb  = size(Greal,3)
    Lreal = size(Greal,5)
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_print_realaxis",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)") reg(fname),"real: write spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"//reg(gf_suffix)
                call splot(reg(suffix),wr,Greal(ispin,ispin,iorb,iorb,:))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)") reg(fname),"real: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_realw"//reg(gf_suffix)
                   call splot(reg(suffix),wr,Greal(ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
          !
       case default                  !print all off-diagonals
          write(*,"(A,1x,A)") reg(fname),"real: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"//reg(gf_suffix)
                      call splot(reg(suffix),wr,Greal(ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_print_realaxis_main


  !> REALAXIS: ineq sites
  subroutine dmft_gf_print_realaxis_ineq(Greal,fname,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Greal
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    logical,optional                             :: itar
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,2)
    Norb  = size(Greal,4)
    Lreal = size(Greal,6)
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_print_realaxis_ineq",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"real: write spin-orbital diagonal elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"//reg(gf_suffix)
                call splot(reg(suffix),wr,Greal(:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_realw"//reg(gf_suffix)
                   call splot(reg(suffix),wr,Greal(:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
          !
       case(3)
          write(*,"(A,1x,A)")reg(fname),"real: write all elements. No Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"//reg(gf_suffix)
                      call splot(reg(suffix),wr,Greal(:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case(4)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"real: write spin-orbital diagonal elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nlat
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call splot(reg(suffix),wr,Greal(ilat,ispin,ispin,iorb,iorb,:))
                enddo
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"
                if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                     pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
             enddo
          enddo
          !
       case(5)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   do ilat=1,Nlat
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                      call splot(reg(suffix),wr,Greal(ilat,ispin,ispin,iorb,jorb,:))
                   enddo
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_realw"
                   if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                        pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"real: write all elements. Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      do ilat=1,Nlat
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                         call splot(reg(suffix),wr,Greal(ilat,ispin,jspin,iorb,jorb,:))
                      enddo
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"
                      if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                           pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_print_realaxis_ineq



  !> REAL AXIS: cluster + ineq sites
#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_gf_print_realaxis_cluster_ineq(Greal,fname,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in) :: Greal
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    logical,optional                             :: itar
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    !
    Nineq = size(Greal,1)
    Nlat  = size(Greal,2)
    Nspin = size(Greal,4)
    Norb  = size(Greal,6)
    Lreal = size(Greal,8)
    call assert_shape(Greal,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_print_realaxis_ineq",reg(fname)//"_mats")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         do iineq=1,Nineq
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//&
                                 "_realw_"//reg(index)//str(iineq,pad)//reg(gf_suffix)
                            call splot(reg(suffix),wr,Greal(iineq,ilat,jlat,ispin,ispin,iorb,jorb,:))
                         enddo
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_realw"
                         if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            do iineq=1,Nineq
                               suffix=reg(fname)//&
                                    "_i"//str(ilat)//str(jlat)//&
                                    "_l"//str(iorb)//str(jorb)//&
                                    "_s"//str(ispin)//str(jspin)//&
                                    "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                               call splot(reg(suffix),wr,Greal(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:))
                            enddo
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_realw"
                            if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                                 pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_print_realaxis_cluster_ineq
#endif

  !> REALAXIS: full GF
  subroutine dmft_gij_print_realaxis(Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Greal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal/Lreal]
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    integer                                        :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb,ilat,jlat
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,3)
    Norb  = size(Greal,5)
    Lreal = size(Greal,7)
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gij_print_realaxis",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(0)
          write(*,"(A,1x,A)")reg(fname),"real: not written on file."
       case(1)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"real: spin-orbital diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"//reg(gf_suffix)
                call splot(reg(suffix),wr,Greal(:,:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       case(2)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"real: write spin-orbital diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      suffix=reg(fname)//&
                           "_i"//str(ilat)//str(jlat)//&
                           "_l"//str(iorb)//str(iorb)//&
                           "_s"//str(ispin)//&
                           "_realw"//str(gf_suffix)
                      call splot(reg(suffix),wr,Greal(ilat,jlat,ispin,ispin,iorb,iorb,:))
                   enddo
                enddo
             enddo
          enddo
       case(3)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"real: spin diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_realw"//reg(gf_suffix)
                   call splot(reg(suffix),wr,Greal(:,:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       case(4)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"real: write spin diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_realw"//str(gf_suffix)
                         call splot(reg(suffix),wr,Greal(ilat,jlat,ispin,ispin,iorb,jorb,:))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       case(5)              !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"real: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"//reg(gf_suffix)
                      call splot(reg(suffix),wr,Greal(:,:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
       case default                  !print all off-diagonals, many files
          write(*,"(A,1x,A)") reg(fname),"rea: write all elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_realw"//reg(gf_suffix)
                            call splot(reg(suffix),wr,Greal(ilat,jlat,ispin,jspin,iorb,jorb,:))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
  end subroutine dmft_gij_print_realaxis





  !> READ GF:
  !> MATSUBARA: single site
  subroutine dmft_gf_read_matsubara_main(Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    Nspin = size(Gmats,1)
    Norb  = size(Gmats,3)
    Lmats = size(Gmats,5)
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_read_matsubara",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"//str(gf_suffix)
                call sread(reg(suffix),wm,Gmats(ispin,ispin,iorb,iorb,:))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_iw"//reg(gf_suffix)
                   call sread(reg(suffix),wm,Gmats(ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
          !
       case default                  !print all off-diagonals
          write(*,"(A,1x,A)") reg(fname),"matsubara: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"//reg(gf_suffix)
                      call sread(reg(suffix),wm,Gmats(ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
#ifdef _MPI    
    if(check_MPI())call Bcast_MPI(MPI_COMM_WORLD,Gmats)
#endif
  end subroutine dmft_gf_read_matsubara_main


  !> MATSUBARA: ineq sites
  subroutine dmft_gf_read_matsubara_ineq(Gmats,fname,iprint,ineq_index,ineq_pad)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Gmats
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,2)
    Norb  = size(Gmats,4)
    Lmats = size(Gmats,6)
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_read_matsubara_ineq",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin-orbital diagonal elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"//reg(gf_suffix)
                call file_gunzip(reg(suffix))
                call sread(reg(suffix),wm,Gmats(:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_iw"//reg(gf_suffix)
                   call file_gunzip(reg(suffix))
                   call sread(reg(suffix),wm,Gmats(:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
          !
       case(3)                  !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. No Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"//reg(gf_suffix)
                      call file_gunzip(reg(suffix))
                      call sread(reg(suffix),wm,Gmats(:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case(4)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin-orbital diagonal elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nlat
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_iw"
                   call file_untargz(tarball=reg(suffix))
                   !
                   !
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call sread(reg(suffix),wm,Gmats(ilat,ispin,ispin,iorb,iorb,:))
                   !
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_iw"
                   call file_targz(tarball=reg(suffix),&
                        pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                enddo
             enddo
          enddo
          !
       case(5)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   do ilat=1,Nlat
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_iw"
                      call file_untargz(tarball=reg(suffix))
                      !
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                      call sread(reg(suffix),wm,Gmats(ilat,ispin,ispin,iorb,jorb,:))
                      !
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_iw"
                      call file_targz(tarball=reg(suffix),&
                           pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      do ilat=1,Nlat
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_iw"
                         call file_untargz(tarball=reg(suffix))
                         !
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                         call sread(reg(suffix),wm,Gmats(ilat,ispin,jspin,iorb,jorb,:))
                         !
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_iw"
                         call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
#ifdef _MPI    
    if(check_MPI())call Bcast_MPI(MPI_COMM_WORLD,Gmats)
#endif
  end subroutine dmft_gf_read_matsubara_ineq



#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_gf_read_matsubara_cluster_ineq(Gmats,fname,iprint,ineq_index,ineq_pad)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in) :: Gmats
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !    
    !
    Nineq = size(Gmats,1)
    Nlat  = size(Gmats,2)
    Nspin = size(Gmats,4)
    Norb  = size(Gmats,6)
    Lmats = size(Gmats,8)
    call assert_shape(Gmats,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_print_matsubara_ineq",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         do iineq=1,Nineq
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//&
                                 "_iw_"//reg(index)//str(iineq,pad)//reg(gf_suffix)
                            call file_gunzip(reg(suffix))
                            call sread(reg(suffix),wm,Gmats(iineq,ilat,jlat,ispin,ispin,iorb,jorb,:))
                         enddo
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_iw"
                         call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            do iineq=1,Nineq
                               suffix=reg(fname)//&
                                    "_i"//str(ilat)//str(jlat)//&
                                    "_l"//str(iorb)//str(jorb)//&
                                    "_s"//str(ispin)//str(jspin)//&
                                    "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                               call file_gunzip(reg(suffix))
                               call sread(reg(suffix),wm,Gmats(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:))
                            enddo
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_iw"
                            call file_targz(tarball=reg(suffix),&
                                 pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_read_matsubara_cluster_ineq
#endif



  !> MATSUBARA: full Gij
  subroutine dmft_gij_read_matsubara(Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Gmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats/Lreal]
    real(8),dimension(size(Gmats,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,3)
    Norb  = size(Gmats,5)
    Lmats = size(Gmats,7)
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gij_read_matsubara",reg(fname)//"_mats")
    !
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    !
    if(mpi_master)then
       select case(iprint)
       case(0)
          write(*,"(A,1x,A)")reg(fname),"matsubara: not written on file."
       case(1)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"matsubara: spin-orbital diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_iw"//reg(gf_suffix)
                call sread(reg(suffix),wm,Gmats(:,:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       case(2)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin-orbital diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      suffix=reg(fname)//&
                           "_i"//str(ilat)//str(jlat)//&
                           "_l"//str(iorb)//str(iorb)//&
                           "_s"//str(ispin)//&
                           "_iw"//str(gf_suffix)
                      call sread(reg(suffix),wm,Gmats(ilat,jlat,ispin,ispin,iorb,iorb,:))
                   enddo
                enddo
             enddo
          enddo
       case(3)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"matsubara: spin diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_iw"//reg(gf_suffix)
                   call sread(reg(suffix),wm,Gmats(:,:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       case(4)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"matsubara: write spin diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_iw"//str(gf_suffix)
                         call sread(reg(suffix),wm,Gmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       case(5)              !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_iw"//reg(gf_suffix)
                      call sread(reg(suffix),wm,Gmats(:,:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
       case default                  !print all off-diagonals, many files
          write(*,"(A,1x,A)") reg(fname),"matsubara: write all elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_iw"//reg(gf_suffix)
                            call sread(reg(suffix),wm,Gmats(ilat,jlat,ispin,jspin,iorb,jorb,:))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
#ifdef _MPI    
    if(check_MPI())call Bcast_MPI(MPI_COMM_WORLD,Gmats)
#endif
  end subroutine dmft_gij_read_matsubara





  !##################################################################
  !##################################################################
  !##################################################################







  !> REALAXIS: single site
  subroutine dmft_gf_read_realaxis_main(Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Greal
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !
    Nspin = size(Greal,1)
    Norb  = size(Greal,3)
    Lreal = size(Greal,5)
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_read_realaxis",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)") reg(fname),"real: write spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"//reg(gf_suffix)
                call sread(reg(suffix),wr,Greal(ispin,ispin,iorb,iorb,:))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)") reg(fname),"real: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_realw"//reg(gf_suffix)
                   call sread(reg(suffix),wr,Greal(ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
          !
       case default                  !print all off-diagonals
          write(*,"(A,1x,A)") reg(fname),"real: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"//reg(gf_suffix)
                      call sread(reg(suffix),wr,Greal(ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
#ifdef _MPI    
    if(check_MPI())call Bcast_MPI(MPI_COMM_WORLD,Greal)
#endif
  end subroutine dmft_gf_read_realaxis_main


  !> REALAXIS: ineq sites
  subroutine dmft_gf_read_realaxis_ineq(Greal,fname,iprint,ineq_index,ineq_pad)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Greal
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,2)
    Norb  = size(Greal,4)
    Lreal = size(Greal,6)
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_read_realaxis_ineq",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"real: write spin-orbital diagonal elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"//reg(gf_suffix)
                call file_gunzip(reg(suffix))
                call sread(reg(suffix),wr,Greal(:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. No Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        "_realw"//reg(gf_suffix)
                   call file_gunzip(reg(suffix))
                   call sread(reg(suffix),wr,Greal(:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
          !
       case(3)
          write(*,"(A,1x,A)")reg(fname),"real: write all elements. No Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"//reg(gf_suffix)
                      call file_gunzip(reg(suffix))
                      call sread(reg(suffix),wr,Greal(:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case(4)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"real: write spin-orbital diagonal elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nlat
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_realw"
                   call file_untargz(tarball=reg(suffix))
                   !
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call sread(reg(suffix),wr,Greal(ilat,ispin,ispin,iorb,iorb,:))
                   !
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_realw"
                   call file_targz(tarball=reg(suffix),&
                        pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                   !
                enddo
             enddo
          enddo
          !
       case(5)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. Split."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   do ilat=1,Nlat
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_realw"
                      call file_untargz(tarball=reg(suffix))
                      !
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                      call sread(reg(suffix),wr,Greal(ilat,ispin,ispin,iorb,jorb,:))
                      !
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           "_realw"
                      call file_targz(tarball=reg(suffix),&
                           pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"real: write all elements. Split."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      do ilat=1,Nlat
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_realw"
                         call file_untargz(tarball=reg(suffix))
                         !
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                         call sread(reg(suffix),wr,Greal(ilat,ispin,jspin,iorb,jorb,:))
                         !
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              "_realw"
                         call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
#ifdef _MPI    
    if(check_MPI())call Bcast_MPI(MPI_COMM_WORLD,Greal)
#endif
  end subroutine dmft_gf_read_realaxis_ineq



#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_gf_read_realaxis_cluster_ineq(Greal,fname,iprint,ineq_index,ineq_pad)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in) :: Greal
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    character(len=*),optional                    :: ineq_index
    integer,optional                             :: ineq_pad
    character(len=:),allocatable                 :: index
    integer                                      :: pad
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(present(ineq_index))then
       index=trim(adjustl(trim(ineq_index)))
    else
       index='indx'
    endif
    pad=6
    if(present(ineq_pad))pad=ineq_pad
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !    
    !
    Nineq = size(Greal,1)
    Nlat  = size(Greal,2)
    Nspin = size(Greal,4)
    Norb  = size(Greal,6)
    Lreal = size(Greal,8)
    call assert_shape(Greal,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_print_realaxis_ineq",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"realaxis: write spin diagonal and all orbitals elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         do iineq=1,Nineq
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//&
                                 "_realw_"//reg(index)//str(iineq,pad)//reg(gf_suffix)
                            call file_gunzip(reg(suffix))
                            call sread(reg(suffix),wm,Greal(iineq,ilat,jlat,ispin,ispin,iorb,jorb,:))
                         enddo
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_realw"
                         call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"realaxis: write all elements. Split."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            do iineq=1,Nineq
                               suffix=reg(fname)//&
                                    "_i"//str(ilat)//str(jlat)//&
                                    "_l"//str(iorb)//str(jorb)//&
                                    "_s"//str(ispin)//str(jspin)//&
                                    "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                               call file_gunzip(reg(suffix))
                               call sread(reg(suffix),wm,Greal(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:))
                            enddo
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_realw"
                            call file_targz(tarball=reg(suffix),&
                                 pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
  end subroutine dmft_gf_read_realaxis_cluster_ineq
#endif




  !> REALAXIS: full GF
  subroutine dmft_gij_read_realaxis(Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Greal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal/Lreal]
    real(8),dimension(size(Greal,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    integer                                        :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb,ilat,jlat
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Retrieve parameters:
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,3)
    Norb  = size(Greal,5)
    Lreal = size(Greal,7)
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gij_read_realaxis",reg(fname)//"_real")
    !
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    if(mpi_master)then
       select case(iprint)
       case(0)
          write(*,"(A,1x,A)")reg(fname),"real: not written on file."
       case(1)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"real: spin-orbital diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     "_realw"//reg(gf_suffix)
                call file_gunzip(reg(suffix))
                call sread(reg(suffix),wr,Greal(:,:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       case(2)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"real: write spin-orbital diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      suffix=reg(fname)//&
                           "_i"//str(ilat)//str(jlat)//&
                           "_l"//str(iorb)//str(iorb)//&
                           "_s"//str(ispin)//&
                           "_realw"//str(gf_suffix)
                      call sread(reg(suffix),wr,Greal(ilat,jlat,ispin,ispin,iorb,iorb,:))
                   enddo
                enddo
             enddo
          enddo
       case(3)                  !print only diagonal elements, single file  
          write(*,"(A,1x,A)")reg(fname),"real: spin diagonal elements on one file."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        "_realw"//reg(gf_suffix)
                   call file_gunzip(reg(suffix))
                   call sread(reg(suffix),wr,Greal(:,:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       case(4)                  !print only diagonal elements, many files
          write(*,"(A,1x,A)") reg(fname),"real: write spin diagonal elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              "_realw"//str(gf_suffix)
                         call sread(reg(suffix),wr,Greal(ilat,jlat,ispin,ispin,iorb,jorb,:))
                      enddo
                   enddo
                enddo
             enddo
          enddo
       case(5)              !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"real: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           "_realw"//reg(gf_suffix)
                      call file_gunzip(reg(suffix))
                      call sread(reg(suffix),wr,Greal(:,:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
       case default                  !print all off-diagonals, many files
          write(*,"(A,1x,A)") reg(fname),"rea: write all elements on many files."
          do ilat=1,Nlat
             do jlat=1,Nlat
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 "_realw"//reg(gf_suffix)
                            call sread(reg(suffix),wr,Greal(ilat,jlat,ispin,jspin,iorb,jorb,:))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
#ifdef _MPI    
    if(check_MPI())call Bcast_MPI(MPI_COMM_WORLD,Greal)
#endif
  end subroutine dmft_gij_read_realaxis




end module DMFT_GFIO
