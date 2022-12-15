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


  !AGNOSTIC INTERFACE:
  interface dmft_write_gf
     module procedure :: dmft_gf_push_zeta
     module procedure :: dmft_gf_print_rank2
     module procedure :: dmft_gf_print_rank4
     module procedure :: dmft_gf_print_rank5
     module procedure :: dmft_gf_print_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_cluster_ineq
#endif
  end interface dmft_write_gf


  integer                          :: Lfreq
  real(8),dimension(:),allocatable :: wfreq

  character(len=128)               :: suffix
  character(len=128)               :: gf_suffix='.dat'
  character(len=8)                 :: w_suffix
  integer                          :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Nineq,Nilso,Ntot
  integer                          :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js,iineq
  !
  real(8)                          :: beta
  real(8)                          :: wini,wfin 
  !
  integer                          :: mpi_ierr
  integer                          :: mpi_rank
  integer                          :: mpi_size
  logical                          :: mpi_master



  public :: set_gf_suffix
  public :: dmft_write_gf


contains



  subroutine dmft_gf_print_rank2(Func,fname,axis,iprint)
    complex(8),dimension(:,:,:),intent(in) :: Func
    character(len=*),intent(in)            :: fname
    character(len=*)                       :: axis
    integer,intent(in),optional            :: iprint
    integer                                :: iprint_
    !
    iprint_=10;if(present(iprint))iprint_=iprint
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !    
    Ntot  = size(Func,1)
    Lfreq = size(Func,3)
    call assert_shape(Func,[Ntot,Ntot,Lfreq],"gf_print_main","Func")    
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)then
       select case(iprint_)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write diagonal elements."
          do io=1,Ntot
             suffix=reg(fname)//"_i"//str(io)//"j"//str(io)//&
                  str(w_suffix)//str(gf_suffix)
             call splot(reg(suffix),wfreq,Func(io,io,:))
          enddo
       case default                  !print all off-diagonals
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write all elements."
          do io=1,Ntot
             do jo=1,Ntot
                suffix=reg(fname)//"_i"//str(io)//"j"//str(jo)//&
                     str(w_suffix)//reg(gf_suffix)
                call splot(reg(suffix),wfreq,Func(io,jo,:))
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wfreq))deallocate(wfreq)
  end subroutine dmft_gf_print_rank2



  subroutine dmft_gf_print_rank4(Func,fname,axis,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Func
    character(len=*),intent(in)                :: fname
    character(len=*)                           :: axis
    integer,intent(in),optional                :: iprint
    integer                                    :: iprint_
    !
    iprint_=10;if(present(iprint))iprint_=iprint
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !    
    Nspin = size(Func,1)
    Norb  = size(Func,3)
    Lfreq = size(Func,5)
    call assert_shape(Func,[Nspin,Nspin,Norb,Norb,Lfreq],"dmft_gf_print_rank4","Func")    
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)then
       select case(iprint_)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)//str(gf_suffix)
                call splot(reg(suffix),wfreq,Func(ispin,ispin,iorb,iorb,:))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(gf_suffix)
                   call splot(reg(suffix),wfreq,Func(ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
          !
       case default                  !print all off-diagonals
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wfreq,Func(ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wfreq))deallocate(wfreq)
  end subroutine dmft_gf_print_rank4





  subroutine dmft_gf_print_rank5(Func,fname,axis,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Func
    character(len=*),intent(in)                  :: fname
    character(len=*)                             :: axis
    integer,intent(in),optional                  :: iprint
    integer                                      :: iprint_
    character(len=*),optional                    :: ineq_index
    character(len=:),allocatable                 :: index
    integer,optional                             :: ineq_pad
    integer                                      :: pad
    logical,optional                             :: itar
    logical                                      :: itar_
    !
    iprint_=10;if(present(iprint))iprint_=iprint
    index='_indx';if(present(ineq_index))index="_"//reg(ineq_index)
    pad=6        ;if(present(ineq_pad))pad=ineq_pad
    itar_=.false.;if(present(itar))itar_=itar
    !
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !
    Nlat  = size(Func,1)
    Nspin = size(Func,2)
    Norb  = size(Func,4)
    Lfreq = size(Func,6)
    call assert_shape(Func,[Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_gf_print_rank5","Func")
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)then
       select case(iprint_)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin-orbital diagonal elements. Single File."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)//reg(gf_suffix)
                call splot(reg(suffix),wfreq,Func(:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin diagonal and all orbitals elements. Single File."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(gf_suffix)
                   call splot(reg(suffix),wfreq,Func(:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
          !
       case(3)                  !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: all elements. Single File."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wfreq,Func(:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case(4)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do ilat=1,Nlat
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(iorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call splot(reg(suffix),wfreq,Func(ilat,ispin,ispin,iorb,iorb,:))
                enddo
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)
                if(itar_)call file_targz(tarball=reg(suffix),&
                     pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
             enddo
          enddo
          !
       case(5)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   do ilat=1,Nlat
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//&
                           str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                      call splot(reg(suffix),wfreq,Func(ilat,ispin,ispin,iorb,jorb,:))
                   enddo
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)
                   if(itar_)call file_targz(tarball=reg(suffix),&
                        pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      do ilat=1,Nlat
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                         call splot(reg(suffix),wfreq,Func(ilat,ispin,jspin,iorb,jorb,:))
                      enddo
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)
                      if(itar_)call file_targz(tarball=reg(suffix),&
                           pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wfreq))deallocate(wfreq)
  end subroutine dmft_gf_print_rank5




  subroutine dmft_gf_print_rank6(Func,fname,axis,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Func ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lfreq/Lreal]
    character(len=*),intent(in)                    :: fname
    character(len=*)                               :: axis
    integer,intent(in),optional                    :: iprint
    integer                                        :: iprint_
    character(len=*),optional                      :: ineq_index
    character(len=:),allocatable                   :: index
    integer,optional                               :: ineq_pad
    integer                                        :: pad
    logical,optional                               :: itar
    logical                                        :: itar_
    !
    iprint_=10;if(present(iprint))iprint_=iprint
    index='_indx';if(present(ineq_index))index="_"//reg(ineq_index)
    pad=6        ;if(present(ineq_pad))pad=ineq_pad
    itar_=.false.;if(present(itar))itar_=itar
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    Nlat  = size(Func,1)
    Nspin = size(Func,3)
    Norb  = size(Func,5)
    Lfreq = size(Func,7)
    call assert_shape(Func,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_gij_print_matsubara","Func")
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)then
       select case(iprint_)
       case(1)!print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin-orbital diagonal elements. Single File."
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)//reg(gf_suffix)
                call splot(reg(suffix),wfreq,Func(:,:,ispin,ispin,iorb,iorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
          !
       case(2)!print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin diagonal and all orbitals elements. Single File."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(gf_suffix)
                   call splot(reg(suffix),wfreq,Func(:,:,ispin,ispin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
          !
       case(3)              !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: all elements. Single File."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wfreq,Func(:,:,ispin,jspin,iorb,jorb,:))
                      call file_gzip(reg(suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       case(4)
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin-orbital diagonal elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                !
                do ilat=1,Nlat
                   do jlat=1,Nlat
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(iorb)//&
                           "_s"//str(ispin)//&
                           str(w_suffix)//reg(index)//&
                           str(ilat,pad)//"_"//str(jlat,pad)//reg(gf_suffix)
                      call splot(reg(suffix),wfreq,Func(ilat,jlat,ispin,ispin,iorb,iorb,:))
                   enddo
                enddo
                suffix=reg(fname)//&
                     "_l"//str(iorb)//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)
                if(itar_)call file_targz(tarball=reg(suffix),&
                     pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
             enddo
          enddo
          !
       case(5)
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   do ilat=1,Nlat
                      do jlat=1,Nlat
                         suffix=reg(fname)//&
                              "_l"//str(iorb)//str(iorb)//&
                              "_s"//str(ispin)//&
                              str(w_suffix)//reg(index)//&
                              str(ilat,pad)//"_"//str(jlat,pad)//reg(gf_suffix)
                         call splot(reg(suffix),wfreq,Func(ilat,jlat,ispin,ispin,iorb,jorb,:))
                      enddo
                   enddo
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)
                   if(itar_)call file_targz(tarball=reg(suffix),&
                        pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                   !
                enddo
             enddo
          enddo
          !
       case default                  !print all off-diagonals, many files
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write all elements."
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      !
                      do ilat=1,Nlat
                         do jlat=1,Nlat
                            suffix=reg(fname)//&
                                 "_l"//str(iorb)//str(iorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 str(w_suffix)//reg(index)//&
                                 str(ilat,pad)//"_"//str(jlat,pad)//reg(gf_suffix)
                            call splot(reg(suffix),wfreq,Func(ilat,jlat,ispin,jspin,iorb,jorb,:))
                         enddo
                      enddo
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//str(iorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)
                      if(itar_)call file_targz(tarball=reg(suffix),&
                           pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
    if(allocated(wfreq))deallocate(wfreq)
  end subroutine dmft_gf_print_rank6


#if __GFORTRAN__ &&  __GNUC__ > 8
  subroutine dmft_gf_print_cluster_ineq(Func,fname,axis,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:,:,:,:,:),intent(in) :: Func
    character(len=*),intent(in)                      :: fname
    integer,intent(in)                               :: iprint
    character(len=*)                                :: axis
    character(len=*),optional                        :: ineq_index
    integer,optional                                 :: ineq_pad
    logical,optional                                 :: itar
    character(len=:),allocatable                     :: index
    integer                                          :: pad
    logical                                          :: itar_
    !
    !MPI setup:
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    index='_indx';if(present(ineq_index))index="_"//reg(ineq_index)
    pad=6        ;if(present(ineq_pad))pad=ineq_pad
    itar_=.false.;if(present(itar))itar_=itar
    !
    !
    Nineq = size(Func,1)
    Nlat  = size(Func,2)
    Nspin = size(Func,4)
    Norb  = size(Func,6)
    Lfreq = size(Func,8)
    call assert_shape(Func,[Nineq,Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lfreq],"dmft_gf_print_cluster_ineq","Func")
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)then
       select case(iprint)
       case(1)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write spin diagonal and all orbitals elements. Split."
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
                                 str(w_suffix)//reg(index)//str(iineq,pad)//reg(gf_suffix)
                            call splot(reg(suffix),wfreq,Func(iineq,ilat,jlat,ispin,ispin,iorb,jorb,:))
                         enddo
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//str(jlat)//&
                              "_l"//str(iorb)//str(jorb)//&
                              "_s"//str(ispin)//&
                              str(w_suffix)
                         if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                              pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write all elements. Split."
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
                                    str(w_suffix)//reg(index)//str(iineq,pad)//reg(gf_suffix)
                               call splot(reg(suffix),wfreq,Func(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:))
                            enddo
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//str(jlat)//&
                                 "_l"//str(iorb)//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 str(w_suffix)
                            if(present(itar).AND.itar)call file_targz(tarball=reg(suffix),&
                                 pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wfreq))deallocate(wfreq)
  end subroutine dmft_gf_print_cluster_ineq
#endif













  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################





  subroutine set_gf_suffix(string)
    character(len=*) :: string
    gf_suffix=reg(string)
  end subroutine set_gf_suffix

  subroutine dmft_gf_push_zeta(zeta)
    real(8),dimension(:) :: zeta
    Lfreq=size(zeta)
    if(allocated(wfreq))deallocate(wfreq)
    allocate(wfreq(Lfreq))
    wfreq=zeta
  end subroutine dmft_gf_push_zeta


  subroutine build_frequency_array(axis)
    character(len=*)              :: axis
    if(allocated(wfreq))then
       if(size(wfreq)/=Lfreq)stop "dmft_gfio_build_frequency_array ERROR: pushed wfreq has wrong size"
    else
       allocate(wfreq(Lfreq))
       select case(axis)
       case default;
          stop "dmft_gfio_build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
       case("matsubara","mats")
          call get_ctrl_var(beta,"BETA")
          wfreq = pi/beta*(2*arange(1,Lfreq)-1)
       case("realaxis","real")
          call get_ctrl_var(wini,"WINI")
          call get_ctrl_var(wfin,"WFIN")
          wfreq = linspace(wini,wfin,Lfreq)
       end select
    endif
    select case(axis)
    case default;stop "dmft_gfio_build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
    case("matsubara","mats");w_suffix="_iw"
    case("realaxis","real") ;w_suffix="_realw"
    end select
    return
  end subroutine build_frequency_array



end module DMFT_GFIO
