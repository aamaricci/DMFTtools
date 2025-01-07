module GF_IO
  USE GF_COMMON
  implicit none
  private


  !DMFT INTERFACE:
  interface dmft_write_gf
     module procedure :: gf_push_zeta
     module procedure :: dmft_gf_print_rank2
     module procedure :: dmft_gf_print_rank3
     module procedure :: dmft_gf_print_rank4
     module procedure :: dmft_gf_print_rank5
     module procedure :: dmft_gf_print_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_cluster_ineq
#endif
  end interface dmft_write_gf

  !AGNOSTIC INTERFACE:
  interface write_gf
     module procedure :: gf_push_zeta
     module procedure :: dmft_gf_print_rank2
     module procedure :: dmft_gf_print_rank3
     module procedure :: dmft_gf_print_rank4
     module procedure :: dmft_gf_print_rank5
     module procedure :: dmft_gf_print_rank6
#if __GFORTRAN__ &&  __GNUC__ > 8
     module procedure :: dmft_gf_print_cluster_ineq
#endif
  end interface write_gf


  public :: dmft_write_gf
  public :: write_gf

contains



  subroutine dmft_gf_print_rank2(Func,fname,axis,iprint,Nvec,labels)
    complex(8),dimension(:,:,:),intent(in)     :: Func
    character(len=*),intent(in)                :: fname
    character(len=*)                           :: axis
    integer,intent(in),optional                :: iprint
    integer                                    :: iprint_
    integer,dimension(:),optional              :: Nvec
    character(len=*),dimension(:),optional     :: labels
    integer                                    :: istride,jstride,Mdim
    integer,dimension(:),allocatable           :: Mvec
    character(len=64),dimension(:),allocatable :: Mlabels
    !
    iprint_=10;if(present(iprint))iprint_=iprint
    !
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
    if(present(Nvec))then
       if(.not.present(labels))stop 'dmft_gfio: present(Nvec)=T but  present(labels)=F'
       if(sum(Nvec)/=Ntot)stop "dmft_gfio: sum(Nvec) != Ntot"
       Mdim = size(Nvec)
       if(size(labels)/=Mdim)stop "dmft_gfio: size(labels) != size(Nvec)"
       allocate(Mvec(Mdim))
       allocate(Mlabels(Mdim))
       Mvec    = Nvec
       do iel=1,Mdim
          Mlabels(iel) = str(labels(iel))
       enddo
    else
       Mdim = 1
       allocate(Mvec(Mdim))
       allocate(Mlabels(Mdim))
       Mvec(1) = Ntot
       Mlabels(1) = str('i')
    endif
    !
    if(mpi_master)then
       select case(iprint_)
       case(1)
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write diagonal elements."
          istride=0
          do iel=1,Mdim
             do io=1,Mvec(iel)
                suffix=reg(fname)//"_"//&
                     str(Mlabels(iel))//str(io)//&
                     str(w_suffix)//str(gf_suffix)
                call splot(reg(suffix),wio,Func(istride+io,istride+io,:))
             enddo
             istride=istride+Mvec(iel)
          enddo
       case default
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write all elements."
          istride=0
          jstride=0
          do iel=1,Mdim
             do io=1,Mvec(iel)
                !
                do jel=1,Mdim
                   do jo=1,Mvec(jel)
                      suffix=reg(fname)//"_"//&
                           str(Mlabels(iel))//str(io)//&
                           str(Mlabels(jel))//str(jo)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wio,Func(istride+io,jstride+jo,:))
                   enddo
                   jstride=jstride+Mvec(jel)
                enddo
                !
             enddo
             istride=istride+Mvec(iel)
          enddo
          !
       end select
    endif
    if(allocated(wio))deallocate(wio)
  end subroutine dmft_gf_print_rank2


  subroutine dmft_gf_print_rank3(Func,fname,axis,iprint,ineq_index,ineq_pad,itar)
    complex(8),dimension(:,:,:,:),intent(in) :: Func
    character(len=*),intent(in)              :: fname
    character(len=*)                         :: axis
    integer,intent(in),optional              :: iprint
    integer                                  :: iprint_
    character(len=*),optional                :: ineq_index
    character(len=:),allocatable             :: index
    integer,optional                         :: ineq_pad
    integer                                  :: pad
    logical,optional                         :: itar
    logical                                  :: itar_
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
    Nso   = size(Func,2)
    Lfreq = size(Func,4)
    call assert_shape(Func,[Nlat,Nso,Nso,Lfreq],"dmft_gf_print_rank3","Func")
    !
    call build_frequency_array(axis)
    !
    if(mpi_master)then
       select case(iprint_)
       case(1)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin-orbital diagonal elements. Single File."
          do io=1,Nso
             suffix=reg(fname)//&
                  "_io"//str(io)//"jo"//str(io)//&
                  str(w_suffix)//reg(gf_suffix)
             call splot(reg(suffix),wio,Func(:,io,jo,:))
             call file_bzip(reg(suffix))
          enddo
          !
       case(2,3)                  !print all off-diagonals
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: all elements. Single File."
          do io=1,Nso
             do jo=1,Nso
                suffix=reg(fname)//&
                     "_io"//str(io)//"jo"//str(jo)//&
                     str(w_suffix)//reg(gf_suffix)
                call splot(reg(suffix),wio,Func(:,io,jo,:))
                call file_bzip(reg(suffix))
             enddo
          enddo
          !
       case(4)                  !print only diagonal elements
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write spin-orbital diagonal elements."
          do io=1,Nso
             do ilat=1,Nlat
                suffix=reg(fname)//&
                     "_io"//str(io)//"jo"//str(io)//&
                     str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                call splot(reg(suffix),wio,Func(ilat,io,io,:))
             enddo
             suffix=reg(fname)//&
                  "_io"//str(io)//"io"//str(io)//&
                  str(w_suffix)
             if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
                  pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
          enddo
          !
       case default
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: write all elements."
          do io=1,Nso
             do jo=1,Nso
                do ilat=1,Nlat
                   suffix=reg(fname)//&
                        "_io"//str(io)//"jo"//str(jo)//&
                        str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call splot(reg(suffix),wio,Func(ilat,io,jo,:))
                enddo
                suffix=reg(fname)//&
                     "_io"//str(io)//"jo"//str(jo)//&
                     str(w_suffix)
                if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
                     pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wio))deallocate(wio)
  end subroutine dmft_gf_print_rank3




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
                     "_l"//str(iorb)//"m"//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)//str(gf_suffix)
                call splot(reg(suffix),wio,Func(ispin,ispin,iorb,iorb,:))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)") reg(fname),"dmft_gfio: write spin diagonal and all orbitals elements."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//"m"//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(gf_suffix)
                   call splot(reg(suffix),wio,Func(ispin,ispin,iorb,jorb,:))
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
                           "_l"//str(iorb)//"m"//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wio,Func(ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wio))deallocate(wio)
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
                     "_l"//str(iorb)//"m"//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)//reg(gf_suffix)
                call splot(reg(suffix),wio,Func(:,ispin,ispin,iorb,iorb,:))
                call file_bzip(reg(suffix))
             enddo
          enddo
          !
       case(2)                  !print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin diagonal and all orbitals elements. Single File."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//"m"//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(gf_suffix)
                   call splot(reg(suffix),wio,Func(:,ispin,ispin,iorb,jorb,:))
                   call file_bzip(reg(suffix))
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
                           "_l"//str(iorb)//"m"//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wio,Func(:,ispin,jspin,iorb,jorb,:))
                      call file_bzip(reg(suffix))
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
                        "_l"//str(iorb)//"m"//str(iorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                   call splot(reg(suffix),wio,Func(ilat,ispin,ispin,iorb,iorb,:))
                enddo
                suffix=reg(fname)//&
                     "_l"//str(iorb)//"m"//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)
                if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
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
                           "_l"//str(iorb)//"m"//str(jorb)//&
                           "_s"//str(ispin)//&
                           str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                      call splot(reg(suffix),wio,Func(ilat,ispin,ispin,iorb,jorb,:))
                   enddo
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//"m"//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)
                   if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
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
                              "_l"//str(iorb)//"m"//str(jorb)//&
                              "_s"//str(ispin)//str(jspin)//&
                              str(w_suffix)//reg(index)//str(ilat,pad)//reg(gf_suffix)
                         call splot(reg(suffix),wio,Func(ilat,ispin,jspin,iorb,jorb,:))
                      enddo
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//"m"//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)
                      if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
                           pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
          !
       end select
    endif
    if(allocated(wio))deallocate(wio)
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
                     "_l"//str(iorb)//"m"//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)//reg(gf_suffix)
                call splot(reg(suffix),wio,Func(:,:,ispin,ispin,iorb,iorb,:))
                call file_bzip(reg(suffix))
             enddo
          enddo
          !
       case(2)!print spin-diagonal, all orbitals 
          write(*,"(A,1x,A)")reg(fname),"dmft_gfio: spin diagonal and all orbitals elements. Single File."
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//"m"//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)//reg(gf_suffix)
                   call splot(reg(suffix),wio,Func(:,:,ispin,ispin,iorb,jorb,:))
                   call file_bzip(reg(suffix))
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
                           "_l"//str(iorb)//"m"//str(jorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)//reg(gf_suffix)
                      call splot(reg(suffix),wio,Func(:,:,ispin,jspin,iorb,jorb,:))
                      call file_bzip(reg(suffix))
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
                           "_l"//str(iorb)//"m"//str(iorb)//&
                           "_s"//str(ispin)//&
                           str(w_suffix)//reg(index)//&
                           str(ilat,pad)//"_"//str(jlat,pad)//reg(gf_suffix)
                      call splot(reg(suffix),wio,Func(ilat,jlat,ispin,ispin,iorb,iorb,:))
                   enddo
                enddo
                suffix=reg(fname)//&
                     "_l"//str(iorb)//"m"//str(iorb)//&
                     "_s"//str(ispin)//&
                     str(w_suffix)
                if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
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
                              "_l"//str(iorb)//"m"//str(iorb)//&
                              "_s"//str(ispin)//&
                              str(w_suffix)//reg(index)//&
                              str(ilat,pad)//"_"//str(jlat,pad)//reg(gf_suffix)
                         call splot(reg(suffix),wio,Func(ilat,jlat,ispin,ispin,iorb,jorb,:))
                      enddo
                   enddo
                   suffix=reg(fname)//&
                        "_l"//str(iorb)//"m"//str(jorb)//&
                        "_s"//str(ispin)//&
                        str(w_suffix)
                   if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
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
                                 "_l"//str(iorb)//"m"//str(iorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 str(w_suffix)//reg(index)//&
                                 str(ilat,pad)//"_"//str(jlat,pad)//reg(gf_suffix)
                            call splot(reg(suffix),wio,Func(ilat,jlat,ispin,jspin,iorb,jorb,:))
                         enddo
                      enddo
                      suffix=reg(fname)//&
                           "_l"//str(iorb)//"m"//str(iorb)//&
                           "_s"//str(ispin)//str(jspin)//&
                           str(w_suffix)
                      if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
                           pattern=reg(suffix)//reg(index)//"*"//reg(gf_suffix))
                   enddo
                enddo
             enddo
          enddo
       end select
    endif
    if(allocated(wio))deallocate(wio)
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
                                 "_i"//str(ilat)//"j"//str(jlat)//&
                                 "_l"//str(iorb)//"m"//str(jorb)//&
                                 "_s"//str(ispin)//&
                                 str(w_suffix)//reg(index)//str(iineq,pad)//reg(gf_suffix)
                            call splot(reg(suffix),wio,Func(iineq,ilat,jlat,ispin,ispin,iorb,jorb,:))
                         enddo
                         suffix=reg(fname)//&
                              "_i"//str(ilat)//"j"//str(jlat)//&
                              "_l"//str(iorb)//"m"//str(jorb)//&
                              "_s"//str(ispin)//&
                              str(w_suffix)
                         if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
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
                                    "_i"//str(ilat)//"j"//str(jlat)//&
                                    "_l"//str(iorb)//"m"//str(jorb)//&
                                    "_s"//str(ispin)//str(jspin)//&
                                    str(w_suffix)//reg(index)//str(iineq,pad)//reg(gf_suffix)
                               call splot(reg(suffix),wio,Func(iineq,ilat,jlat,ispin,jspin,iorb,jorb,:))
                            enddo
                            suffix=reg(fname)//&
                                 "_i"//str(ilat)//"j"//str(jlat)//&
                                 "_l"//str(iorb)//"m"//str(jorb)//&
                                 "_s"//str(ispin)//str(jspin)//&
                                 str(w_suffix)
                            if(itar_)call file_targz(tarball="tar_"//reg(suffix),&
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
    if(allocated(wio))deallocate(wio)
  end subroutine dmft_gf_print_cluster_ineq
#endif













  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################


  ! subroutine set_gf_suffix(string)
  !   character(len=*) :: string
  !   gf_suffix=reg(string)
  ! end subroutine set_gf_suffix



  ! subroutine dmft_gf_push_zeta(zeta)
  !   real(8),dimension(:) :: zeta
  !   Lfreq=size(zeta)
  !   if(allocated(wio))deallocate(wio)
  !   allocate(wio(Lfreq))
  !   wio=zeta
  ! end subroutine dmft_gf_push_zeta


  ! subroutine build_frequency_array(axis)
  !   character(len=*)              :: axis
  !   if(allocated(wio))then
  !      if(size(wio)/=Lfreq)stop "dmft_gfio_build_frequency_array ERROR: pushed wio has wrong size"
  !   else
  !      allocate(wio(Lfreq))
  !      select case(axis)
  !      case default;
  !         stop "dmft_gfio_build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
  !      case("matsubara","mats","m")
  !         call get_ctrl_var(beta,"BETA")
  !         wio = pi/beta*(2*arange(1,Lfreq)-1)
  !      case("realaxis","real","r")
  !         call get_ctrl_var(wini,"WINI")
  !         call get_ctrl_var(wfin,"WFIN")
  !         wio = linspace(wini,wfin,Lfreq)
  !      end select
  !   endif
  !   select case(axis)
  !   case default;stop "dmft_gfio_build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
  !   case("matsubara","mats");w_suffix="_iw"
  !   case("realaxis","real") ;w_suffix="_realw"
  !   end select
  !   return
  ! end subroutine build_frequency_array



end module GF_IO
