module TB_WANNIER90
  USE TB_COMMON
  USE TB_BASIS
  USE TB_IO
  USE TB_EFERMI
  implicit none



  type :: w90_structure
     character(len=:),allocatable               :: w90_file   !Name of the W90 file 
     integer                                    :: Num_wann=0 !Number of Wannier Orbitals
     integer                                    :: Nrpts=0    !Number of Wigner-Seitz cells
     character(len=20)                          :: DegFmt='(15I5)' !Ndegen format
     integer,allocatable,dimension(:)           :: Nlat            !Number of sites per element in the UC
     integer,allocatable,dimension(:)           :: Norb            !Number of orbitals per element in the UC
     integer,allocatable,dimension(:)           :: Nspin           !Number of spin per element in the UC
     integer                                    :: Ncomp=0         !Number of elements in the Unit Cell
     integer                                    :: Ntot=0          !Total number of lattice-spin-orbitals in the UC
     integer,allocatable,dimension(:)           :: Ndegen          !
     integer,allocatable,dimension(:,:)         :: Rvec   !Store the direct lattice basis
     real(8),allocatable,dimension(:,:)         :: Rgrid  !The real lattice grid 
     complex(8),allocatable,dimension(:,:,:)    :: Hij    !The H(Ri,Rj) Hamiltonian
     complex(8),allocatable,dimension(:,:)      :: Hloc   !The local part H(R_0,R_0)
     real(8),allocatable,dimension(:,:)         :: Kgrid  !The grid in the reciprocal wave vector lattice 
     real(8)                                    :: Efermi !The Fermi energy
     complex(8),allocatable,dimension(:,:)      :: Zeta   !The renormalization factor matrix H*=Z.(H+Self)
     real(8),allocatable,dimension(:,:)         :: Self   !The self-energy shift H*=Z.(H+Self)
     real(8),dimension(3)                       :: BZorigin !The origin of the BZ
     logical                                    :: Spinor =.false. !Spinor flag, spin included in W90 or not
     logical                                    :: Header =.true.  !Header flag, T=present in the file, F=not
     logical                                    :: iRenorm=.false. !H--> H*
     logical                                    :: iFermi =.false. !Fermi level evaluated
     logical                                    :: verbose=.false. !Verbosity level
     logical                                    :: hcheck =.true.  !Check Hermiticity upong getting H(k)
     logical                                    :: status =.false. !Allocated
  end type w90_structure


  type(w90_structure) :: TB_w90


contains




  !< Setup the default_w90 structure with information coming from specified w90_file
  subroutine setup_w90(w90_file,nlat,norb,nspin,origin,spinor,header,verbose,Hcheck)
    character(len=*),intent(in) :: w90_file
    integer                     :: Nlat(:)
    integer                     :: Norb(size(Nlat))
    integer,optional            :: Nspin
    real(8),optional            :: origin(:)
    logical,optional            :: spinor,verbose,Hcheck,header
    logical                     :: spinor_,verbose_,hcheck_,header_
    logical                     :: spin_blocks
    integer                     :: unitIO
    integer                     :: Num_wann
    integer                     :: Nrpts,Len,Nspin_,Ncomp
    integer                     :: i,j,ir,a,b
    integer                     :: rx,ry,rz
    real(8)                     :: re,im,origin_(3)
    !
    nspin_   = 1       ;if(present(nspin))nspin_    = nspin
    origin_  = 0d0     ;if(present(origin))origin_(:size(origin))=origin
    spinor_  = .false. ;if(present(spinor))spinor_  = spinor
    header_ = .true.   ;if(present(header))header_   = header
    verbose_ = .false. ;if(present(verbose))verbose_= verbose
    hcheck_  = .true.  ;if(present(hcheck))hcheck_  = hcheck
    !
    spin_blocks = (.not.spinor_).AND.(nspin_==2) !if Nspin=2 but spinor=F build H(k) with two independent spin-blocks
    !
    if(spinor_ .AND. Nspin_==1)then
       write(*,*)"WARNING setup_w90: Spinor=T AND Nspin==1 is not allowed. Setting Nspin=2"
       call sleep(2)
       Nspin_=2
    endif
    !
    if(spin_blocks)then
       write(*,*)"WARNING setup_w90: Spinor=F AND Nspin=2. Building two spin blocks: Num_wann**2 x SpinUp - Num_wann**2 x SpinDw"
       call sleep(2)
    endif
    !
    if(.not.set_eivec)then
       write(*,*)"WARNING setup_w90: set_eivec=false. Using default basis."
       call sleep(2)
    endif
    !
    call delete_w90
    !
    TB_w90%w90_file = str(w90_file)
    !
    !
    !Read from file initial info OR reconstruct them is Header is missing
    open(free_unit(unitIO),&
         file=TB_w90%w90_file,&
         status="old",&
         action="read")
    if(header_)then
       read(unitIO,*)          !skip first line
       read(unitIO,*) Num_wann !Number of Wannier orbitals
       read(unitIO,*) Nrpts    !Number of Wigner-Seitz vectors
    else
       Len      = file_length(TB_w90%w90_file)       
       Num_wann = sum(Nlat*Norb) !Nlat(1)*Norb(1) + ... + Nlat(Nineq)*Norb(Nineq)
       if(spinor_)Num_wann=Num_wann*Nspin_
       if( mod(len,Num_wann*Num_wann)/=0 )stop "ERROR setup_w90: no Header. Wrong Num_wann."
       Nrpts    = Len/Num_wann/Num_wann
    endif
    !
    !
    Ncomp = size(Nlat)
    TB_w90%Ncomp=Ncomp
    allocate(TB_w90%Nlat(Ncomp))
    allocate(TB_w90%Norb(Ncomp))
    allocate(TB_w90%Nspin(Ncomp))
    !
    !Setting structure data:
    !Spinor=F AND Nspin=2: build H(k) in the two spin-blocks. We assume Num_wann=Nlat*Norb as spinor=F. 
    !Spinor=T: Nspin must be 2. H(k) contains spin degrees of freedom, SOC terms or in-plane terms. Num_wann=Nlat*Norb*Nspin
    TB_w90%Nlat     = Nlat
    TB_w90%Norb     = Norb
    TB_w90%Nspin    = Nspin_
    TB_w90%Num_wann = Num_wann
    if(spinor_)then
       if(Num_wann /= sum(TB_w90%Nlat*TB_w90%Norb*TB_w90%Nspin))stop "setup_w90: Num_wann != sum(Nlat*Norb)*Nspin"
    else
       if(Num_wann /= sum(TB_w90%Nlat*TB_w90%Norb))stop "setup_w90: Num_wann != sum(Nlat*Norb)"
    endif
    TB_w90%Nrpts    = Nrpts
    TB_w90%Ntot     = sum(TB_w90%Nlat*TB_w90%Norb*TB_w90%Nspin)
    allocate(TB_w90%Ndegen(Nrpts))
    allocate(TB_w90%Rvec(Nrpts,3))
    allocate(TB_w90%Rgrid(Nrpts,3))
    allocate(TB_w90%Hij(  TB_w90%Ntot, TB_w90%Ntot, Nrpts))
    allocate(TB_w90%Hloc( TB_w90%Ntot, TB_w90%Ntot ))
    allocate(TB_w90%Zeta( TB_w90%Ntot, TB_w90%Ntot ))
    allocate(TB_w90%Self( TB_w90%Ntot, TB_w90%Ntot ))
    TB_w90%Ndegen   = 0
    TB_w90%Rvec     = 0
    TB_w90%Hij      = zero
    TB_w90%Hloc     = zero
    TB_w90%Zeta     = eye(TB_w90%Ntot )
    TB_w90%Self     = 0d0
    TB_w90%Efermi   = 0d0
    TB_w90%verbose  = Verbose_
    TB_w90%hcheck   = Hcheck_
    TB_w90%BZorigin = Origin_
    TB_w90%Spinor   = Spinor_
    TB_w90%Header   = Header_
    TB_w90%status   =.true.
    !
    !Read Ndegen from file: (as in original W90)
    if(TB_w90%Header)read(unitIO,str(TB_w90%DegFmt))(TB_w90%Ndegen(ir),ir=1,Nrpts)
    !
    !Read w90 TB Hamiltonian
    do ir=1,Nrpts
       do i=1,Num_wann
          do j=1,Num_wann
             !
             read(unitIO,*)rx,ry,rz,a,b,re,im
             !
             TB_w90%Rvec(ir,:)  = [rx,ry,rz]
             TB_w90%Rgrid(ir,:) = rx*ei_x + ry*ei_y + rz*ei_z
             !
             TB_w90%Hij(a,b,ir)=dcmplx(re,im)
             if(spin_blocks)TB_w90%Hij(a+num_wann, b+num_wann, ir)=dcmplx(re,im)

             !
             if( all(TB_w90%Rvec(ir,:)==0) )then
                TB_w90%Hloc(a,b)=dcmplx(re,im)
                if(spin_blocks)TB_w90%Hloc(a+num_wann, b+num_wann)=dcmplx(re,im)
             endif
          enddo
       enddo
    enddo
    close(unitIO)
    !
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    if(mpi_master)then
       write(*,*)
       write(*,'(1A)')         "-------------- H_LDA --------------"
       write(*,'(A,I6)')      "  Number of Wannier functions:   ",TB_w90%num_wann
       write(*,'(A,I6)')      "  Number of Wigner-Seitz vectors:",TB_w90%nrpts
       write(*,'(A,I6)')      "  Consider re-ordering Hij/Hk with TB_reshape:"
    endif
    !
  end subroutine setup_w90


  subroutine delete_w90()
    if(.not.TB_w90%status)return
    TB_w90%Num_wann =0
    TB_w90%Nrpts    =0
    TB_w90%DegFmt   =''
    TB_w90%Nspin    =0
    TB_w90%Ntot     =0
    TB_w90%Efermi   =0d0
    TB_w90%verbose  =.false.
    TB_w90%spinor   =.false.
    TB_w90%header   =.false.
    deallocate(TB_w90%Nlat)
    deallocate(TB_w90%Norb)
    deallocate(TB_w90%w90_file)
    deallocate(TB_w90%Ndegen)
    deallocate(TB_w90%Rvec)
    deallocate(TB_w90%Rgrid)
    deallocate(TB_w90%Hij)
    deallocate(TB_w90%Hloc)
    TB_w90%status=.false.
  end subroutine delete_w90




  subroutine transform_Hij_w90_1R(transform_Hij)
    interface 
       function transform_Hij(Hij,N)
         integer                   :: N
         complex(8),dimension(N,N) :: Hij
         complex(8),dimension(N,N) :: transform_Hij
       end function transform_Hij
    end interface
    complex(8),dimension(:,:),allocatable :: Htmp
    integer                               :: i,j,ir,a,b
    !
    !Read w90 TB Hamiltonian
    allocate(Htmp(TB_w90%Ntot,TB_w90%Ntot))
    !
    do ir=1,TB_w90%Nrpts
       !
       !Copy H(R) into a tmp array
       Htmp               = TB_w90%Hij(:,:,ir)
       !Apply transformation on the single \bar{H}(R) = T(H(R))
       TB_w90%Hij(:,:,ir) = transform_Hij(Htmp,TB_w90%Ntot)
       !
       if( all(TB_w90%Rvec(ir,:)==0) )then
          Htmp        = TB_w90%Hloc
          TB_w90%Hloc = transform_Hij(Htmp,TB_w90%Ntot)
       endif
    enddo
    !
    deallocate(Htmp)
  end subroutine transform_Hij_w90_1R


  subroutine transform_Hij_w90_allR(transform_Hij)
    interface 
       function transform_Hij(Hij,N,Nr)
         integer                      :: N,Nr
         complex(8),dimension(N,N,Nr) :: Hij
         complex(8),dimension(N,N,Nr) :: transform_Hij
       end function transform_Hij
    end interface
    complex(8),dimension(:,:,:),allocatable :: Htmp
    integer                                 :: i,j,ir,a,b
    !
    !Read w90 TB Hamiltonian
    allocate(Htmp(TB_w90%Ntot,TB_w90%Ntot,TB_w90%Nrpts))
    !
    !Copy H into a tmp array
    Htmp               = TB_w90%Hij
    !Apply transformation on the single \bar{H}(R) = T(H(R))
    TB_w90%Hij = transform_Hij(Htmp,TB_w90%Ntot,TB_w90%Nrpts)
    !
    deallocate(Htmp)
  end subroutine transform_Hij_w90_allR



  !< generate an internal function to be called in the building procedure
  function w90_hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: Hk,Hk_f
    complex(8)                :: kExp
    integer                   :: nDim
    integer                   :: i,j,ir
    real(8)                   :: rvec(size(kvec)),rdotk,WSdeg
    !
    if(.not.TB_w90%status)stop "w90_hk_model: TB_w90 was not setup"
    !
    nDim = size(kvec)
    if(TB_w90%Ntot /= N )stop "w90_hk_model: Nlso != N"
    Hk = zero
    do ir=1,TB_w90%Nrpts
       Rvec  = TB_w90%Rgrid(ir,:nDim)
       rdotk = dot_product(Kvec,Rvec)
       kExp  = dcmplx(cos(rdotk),-sin(rdotk))
       WSdeg = TB_w90%Ndegen(ir)
       !
       Hk = Hk+TB_w90%Hij(:,:,ir)*kExp/WSdeg
    enddo
    !
    where(abs(Hk)<1d-9)Hk=zero
    !
    if(TB_w90%hcheck)call herm_check(Hk)
    if(TB_w90%Ifermi)Hk = Hk - TB_w90%Efermi*eye(N)
    !
    if(TB_w90%Irenorm)then
       Hk   = Hk - TB_w90%Hloc
       Hk_f = (TB_w90%Zeta .x. Hk) .x. TB_w90%Zeta
       Hk   = Hk_f + TB_w90%Hloc + TB_w90%Self
    endif
    !
  end function w90_hk_model




  subroutine build_hk_w90(Hk,N,Nkvec,Kpts_grid)
    integer                                                :: N
    integer,dimension(:),intent(in)                        :: Nkvec
    complex(8),dimension(N,N,product(Nkvec))               :: Hk,haux
    real(8),dimension(product(Nkvec),size(Nkvec)),optional :: Kpts_grid ![Nk][Ndim]
    !
    real(8),dimension(product(Nkvec),size(Nkvec))          :: kgrid ![Nk][Ndim]
    integer                                                :: Nk
    logical                                                :: IOfile
    integer                                                :: i,j,ik
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
    Nk =product(Nkvec)
    !
    if(.not.TB_w90%status)stop "ERROR build_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
    if(N/=TB_w90%Ntot)stop "ERROR build_hk_w90: incorrect N."
    !
    if(allocated(TB_w90%Kgrid))deallocate(TB_w90%Kgrid)
    allocate(TB_w90%Kgrid(Nk,size(Nkvec)))
    call build_kgrid(Nkvec,Kgrid,TB_w90%BZorigin) !check bk_1,2,3 vectors have been set
    TB_w90%Kgrid=Kgrid
    if(present(Kpts_grid))Kpts_grid=Kgrid
    !
    Haux=zero
    do ik=1+mpi_rank,Nk,mpi_size
       haux(:,:,ik) = w90_hk_model(Kgrid(ik,:),N)
    enddo
#ifdef _MPI
    if(check_MPI())then
       Hk = zero
       call AllReduce_MPI(MPI_COMM_WORLD,Haux,Hk)
    else
       Hk = Haux
    endif
#else
    Hk = Haux
#endif
    !
  end subroutine build_hk_w90




  subroutine build_hk_w90_path(Hk,N,Kpath,Nkpath)
    integer                                                   :: N
    real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
    integer                                                   :: Nkpath
    integer                                                   :: Npts,Nktot
    complex(8),dimension(N,N,(size(kpath,1)-1)*Nkpath)        :: Hk,haux
    !
    real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
    integer                                                   :: Nk
    logical                                                   :: IOfile
    integer                                                   :: i,j,ik
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
    if(.not.TB_w90%status)stop "ERROR build_hk_w90_path: TB_w90 structure not allocated. Call setup_w90 first."
    if(N/=TB_w90%Ntot)stop "ERROR build_hk_w90_path: incorrect N."
    !
    Npts  =  size(kpath,1)          !# of k-points along the path
    Nktot = (Npts-1)*Nkpath
    !
    call kgrid_from_path_grid(kpath,Nkpath,kgrid)
    !
    Haux=zero
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = w90_hk_model(Kgrid(ik,:),N)
    enddo
#ifdef _MPI
    if(check_MPI())then
       Hk = zero
       call AllReduce_MPI(MPI_COMM_WORLD,Haux,Hk)
    else
       Hk = Haux
    endif
#else
    Hk = Haux
#endif
    !
  end subroutine build_hk_w90_path



  subroutine Hloc_w90(Hloc)
    complex(8),dimension(:,:) :: Hloc
    call assert_shape(Hloc,[TB_w90%Ntot,TB_w90%Ntot],"Hloc_w90","Hloc")
    Hloc = TB_w90%Hloc
    if(TB_w90%Ifermi)Hloc = Hloc - TB_w90%Efermi*eye(TB_w90%Ntot)
    !
  end subroutine Hloc_w90

  subroutine FermiLevel_w90(Nkvec,filling,Ef)
    integer,dimension(:),intent(in)          :: Nkvec
    real(8)                                  :: filling,Efermi
    real(8),optional                         :: Ef
    complex(8),dimension(:,:,:),allocatable  :: Hk
    integer                                  :: Nk
    if(TB_w90%Ifermi)return
    Nk   = product(Nkvec)
    allocate(Hk(TB_w90%Ntot,TB_w90%Ntot,Nk))
    call build_hk_w90(Hk,TB_w90%Ntot,Nkvec)
    call TB_FermiLevel(Hk,filling,Efermi,TB_w90%Nspin(1),TB_w90%verbose)
    TB_w90%Efermi = Efermi
    if(present(Ef))Ef=Efermi
    TB_w90%Ifermi=.true.
  end subroutine FermiLevel_w90

  subroutine Zeta_w90_vector(zeta)
    real(8),dimension(:)          :: zeta
    real(8),dimension(size(zeta)) :: sq_zeta
    call assert_shape(zeta,[TB_w90%Ntot],"Zeta_w90","zeta")
    sq_zeta     = sqrt(zeta)
    TB_w90%zeta = diag(sq_zeta)
    TB_w90%Irenorm=.true.
  end subroutine Zeta_w90_vector

  subroutine Zeta_w90_matrix(zeta)
    real(8),dimension(:,:) :: zeta
    call assert_shape(zeta,[TB_w90%Ntot,TB_w90%Ntot],"Zeta_w90","zeta")  
    TB_w90%zeta = sqrt(zeta)
    TB_w90%Irenorm=.true.
  end subroutine Zeta_w90_matrix


  subroutine Self0_w90(Self0)
    real(8),dimension(:,:) :: self0
    call assert_shape(self0,[TB_w90%Ntot,TB_w90%Ntot],"Self0_w90","Self0")  
    TB_w90%Self = Self0
    TB_w90%Irenorm=.true.
  end subroutine Self0_w90










  subroutine write_hk_w90_func(file,Nkvec)
    character(len=*)                              :: file
    integer                                       :: Nkvec(:)
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
    real(8),dimension(3)                          :: kvec
    integer                                       :: Nktot,unit,Dim
    integer                                       :: i,ik,io,jo
    complex(8),dimension(:,:),allocatable         :: hk
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(.not.TB_w90%status)stop "read_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
    !
    if(.not.allocated(TB_w90%Kgrid))then
       call build_kgrid(Nkvec,Kgrid,TB_w90%BZorigin) !check bk_1,2,3 vectors have been set
    else
       call assert_shape(TB_w90%Kgrid,[product(Nkvec),size(Nkvec)],"write_hk_w90","TB_w90%Kgrid")
       Kgrid = TB_w90%Kgrid
    endif
    !
    Dim   = size(Nkvec)
    Nktot = product(Nkvec)
    !
    if(mpi_master)then
       open(free_unit(unit),file=reg(file))
       write(unit,'('//str(1+3*TB_w90%Ncomp)//'(I10,1x),F15.9)')Nktot,TB_w90%Nlat,TB_w90%Nspin,TB_w90%Norb,TB_w90%Efermi
       write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,Dim)
       allocate(Hk(TB_w90%Ntot,TB_w90%Ntot))
       do ik=1,Nktot
          Hk(:,:) = w90_hk_model(Kgrid(ik,:),TB_w90%Ntot)
          kvec=0d0 ; kvec(:Dim) = kgrid(ik,:)
          write(unit,"(3(F15.9,1x))")(kvec(i),i=1,3) 
          do io=1,TB_w90%Ntot
             write(unit,"(1000(F15.9))")(dreal(Hk(io,jo)),jo=1,TB_w90%Ntot)
          enddo
          do io=1,TB_w90%Ntot
             write(unit,"(1000(F15.9))")(dimag(Hk(io,jo)),jo=1,TB_w90%Ntot)
          enddo
       enddo
       close(unit)
    endif
    !
  end subroutine write_hk_w90_func

  subroutine write_hk_w90_array(Hk,file,Nkvec)
    complex(8),dimension(:,:,:)                   :: Hk
    character(len=*)                              :: file
    integer                                       :: Nkvec(:)
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
    real(8),dimension(3)                          :: kvec
    integer                                       :: Nktot,unit,Dim
    integer                                       :: i,ik,io,jo  

    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(.not.TB_w90%status)stop "read_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
    !
    if(.not.allocated(TB_w90%Kgrid))then
       call build_kgrid(Nkvec,Kgrid,TB_w90%BZorigin) !check bk_1,2,3 vectors have been set
    else
       call assert_shape(TB_w90%Kgrid,[product(Nkvec),size(Nkvec)],"write_hk_w90","TB_w90%Kgrid")
       Kgrid = TB_w90%Kgrid
    endif
    !
    Dim   = size(Nkvec)
    Nktot = product(Nkvec)
    !
    call assert_shape(Hk,[TB_w90%Ntot,TB_w90%Ntot,Nktot],"TB_WANNIER/write_hk_w90_array","Hk")
    !
    if(mpi_master)then
       open(free_unit(unit),file=reg(file))
       write(unit,'('//str(1+3*TB_w90%Ncomp)//'(I10,1x),F15.9)')Nktot,TB_w90%Nlat,TB_w90%Nspin,TB_w90%Norb,TB_w90%Efermi
       write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,Dim)
       do ik=1,Nktot
          kvec=0d0 ; kvec(:Dim) = kgrid(ik,:)
          write(unit,"(3(F15.9,1x))")(kvec(i),i=1,3) 
          do io=1,TB_w90%Ntot
             write(unit,"(1000(F15.9))")(dreal(Hk(io,jo,ik)),jo=1,TB_w90%Ntot)
          enddo
          do io=1,TB_w90%Ntot
             write(unit,"(1000(F15.9))")(dimag(Hk(io,jo,ik)),jo=1,TB_w90%Ntot)
          enddo
       enddo
       close(unit)
    endif
    !
  end subroutine write_hk_w90_array




  subroutine read_hk_w90_array(Hk,file,Nkvec)
    complex(8),dimension(:,:,:),allocatable :: Hk
    character(len=*)                        :: file
    integer,intent(inout)                   :: Nkvec(:)
    real(8),dimension(:,:),allocatable      :: kgrid
    real(8),dimension(3)                    :: kvec
    integer                                 :: Nktot,unit,Nk(3),Dim
    integer                                 :: i,ik,ix,iy,iz,io,jo
    real(8)                                 :: kx,ky,kz,Ef
    logical                                 :: ioexist
    character(len=1)                        :: achar
    real(8),dimension(:,:),allocatable      :: reH,imH
    !
    if(.not.TB_w90%status)stop "read_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
    !
    inquire(file=reg(file),exist=ioexist)
    if(.not.ioexist)then
       write(*,*)"can not find file:"//reg(file)
       stop
    endif
    !
    open(free_unit(unit),file=reg(file))
    read(unit,'('//str(1+3*TB_w90%Ncomp)//'(I10,1x),F15.9)')Nktot,TB_w90%Nlat,TB_w90%Nspin,TB_w90%Norb,TB_w90%Efermi
    read(unit,'(1A1,3(I12,1x))')achar,( Nk(ik),ik=1,3 )
    TB_w90%Ifermi =.true.
    !
    Dim    = size(Nkvec)
    Nkvec  = Nk(:Dim)
    if(Nktot  /= product(Nk))stop "read_hk_w90: Nktot != product(Nk)"
    Nktot  = product(Nk)
    !
    allocate(Kgrid(Nktot,Dim))
    allocate(Hk(TB_w90%Ntot,TB_w90%Ntot,Nktot))
    allocate(reH(TB_w90%Ntot,TB_w90%Ntot),imH(TB_w90%Ntot,TB_w90%Ntot))
    !
    ik=0
    do iz=1,Nk(3)
       do iy=1,Nk(2)
          do ix=1,Nk(1)
             ik = ik+1
             read(unit,"(3(F15.9,1x))")kx,ky,kz
             kvec = [kx,ky,kz]
             kgrid(ik,:) = kvec(:Dim)
             do io=1,TB_w90%Ntot
                read(unit,"(1000(F15.9))")(reH(io,jo),jo=1,TB_w90%Ntot)
             enddo
             do io=1,TB_w90%Ntot
                read(unit,"(1000(F15.9))")(imH(io,jo),jo=1,TB_w90%Ntot)
             enddo
             Hk(:,:,ik) = dcmplx(reH,imH)
          enddo
       enddo
    enddo
    close(unit)
    !
  end subroutine read_hk_w90_array






















  subroutine fix_w90file(file,nrpts,num_wann)
    character(len=*)             :: file
    character(len=:),allocatable :: file_bkp
    integer                      :: nrpts,num_wann,len
    integer                      :: ubkp,unit
    integer                      :: ir,i,j
    logical                      :: bool
    integer,allocatable          :: ndegen(:)
    integer                      :: rx,ry,rz,a,b
    real(8)                      :: re,im
    character(len=33)           :: header
    character(len=9)             :: cdate, ctime
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    if(mpi_master)then
       inquire(file=reg(file),exist=bool)
       if(.not.bool)stop "TB_fix_w90: file not found"
       !
       len = file_length(reg(file))
       if(Nrpts/=len/Num_wann/Num_wann)stop "TB_fix_w90: wrong file length"
       !
       allocate(ndegen(Nrpts));Ndegen=1
       !
       file_bkp = reg(file)//".backup"
       !
       call system("cp -vf "//reg(file)//" "//reg(file_bkp))
       open(free_unit(ubkp),file=reg(file_bkp))
       open(free_unit(unit),file=reg(file))
       !
       call io_date(cdate,ctime)
       header = 'fixed on '//cdate//' at '//ctime
       write(unit,*)header ! Date and time
       write(unit,*)num_wann
       write(unit,*)nrpts
       write(unit,"(15I5)")(ndegen(i), i=1, nrpts)
       !
       do ir=1,Nrpts
          do i=1,Num_wann
             do j=1,Num_wann
                !
                read(ubkp,*)rx,ry,rz,a,b,re,im
                write(unit,"(5I5,2F12.6)")rx,ry,rz,a,b,re,im
             enddo
          enddo
       enddo
       close(ubkp)
       close(unit)
    endif
#ifdef _MPI    
    if(check_MPI())call MPI_Barrier(MPI_COMM_WORLD,mpi_ierr)
#endif
  contains
    !From Wannier90
    subroutine io_date(cdate, ctime)
      character(len=9), intent(out)   :: cdate
      character(len=9), intent(out)   :: ctime
      character(len=3), dimension(12) :: months
      data months/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
           'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
      integer                         :: date_time(8)
      !
      call date_and_time(values=date_time)
      !
      write (cdate, '(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
      write (ctime, '(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)
    end subroutine io_date
  end subroutine fix_w90file






  !OLD STUFF:
  !< read the real space hopping matrix from Wannier90 output and create H(k)
  include "w90hr/tight_binding_build_hk_from_w90hr.f90"
#ifdef _MPI
  include "w90hr/tight_binding_build_hk_from_w90hr_mpi.f90"
#endif
  !< read the real space lattice position and compute the space integrals
  !          <t2g| [x,y,z] |t2g> using atomic orbital as a subset.
  include "w90hr/dipole_w90hr.f90"
#ifdef _MPI
  include "w90hr/dipole_w90hr_mpi.f90"
#endif

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




  subroutine herm_check(A)
    complex(8),intent(in) ::   A(:,:)
    integer               ::   row,col,i,j
    row=size(A,1)
    col=size(A,2)
    do i=1,col
       do j=1,row
          if(abs(A(i,j))-abs(A(j,i)) .gt. 1d-6)then
             write(*,'(1A)') "--> NON HERMITIAN MATRIX <--"
             write(*,'(2(1A7,I3))') " row: ",i," col: ",j
             write(*,'(1A)') "  A(i,j)"
             write(*,'(2F22.18)') real(A(i,j)),aimag(A(i,j))
             write(*,'(1A)') "  A(j,i)"
             write(*,'(2F22.18)') real(A(j,i)),aimag(A(j,i))
             stop
          endif
       enddo
    enddo
  end subroutine herm_check



END MODULE TB_WANNIER90
















