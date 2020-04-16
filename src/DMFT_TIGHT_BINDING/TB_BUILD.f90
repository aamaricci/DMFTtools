module TB_BUILD
  USE TB_COMMON
  USE TB_IO
  USE TB_BASIS
  USE TB_WANNIER90
  implicit none


  interface TB_build_model
     module procedure :: build_hk_model_kgrid
     module procedure :: build_hk_model_nkvec
     module procedure :: build_hkR_model_kgrid
     module procedure :: build_hkR_model_nkvec
     module procedure :: build_hk_path
     module procedure :: build_hkR_path
     module procedure :: build_hk_w90
     !
     module procedure :: build_Hij_Nrvec
  end interface TB_build_model


  interface TB_w90_FermiLevel
     module procedure :: FermiLevel_w90
  end interface TB_w90_FermiLevel



  interface TB_hr_to_hk
     module procedure :: hk_from_w90_hr
     module procedure :: hloct_from_w90_hr
#ifdef _MPI
     module procedure :: hk_from_w90_hr_mpi
     module procedure :: hkt_from_w90_hr_mpi
#endif
  end interface TB_hr_to_hk


  interface TB_dipole
     module procedure :: dipole_t2g_LDA
#ifdef _MPI
     module procedure :: dipole_t2g_LDA_mpi
#endif
  end interface TB_dipole



contains




  !< build the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) Hamiltonian matrix
  ! from the function user defined hk_model procedure.
  subroutine build_hk_model_kgrid(Hk,hk_model,Norb,kgrid,wdos)
    integer                                       :: Norb
    integer                                       :: Nktot
    real(8),dimension(:,:)                        :: kgrid ![Nktot][Ndim]
    real(8)                                       :: kvec(size(kgrid,2))
    integer                                       :: ik
    complex(8),dimension(Norb,Norb,size(kgrid,1)) :: hk,haux
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)                     :: kpoint
         integer                                  :: N
         complex(8),dimension(N,N)                :: hk_model
       end function hk_model
    end interface
    logical,optional                              :: wdos
    logical                                       :: wdos_
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
    wdos_=.false.;if(present(wdos))wdos_=wdos
    !
    Nktot  = size(kgrid,1)
    !
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = hk_model(kgrid(ik,:),Norb)
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
    if(wdos_)then
       allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
       allocate(dos_wtk(Nktot))
       dos_wtk=1d0/Nktot
#ifdef _MPI
       if(check_MPI())then
          call dmft_gloc_realaxis(MPI_COMM_WORLD,Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
       else
          call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
       endif
#else
       call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
#endif
       if(mpi_master)call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
    endif
  end subroutine build_hk_model_kgrid



  subroutine build_hk_model_Nkvec(Hk,hk_model,Norb,Nkvec,wdos)
    integer,dimension(:),intent(in)                :: Nkvec
    integer                                        :: Norb
    real(8),dimension(product(Nkvec),size(Nkvec))  :: kgrid ![Nk][Ndim]
    integer                                        :: Nktot
    integer                                        :: ik
    complex(8),dimension(Norb,Norb,product(Nkvec)) :: hk,haux
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    logical,optional                           :: wdos
    logical                                    :: wdos_
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
    wdos_=.false.;if(present(wdos))wdos_=wdos
    !
    call TB_build_kgrid(Nkvec,kgrid,.true.)
    !
    Nktot  = product(Nkvec)
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = hk_model(kgrid(ik,:),Norb)
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
    if(wdos_)then
       allocate(dos_Greal(1,1,1,Norb,Norb,dos_Lreal))
       allocate(dos_wtk(Nktot))
       dos_wtk=1d0/Nktot
#ifdef _MPI
       if(check_MPI())then
          call dmft_gloc_realaxis(MPI_COMM_WORLD,Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
       else
          call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
       endif
#else
       call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(1,1,1,Norb,Norb,dos_Lreal))
#endif
       if(mpi_master)call dmft_print_gf_realaxis(dos_Greal(1,:,:,:,:,:),trim(dos_file),iprint=1)
    endif
  end subroutine build_hk_model_Nkvec









  subroutine build_hkr_model_kgrid(hk,hkr_model,Nlat,Norb,kgrid,pbc,wdos)
    integer                                                 :: Nlat,Norb
    logical                                                 :: pbc
    real(8),dimension(:,:)                                  :: kgrid ![Nktot][Ndim]
    integer                                                 :: Nktot
    integer                                                 :: ik
    complex(8),dimension(Nlat*Norb,Nlat*Norb,size(kgrid,1)) :: hk,haux
    interface 
       function hkr_model(kpoint,Nlat,Norb,pbc)
         real(8),dimension(:)                               :: kpoint
         integer                                            :: Nlat,Norb
         logical                                            :: pbc
         complex(8),dimension(Nlat*Norb,Nlat*Norb)          :: hkr_model
       end function hkr_model
    end interface
    logical,optional                           :: wdos
    logical                                    :: wdos_
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
    wdos_=.false.;if(present(wdos))wdos_=wdos
    !
    Nktot  = size(kgrid,1)
    !
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
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
    if(wdos_)then
       allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
       allocate(dos_wtk(Nktot))
       dos_wtk=1d0/Nktot
#ifdef _MPI
       if(check_MPI())then
          call dmft_gloc_realaxis(MPI_COMM_WORLD,Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
       else
          call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
       endif
#else
       call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
#endif
       if(mpi_master)call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
    endif
  end subroutine build_hkr_model_kgrid




  subroutine build_hkr_model_nkvec(hk,hkr_model,Nlat,Norb,Nkvec,pbc,wdos)
    integer                                                  :: Nlat,Norb
    logical                                                  :: pbc
    integer,dimension(:),intent(in)                          :: Nkvec
    real(8),dimension(product(Nkvec),size(Nkvec))            :: kgrid ![Nk][Ndim]
    integer                                                  :: Nktot
    integer                                                  :: ik
    complex(8),dimension(Nlat*Norb,Nlat*Norb,product(Nkvec)) :: hk,haux
    interface 
       function hkr_model(kpoint,Nlat,Norb,pbc)
         real(8),dimension(:)                                :: kpoint
         integer                                             :: Nlat,Norb
         logical                                             :: pbc
         complex(8),dimension(Nlat*Norb,Nlat*Norb)           :: hkr_model
       end function hkr_model
    end interface
    logical,optional                                         :: wdos
    logical                                                  :: wdos_
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
    wdos_=.false.;if(present(wdos))wdos_=wdos
    !
    call TB_build_kgrid(Nkvec,kgrid,.true.)
    !
    Nktot  = product(Nkvec)
    !
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
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
    if(wdos_)then
       allocate(dos_Greal(Nlat,1,1,Norb,Norb,dos_Lreal))
       allocate(dos_wtk(Nktot))
       dos_wtk=1d0/Nktot
#ifdef _MPI
       if(check_MPI())then
          call dmft_gloc_realaxis(MPI_COMM_WORLD,Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
       else
          call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
       endif
#else
       call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(Nlat,1,1,Norb,Norb,dos_Lreal))
#endif
       if(mpi_master)call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:,:),trim(dos_file),iprint=1)
    endif
  end subroutine build_hkr_model_nkvec




  subroutine build_hk_path(hk,hk_model,Norb,kpath,Nkpath)
    integer                                                   :: Norb
    real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
    integer                                                   :: Nkpath
    integer                                                   :: Npts,Nktot
    integer                                                   :: ik,i
    real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid
    complex(8),dimension(Norb,Norb,(size(kpath,1)-1)*Nkpath)  :: hk,haux
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
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
    Npts  =  size(kpath,1)          !# of k-points along the path
    Nktot = (Npts-1)*Nkpath
    !
    call TB_build_kgrid(kpath,Nkpath,kgrid)
    !
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = hk_model(kgrid(ik,:) , Norb)
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
  end subroutine build_hk_path














  subroutine build_hkR_path(hk,hkr_model,Nlat,Norb,kpath,Nkpath,pbc)
    integer                                                             :: Nlat,Norb
    logical                                                             :: pbc
    real(8),dimension(:,:)                                              :: kpath ![Npts][Ndim]
    integer                                                             :: Nkpath
    integer                                                             :: Npts,Nktot
    integer                                                             :: ik,i
    real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2))           :: kgrid
    complex(8),dimension(Nlat*Norb,Nlat*Norb,(size(kpath,1)-1)*Nkpath)  :: hk,haux
    interface 
       function hkr_model(kpoint,Nlat,Norb,pbc)
         real(8),dimension(:)                      :: kpoint
         integer                                   :: Nlat,Norb
         logical                                   :: pbc
         complex(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
       end function hkr_model
    end interface
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
    Npts  =  size(kpath,1)          !# of k-points along the path
    Nktot = (Npts-1)*Nkpath
    !
    call TB_build_kgrid(kpath,Nkpath,kgrid)
    !
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = hkr_model(kgrid(ik,:),Nlat,Norb,pbc)
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
  end subroutine build_hkR_path





  subroutine build_hk_w90(Hk,Nlso,Nkvec,Kpts_grid,wdos)
    integer                                                :: Nlso
    integer,dimension(:),intent(in)                        :: Nkvec
    complex(8),dimension(Nlso,Nlso,product(Nkvec))         :: Hk,haux
    !
    real(8),dimension(product(Nkvec),size(Nkvec)),optional :: Kpts_grid ![Nk][Ndim]
    logical,optional                                       :: wdos
    logical                                                :: wdos_
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
    wdos_=.false.;if(present(wdos))wdos_=wdos
    !
    Nk =product(Nkvec)
    !
    if(.not.TB_w90%status)stop "build_hk_w90: TB_w90 structure not allocated. Call setup_w90 first."
    !
    call TB_build_kgrid(Nkvec,Kgrid,.true.) !check bk_1,2,3 vectors have been set
    if(present(Kpts_grid))Kpts_grid=Kgrid
    !
    do ik=1+mpi_rank,Nk,mpi_size
       haux(:,:,ik) = w90_hk_model(Kgrid(ik,:),Nlso)
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
    if(wdos_)then
       allocate(dos_Greal(TB_w90%Nlat,TB_w90%Nspin,TB_w90%Nspin,TB_w90%Norb,TB_w90%Norb,dos_Lreal))
       allocate(dos_wtk(Nk))
       dos_wtk=1d0/Nk
       call dmft_gloc_realaxis(Hk,dos_wtk,dos_Greal,zeros(TB_w90%Nlat,TB_w90%Nspin,TB_w90%Nspin,TB_w90%Norb,TB_w90%Norb,dos_Lreal))
       call dmft_print_gf_realaxis(dos_Greal,trim(dos_file),iprint=1)
    endif
  end subroutine build_hk_w90




  subroutine build_hk_w90_path(Hk,Nlso,Kpath,Nkpath)
    integer                                                   :: Nlso
    real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
    integer                                                   :: Nkpath
    integer                                                   :: Npts,Nktot
    complex(8),dimension(Nlso,Nlso,(size(kpath,1)-1)*Nkpath)  :: Hk,haux
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
    if(.not.TB_w90%status)stop "build_hk_w90_path: TB_w90 structure not allocated. Call setup_w90 first."
    !
    Npts  =  size(kpath,1)          !# of k-points along the path
    Nktot = (Npts-1)*Nkpath
    !
    call TB_build_kgrid(kpath,Nkpath,kgrid)
    !
    do ik=1+mpi_rank,Nktot,mpi_size
       haux(:,:,ik) = w90_hk_model(Kgrid(ik,:),Nlso)
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






  subroutine build_Hij_Nrvec(Hij,ts_model,Nso,Nrvec,Links,pbc,wdos)
    integer                                                     :: Nso
    integer,dimension(:),intent(in)                             :: Nrvec
    integer,dimension(:,:),intent(in)                           :: Links ![Nlink][dim]
    logical,optional                                            :: pbc
    logical                                                     :: pbc_
    integer,dimension(product(Nrvec),size(Nrvec))               :: RNgrid
    integer                                                     :: Nlat,Nlink
    integer                                                     :: ilat,jlat
    integer,dimension(size(Nrvec))                              :: Ri,Rj
    integer                                                     :: i,ilink
    complex(8),dimension(Nso,Nso,product(Nrvec),product(Nrvec)) :: Hij
    integer                                                     :: ir,ix,iy,iz,Nr(3)
    !
    interface
       function ts_model(link,Nso)
         integer                       :: link
         integer                       :: Nso
         complex(8),dimension(Nso,Nso) :: ts_model
       end function ts_model
    end interface
    logical,optional                           :: wdos
    logical                                    :: wdos_
    wdos_=.false.;if(present(wdos))wdos_=wdos
    !
    pbc_ = .true. ; if(present(pbc))pbc_=pbc
    !
    if(size(Nrvec)/=size(Links,2))stop "TB_build_Hij ERROR: size(Nvec) != size(Links,2)"
    !
    Nlat  = product(Nrvec)
    Nlink = size(Links,1)
    !
    call TB_build_CoordGrid(Nrvec,RNgrid)
    !
    Hij = zero
    !
    do_lattice: do ilat = 1,Nlat
       Hij(:,:,ilat,ilat) = Hij(:,:,ilat,ilat) + ts_model(0,Nso)
       Ri = RNgrid(ilat,:)
       do_links: do ilink=1,Nlink
          Rj = Ri + Links(ilink,:)
          if(pbc_)then
             do i=1,size(Nrvec)
                if( Rj(i)==0 )Rj(i)=Nrvec(i)
                if( Rj(i)==Nrvec(i)+1)Rj(i)=1
             enddo
          endif
          jlat = TB_find_IndxCoord(Rj,Nrvec)
          if(jlat==0)cycle do_links
          !
          !>Build Hij
          Hij(:,:,ilat,jlat) = Hij(:,:,ilat,jlat) + ts_model(ilink,Nso)
       enddo do_links
    enddo do_lattice
    return
  end subroutine build_Hij_Nrvec



  subroutine FermiLevel_w90(Nkvec,filling,Ef)
    integer,dimension(:),intent(in)          :: Nkvec
    real(8)                                  :: filling,Efermi
    real(8),optional                         :: Ef
    complex(8),dimension(:,:,:),allocatable  :: Hk
    integer                                  :: Nlso,Nk
    if(TB_w90%Ifermi)return
    Nlso = TB_w90%Nspin*TB_w90%Num_Wann
    Nk   = product(Nkvec)
    allocate(Hk(Nlso,Nlso,Nk))
    call TB_build_model(Hk,Nlso,Nkvec)
    call TB_get_FermiLevel(Hk,filling,Efermi)
    if(TB_w90%verbose)write(*,*)"w90 Fermi Level: ",Efermi
    TB_w90%Efermi = Efermi
    if(present(Ef))Ef=Efermi
    TB_w90%Ifermi=.true.
  end subroutine FermiLevel_w90





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




END MODULE TB_BUILD











