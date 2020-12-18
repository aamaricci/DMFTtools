module TB_BUILD
  USE TB_COMMON
  USE TB_IO
  USE TB_BASIS
  implicit none





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
    wdos_=.false.
    !
    Nktot  = size(kgrid,1)
    !
    Haux   = zero
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
    wdos_=.false.
    !
    call build_kgrid(Nkvec,kgrid)
    !
    Nktot  = product(Nkvec)
    Haux   = zero
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
    wdos_=.false.
    !
    Nktot  = size(kgrid,1)
    Haux   = zero
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
    wdos_=.false.
    !
    call build_kgrid(Nkvec,kgrid)
    !
    Nktot  = product(Nkvec)
    Haux   = zero
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
    call kgrid_from_path_grid(kpath,Nkpath,kgrid)
    !
    Haux  = zero
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
    call kgrid_from_path_grid(kpath,Nkpath,kgrid)
    !
    Haux  = zero
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






END MODULE TB_BUILD











