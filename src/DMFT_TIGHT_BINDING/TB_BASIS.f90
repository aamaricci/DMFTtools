module TB_BASIS
  USE TB_COMMON
  implicit none



contains


  subroutine TB_reset_ei
    ei_x=[1d0,0d0,0d0]
    ei_y=[0d0,1d0,0d0]
    ei_z=[0d0,0d0,1d0]
    set_eivec=.false.
  end subroutine TB_reset_ei


  subroutine TB_reset_bk
    bk_x=[1d0,0d0,0d0]*pi2
    bk_y=[0d0,1d0,0d0]*pi2
    bk_z=[0d0,0d0,1d0]*pi2
    set_bkvec=.false.
  end subroutine TB_reset_bk


  subroutine TB_set_ei(eix,eiy,eiz)
    real(8),dimension(:),intent(in)                  :: eix
    real(8),dimension(size(eix)),intent(in),optional :: eiy
    real(8),dimension(size(eix)),intent(in),optional :: eiz
    ei_x = 0d0
    ei_x(1:size(eix)) = eix
    !
    if(present(eiy))then
       ei_y = 0d0
       ei_y(1:size(eiy)) = eiy
    endif
    !
    if(present(eiz))then
       ei_z = 0d0
       ei_z(1:size(eiz)) = eiz
    endif
    !
    set_eivec=.true.
  end subroutine TB_set_ei


  subroutine TB_set_bk(bkx,bky,bkz)
    real(8),dimension(:),intent(in)                  :: bkx
    real(8),dimension(size(bkx)),intent(in),optional :: bky
    real(8),dimension(size(bkx)),intent(in),optional :: bkz
    bk_x = 0d0
    bk_x(1:size(bkx)) = bkx
    !
    if(present(bky))then
       bk_y = 0d0
       bk_y(1:size(bky)) = bky
    endif
    !
    if(present(bkz))then
       bk_z = 0d0
       bk_z(1:size(bkz)) = bkz
    endif
    !
    set_bkvec=.true.
  end subroutine TB_set_bk





  subroutine TB_get_bk(bkx,bky,bkz)
    real(8),dimension(:),intent(inout)                  :: bkx
    real(8),dimension(size(bkx)),intent(inout),optional :: bky
    real(8),dimension(size(bkx)),intent(inout),optional :: bkz
    real(8),dimension(3)                                :: b1,b2,b3
    !
    call TB_reciprocal_basis(a1=ei_x,a2=ei_y,a3=ei_z,b1=b1,b2=b2,b3=b3)
    !
    bkx = b1(1:size(bkx))
    if(present(bky))bky = b2(1:size(bky))
    if(present(bkz))bkz = b3(1:size(bkz))
    !
  end subroutine TB_get_bk

  subroutine TB_get_ei(eix,eiy,eiz)
    real(8),dimension(:),intent(inout)                  :: eix
    real(8),dimension(size(eix)),intent(inout),optional :: eiy
    real(8),dimension(size(eix)),intent(inout),optional :: eiz
    real(8),dimension(3)                                :: a1,a2,a3
    !
    call TB_reciprocal_basis(a1=bk_x,a2=bk_y,a3=bk_z,b1=a1,b2=a2,b3=a3)
    !
    eix = a1(1:size(eix))
    if(present(eiy))eiy = a2(1:size(eiy))
    if(present(eiz))eiz = a3(1:size(eiz))
    !
  end subroutine TB_get_ei


  subroutine TB_ei_length(len)
    real(8),dimension(:) :: len
    integer              :: i,n
    n=size(len)
    if(n>0)len(1) = sqrt(dot_product(ei_x,ei_x))
    if(n>1)len(2) = sqrt(dot_product(ei_y,ei_y))
    if(n>2)len(3) = sqrt(dot_product(ei_z,ei_z))
  end subroutine TB_ei_length


  subroutine TB_bk_length(len)
    real(8),dimension(:) :: len
    integer              :: i,n
    n=size(len)
    if(n>0)len(1) = sqrt(dot_product(bk_x,bk_x))
    if(n>1)len(2) = sqrt(dot_product(bk_y,bk_y))
    if(n>2)len(3) = sqrt(dot_product(bk_z,bk_z))
  end subroutine TB_bk_length



  subroutine TB_build_bk(verbose)
    logical,optional :: verbose
    logical          :: verbose_
    if(.not.set_eivec)stop "TB_build_bk ERROR: Direct basis not set, set_eivec=F"
    verbose_=.false.;if(present(verbose))verbose_=verbose
    call TB_reciprocal_basis(&
         a1=ei_x,a2=ei_y,a3=ei_z,&
         b1=bk_x,b2=bk_y,b3=bk_z)
    set_bkvec=.true.
    if(verbose_)call print_bk
  end subroutine TB_build_bk

  subroutine TB_build_ei(verbose)
    logical,optional :: verbose
    logical          :: verbose_
    if(.not.set_bkvec)stop "TB_build_ei ERROR: Reciprocal basis not set, set_bkvec=F"
    verbose_=.false.;if(present(verbose))verbose_=verbose
    !note that we exchange the direct basis with the reciprocal
    call TB_reciprocal_basis(&
         b1=ei_x,b2=ei_y,b3=ei_z,&
         a1=bk_x,a2=bk_y,a3=bk_z)
    set_eivec=.true.
    if(verbose_)call print_ei
  end subroutine TB_build_ei





  subroutine print_ei(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    if(io_eivec)return
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif    
    if(mpi_master)then
       unit=6
       if(present(pfile))open(free_unit(unit),file=reg(pfile))
       write(unit,"(A)")"Using Direct Lattice vectors:"
       write(unit,"(A,3F8.4,A1)")"ei_x = [",(ei_x(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"ei_y = [",(ei_y(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"ei_z = [",(ei_z(i),i=1,3),"]"
       if(present(pfile))close(unit)
    endif
    io_eivec=.true.
  end subroutine print_ei

  subroutine print_bk(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    if(io_bkvec)return
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    if(mpi_master)then
       unit=6
       if(present(pfile))open(free_unit(unit),file=reg(pfile))
       write(unit,"(A)")"Using Reciprocal Lattice vectors:"
       write(unit,"(A,3F8.4,A1)")"bk_x = [",(bk_x(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"bk_y = [",(bk_y(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"bk_z = [",(bk_z(i),i=1,3),"]"
       if(present(pfile))close(unit)
    endif
    io_bkvec=.true.
  end subroutine print_bk





  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Build the reciprocal space {b1[,b2,b3]} given the direct basis {a1[,a2,a3]}
  !  in dimensions up to 3
  !-------------------------------------------------------------------------------------------
  subroutine TB_reciprocal_basis(a1,a2,a3, b1,b2,b3)
    !This routine generates the reciprocal lattice vectors b1,b2,b3
    !given the space vectors a1,a2,a3.
    !the vectors are in units of the lattice constant (a=1)
    real(8),dimension(:),intent(in)                 :: a1
    real(8),dimension(size(a1)),intent(in),optional :: a2
    real(8),dimension(size(a1)),intent(in),optional :: a3
    real(8),dimension(size(a1))                     :: b1
    real(8),dimension(size(a1)),optional            :: b2
    real(8),dimension(size(a1)),optional            :: b3
    !
    real(8),dimension(3)                            :: ar1,ar2,ar3
    real(8),dimension(3)                            :: bk1,bk2,bk3
    real(8)                                         :: den, s
    integer                                         :: iperm, i, j, k, l, ipol
    integer                                         :: N
    real(8),dimension(3,3)                          :: Mat
    !
    N = size(a1)
    !
    ar1=[1d0,0d0,0d0]
    ar2=[0d0,1d0,0d0]
    ar3=[0d0,0d0,1d0]
    !
    ar1(:N)=a1
    if(present(a2))ar2(:N)=a2
    if(present(a3))ar3(:N)=a3
    !
    den = det(transpose(reshape([ar1,ar2,ar3],shape=[3,3])))
    !
    !    here we compute the reciprocal vectors
    i = 1
    j = 2
    k = 3
    do ipol = 1, 3
       bk1(ipol) = (ar2(j)*ar3(k) - ar2(k)*ar3 (j) )/den*pi2
       bk2(ipol) = (ar3(j)*ar1(k) - ar3(k)*ar1 (j) )/den*pi2
       bk3(ipol) = (ar1(j)*ar2(k) - ar1(k)*ar2 (j) )/den*pi2
       l = i
       i = j
       j = k
       k = l
    enddo
    b1=bk1(:N)
    if(present(b2))b2=bk2(:N)
    if(present(b3))b3=bk3(:N)
    return
  end subroutine TB_reciprocal_basis













  subroutine build_kgrid(Nkvec,kgrid,check_bk,origin)
    integer,dimension(:)                    :: Nkvec
    real(8),dimension(:,:)                  :: kgrid ![Nk][Ndim]
    logical,intent(in),optional             :: check_bk
    real(8),dimension(size(Nkvec)),optional :: origin
    !
    logical                                 :: check_bk_    
    real(8),dimension(size(Nkvec))          :: kvec
    real(8),dimension(:),allocatable        :: grid_x,grid_y,grid_z
    integer                                 :: ik,Ivec(3),Nk(3),ndim,Nktot,i
    real(8)                                 :: Lb(3)
    real(8),dimension(3,3)                  :: bk_grid
    real(8),dimension(3)                    :: ktmp
    !
    Nktot = product(Nkvec)
    Ndim  = size(Nkvec)          !dimension of the grid to be built
    call assert_shape(kgrid,[Nktot,Ndim],"build_kgrid","kgrid")
    !
    if(present(origin))BZ_origin(:Ndim)=origin
    !
    Nk=1
    do ik=1,Ndim
       Nk(ik)=Nkvec(ik)
    enddo
    if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
    !
    call print_bk()
    if(check_bk_.AND..not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
    !
    call TB_get_bk(bk_grid(1,:),bk_grid(2,:),bk_grid(3,:))
    !    
    allocate(grid_x(Nk(1)))
    allocate(grid_y(Nk(2)))   
    allocate(grid_z(Nk(3)))
    !
    grid_x = linspace(0d0,1d0,Nk(1),iend=.false.) + BZ_origin(1)
    grid_y = linspace(0d0,1d0,Nk(2),iend=.false.) + BZ_origin(2)
    grid_z = linspace(0d0,1d0,Nk(3),iend=.false.) + BZ_origin(3)
    !
    do ik=1,Nktot
       ivec = i2indices(ik,Nk)       
       ktmp = [grid_x(ivec(1)), grid_y(ivec(2)), grid_z(ivec(3))]
       forall(i=1:Ndim)kvec(i) = dot_product(ktmp,bk_grid(:,i))
       kgrid(ik,:)=kvec
    end do
  end subroutine build_kgrid


  subroutine build_kgrid_generic(Nkvec,kgrid_x,kgrid_y,kgrid_z,check_bk)
    integer,dimension(:)                          :: Nkvec   ! Nk=product(Nkvec);Ndim=size(Nkvec)
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_x ![Nk][Ndim]
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_y ![Nk][Ndim]
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_z ![Nk][Ndim]
    logical,intent(in),optional                   :: check_bk
    logical                                       :: check_bk_
    real(8),dimension(size(Nkvec))                :: kvec    ![Ndim]
    integer                                       :: ik,ix,iy,iz,Nk(3),ndim
    real(8)                                       :: kx,ky,kz
    !
    check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
    !
    ndim = size(Nkvec)          !dimension of the grid to be built
    !
    Nk=1
    do ik=1,size(Nkvec)
       Nk(ik)=Nkvec(ik)
    enddo
    if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
    !
    call print_bk()
    if(check_bk_.AND..not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
    !
    ik=0
    do iz=1,Nk(3)
       kz = dble(iz-1)/Nk(3)
       do iy=1,Nk(2)
          ky = dble(iy-1)/Nk(2)
          do ix=1,Nk(1)
             kx = dble(ix-1)/Nk(1)
             ik=ik+1
             kgrid_x(ik,:)=kx*bk_x(:ndim)
             kgrid_y(ik,:)=ky*bk_y(:ndim)
             kgrid_z(ik,:)=kz*bk_z(:ndim)
          enddo
       enddo
    enddo
  end subroutine build_kgrid_generic


  subroutine build_rgrid(Nrvec,Rgrid,check_ei)
    integer,dimension(:)                          :: Nrvec ! Nr=product(Nrvec);Ndim=size(Nrvec)
    real(8),dimension(product(Nrvec),size(Nrvec)) :: Rgrid ![Nk][Ndim]
    logical,intent(in),optional                   :: check_ei
    logical                                       :: check_ei_
    real(8),dimension(size(Nrvec))                :: Rvec  ![Ndim]
    integer                                       :: ir,ix,iy,iz,Nr(3)
    !
    check_ei_=.false.;if(present(check_ei))check_ei_=check_ei
    !
    Nr=1
    do ir=1,size(Nrvec)
       Nr(ir)=Nrvec(ir)
    enddo
    if(product(Nr)/=product(Nrvec))stop "TB_build_Rgrid ERROR: product(Nrvec) != product(Nr)"
    !
    call print_ei()
    if(check_ei_.AND..not.set_eivec)stop "TB_build_Rgrid ERROR: Ei vectors not set"
    !
    ir=0
    do iz=1,Nr(3)
       do iy=1,Nr(2)
          do ix=1,Nr(1)
             Rvec = ix*ei_x + iy*ei_y + iz*ei_z
             ir=ir+1
             Rgrid(ir,:)=Rvec
          enddo
       enddo
    enddo
  end subroutine build_rgrid






  subroutine TB_refine_kgrid(Nkvec,kgrid_refine,kcenters,DeltaK)
    integer,dimension(:)                    :: Nkvec
    real(8),dimension(:,:),allocatable      :: kgrid_refine
    real(8),dimension(:,:),intent(in)       :: kcenters
    real(8),dimension(size(Nkvec))          :: DeltaK
    !
    real(8),dimension(size(Nkvec))          :: kvec
    real(8),dimension(:),allocatable        :: grid_x,grid_y,grid_z
    integer                                 :: Ndim,Nktot,Ncntr
    integer                                 :: i,ik,Ivec(3),Nk(3),icntr,jcntr,ic
    real(8),dimension(3,3)                  :: bk_grid
    real(8),dimension(size(kcenters,1),3)   :: kstart
    real(8),dimension(3)                    :: Lb,Dk,Kpt
    real(8),dimension(3)                    :: Kshift_
    logical                                 :: boolBZ,boolOlap
    logical,dimension(size(kcenters,1)-1)   :: cond_Lvec
    !
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    Ncntr = size(kcenters,1)
    Nktot = product(Nkvec)
    Ndim  = size(Nkvec)          !dimension of the grid to be built
    call assert_shape(kcenters,[Ncntr,Ndim])
    !
    !
    Nk=1
    do ik=1,Ndim
       Nk(ik)=Nkvec(ik)
    enddo
    if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
    !
    allocate(grid_x(Nk(1)))
    allocate(grid_y(Nk(2)))   
    allocate(grid_z(Nk(3)))
    !
    !< get a local copy of the basis vectors:
    call TB_get_bk(bk_grid(1,:),bk_grid(2,:),bk_grid(3,:))
    !
    do i=1,3
       Lb(i) = sqrt(dot_product(bk_grid(i,:),bk_grid(i,:)))
    enddo
    !
    Dk        = 0d0
    Dk(:Ndim) = DeltaK(:Ndim)
    !
    kstart = 0d0
    do icntr=1,Ncntr
       kstart(icntr,:Ndim) = kcenters(icntr,:Ndim)/Lb(:Ndim) - Dk(:Ndim)/2
    enddo
    !
    !
    if(mpi_master)call start_timer()
    do icntr=1,Ncntr
       if(mpi_master)call eta(icntr,Ncntr)
       !
       grid_x = linspace(0d0,Dk(1),Nk(1),iend=.false.) + kstart(icntr,1)
       grid_y = linspace(0d0,Dk(2),Nk(2),iend=.false.) + kstart(icntr,2)
       grid_z = linspace(0d0,Dk(3),Nk(3),iend=.false.) + kstart(icntr,3)
       !
       do ik=1,Nktot
          ivec = i2indices(ik,Nk)
          Kpt  = [grid_x(ivec(1)), grid_y(ivec(2)), grid_z(ivec(3))]
          !
          !if the point is not in the BZ cycle
          boolBZ = in_rectangle(BZ_origin(:Ndim),dble(ones(Ndim)),Kpt(:Ndim),.true.)
          if(.not.boolBZ)cycle
          !
          !if the point is in any other patch: cycle
          if(icntr>1)then
             forall(jcntr=icntr-1:1:-1)&
                  cond_Lvec(jcntr) = in_rectangle(kstart(jcntr,:Ndim),Dk(:Ndim),Kpt(:Ndim),.false.)
             if(any(cond_Lvec(icntr-1:1:-1)))cycle
          endif
          !
          forall(i=1:Ndim)kvec(i) = dot_product(Kpt,bk_grid(:,i))
          call add_to(kgrid_refine,kvec)
       end do
    enddo
    if(mpi_master)call stop_timer("TB_refine_kgrid")
    !
  contains
    !
    pure function in_rectangle(o,L,p,equal) result(bool)
      real(8),dimension(:),intent(in)       :: o
      real(8),dimension(size(o)),intent(in) :: L
      real(8),dimension(size(o)),intent(in) :: p
      logical,intent(in)                    :: equal
      logical                               :: bool
      integer                               :: idim
      if(equal)then
         bool = p(1)>=o(1) .AND. p(1)<=o(1)+L(1)
         do idim=2,size(o)
            bool = bool .AND. p(idim)>=o(idim) .AND. p(idim)<=o(idim)+L(idim)
         enddo
      else
         bool = p(1)>o(1) .AND. p(1)<o(1)+L(1)
         do idim=2,size(o)
            bool = bool .AND. p(idim)>o(idim) .AND. p(idim)<o(idim)+L(idim)
         enddo
      endif
    end function in_rectangle
    !
    pure function in_circle(c,r,x) result(bool)
      real(8),dimension(:),intent(in)       :: c
      real(8),intent(in)                    :: r
      real(8),dimension(size(c)),intent(in) :: x
      logical                               :: bool
      bool = sum((x(:)-c(:))**2) < r**2
    end function in_circle
    !
  end subroutine TB_refine_kgrid




  
  subroutine TB_write_grid(Grid,file_grid)
    real(8),dimension(:,:) :: Grid ![Nk/Nr,Ndim]
    character(len=*)       :: file_grid
    integer                :: i,j,unit
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    if(mpi_master)then
       open(free_unit(unit),file=reg(file_grid),status="unknown",action="write",position="rewind")
       do i=1,size(Grid,1)
          write(unit,'(3F15.7)') (Grid(i,j),j=1,size(Grid,2))
       enddo
       close(unit)
    endif
  end subroutine TB_write_grid




  subroutine TB_build_CoordGrid(Nrvec,RNgrid)
    integer,dimension(:)                          :: Nrvec  ! Nr=product(Nrvec);Ndim=size(Nrvec)
    integer,dimension(product(Nrvec),size(Nrvec)) :: RNgrid ![Nr][Ndim]
    integer                                       :: ir,ix,iy,iz,Nr(3),Rvec(3)
    !
    Nr=1
    do ir=1,size(Nrvec)
       Nr(ir)=Nrvec(ir)
    enddo
    !
    ir=0
    do iz=1,Nr(3)
       do iy=1,Nr(2)
          do ix=1,Nr(1)
             ir=ir+1
             Rvec = [ix,iy,iz]
             RNgrid(ir,:) = Rvec(1:size(Nrvec))
          enddo
       enddo
    enddo
    ir=0
  end subroutine TB_build_CoordGrid


  function TB_find_IndxCoord(RNvec,Nrvec) result(ir)
    integer,dimension(:)           :: Nrvec
    integer,dimension(size(Nrvec)) :: RNvec
    integer                        :: Nr(3),Rvec(3)
    integer                        :: ir,ix,iy,iz
    logical :: bool
    !
    Nr=1
    Rvec=1
    do ir=1,size(Nrvec)
       Nr(ir)=Nrvec(ir)
       Rvec(ir)=RNvec(ir)
    enddo
    !
    ir=0
    do iz=1,Nr(3)
       do iy=1,Nr(2)
          do ix=1,Nr(1)
             ir=ir+1
             bool = Rvec(1)==ix.AND.Rvec(2)==iy.AND.Rvec(3)==iz
             if(bool)return
          enddo
       enddo
    enddo
    ir=0
  end function TB_find_IndxCoord




  subroutine kgrid_from_path_grid(kpath,Nkpath,kgrid)
    real(8),dimension(:,:)                                    :: kpath ![Npts][Ndim]
    integer                                                   :: Nkpath
    real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2)) :: kgrid ![(Npts-1)*Nkpath][Ndim]
    real(8),dimension(size(kpath,2))                          :: kstart,kstop,kpoint,kdiff
    integer                                                   :: ipts,ik,ic,dim,Npts
    Npts = size(kpath,1)
    ic=0
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nkpath)
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          kgrid(ic,:)=kpoint
       enddo
    enddo
  end subroutine kgrid_from_path_grid

  subroutine kgrid_from_path_dim(kpath,Nkpath,dim,kxpath)
    real(8),dimension(:,:)                      :: kpath ![Npts][Ndim]
    integer                                     :: Nkpath
    integer                                     :: dim
    real(8),dimension((size(kpath,1)-1)*Nkpath) :: kxpath ![(Npts-1)*Nkpath][Ndim]
    real(8),dimension(size(kpath,2))            :: kstart,kstop,kpoint,kdiff
    integer                                     :: ipts,ik,ic,Npts
    Npts = size(kpath,1)
    if(dim>size(kpath,2))stop "TB_build_grid ERROR: dim > Ndim"
    ic=0
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nkpath)
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          kxpath(ic)=kpoint(dim)
       enddo
    enddo
  end subroutine kgrid_from_path_dim




END MODULE TB_BASIS












! function refine_kgrid(Nkvec,kgrid,kcenter) result(kgrid_refine)
!   integer,dimension(:)                               :: Nkvec
!   real(8),dimension(:,:),intent(in)                  :: kgrid    ![Nk][Ndim]
!   real(8),dimension(size(Nkvec)),intent(in),optional :: kcenter  ![Ndim]
!   !
!   real(8),dimension(size(kgrid,1),size(kgrid,2))     :: kgrid_refine ![Nk][Ndim]
!   !
!   real(8),dimension(size(Nkvec))                     :: kvec  ![Ndim]
!   real(8),dimension(:),allocatable                   :: grid_x,grid_y,grid_z
!   integer                                            :: i,ik,ix,iy,iz,Nk(3),ndim,Nktot
!   real(8)                                            :: kx,ky,kz,Lb(3),Dk(3)
!   real(8),dimension(3)                               :: kstart_
!   !
!   Nktot = product(Nkvec)
!   Ndim  = size(Nkvec)          !dimension of the grid to be built
!   !
!   call assert_shape(kgrid,[Nktot,Ndim],"build_kgrid","kgrid")
!   !
!   Nk=1
!   do ik=1,Ndim
!      Nk(ik)=Nkvec(ik)
!   enddo
!   if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
!   !
!   Lb(1) = sqrt(dot_product(bk_x,bk_x))
!   Lb(2) = sqrt(dot_product(bk_y,bk_y))
!   Lb(3) = sqrt(dot_product(bk_z,bk_z))
!   Dk    = 0d0
!   ik = 2
!   do i=1,ndim
!      ik = ik + (i-1)*product(Nk(1:i-1))
!      Dk(i) = kgrid(ik,i)-kgrid(1,i)
!   enddo
!   !
!   kstart_ = 0d0 ; if(present(kcenter))kstart_(:Ndim) = kcenter - 2*Dk(:Ndim)
!   ! do i=1,ndim       
!   !    if(kstart_(i)< minval(kgrid(:,i)) )kstart_(i)=minval(kgrid(:,i))
!   !    if(kstart_(i)> maxval(kgrid(:,i)) )kstart_(i)=maxval(kgrid(:,i))
!   ! enddo
!   kstart_ = kstart_/Lb
!   Dk      = Dk/Lb
!   !    
!   allocate(grid_x(Nk(1)))
!   allocate(grid_y(Nk(2)))   
!   allocate(grid_z(Nk(3)))
!   !
!   grid_x = linspace(0d0,1d0,Nk(1),iend=.false.)*4*Dk(1) + kstart_(1)
!   grid_y = linspace(0d0,1d0,Nk(2),iend=.false.)*4*Dk(2) + kstart_(2)
!   grid_z = linspace(0d0,1d0,Nk(3),iend=.false.)*4*Dk(3) + kstart_(3)
!   !
!   do iz=1,Nk(3)
!      do iy=1,Nk(2)
!         do ix=1,Nk(1)
!            kx   = grid_x(ix)
!            ky   = grid_y(iy)
!            kz   = grid_z(iz)
!            ik   = indices2i([ix,iy,iz],Nk)
!            kvec = kx*bk_x(:ndim) + ky*bk_y(:ndim) + kz*bk_z(:ndim)
!            kgrid_refine(ik,:)=kvec
!         end do
!      end do
!   end do
! end function refine_kgrid

! subroutine build_kgrid(Nkvec,kgrid,check_bk)
!   integer,dimension(:)                          :: Nkvec ! Nk=product(Nkvec);Ndim=size(Nkvec)
!   real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
!   logical,intent(in),optional                   :: check_bk
!   logical                                       :: check_bk_
!   real(8),dimension(size(Nkvec))                :: kvec  ![Ndim]
!   integer                                       :: ik,ix,iy,iz,Nk(3),ndim
!   real(8)                                       :: kx,ky,kz
!   !
!   check_bk_=.false.;if(present(check_bk))check_bk_=check_bk
!   !
!   ndim = size(Nkvec)          !dimension of the grid to be built
!   !
!   Nk=1
!   do ik=1,size(Nkvec)
!      Nk(ik)=Nkvec(ik)
!   enddo
!   if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
!   !
!   call print_bk()
!   if(check_bk_.AND..not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
!   !
!   ik=0
!   do iz=1,Nk(3)
!      kz = dble(iz-1)/Nk(3)
!      do iy=1,Nk(2)
!         ky = dble(iy-1)/Nk(2)
!         do ix=1,Nk(1)
!            kx = dble(ix-1)/Nk(1)
!            kvec = kx*bk_x(:ndim) + ky*bk_y(:ndim) + kz*bk_z(:ndim)
!            ik=ik+1
!            kgrid(ik,:)=kvec
!         enddo
!      enddo
!   enddo
! end subroutine build_kgrid
