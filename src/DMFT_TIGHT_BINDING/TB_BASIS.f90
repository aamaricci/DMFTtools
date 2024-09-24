!There are different ways to set the direct-reciprocal lattice basis ei-bk:
! 1) user calls TB_set_ei/bk (user knows the vectors or uses default ones)
!
! 2) user calls either TB_set_ei OR TB_set_bk. IF set_bkvec/set_eivec=F  the TB_
! procedures will build the requied basis thru TB_build_bk/ei with a warning to the user.
! IF no basis is set the grid creation raises error  
!
! 3) user calls either TB_set_ei OR TB_set_bk. Then user calls TB_build_bk/TB_build_ei
!    to implicitly construct the missing basis without retrieving it. This will be used by the TB_
!
! 4) user calls either TB_set_ei OR TB_set_bk. Then user constructs and retrieve reciprocal basis
! explicitly then pushes them using again TB_set_bk OR TB_set_ei.


module TB_BASIS
  USE TB_COMMON
  implicit none


contains


  !This check that at least one basis is set
  !F = F.or.F, T otherwise
  function TB_basis_check() result(bool)
    logical :: bool
    bool = set_eivec.OR.set_bkvec    
  end function TB_basis_check


  subroutine TB_reset_ei
    ei_1=[1d0,0d0,0d0]
    ei_2=[0d0,1d0,0d0]
    ei_3=[0d0,0d0,1d0]
    set_eivec=.false.
  end subroutine TB_reset_ei


  subroutine TB_reset_bk
    bk_1=[1d0,0d0,0d0]*pi2
    bk_2=[0d0,1d0,0d0]*pi2
    bk_3=[0d0,0d0,1d0]*pi2
    set_bkvec=.false.
  end subroutine TB_reset_bk


  !Called with no argument will just use default basis vectors
  subroutine TB_set_ei(eix,eiy,eiz)
    real(8),dimension(:),intent(in),optional :: eix
    real(8),dimension(:),intent(in),optional :: eiy
    real(8),dimension(:),intent(in),optional :: eiz
    logical :: bool
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !Reset to default
    call TB_reset_ei()
    bool=.true.
    !
    if(present(eix))then
       ei_1 = 0d0
       ei_1(1:size(eix)) = eix
       bool=.false.
    endif
    !
    if(present(eiy))then
       ei_2 = 0d0
       ei_2(1:size(eiy)) = eiy
       bool=.false.
    endif
    !
    if(present(eiz))then
       ei_3 = 0d0
       ei_3(1:size(eiz)) = eiz
       bool=.false.
    endif
    !
    if(bool.and.mpi_master)write(*,"(A)")"TB_set_ei: using default values"
    set_eivec=.true.
    !if(mpi_master)write(*,*) set_eivec
  end subroutine TB_set_ei

  subroutine TB_set_bk(bkx,bky,bkz)
    real(8),dimension(:),intent(in),optional :: bkx
    real(8),dimension(:),intent(in),optional :: bky
    real(8),dimension(:),intent(in),optional :: bkz
    logical :: bool
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif    
    !
    !Reset to default
    call TB_reset_bk()
    bool=.true.
    !
    if(present(bkx))then
       bk_1 = 0d0
       bk_1(1:size(bkx)) = bkx
       bool=.false.
    endif
    !
    if(present(bky))then
       bk_2 = 0d0
       bk_2(1:size(bky)) = bky
       bool=.false.
    endif
    !
    if(present(bkz))then
       bk_3 = 0d0
       bk_3(1:size(bkz)) = bkz
       bool=.false.
    endif
    !
    if(bool.AND.mpi_master)write(*,"(A)")"TB_set_bk: using default values"
    set_bkvec=.true.
  end subroutine TB_set_bk




  subroutine TB_build_bk(verbose)
    logical,optional :: verbose
    logical          :: verbose_
    if(.not.set_eivec)stop "TB_build_bk ERROR: Direct basis not set, set_eivec=F"
    verbose_=.false.;if(present(verbose))verbose_=verbose
    call TB_reciprocal_basis(&
         a1=ei_1,a2=ei_2,a3=ei_3,&
         b1=bk_1,b2=bk_2,b3=bk_3)
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
         b1=ei_1,b2=ei_2,b3=ei_3,&
         a1=bk_1,a2=bk_2,a3=bk_3)
    set_eivec=.true.
    if(verbose_)call print_ei
  end subroutine TB_build_ei




  subroutine TB_get_bk(b1,b2,b3)
    real(8),dimension(:),intent(inout)                 :: b1
    real(8),dimension(size(b1)),intent(inout),optional :: b2
    real(8),dimension(size(b1)),intent(inout),optional :: b3
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(.not.tb_basis_check())stop "TB_get_bk ERROR: neiter ei nor bk basis are set."
    !
    if(.not.set_bkvec)then !set_eivec==T 
       if(mpi_master)write(*,"(A)")"Building bk from ei:"
       call TB_build_bk(.true.)
    endif
    !
    b1 = bk_1(1:size(b1))
    if(present(b2))b2 = bk_2(1:size(b2))
    if(present(b3))b3 = bk_3(1:size(b3))
    !
  end subroutine TB_get_bk

  subroutine TB_get_ei(e1,e2,e3)
    real(8),dimension(:),intent(inout)                 :: e1
    real(8),dimension(size(e1)),intent(inout),optional :: e2
    real(8),dimension(size(e1)),intent(inout),optional :: e3
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(.not.tb_basis_check())stop "TB_get_ei ERROR: neiter ei nor bk basis are set."
    !
    if(.not.set_eivec)then  !set_bkvec==T 
       if(mpi_master)write(*,"(A)")"Building ei from bk:"
       call TB_build_ei(.true.)
    endif
    !
    e1 = ei_1(1:size(e1))
    if(present(e2))e2 = ei_2(1:size(e2))
    if(present(e3))e3 = ei_3(1:size(e3))
    !
  end subroutine TB_get_ei


  subroutine TB_ei_length(len)
    real(8),dimension(:) :: len
    integer              :: i,n
    if(.not.set_eivec)stop "TB_ei_length error: ei basis not set"
    n=size(len)
    if(n>0)len(1) = sqrt(dot_product(ei_1,ei_1))
    if(n>1)len(2) = sqrt(dot_product(ei_2,ei_2))
    if(n>2)len(3) = sqrt(dot_product(ei_3,ei_3))
  end subroutine TB_ei_length


  subroutine TB_bk_length(len)
    real(8),dimension(:) :: len
    integer              :: i,n
    if(.not.set_bkvec)stop "TB_bk_length error: bk basis not set"
    n=size(len)
    if(n>0)len(1) = sqrt(dot_product(bk_1,bk_1))
    if(n>1)len(2) = sqrt(dot_product(bk_2,bk_2))
    if(n>2)len(3) = sqrt(dot_product(bk_3,bk_3))
  end subroutine TB_bk_length








  subroutine print_ei(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    ! if(io_eivec)return
    if(.not.set_eivec)stop "TB_print_ei error: ei basis not set"
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif    
    if(mpi_master)then
       unit=6
       if(present(pfile))open(free_unit(unit),file=reg(pfile))
       write(unit,"(A)")"Using Direct Lattice vectors:"
       write(unit,"(A,3F8.4,A1)")"ei_1 = [",(ei_1(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"ei_2 = [",(ei_2(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"ei_3 = [",(ei_3(i),i=1,3),"]"
       if(present(pfile))close(unit)
    endif
    ! io_eivec=.true.
  end subroutine print_ei

  subroutine print_bk(pfile)
    character(len=*),optional :: pfile
    integer                   :: unit,i
    ! if(io_bkvec)return
    if(.not.set_bkvec)stop "TB_bk_length error: bk basis not set"
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    if(mpi_master)then
       unit=6
       if(present(pfile))open(free_unit(unit),file=reg(pfile))
       write(unit,"(A)")"Using Reciprocal Lattice vectors:"
       write(unit,"(A,3F8.4,A1)")"bk_1 = [",(bk_1(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"bk_2 = [",(bk_2(i),i=1,3),"]"
       write(unit,"(A,3F8.4,A1)")"bk_3 = [",(bk_3(i),i=1,3),"]"
       if(present(pfile))close(unit)
    endif
    ! io_bkvec=.true.
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













  subroutine build_kgrid(Nkvec,kgrid,origin,project,width)
    integer,dimension(:)                    :: Nkvec
    real(8),dimension(:,:)                  :: kgrid ![Nk][Ndim]
    logical,optional :: project
    real(8),dimension(size(Nkvec)),optional :: origin
    real(8),dimension(size(Nkvec)),optional :: width
    real(8),dimension(size(Nkvec))          :: kvec
    real(8),dimension(:),allocatable        :: grid_x,grid_y,grid_z
    integer                                 :: ik,Ivec(3),Nk(3),ndim,Nktot,i
    real(8),dimension(3)                    :: ktmp,wtmp
    logical :: project_
    !
    if(.not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
    !
    Nktot = product(Nkvec)
    Ndim  = size(Nkvec)          !dimension of the grid to be built
    wtmp      = 1d0
    BZ_origin = 0d0
    call assert_shape(kgrid,[Nktot,Ndim],"build_kgrid","kgrid")
    !
    project_=.true.;if(present(project))project_=project
    if(present(origin))BZ_origin(:Ndim)=origin
    if(present(width))wtmp(1:ndim)=width
    Nk=1
    do ik=1,Ndim
       Nk(ik)=Nkvec(ik)
    enddo
    if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
    !
    call print_bk()
    !
    allocate(grid_x(Nk(1)))
    allocate(grid_y(Nk(2)))   
    allocate(grid_z(Nk(3)))
    !
    grid_x = linspace(0d0,wtmp(1),Nk(1),iend=.false.) + BZ_origin(1)
    grid_y = linspace(0d0,wtmp(2),Nk(2),iend=.false.) + BZ_origin(2)
    grid_z = linspace(0d0,wtmp(3),Nk(3),iend=.false.) + BZ_origin(3)
    !
    do ik=1,Nktot
       ivec = i2indices(ik,Nk)
       ktmp = [grid_x(ivec(1)), grid_y(ivec(2)), grid_z(ivec(3))]
       kgrid(ik,:) = ktmp(1:ndim)
       if(project_)&
            kgrid(ik,:) = ktmp(1)*bk_1(:ndim) + ktmp(2)*bk_2(:ndim) + ktmp(3)*bk_3(:ndim)
    end do
  end subroutine build_kgrid


  subroutine refine_kgrid(kgridIN,Nkvec,kcenters,lambda,KgridOut,WkOut)
    real(8),dimension(:,:)                :: kgridIN
    integer,dimension(:)                  :: Nkvec
    real(8),dimension(:,:)                :: kcenters
    real(8),dimension(size(Nkvec))        :: lambda
    real(8),dimension(:,:),allocatable    :: KgridOut
    real(8),dimension(:),allocatable      :: WkOut
    !
    real(8),dimension(size(Nkvec))        :: kvec
    real(8),dimension(:),allocatable      :: grid_x,grid_y,grid_z
    integer                               :: ik,Ivec(3),Nk(3),ndim,Nktot,i,Ncntr,NkIN,icntr,jcntr,NkClean,NkRefined,NkOut
    real(8),dimension(3)                  :: ktmp,wtmp
    real(8),dimension(3,3)                :: bk_grid
    real(8),dimension(size(kcenters,1),3) :: kstart
    real(8),dimension(3)                  :: Lb,Dk,Kpt
    logical                               :: boolBZ
    logical,dimension(size(kcenters,1))   :: cond_Lvec
    !
    if(.not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    !
    NkIN  = size(kgridIN,1)
    Ndim  = size(kgridIN,2)
    Ncntr = size(kcenters,1)
    if(Ndim /= size(Nkvec))stop "refine_kgrid error: Nkvec has bad dimension"
    call assert_shape(kcenters,[Ncntr,Ndim],"refine_kgrid","kcenters")
    !
    Nk=1
    do ik=1,Ndim
       Nk(ik)=Nkvec(ik)
    enddo
    if(product(Nk)/=product(Nkvec))stop "refine_grid error: product(Nkvec) != product(Nk)"
    Nktot = product(Nk)
    !
    !< get a local copy of the basis vectors:
    call TB_get_bk(bk_grid(1,:),bk_grid(2,:),bk_grid(3,:))
    do i=1,3
       Lb(i) = sqrt(dot_product(bk_grid(i,:),bk_grid(i,:)))
    enddo
    !
    Dk     = 0d0
    kstart = 0d0
    !
    Dk(:Ndim) = lambda
    do icntr=1,Ncntr
       kstart(icntr,:Ndim) = kcenters(icntr,:Ndim)/Lb(:Ndim) - Dk(:Ndim)/2
    enddo
    !
    !Check if, for every center, Ktmp is in the patch. If not add it to the actual grid
    do ik=1,size(kgridIN,1)
       ktmp=0d0;ktmp(:Ndim) = KgridIN(ik,:)/Lb(:Ndim)
       do icntr=1,Ncntr
          cond_Lvec(icntr) = in_rectangle(kstart(icntr,:Ndim),Dk(:Ndim),Ktmp(:Ndim),.true.)
       enddo
       if(any(cond_Lvec))cycle
       call add_to(kgridOut,kgridIN(ik,:))
    enddo
    NkClean = size(KgridOut,1)
    !
    !Now refine the grid:
    do icntr=1,Ncntr
       !
       grid_x = linspace(0d0,Dk(1),Nk(1),iend=.false.) + kstart(icntr,1)
       grid_y = linspace(0d0,Dk(2),Nk(2),iend=.false.) + kstart(icntr,2)
       grid_z = linspace(0d0,Dk(3),Nk(3),iend=.false.) + kstart(icntr,3)
       !
       do ik=1,Nktot
          ivec = i2indices(ik,Nk)
          Kpt  = [grid_x(ivec(1)), grid_y(ivec(2)), grid_z(ivec(3))]
          where(Kpt(:Ndim)<0d0)Kpt(:Ndim)=Kpt(:Ndim)+1d0
          !
          !if the point is not in the BZ cycle
          boolBZ = in_rectangle(BZ_origin(:Ndim),dble(ones(Ndim)),Kpt(:Ndim),.true.)
          !
          !if the point is in any other patch: cycle
          if(icntr>1)then
             do jcntr=icntr-1,1,-1
                cond_Lvec(jcntr) = in_rectangle(kstart(jcntr,:Ndim),Dk(:Ndim),Kpt(:Ndim),.false.)
             enddo
             if(any(cond_Lvec(icntr-1:1:-1)))cycle
          endif
          !
          forall(i=1:Ndim)kvec(i) = dot_product(Kpt,bk_grid(:,i))
          call add_to(KgridOut,kvec)
       end do
    enddo
    NkOut     = size(KgridOut,1)
    NkRefined = NkOut-NkClean
    allocate(WkOut(Nkout))
    WkOut(1:NkClean) = (1d0-product(lambda))/NkClean
    WkOut(NkClean+1:)= product(lambda)/NkRefined
    !
  contains
    !
    function in_rectangle(o,L,p,equal) result(bool)
      real(8),dimension(:),intent(in)       :: o
      real(8),dimension(size(o)),intent(in) :: L
      real(8),dimension(size(o)),intent(in) :: p
      real(8),dimension(size(o))            :: a,b
      logical,intent(in)                    :: equal
      logical                               :: bool
      integer                               :: idim,icase
      if(any(L>1d0))stop "in_rectangle error: L>1 c`mon! "
      !Ensure a < b
      a = o
      b = o+L
      !check point p is in [0:1)
      if(any(p<0d0))stop "in_rectangle error: any(p<0d0)"
      if(any(p>1d0))stop "in_rectangle error: any(p>1d0)"
      !check different case with respect to [0:1) periodicity:
      !Case 1: 0<a_i<1, 0<b_i<1 DEFAULT
      icase = 1
      !Case 2: a_i < 0 + b=a+L < 1
      if(any(a<0d0).AND.any(b<1d0))icase = 2
      !Case 3: b_i > 1 + a=b-L < 1 
      if(any(b>1d0).AND.any(a<1d0))icase = 3
      !Case 4: a < b < 0
      if(any(a<0d0).AND.any(b<0d0))icase = 4
      !Case 5: b > a > 1
      if(any(a>1d0).AND.any(b>1d0))icase = 5
      !
      !p in [0:1)
      if(equal)then
         select case(icase)
         case default           ! 0 < a_i <= p_i <= b_i < 1
            bool = p(1)>=a(1) .AND. p(1)<=b(1)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>=a(idim) .AND. p(idim)<=b(idim))
            enddo
         case(2)                ! 0 < a_i+1 <= p_i || p_i <= b_i < 1
            bool = p(1)>=(a(1)+1d0) .OR. p(1)<=b(1)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>=(a(idim)+1d0) .OR. p(idim)<=b(idim))
            enddo
         case(3)                ! 0 < a_i <= p_i || p_i < b_i-1 < 1
            bool = p(1)>=a(1) .OR. p(1)<=(b(1)-1d0)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>=a(idim) .OR. p(idim)<=(b(idim)-1d0))
            enddo
         case(4)                ! 0 < a_i+1 < p_i && p_i < b_i+1 < 1
            bool = p(1)>=(a(1)+1d0) .AND. p(1)<=(b(1)+1d0)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>=(a(idim)+1d0) .AND. p(idim)<=(b(idim)+1d0))
            enddo
         case(5)                ! 0 < a_i-1 < p_i && p_i < b_i-1 < 1
            bool = p(1)>=(a(1)-1d0) .AND. p(1)<=(b(1)-1d0)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>=(a(idim)-1d0) .AND. p(idim)<=(b(idim)-1d0))
            enddo
         end select
      else
         select case(icase)
         case default           ! 0 < a_i < p_i < b_i < 1
            bool = p(1)>a(1) .AND. p(1)<b(1)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>a(idim) .AND. p(idim)<b(idim))
            enddo
         case(2)                ! 0 < a_i+1 < p_i || p_i < b_i < 1
            bool = p(1)>(a(1)+1d0) .OR. p(1)<b(1)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>(a(idim)+1d0) .OR. p(idim)<b(idim))
            enddo
         case(3)                ! 0 < a_i < p_i || p_i < b_i-1 < 1
            bool = p(1)>a(1) .OR. p(1)<(b(1)-1d0)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>a(idim) .OR. p(idim)<(b(idim)-1d0))
            enddo
         case(4)                ! 0 < a_i+1 < p_i && p_i < b_i+1 < 1
            bool = p(1)>(a(1)+1d0) .AND. p(1)<(b(1)+1d0)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>(a(idim)+1d0) .AND. p(idim)<(b(idim)+1d0))
            enddo
         case(5)                ! 0 < a_i-1 < p_i && p_i < b_i-1 < 1
            bool = p(1)>(a(1)-1d0) .AND. p(1)<(b(1)-1d0)
            do idim=2,size(o)
               bool = bool .AND. (p(idim)>(a(idim)-1d0) .AND. p(idim)<(b(idim)-1d0))
            enddo
         end select
      endif
    end function in_rectangle
    !
  end subroutine refine_kgrid


  subroutine build_kgrid_generic(Nkvec,kgrid_x,kgrid_y,kgrid_z)
    integer,dimension(:)                          :: Nkvec   ! Nk=product(Nkvec);Ndim=size(Nkvec)
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_x ![Nk][Ndim]
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_y ![Nk][Ndim]
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid_z ![Nk][Ndim]
    real(8),dimension(size(Nkvec))                :: kvec    ![Ndim]
    integer                                       :: ik,ix,iy,iz,Nk(3),ndim
    real(8)                                       :: kx,ky,kz
    !
    if(.not.set_bkvec)stop "TB_build_grid ERROR: bk vectors not set"
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
    !
    ik=0
    do iz=1,Nk(3)
       kz = dble(iz-1)/Nk(3)
       do iy=1,Nk(2)
          ky = dble(iy-1)/Nk(2)
          do ix=1,Nk(1)
             kx = dble(ix-1)/Nk(1)
             ik=ik+1
             kgrid_x(ik,:)=kx*bk_1(:ndim)
             kgrid_y(ik,:)=ky*bk_2(:ndim)
             kgrid_z(ik,:)=kz*bk_3(:ndim)
          enddo
       enddo
    enddo
  end subroutine build_kgrid_generic


  subroutine build_rgrid(Nrvec,Rgrid)
    integer,dimension(:)                          :: Nrvec ! Nr=product(Nrvec);Ndim=size(Nrvec)
    real(8),dimension(product(Nrvec),size(Nrvec)) :: Rgrid ![Nk][Ndim]
    real(8),dimension(size(Nrvec))                :: Rvec  ![Ndim]
    integer                                       :: ir,ix,iy,iz,Nr(3),D
    !
    if(.not.set_eivec)stop "TB_build_Rgrid ERROR: Ei vectors not set"
    !
    D = size(Nrvec)
    Nr=1
    do ir=1,D
       Nr(ir)=Nrvec(ir)
    enddo
    if(product(Nr)/=product(Nrvec))stop "TB_build_Rgrid ERROR: product(Nrvec) != product(Nr)"
    !
    call print_ei()
    !
    ir=0
    do iz=1,Nr(3)
       do iy=1,Nr(2)
          do ix=1,Nr(1)
             Rvec = ix*ei_1(1:D) + iy*ei_2(1:D) + iz*ei_3(1:D)
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
    if(.not.set_bkvec)stop "TB_regine_kgrid ERROR: bk vectors not set"
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

  subroutine klen_from_path(kpath,Nkpath,karray)
    real(8),dimension(:,:)                      :: kpath ![Npts][Ndim]
    integer                                     :: Nkpath
    real(8),dimension((size(kpath,1)-1)*Nkpath) :: karray ![(Npts-1)*Nkpath]
    real(8),dimension(size(kpath,2))            :: kstart,kstop,kdiff
    real(8) :: klen
    integer                                     :: ipts,ik,ic,dim,Npts
    Npts = size(kpath,1)
    ic=0
    klen=0d0  
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nkpath)
       do ik=1,Nkpath
          ic=ic+1
          karray(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))          
       enddo
    enddo
  end subroutine klen_from_path



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
!   Lb(1) = sqrt(dot_product(bk_1,bk_1))
!   Lb(2) = sqrt(dot_product(bk_2,bk_2))
!   Lb(3) = sqrt(dot_product(bk_3,bk_3))
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
!            kvec = kx*bk_1(:ndim) + ky*bk_2(:ndim) + kz*bk_3(:ndim)
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
!            kvec = kx*bk_1(:ndim) + ky*bk_2(:ndim) + kz*bk_3(:ndim)
!            ik=ik+1
!            kgrid(ik,:)=kvec
!         enddo
!      enddo
!   enddo
! end subroutine build_kgrid
