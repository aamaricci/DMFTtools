module TB_FSURF
  USE TB_COMMON
  USE TB_BASIS
  USE TB_IO
  USE TB_WANNIER90
  USE TB_BUILD
  implicit none
  private


  interface add_to
     module procedure :: add_to_A1
     module procedure :: add_to_A2
     module procedure :: add_to_A3
  end interface add_to


  interface TB_fsurface
     module procedure :: TB_fsurf_nkvec
     module procedure :: TB_fsurf_w90_nkvec
  end interface TB_fsurface
  public :: TB_fsurface


contains




  subroutine TB_fsurf_nkvec(hk_model,Nlso,Ef,Nkvec,colors_name,file,cutoff,pi_shift,max_order,deltak,iwrite)
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    integer                                       :: Nlso
    real(8)                                       :: Ef
    integer,dimension(:),intent(in)               :: Nkvec
    !
    type(rgb_color),dimension(Nlso),optional      :: colors_name
    real(8),optional                              :: cutoff
    character(len=*),optional                     :: file
    logical,optional                              :: pi_shift,iwrite
    integer,optional                              :: max_order
    real(8),optional                              :: deltak
    !
    real(8)                                       :: cutoff_,deltak_
    logical                                       :: pi_shift_,iwrite_
    integer                                       :: max_order_
    !
    character(len=256)                            :: file_
    type(rgb_color),dimension(Nlso)               :: colors_name_
    !
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid
    !
    real(8),dimension(:,:),allocatable            :: rfd_kgrid
    real(8),dimension(size(Nkvec))                :: kpoint,kcenter,kstart
    real(8),dimension(product(Nkvec),size(Nkvec)) :: tmp_kgrid
    real(8),dimension(:,:),allocatable            :: kpts
    !
    !
    complex(8),dimension(Nlso,Nlso)               :: evec
    real(8),dimension(Nlso)                       :: eval
    real(8),dimension(Nlso)                       :: coeff
    type(rgb_color),dimension(Nlso)               :: corb,c
    real(8),dimension(size(Nkvec))                :: bk_len
    real(8)                                       :: min_xrange,min_yrange
    real(8)                                       :: max_xrange,max_yrange
    integer                                       :: ik,io,jo,Nktot,Ndim,unit,ipt,ic,Ncntr,iorder
    !
    file_       = "FSurf"  ;if(present(file))file_=file
    cutoff_     = 1d-1     ;if(present(cutoff))cutoff_=cutoff
    colors_name_=black     ;if(present(colors_name))colors_name_=colors_name
    pi_shift_   = .false.  ;if(present(pi_shift))pi_shift_=pi_shift
    max_order_  = 4        ;if(present(max_order))max_order_=max_order
    deltak_     = 0.13d0    ;if(present(deltak))deltak_=deltak
    iwrite_     = .false.  ;if(present(iwrite))iwrite_=iwrite
    !
    !< Get colors
    do io=1,Nlso
       corb(io) = colors_name_(io)
    enddo
    !
    !< Build K-grid
    ! if(.not.set_bkvec)stop "TB_fsurf_nkvec ERROR: bk vectors not set!"
    if(pi_shift_)then
       kstart = -0.5d0
       call TB_build_kgrid(Nkvec,kgrid,.true.,kstart=kstart)
    else
       kstart = 0d0
       call TB_build_kgrid(Nkvec,kgrid,.true.)
    endif
    !
    !< Count points in the actual grid:
    Nktot  = product(Nkvec)
    Ndim   = size(Nkvec)
    !
    call TB_bk_length(bk_len)
    min_xrange = (0d0+kstart(1))*bk_len(1) !!minval(kgrid(:,1))
    max_xrange = (1d0+kstart(1))*bk_len(1) !!maxval(kgrid(:,1))
    if(Ndim>1)then
       min_yrange = (0d0+kstart(2))*bk_len(2)!!minval(kgrid(:,2))
       max_yrange = (1d0+kstart(2))*bk_len(2)!!maxval(kgrid(:,2))
    endif
    !
    open(unit,file=reg(file_)//".gp")
    write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
    write(unit,*)"#set out '"//reg(file_)//".png'"
    write(unit,*)""
    write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
    write(unit,*)"#set out '"//reg(file_)//".svg'"
    write(unit,*)""
    write(unit,*)"set term postscript eps enhanced color 'Times'"
    write(unit,*)"set output '|ps2pdf -dEPSCrop - "//reg(file_)//".pdf'"
    !
    write(unit,*)"set font 'Times,20'"  
    write(unit,*)"set size square"
    write(unit,*)"set xlabel 'k_x' font 'Times Italic,20'"
    write(unit,*)"set ylabel 'k_y' font 'Times Italic,20'"
    write(unit,*)"unset key"
    !
    write(unit,*)"set tics font ',20'"
    write(unit,*)"set xtics ('-{/Symbol p}' -"//str(max_xrange/2)//&
         ", '0' 0, '{/Symbol p}' "//str(max_xrange/2)//", '2{/Symbol p}' "//str(max_xrange)//")"
    write(unit,*)"set ytics ('-{/Symbol p}' -"//str(max_yrange/2)//&
         ", '0' 0, '{/Symbol p}' "//str(max_yrange/2)//", '2{/Symbol p}' "//str(max_yrange)//")"
    !
    write(unit,*)"#set label '{/Symbol G}' at 0,0 right font ',24'"
    write(unit,*)"#set label '' at 0,0 point pt 4"
    write(unit,*)"#set label 'X' at pi,0.25 center font ',24'"
    write(unit,*)"#set label '' at pi,0 point pt 4"
    write(unit,*)"#set label 'Y' at 0.25,pi center font ',24'"
    write(unit,*)"#set label '' at 0,pi point pt 4"
    write(unit,*)"#set label 'M' at pi+0.2,pi+0.2 center font ',24'"
    write(unit,*)"#set label '' at pi,pi point pt 4"
    !
    write(unit,*)"set xrange ["//str(min_xrange)//":"//str(max_xrange)//"]"
    write(unit,*)"set yrange ["//str(min_yrange)//":"//str(max_yrange)//"]"
    !
    write(unit,*)"plot '"//reg(file_)//"' u 2:3:(1):1 w p pt 7 ps 0.4 lc rgb variable notit"
    !
    close(unit)
    call system("chmod +x "//reg(file_)//".gp")
    !
    !
    !
    allocate(kpts(0,Ndim))
    !
    open(free_unit(unit),file=reg(file_))
    call start_timer()
    do ik=1,Nktot
       kpoint = kgrid(ik,:)
       Evec = hk_model(kpoint,Nlso)
       call eigh(Evec,Eval)
       if( all(abs(Eval-Ef)>cutoff_) )cycle
       do io=1,Nlso
          coeff(:)=abs(Evec(:,io))**2
          c(io)   = coeff.dot.corb
          if(abs(Eval(io)-Ef) < cutoff_)then
             call add_to(kpts,kpoint)
             write(unit,*)rgb(c(io)),(kpoint(jo),jo=1,Ndim),Eval(io)
          endif
       enddo
    enddo
    close(unit)
    Ncntr=size(kpts,1)
    cutoff_=cutoff_/3d0
    if(iwrite_)call write_grid(kgrid,"kgrid_0")
    call stop_timer("TB_FSurface: 1st-order - kpts: "//str(Ncntr))
    !
    !
    do iorder=2,max_order_
       if(Ncntr==0)exit
       !
       open(free_unit(unit),file=reg(file_)//".1")
       rewind(unit)
       !
       allocate(rfd_kgrid(0,Ndim))
       call TB_refine_kgrid(Nkvec,rfd_kgrid,kpts,deltak_)
       !
       call start_timer()
       deallocate(kpts) ; allocate(kpts(0,Ndim))
       do ik=1,size(rfd_kgrid,1)
          kpoint = rfd_kgrid(ik,:)
          Evec = hk_model(kpoint,Nlso)
          call eigh(Evec,Eval)
          if( all(abs(Eval-Ef)>cutoff_) )cycle
          do io=1,Nlso
             coeff(:)=abs(Evec(:,io))**2
             c(io)   = coeff.dot.corb
             if(abs(Eval(io)-Ef) < cutoff_)then
                call add_to(kpts,kpoint)
                write(unit,*)rgb(c(io)),(kpoint(jo),jo=1,Ndim),Eval(io)
             endif
          enddo
       enddo
       Ncntr = size(kpts,1)
       call stop_timer("TB_FSurface "//str(iorder)//"th-order - kpts: "//str(Ncntr))
       close(unit)
       call system("mv -f "//reg(file_)//".1 "//reg(file_))
       if(iwrite_)call write_grid(rfd_kgrid,"rfd_kgrid_"//str(iorder))
       deallocate(rfd_kgrid)
       !
       cutoff_=cutoff_/3d0
       deltak_=deltak_/2d0
    enddo
    !
  end subroutine TB_fsurf_nkvec






  subroutine TB_fsurf_w90_nkvec(Nlso,Ef,Nkvec,colors_name,file,cutoff,pi_shift,max_order,deltak)
    integer                                  :: Nlso
    real(8)                                  :: Ef
    integer,dimension(:),intent(in)          :: Nkvec
    !
    type(rgb_color),dimension(Nlso),optional :: colors_name
    real(8),optional                         :: cutoff
    character(len=*),optional                :: file
    logical,optional                         :: pi_shift
    integer,optional                         :: max_order
    real(8),optional                         :: deltak
    !
    real(8)                                  :: cutoff_,deltak_
    logical                                  :: pi_shift_
    integer                                  :: max_order_
    !
    character(len=256)                       :: file_
    type(rgb_color),dimension(Nlso)          :: colors_name_
    !
    !
    if(.not.TB_w90%status)stop "TB_fsurf_w90_nkvec: TB_w90 structure not allocated. Call setup_w90 first."
    !
    file_       = "w90FSurf"  ;if(present(file))file_=file
    cutoff_     = 1d-1        ;if(present(cutoff))cutoff_=cutoff
    colors_name_=black        ;if(present(colors_name))colors_name_=colors_name
    pi_shift_   = .false.     ;if(present(pi_shift))pi_shift_=pi_shift
    max_order_  = 3           ;if(present(max_order))max_order_=max_order
    deltak_     = 0.13d0      ;if(present(deltak))deltak_=deltak
    !
    call TB_fsurf_nkvec(w90_hk_model,Nlso,Ef,Nkvec,colors_name_,file_,cutoff_,pi_shift_,max_order_,deltak_)
  end subroutine TB_fsurf_w90_nkvec






  subroutine TB_refine_kgrid(Nkvec,kgrid_refine,kcenters,DeltaK)
    integer,dimension(:)                                 :: Nkvec
    real(8),dimension(:,:),allocatable                   :: kgrid_refine
    real(8),dimension(:,:),intent(in)                    :: kcenters
    real(8)                                              :: DeltaK
    !
    real(8),dimension(size(Nkvec))                       :: kvec
    real(8),dimension(:),allocatable                     :: grid_x,grid_y,grid_z
    integer                                              :: Ndim,Nktot,Ncntr
    integer                                              :: i,ik,ix,iy,iz,Nk(3),icntr,jcntr,ic
    real(8)                                              :: kx,ky,kz
    real(8),dimension(3)                                 :: bk_x
    real(8),dimension(3)                                 :: bk_y
    real(8),dimension(3)                                 :: bk_z
    real(8),dimension(size(kcenters,1),3)                :: kstart
    real(8),dimension(3)                                 :: Lb,Dk,Origin,Kpt,Unity
    logical                                              :: boolBZ,boolOlap
    logical,dimension(size(kcenters,1)-1)                :: cond_Lvec
    !
    ! DeltaK = 1d0/sqrt(99d0)        !use an input for this
    !
    !
    Ncntr = size(kcenters,1)
    Nktot = product(Nkvec)
    Ndim  = size(Nkvec)          !dimension of the grid to be built
    call assert_shape(kcenters,[Ncntr,Ndim])
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
    Origin = 0d0
    Unity  = 1d0
    !
    call TB_get_bk(bk_x,bk_y,bk_z)
    !
    Lb(1) = sqrt(dot_product(bk_x,bk_x))
    Lb(2) = sqrt(dot_product(bk_y,bk_y))
    Lb(3) = sqrt(dot_product(bk_z,bk_z))
    !
    Dk        = 0d0
    Dk(:Ndim) = DeltaK  !half way back in each direction
    !
    kstart = 0d0
    do icntr=1,Ncntr
       kstart(icntr,:Ndim) = kcenters(icntr,:Ndim)/Lb(:Ndim) - Dk(:Ndim)/2
    enddo
    !
    !
    call start_timer()
    do icntr=1,Ncntr
       call eta(icntr,Ncntr)
       !
       grid_x = linspace(0d0,Dk(1),Nk(1),iend=.false.) + kstart(icntr,1)
       grid_y = linspace(0d0,Dk(2),Nk(2),iend=.false.) + kstart(icntr,2)
       grid_z = linspace(0d0,Dk(3),Nk(3),iend=.false.) + kstart(icntr,3)
       !
       do iz=1,Nk(3)
          do iy=1,Nk(2)
             do ix=1,Nk(1)
                kx   = grid_x(ix)
                ky   = grid_y(iy)
                kz   = grid_z(iz)
                ik   = indices2i([ix,iy,iz],Nk)
                Kpt  = [kx,ky,kz]
                !if the point is not in the BZ cycle
                boolBZ = in_rectangle(Origin(:Ndim),Unity(:Ndim),Kpt(:Ndim))
                if(.not.boolBZ)cycle
                !
                !if the point is in any other patch: cycle
                if(icntr>1)then
                   forall(jcntr=icntr-1:1:-1)&
                        cond_Lvec(jcntr) = in_rectangle(kstart(jcntr,:Ndim),Dk(:Ndim),Kpt(:Ndim))
                   if(any(cond_Lvec(icntr-1:1:-1)))cycle
                endif
                !
                kvec = kx*bk_x(:ndim) + ky*bk_y(:ndim) + kz*bk_z(:ndim)
                call add_to(kgrid_refine,kvec)
             end do
          end do
       end do
    enddo
    call stop_timer("TB_refine_kgrid")
  end subroutine TB_refine_kgrid

  pure function in_rectangle(origin,side,point) result(bool)
    real(8),dimension(:),intent(in)            :: origin
    real(8),dimension(size(origin)),intent(in) :: side
    real(8),dimension(size(origin)),intent(in) :: point
    logical                                    :: bool
    integer                                    :: idim
    bool = point(1)>=origin(1) .AND. point(1)<=origin(1)+side(1)
    do idim=2,size(origin)
       bool = bool .AND. point(idim)>=origin(idim) .AND. point(idim)<=origin(idim)+side(idim)
    enddo
  end function in_rectangle



  subroutine write_grid(Grid,file_grid)
    real(8),dimension(:,:) :: Grid ![Nk/Nr,Ndim]
    character(len=*)       :: file_grid
    integer                :: i,j,unit
    open(free_unit(unit),file=file_grid)!,status="unknown",action="write",position="rewind")
    do i=1,size(Grid,1)
       write(unit,*) Grid(i,:)
    enddo
    close(unit)
  end subroutine write_grid




  subroutine add_to_A1(vec,val)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),intent(in)                             :: val  
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_A1



  subroutine add_to_A2(vec,val)
    real(8),dimension(:,:),allocatable,intent(inout) :: vec
    real(8),intent(in),dimension(size(vec,2))        :: val  
    real(8),dimension(:,:),allocatable               :: tmp
    integer                                          :: n,ndim
    !
    ndim = size(vec,2)
    !
    if (allocated(vec)) then
       n = size(vec,1)
       allocate(tmp(n+1,ndim))
       tmp(1:n,:) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n,ndim))
    end if
    !
    !Put val as last entry:
    vec(n,:) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_A2


  subroutine add_to_A3(vec,val)
    real(8),dimension(:,:,:),allocatable,intent(inout)    :: vec
    real(8),dimension(size(vec,2),size(vec,3)),intent(in) :: val  
    real(8),dimension(:,:,:),allocatable                  :: tmp
    integer                                               :: n,n2,n3
    !
    n2 = size(vec,2)
    n3 = size(vec,3)
    !
    if (allocated(vec)) then
       n = size(vec,1)
       allocate(tmp(n+1,n2,n3))
       tmp(1:n,:,:) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n,n2,n3))
    end if
    !
    !Put val as last entry:
    vec(n,:,:) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_A3



END MODULE TB_FSURF











