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



  subroutine TB_fsurf_nkvec(hk_model,Nlso,Ef,Nkvec,colors_name,file,cutoff,pi_shift)
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
    logical,optional                              :: pi_shift
    !
    real(8)                                       :: cutoff_
    logical                                       :: pi_shift_
    !
    character(len=256)                            :: file_
    type(rgb_color),dimension(Nlso)               :: colors_name_
    !
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
    real(8),dimension(size(Nkvec))                :: kpoint,kcenter,kstart
    real(8),dimension(product(Nkvec),size(Nkvec)) :: tmp_kgrid
    real(8),dimension(:,:),allocatable            :: kpts0
    real(8),dimension(:,:),allocatable            :: kpts1
    real(8),dimension(:,:),allocatable            :: kpts2
    !
    real(8),dimension(:,:,:),allocatable          :: rfd_kgrid1
    real(8),dimension(:,:,:),allocatable          :: rfd_kgrid2
    !
    complex(8),dimension(Nlso,Nlso)               :: evec
    real(8),dimension(Nlso)                       :: eval
    real(8),dimension(Nlso)                       :: coeff
    type(rgb_color),dimension(Nlso)               :: corb,c
    real(8),dimension(size(Nkvec))                :: bk_len
    real(8)                                       :: min_xrange,min_yrange
    real(8)                                       :: max_xrange,max_yrange
    integer                                       :: ik,io,jo,Nktot,Ndim,unit,ipt,ic

    file_       = "FSurf"  ;if(present(file))file_=file
    cutoff_     = 1d-1     ;if(present(cutoff))cutoff_=cutoff
    colors_name_=black     ;if(present(colors_name))colors_name_=colors_name
    pi_shift_   = .false.  ;if(present(pi_shift))pi_shift_=pi_shift

    !< Get colors
    do io=1,Nlso
       corb(io) = colors_name_(io)
    enddo
    !
    !< Build K-grid
    if(.not.set_bkvec)stop "TB_fsurf_nkvec ERROR: bk vectors not set!"
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
    allocate(kpts0(0,Ndim))
    open(free_unit(unit),file=reg(file_))
    call start_timer()
    do ik=1,Nktot
       call eta(ik,Nktot)
       kpoint = kgrid(ik,:)
       Evec = hk_model(kpoint,Nlso)
       call eigh(Evec,Eval)
       if( all(abs(Eval-Ef)>cutoff_) )cycle
       do io=1,Nlso
          coeff(:)=abs(Evec(:,io))**2
          c(io)   = coeff.dot.corb
          if(abs(Eval(io)-Ef) < cutoff_)then
             call add_to(kpts0,kpoint)
             write(unit,*)rgb(c(io)),(kpoint(jo),jo=1,Ndim),Eval(io)
          endif
       enddo
    enddo
    close(unit)
    call stop_timer("TB_FSurface: 0-th order")
    write(*,*)"# points:",size(kpts0,1)

    !Refine 1 
    if(size(kpts0,1)>0)then
       cutoff_=cutoff_/10
       open(free_unit(unit),file=reg(file_))
       rewind(unit)       
       allocate(rfd_kgrid1(0,Nktot,Ndim))
       allocate(kpts1(0,Ndim))
       call start_timer()
       do ipt=1,size(kpts0,1)
          kcenter = kpts0(ipt,:)
          tmp_kgrid = TB_refine_kgrid(Nkvec,Kgrid,Kcenter)
          !
          do ik=1,Nktot
             call eta(ik + (ipt-1)*Nktot,Nktot*size(kpts0,1))
             kpoint = tmp_kgrid(ik,:)
             Evec = hk_model(kpoint,Nlso)
             call eigh(Evec,Eval)
             if( all(abs(Eval-Ef)>cutoff_) )cycle
             do io=1,Nlso
                coeff(:)=abs(Evec(:,io))**2
                c(io)   = coeff.dot.corb
                if(abs(Eval(io)-Ef) < cutoff_)then
                   call add_to(kpts1,kpoint)
                   write(unit,*)rgb(c(io)),(kpoint(jo),jo=1,Ndim),Eval(io)
                endif
             enddo
          enddo
          if(size(kpts1,1)>0)call add_to(rfd_kgrid1,tmp_kgrid)
       enddo
       close(unit)
       call stop_timer("TB_FSurface: 1-st order")
       write(*,*)"# points:",size(kpts1,1)
    endif


    !Refine 2 
    if(size(kpts1,1)>0)then
       cutoff_=cutoff_/10
       open(free_unit(unit),file=reg(file_))
       rewind(unit)
       allocate(rfd_kgrid2(0,Nktot,Ndim))
       allocate(kpts2(0,Ndim))
       call start_timer()
       do ipt=1,size(kpts1,1)
          kcenter = kpts1(ipt,:)
          tmp_kgrid = TB_refine_kgrid(Nkvec,rfd_Kgrid1(ipt,:,:),Kcenter)
          !
          do ik=1,Nktot
             call eta(ik + (ipt-1)*Nktot,Nktot*size(kpts1,1))
             kpoint = tmp_kgrid(ik,:)
             Evec = hk_model(kpoint,Nlso)
             call eigh(Evec,Eval)
             if( all(abs(Eval-Ef)>cutoff_) )cycle
             do io=1,Nlso
                coeff(:)=abs(Evec(:,io))**2
                c(io)   = coeff.dot.corb
                if(abs(Eval(io)-Ef) < cutoff_)then
                   call add_to(kpts2,kpoint)
                   write(unit,*)rgb(c(io)),(kpoint(jo),jo=1,Ndim),Eval(io)
                endif
             enddo
          enddo
          if(size(kpts2,1)>0)call add_to(rfd_kgrid2,tmp_kgrid)
       enddo
       close(unit)
       call stop_timer("TB_FSurface: 2-nd order")
       write(*,*)"# points:",size(kpts2,1)
    endif


    call TB_bk_length(bk_len)
    min_xrange = (0d0+kstart(1))*bk_len(1) !!minval(kgrid(:,1))
    max_xrange = (1d0+kstart(1))*bk_len(1) !!maxval(kgrid(:,1))
    if(Ndim>1)then
       min_yrange = (0d0+kstart(2))*bk_len(2)!!minval(kgrid(:,2))
       max_yrange = (1d0+kstart(2))*bk_len(2)!!maxval(kgrid(:,2))
    endif

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
    write(unit,*)"set xtics ('-{/Symbol p}' -pi, '0' 0, '{/Symbol p}' pi, '2{/Symbol p}' 2*pi)"
    write(unit,*)"set ytics ('-{/Symbol p}' -pi, '0' 0, '{/Symbol p}' pi, '2{/Symbol p}' 2*pi)"
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
    write(unit,*)"plot '"//reg(file_)//"' u 2:3:(1):1 w p pt 7 ps 0.3 lc rgb variable notit"
    !
    close(unit)
    call system("chmod +x "//reg(file_)//".gp")

  end subroutine TB_fsurf_nkvec









  subroutine TB_fsurf_w90_nkvec(Nlso,Ef,Nkvec,colors_name,file,cutoff,pi_shift)
    integer                                  :: Nlso
    real(8)                                  :: Ef
    integer,dimension(:),intent(in)          :: Nkvec
    !
    type(rgb_color),dimension(Nlso),optional :: colors_name
    real(8),optional                         :: cutoff
    character(len=*),optional                :: file
    logical,optional                         :: pi_shift
    !
    real(8)                                  :: cutoff_
    logical                                  :: pi_shift_
    character(len=256)                       :: file_
    type(rgb_color),dimension(Nlso)          :: colors_name_
    !
    if(.not.TB_w90%status)stop "TB_fsurf_w90_nkvec: TB_w90 structure not allocated. Call setup_w90 first."
    !
    file_       = "FSurf"  ;if(present(file))file_=file
    cutoff_     = 1d-1     ;if(present(cutoff))cutoff_=cutoff
    colors_name_=black     ;if(present(colors_name))colors_name_=colors_name
    pi_shift_   = .false.  ;if(present(pi_shift))pi_shift_=pi_shift
    !
    call TB_fsurf_nkvec(w90_hk_model,Nlso,Ef,Nkvec,colors_name_,file_,cutoff_,pi_shift_)
  end subroutine TB_fsurf_w90_nkvec







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











