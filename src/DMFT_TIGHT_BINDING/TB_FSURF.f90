module TB_FSURF
  USE TB_COMMON
  USE TB_BASIS
  USE TB_IO
  USE TB_WANNIER90
  USE TB_BUILD
  implicit none
  private


  interface TB_fsurface
     module procedure :: TB_fsurf_nkvec
     module procedure :: TB_fsurf_w90_nkvec
  end interface TB_fsurface
  public :: TB_fsurface


contains




  subroutine TB_fsurf_nkvec(hk_model,Nlso,Ef,Nkvec,colors_name,file,cutoff,max_order,deltak,BZorigin,iwrite)
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
    logical,optional                              :: iwrite
    integer,optional                              :: max_order
    real(8),optional                              :: deltak
    real(8),dimension(size(Nkvec)),optional       :: BZorigin
    !
    real(8)                                       :: cutoff_,deltak_
    logical                                       :: iwrite_
    integer                                       :: max_order_
    !
    character(len=256)                            :: file_
    type(rgb_color),dimension(Nlso)               :: colors_name_
    !
    real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid
    real(8),dimension(:,:),allocatable            :: rfd_kgrid
    real(8),dimension(size(Nkvec))                :: kpoint,kcenter,kstart
    real(8),dimension(product(Nkvec),size(Nkvec)) :: tmp_kgrid
    real(8),dimension(:,:),allocatable            :: kpts
    complex(8),dimension(Nlso,Nlso)               :: evec
    real(8),dimension(Nlso)                       :: eval
    real(8),dimension(Nlso)                       :: coeff
    type(rgb_color),dimension(Nlso)               :: corb,c
    real(8),dimension(size(Nkvec))                :: bk_len
    real(8)                                       :: min_range(2)
    real(8)                                       :: max_range(2)
    integer                                       :: ik,io,jo,Nktot,Ndim,unit,ipt,ic,Ncntr,iorder
    real(8),dimension(2)                          :: Gp=[0d0,0d0]
    real(8),dimension(2)                          :: Xp=[0.5d0,0d0]
    real(8),dimension(2)                          :: Yp=[0d0,0.5d0]
    real(8),dimension(2)                          :: Mp=[0.5d0,0.5d0]
    !
    file_       = "FSurf"  ;if(present(file))file_=file
    cutoff_     = 1d-1     ;if(present(cutoff))cutoff_=cutoff
    colors_name_=black     ;if(present(colors_name))colors_name_=colors_name
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
    if(.not.set_bkvec)stop "TB_fsurf_nkvec ERROR: bk vectors not set!"
    !
    !< Count points in the actual grid:
    Nktot  = product(Nkvec)
    Ndim   = size(Nkvec)
    if(present(BZorigin))BZ_origin(:Ndim)=BZorigin
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
    write(unit,*)"set tics font ',20'"
    write(unit,*)"set xtics pi/4"
    write(unit,*)"set ytics pi/4"
    write(unit,*)"set format x ''"
    write(unit,*)"set format y ''"
    !
    call TB_bk_length(bk_len)
    Gp = Gp*bk_len
    Xp = Xp*bk_len
    Yp = Yp*bk_len
    Mp = Mp*bk_len
    write(unit,*)"set label '{/Symbol G}' at "&
         //str(Gp(1)-0.03*bk_len(1))//","//str(Gp(2)-0.03*bk_len(2))//&
         " right font ',24'"
    !
    write(unit,*)"set label 'X' at "&
         //str(Xp(1)-0.03*bk_len(1))//","//str(Xp(2)-0.03*bk_len(2))//&
         " center font ',24'"
    !
    write(unit,*)"set label 'Y' at "&
         //str(Yp(1)-0.03*bk_len(1))//","//str(Yp(2)-0.03*bk_len(1))//&
         " center font ',24'"
    !
    write(unit,*)"set label 'M' at "&
         //str(Mp(1)-0.03*bk_len(1))//","//str(Mp(2)-0.03*bk_len(2))//&
         " center font ',24'"
    !
    write(unit,*)"set label '' at "//str(Gp(1))//","//str(Gp(2))//" point pt 4"
    write(unit,*)"set label '' at "//str(Xp(1))//","//str(Xp(2))//" point pt 4"
    write(unit,*)"set label '' at "//str(Yp(1))//","//str(Yp(2))//" point pt 4"
    write(unit,*)"set label '' at "//str(Mp(1))//","//str(Mp(2))//" point pt 4"
    !
    min_range = (0d0+BZ_origin(1))*bk_len(:)
    max_range = (1d0+BZ_origin(1))*bk_len(:)
    write(unit,*)"set xrange ["//str(min_range(1))//":"//str(max_range(1))//"]"
    write(unit,*)"set yrange ["//str(min_range(2))//":"//str(max_range(2))//"]"
    !
    write(unit,*)"plot '"//reg(file_)//"' u 2:3:(1):1 w p pt 7 ps 0.4 lc rgb variable notit"
    !
    close(unit)
    call system("chmod +x "//reg(file_)//".gp")
    !
    !
    call TB_build_kgrid(Nkvec,kgrid,.true.,BZ_origin)
    if(iwrite_)call TB_write_grid(kgrid,"rfd_kgrid_1")
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
       if(iwrite_)call TB_write_grid(rfd_kgrid,"rfd_kgrid_"//str(iorder))
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
       deallocate(rfd_kgrid)
       !
       cutoff_=cutoff_/3d0
       deltak_=deltak_/2d0
    enddo
    !
  end subroutine TB_fsurf_nkvec






  subroutine TB_fsurf_w90_nkvec(Nlso,Ef,Nkvec,colors_name,file,cutoff,max_order,deltak)
    integer                                  :: Nlso
    real(8)                                  :: Ef
    integer,dimension(:),intent(in)          :: Nkvec
    !
    type(rgb_color),dimension(Nlso),optional :: colors_name
    real(8),optional                         :: cutoff
    character(len=*),optional                :: file
    integer,optional                         :: max_order
    real(8),optional                         :: deltak
    !
    real(8)                                  :: cutoff_,deltak_
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
    max_order_  = 3           ;if(present(max_order))max_order_=max_order
    deltak_     = 0.13d0      ;if(present(deltak))deltak_=deltak
    !
    call TB_fsurf_nkvec(w90_hk_model,Nlso,Ef,Nkvec,&
         colors_name_,&
         file_,&
         cutoff_,&
         max_order_,&
         deltak_,&
         BZorigin=TB_w90%BZorigin,&
         iwrite=TB_w90%verbose)
  end subroutine TB_fsurf_w90_nkvec



END MODULE TB_FSURF











