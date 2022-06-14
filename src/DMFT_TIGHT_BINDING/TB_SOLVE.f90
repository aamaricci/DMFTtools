module TB_SOLVE
  USE TB_COMMON
  USE TB_BASIS
  USE TB_IO
  USE TB_WANNIER90
  USE TB_BUILD
  implicit none



contains



  !<  solve the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) along a given linear
  ! path in the Brillouin Zone. A GNUPLOT script to plot the bands together with their
  ! character is generated.
  subroutine solve_Hk_along_BZpath(hk_model,Nlso,kpath,Nkpath,colors_name,points_name,file,iproject)
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    integer                                   :: Nlso
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath
    type(rgb_color),dimension(Nlso)           :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: iproject
    character(len=256)                        :: file_
    logical                                   :: iproject_
    character(len=256)                        :: xtics
    integer                                   :: Npts,Ndim,Nktot
    integer                                   :: ipts,ik,ic,unit,iorb
    real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
    real(8)                                   :: eval(Nlso),coeff(Nlso),klen,ktics(size(Kpath,1))
    complex(8)                                :: h(Nlso,Nlso)
    type(rgb_color)                           :: corb(Nlso),c(Nlso)
    character(len=10)                         :: chpoint
    character(len=32)                         :: fmt
    real(8),allocatable                       :: kseg(:),Ekval(:,:)
    integer,allocatable                       :: Ekcol(:,:)
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    file_    = "Eigenbands.tb";if(present(file))file_=file
    iproject_= TB_w90%status
    if(TB_w90%status)write(*,*)"Using iproject=.TRUE. in W90 interface. Disable it explicitly using iproject=.false. "
    if(present(iproject))iproject_=iproject
    !
    Npts = size(kpath,1)
    Ndim = size(kpath,2)
    Nktot= (Npts-1)*Nkpath
    do iorb=1,Nlso
       corb(iorb) = colors_name(iorb)
    enddo
    !
    if(.not.set_bkvec)stop "solve_w90hk_along_BZpath ERROR: bk vectors not set!"
    !  
    if(iproject_)then
       select case(Ndim)
       case (1)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
       case(2)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
       case (3)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
       end select
    endif
    !
    !
    if(mpi_master)then
       write(*,*)"solving model along the path:"
       do ipts=1,Npts
          write(*,*)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic = 0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nlso))
    allocate(ekcol(Nktot,Nlso))
    klen=0d0  
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nkpath)
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hk_model(kpoint,Nlso)
          call eigh(h,Eval)
          do iorb=1,Nlso
             coeff(:)=h(:,iorb)*conjg(h(:,iorb))
             c(iorb) = coeff.dot.corb
             Ekval(ic,iorb) = Eval(iorb)
             Ekcol(ic,iorb) = rgb(c(iorb))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = kseg(ic-1)
    !
    if(mpi_master)then
       unit=free_unit()
       open(unit,file=reg(file_))
       do iorb=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval(ic,iorb),Ekcol(ic,iorb)
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
       do ipts=2,Npts-1
          xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
       enddo
       xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
       !
       open(unit,file=reg(file_)//".gp")
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf -dEPSCrop - "//reg(file_)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid noytics xtics"
       !
       do iorb=1,Nlso
          chpoint=str(0.95d0-(iorb-1)*0.05d0)
          write(unit,"(A)")str("#set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
               " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
       enddo
       !
       write(unit,*)"plot '"//reg(file_)//"' every :::0 u 1:2:3 w l lw 3 lc rgb variable"
       write(unit,*)"# to print from the i-th to the j-th block use: every :::i::j"
       !
       close(unit)
       call system("chmod +x "//reg(file_)//".gp")
    endif
  end subroutine solve_Hk_along_BZpath

  subroutine solve_w90Hk_along_BZpath(Nlso,kpath,Nkpath,colors_name,points_name,file,iproject)
    integer                                   :: Nlso
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath
    type(rgb_color),dimension(Nlso)           :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: iproject
    character(len=256)                        :: file_
    logical                                   :: iproject_
    character(len=256)                        :: xtics
    !
    if(.not.TB_w90%status)stop "solve_w90hk_along_BZpath: TB_w90 structure not allocated. Call setup_w90 first."
    !
    file_    = "Eigenbands.tb";if(present(file))file_=file
    iproject_= TB_w90%status
    if(TB_w90%status)write(*,*)"Using iproject=.TRUE. in W90 interface. Disable it explicitly using iproject=.false. "
    if(present(iproject))iproject_=iproject
    !    character(len=256)                        :: xtics

    call solve_Hk_along_BZpath(w90_hk_model,Nlso,kpath,Nkpath,colors_name,points_name,file_,iproject_)
    !
  end subroutine solve_w90Hk_along_BZpath



  subroutine solve_HkR_along_BZpath(hkr_model,Nlat,Nso,kpath,Nkpath,colors_name,points_name,file,pbc,iproject)
    interface 
       function hkr_model(kpoint,Nlat,Nso,pbc)
         real(8),dimension(:)                    :: kpoint
         integer                                 :: Nlat,Nso
         logical                                 :: pbc
         complex(8),dimension(Nlat*Nso,Nlat*Nso) :: hkr_model
       end function hkr_model
    end interface
    integer                                   :: Nlat,Nso,Nlso
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath,Nktot
    type(rgb_color),dimension(Nlat,Nso)       :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: pbc,iproject
    character(len=256)                        :: file_
    logical                                   :: pbc_,iproject_
    character(len=256)                        :: xtics
    integer                                   :: Npts,Ndim
    integer                                   :: ipts,ik,ic,unit,iorb,ilat,io,nrot,u1,u2
    real(8)                                   :: coeff(Nlat*Nso),klen,ktics(size(Kpath,1))
    type(rgb_color)                           :: corb(Nlat*Nso),c(Nlat*Nso)
    real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
    complex(8),dimension(Nlat*Nso,Nlat*Nso)   :: h
    real(8),dimension(Nlat*Nso)               :: Eval
    real(8),allocatable                       :: kseg(:),Ekval(:,:)
    integer,allocatable                       :: Ekcol(:,:)
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    file_    = "Eigenbands.tb";if(present(file))file_=file
    iproject_= .false.        ;if(present(iproject))iproject_=iproject
    pbc_     = .true.         ;if(present(pbc))pbc_=pbc
    !
    Nlso  = Nlat*Nso
    Npts  = size(kpath,1)
    Ndim  = size(kpath,2)
    Nktot = (Npts-1)*Nkpath
    !
    do ilat=1,Nlat
       do io=1,Nso
          corb(io + (ilat-1)*Nso) = colors_name(ilat,io)
       enddo
    enddo
    !
    if(.not.set_bkvec)stop "solve_w90hk_along_BZpath ERROR: bk vectors not set!"
    !
    if(iproject_)then
       select case(Ndim)
       case (1)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
       case(2)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
       case (3)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
       end select
    endif
    !
    if(mpi_master)then
       write(*,*)"Solving model along the path:"
       do ipts=1,Npts
          write(*,"(A,10(A,1x),A1)")"Point"//str(ipts)//": [",(str(kpath(ipts,ic)),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic=0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nlso))
    allocate(ekcol(Nktot,Nlso))
    klen = 0d0
    if(mpi_master)call start_timer()
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/Nkpath
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hkr_model(kpoint,Nlat,Nso,pbc)
          call eigh(h,Eval)
          if(mpi_master)call eta(ic,Nktot)
          do io=1,Nlso
             coeff(:)=h(:,io)*conjg(h(:,io))
             c(io) = coeff.dot.corb
             Ekval(ic,io) = Eval(io)
             Ekcol(ic,io) = rgb(c(io))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = Kseg(ic-1)
    if(mpi_master)call stop_timer()
    !
    if(mpi_master)then
       open(free_unit(unit),file=str(file_))
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval(ic,io),Ekcol(ic,io)
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       xtics=""
       xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
       do ipts=2,Npts-1
          xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
       enddo
       xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
       !
       open(unit,file=reg(file_)//".gp")
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid ytics xtics"
       !
       write(unit,*)"plot '"//reg(file_)//"' every :::0 u 1:2:3 w l lw 3 lc rgb variable"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_)//".gp")
    endif
  end subroutine solve_HkR_along_BZpath



  include "w90hr/read_Hr_w90_solve_Hk_along_BZpath.f90"



END MODULE TB_SOLVE











