subroutine solve_Hk_along_BZpath(hk_model,Nlso,kpath,Nk,colors_name,points_name,file)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  integer                                   :: Nlso
  real(8),dimension(:,:)                    :: kpath
  integer                                   :: Nk
  type(rgb_color),dimension(Nlso)           :: colors_name
  character(len=*),dimension(size(kpath,1)) :: points_name
  character(len=*),optional                 :: file
  character(len=256)                        :: file_,xtics
  integer                                   :: Npts,Nrot,Ndim
  integer                                   :: ipts,ik,ic,unit,u1,u2,iorb
  real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
  real(8)                                   :: eval(Nlso),coeff(Nlso),klen,ktics(size(Kpath,1))
  complex(8)                                :: h(Nlso,Nlso),u(Nlso,Nlso)
  type(rgb_color)                           :: corb(Nlso),c(Nlso)
  character(len=10)                         :: chpoint
  character(len=32)                         :: fmt
  real(8),allocatable                       :: kseg(:),Ekval(:,:)
  integer,allocatable                       :: Ekcol(:,:)

  file_="Eigenbands.tb";if(present(file))file_=file
  Npts=size(kpath,1)
  Ndim=size(kpath,2)
  do iorb=1,Nlso
     corb(iorb) = colors_name(iorb)
  enddo
  !
  if(.not.set_bkvec)stop "solve_w90hk_along_BZpath ERROR: bk vectors not set!"
  select case(Ndim)
  case (1)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
  case(2)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
  case (3)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
  end select
  !
  !
  write(*,*)"Solving model along the path:"
  write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
  do ipts=1,Npts
     write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
  enddo
  !
  ic = 0
  allocate(kseg(Nk*(Npts-1)))
  allocate(ekval(Nk*(Npts-1),Nlso))
  allocate(ekcol(Nk*(Npts-1),Nlso))
  ! unit=free_unit()
  ! open(unit,file=reg(file_))
  klen=0d0  
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/dble(Nk)
     ktics(ipts)  = klen
     do ik=1,Nk
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
        ! write(unit,'(F18.12,100(F18.12,I18))')klen,(Eval(iorb),rgb(c(iorb)),iorb=1,Nlso)
        klen = klen + sqrt(dot_product(kdiff,kdiff))
     enddo
  enddo
  ktics(Npts) = kseg(ic-1)
  ! close(unit)
  !
  unit=free_unit()
  open(unit,file=reg(file_))
  do iorb=1,Nlso
     do ic=1,Nk*(Npts-1)
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
  if(Nlso<10)then !a safe compromise to have all orbitals key readable.
     do iorb=1,Nlso
        chpoint=str(0.95d0-(iorb-1)*0.05d0)
        write(unit,"(A)")str("set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
             " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
     enddo
  else
     do iorb=1,Nlso
        chpoint=str(0.95d0-(iorb-1)*0.05d0)
        write(unit,"(A)")str("#set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
             " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
     enddo
  endif
  !
  write(unit,*)"plot '"//reg(file_)//"' u 1:2:3 w l lw 3 lc rgb variable" !,\"
  ! do iorb=2,Nlso-1
  !    u1=2+(iorb-1)*2
  !    u2=3+(iorb-1)*2
  !    write(unit,*)"'"//reg(file_)//"' u 1:"&
  !         //reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable,\"
  ! enddo
  ! u1=2+(Nlso-1)*2
  ! u2=3+(Nlso-1)*2
  ! write(unit,*)"'"//reg(file_)//"' u 1:"&
  !      //reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable"
  ! !
  close(unit)
  call system("chmod +x "//reg(file_)//".gp")
end subroutine solve_Hk_along_BZpath



subroutine solve_w90Hk_along_BZpath(Nlso,kpath,Nk,colors_name,points_name,file)
  integer                                   :: Nlso
  real(8),dimension(:,:)                    :: kpath
  integer                                   :: Nk
  type(rgb_color),dimension(Nlso)           :: colors_name
  character(len=*),dimension(size(kpath,1)) :: points_name
  character(len=*),optional                 :: file
  character(len=256)                        :: file_,xtics
  integer                                   :: Npts,Nrot,Ndim
  integer                                   :: ipts,ik,ic,unit,u1,u2,iorb
  real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
  real(8)                                   :: eval(Nlso),coeff(Nlso),klen,ktics(size(Kpath,1))
  complex(8)                                :: h(Nlso,Nlso),u(Nlso,Nlso)
  type(rgb_color)                           :: corb(Nlso),c(Nlso)
  character(len=10)                         :: chpoint
  character(len=32)                         :: fmt
  real(8),allocatable                       :: kseg(:),Ekval(:,:)
  integer,allocatable                       :: Ekcol(:,:)
  !
  if(.not.TB_w90%status)stop "solve_w90hk_along_BZpath: TB_w90 structure not allocated. Call setup_w90 first."
  !
  file_="Eigenbands.tb";if(present(file))file_=file
  Npts=size(kpath,1)
  Ndim=size(kpath,2)
  do iorb=1,Nlso
     corb(iorb) = colors_name(iorb)
  enddo
  !
  if(.not.set_bkvec)stop "solve_w90hk_along_BZpath ERROR: bk vectors not set!"
  select case(Ndim)
  case (1)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
  case(2)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
  case (3)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
  end select
  !
  !
  write(*,*)"Solving model along the path:"
  write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
  do ipts=1,Npts
     write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
  enddo
  !
  ic = 0
  allocate(kseg(Nk*(Npts-1)))
  allocate(ekval(Nk*(Npts-1),Nlso))
  allocate(ekcol(Nk*(Npts-1),Nlso))
  ! unit=free_unit()
  ! open(unit,file=reg(file_))
  klen=0d0  
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/dble(Nk)
     ktics(ipts)  = klen
     do ik=1,Nk
        ic=ic+1
        kpoint = kstart + (ik-1)*kdiff
        !
        h = w90_hk_model(kpoint,Nlso)
        !
        call eigh(h,Eval)
        do iorb=1,Nlso
           coeff(:)=h(:,iorb)*conjg(h(:,iorb))
           c(iorb) = coeff.dot.corb
           Ekval(ic,iorb) = Eval(iorb)
           Ekcol(ic,iorb) = rgb(c(iorb))
        enddo
        kseg(ic) = klen
        ! write(unit,'(F18.12,100(F18.12,I18))')klen,(Eval(iorb),rgb(c(iorb)),iorb=1,Nlso)
        klen = klen + sqrt(dot_product(kdiff,kdiff))
     enddo
  enddo
  ktics(Npts) = kseg(ic-1)
  ! close(unit)
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
  if(Nlso<10)then !a safe compromise to have all orbitals key readable.
     do iorb=1,Nlso
        chpoint=str(0.95d0-(iorb-1)*0.05d0)
        write(unit,"(A)")str("set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
             " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
     enddo
  else
     do iorb=1,Nlso
        chpoint=str(0.95d0-(iorb-1)*0.05d0)
        write(unit,"(A)")str("#set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
             " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
     enddo
  endif
  !
  write(unit,*)"plot '"//reg(file_)//"' u 1:2:3 w l lw 3 lc rgb variable" !,\"
  ! do iorb=2,Nlso-1
  !    u1=2+(iorb-1)*2
  !    u2=3+(iorb-1)*2
  !    write(unit,*)"'"//reg(file_)//"' u 1:"&
  !         //reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable,\"
  ! enddo
  ! u1=2+(Nlso-1)*2
  ! u2=3+(Nlso-1)*2
  ! write(unit,*)"'"//reg(file_)//"' u 1:"&
  !      //reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable"
  !
  close(unit)
  call system("chmod +x "//reg(file_)//".gp")
end subroutine solve_w90Hk_along_BZpath





subroutine solve_HkR_along_BZpath(hkr_model,Nlat,Nso,kpath,Nkpath,colors_name,points_name,file,pbc)
  interface 
     function hkr_model(kpoint,Nlat,Nso,pbc)
       real(8),dimension(:)                    :: kpoint
       integer                                 :: Nlat,Nso
       logical                                 :: pbc
       complex(8),dimension(Nlat*Nso,Nlat*Nso) :: hkr_model
     end function hkr_model
  end interface
  integer                                      :: Nlat,Nso,Nlso
  real(8),dimension(:,:)                       :: kpath
  integer                                      :: Nkpath,Nktot
  type(rgb_color),dimension(Nlat,Nso)          :: colors_name
  character(len=*),dimension(size(kpath,1))    :: points_name
  character(len=*),optional                    :: file
  logical                                      :: pbc
  character(len=256)                           :: file_,xtics
  integer                                      :: Npts,Ndim!,units(Nlat*Nso)
  integer                                      :: ipts,ik,ic,unit,unit_,iorb,ilat,io,nrot,u1,u2
  real(8)                                      :: coeff(Nlat*Nso),klen,ktics(size(Kpath,1))
  type(rgb_color)                              :: corb(Nlat*Nso),c(Nlat*Nso)
  real(8),dimension(size(kpath,2))             :: kstart,kstop,kpoint,kdiff
  complex(8),dimension(Nlat*Nso,Nlat*Nso)      :: h
  real(8),dimension(Nlat*Nso)                  :: Eval
  character(len=64)                            :: fmt
  !
  Nlso  = Nlat*Nso
  write(fmt,"(A,I0,A)")"(I12,",Nlso,"(F18.12,I18))"
  file_ = "Eigenbands.tb";if(present(file))file_=file
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
  select case(Ndim)
  case (1)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
  case(2)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
  case (3)
     forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
  end select
  !
  write(*,*)"Solving model along the path:"
  do ipts=1,Npts
     write(*,"(A,10(A,1x),A1)")"Point"//str(ipts)//": [",(str(kpath(ipts,ic)),ic=1,size(kpath,2)),"]"
  enddo
  !
  do io=1,Nlso
     open(free_unit(unit_),file="site_"//str(io,4)//"_"//str(file_))
     rewind(unit_)
     close(unit_)
  enddo
  !
  ic=0
  klen = 0d0
  call start_timer()
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
        call eta(ic,Nktot)       
        do io=1,Nlso
           coeff(:)=h(:,io)*conjg(h(:,io))
           c(io) = coeff.dot.corb
           open(free_unit(unit_),file="site_"//str(io,4)//"_"//str(file_),position="append")
           write(unit_,"(F18.12,F18.12,I18)")klen,Eval(io),rgb(c(io))
           close(unit_)
        enddo
        klen = klen + sqrt(dot_product(kdiff,kdiff))
     enddo
  enddo
  ktics(Npts) = klen
  call stop_timer()
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
  io=1
  write(unit,*)"plot 'site_"//reg(txtfy(io,4))//"_"//reg(file_)//"' u 1:2:3 w l lw 2 lc rgb variable,\"
  do io=2,Nlso-1
     write(unit,*)"'site_"//reg(txtfy(io,4))//"_"//reg(file_)//"' u 1:2:3 w l lw 2 lc rgb variable,\"
  enddo
  write(unit,*)"'site_"//reg(txtfy(Nlso,4))//"_"//reg(file_)//"' u 1:2:3 w l lw 2 lc rgb variable"
  close(unit)
  call system("chmod +x "//reg(file_)//".gp")
end subroutine solve_HkR_along_BZpath

