subroutine solve_Hk_along_BZpath(hk_model,Norb,kpath,Nk,colors_name,points_name,file)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  integer                                   :: Norb
  real(8),dimension(:,:)                    :: kpath
  integer                                   :: Nk
  type(rgb_color),dimension(Norb)           :: colors_name
  character(len=*),dimension(size(kpath,1)) :: points_name
  character(len=*),optional                 :: file
  character(len=256)                        :: file_,xtics
  integer                                   :: Npts,Nrot
  integer                                   :: ipts,ik,ic,unit,u1,u2,iorb
  real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
  real(8)                                   :: eval(Norb),coeff(Norb)
  complex(8)                                :: h(Norb,Norb),u(Norb,Norb)
  type(rgb_color)                           :: corb(Norb),c(Norb)
  !

  !
  file_="Eigenbands.tb";if(present(file))file_=file
  Npts=size(kpath,1)
  do iorb=1,Norb
     corb(iorb) = colors_name(iorb)
  enddo
  unit=free_unit()
  open(unit,file=reg(file_))
  !
  write(*,*)"Solving model along the path:"
  do ipts=1,Npts
     write(*,"(A,10(A6,2x),A1)")"Point"//reg(str(ipts))//": [",&
          (reg(str(kpath(ipts,ic))),ic=1,size(kpath,2)-1),reg(str(kpath(ipts,size(kpath,2))))//"]"
  enddo
  !
  ic = 0
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/dble(Nk)
     do ik=1,Nk
        ic=ic+1
        kpoint = kstart + (ik-1)*kdiff
        h = hk_model(kpoint,Norb)
        call eigh(h,Eval)
        do iorb=1,Norb
           coeff(:)=h(:,iorb)*conjg(h(:,iorb))
           c(iorb) = coeff.dot.corb
        enddo
        write(unit,'(I12,100(F18.12,I18))')ic,(Eval(iorb),rgb(c(iorb)),iorb=1,Norb)
     enddo
  enddo
  close(unit)
  xtics="'"//reg(points_name(1))//"'1,"
  do ipts=2,Npts-1
     xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//reg(txtfy((ipts-1)*Nk+1))//","
  enddo
  xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//reg(txtfy((Npts-1)*Nk))//""
  open(unit,file=reg(file_)//".gp")
  write(unit,*)"set term wxt"
  write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
  write(unit,*)"#set out '"//reg(file_)//".png'"
  write(unit,*)""
  write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
  write(unit,*)"#set out '"//reg(file_)//".svg'"
  write(unit,*)""
  write(unit,*)"#set term postscript eps enhanced color 'Times'"
  write(unit,*)"#set output '|ps2pdf - "//reg(file_)//".pdf'"
  write(unit,*)"unset key"
  write(unit,*)"set xtics ("//reg(xtics)//")"
  write(unit,*)"set grid noytics xtics"
  !
  do iorb=1,Norb
     write(unit,"(A)")reg("set label 'Orb "//reg(str(iorb))//"' tc rgb "//reg(str(rgb(corb(iorb))))//&
          " at graph 0.9,"//reg(str(0.95d0-(iorb-1)*0.05d0))//" font 'Times-Italic,11'")
  enddo
  !
  write(unit,*)"plot '"//reg(file_)//"' u 1:2:3 w l lw 3 lc rgb variable,\"
  do iorb=2,Norb-1
     u1=2+(iorb-1)*2
     u2=3+(iorb-1)*2
     write(unit,*)"'"//reg(file_)//"' u 1:"&
          //reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable,\"
  enddo
  u1=2+(Norb-1)*2
  u2=3+(Norb-1)*2
  write(unit,*)"'"//reg(file_)//"' u 1:"&
       //reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable"
  !
  close(unit)
  call system("chmod +x "//reg(file_)//".gp")
end subroutine solve_Hk_along_BZpath




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
  integer                                      :: Npts,units(Nlat*Nso)
  integer                                      :: ipts,ik,ic,unit,iorb,ilat,io,nrot,u1,u2
  real(8)                                      :: coeff(Nlat*Nso)
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
  Nktot = (Npts-1)*Nkpath
  !
  do ilat=1,Nlat
     do io=1,Nso
        corb(io + (ilat-1)*Nso) = colors_name(ilat,io)
     enddo
  enddo
  !
  write(*,*)"Solving model along the path:"
  do ipts=1,Npts
     write(*,"(A,10(A,1x),A1)")"Point"//reg(str(ipts))//": [",(reg(str(kpath(ipts,ic))),ic=1,size(kpath,2)),"]"
  enddo
  !
  units=free_units(Nlso)
  do io=1,Nlso
     open(units(io),file="site_"//reg(txtfy(io,4))//"_"//reg(file_))
  enddo
  open(free_unit(unit),file=reg(file_))
  !
  ic=0
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/Nkpath
     do ik=1,Nkpath
        ic=ic+1
        kpoint = kstart + (ik-1)*kdiff
        h = hkr_model(kpoint,Nlat,Nso,pbc)
        call eigh(h,Eval)
        do io=1,Nlso
           coeff(:)=h(:,io)*conjg(h(:,io))
           c(io) = coeff.dot.corb
           write(units(io),"(I12,F18.12,I18)")ic,Eval(io),rgb(c(io))
        enddo
        write(unit,fmt)ic,(Eval(io),rgb(c(io)),io=1,Nlso)
     enddo
  enddo
  close(unit)
  do io=1,Nlso
     close(units(io))
  enddo
  xtics="'"//reg(points_name(1))//"'1,"
  do ipts=2,Npts-1
     xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//reg(txtfy((ipts-1)*Nkpath+1))//","
  enddo
  xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//reg(txtfy((Npts-1)*Nkpath))//""
  open(unit,file=reg(file_)//".gp")
  write(unit,*)"set term wxt"
  write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
  write(unit,*)"#set out '"//reg(file_)//".png'"
  write(unit,*)""
  write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
  write(unit,*)"#set out '"//reg(file_)//".svg'"
  write(unit,*)""
  write(unit,*)"#set term postscript eps enhanced color 'Times'"
  write(unit,*)"#set output '|ps2pdf - "//reg(file_)//".pdf'"
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

!>OLD VERSION:
! subroutine solve_HkR_along_BZpath(hkr_model,Nlat,Norb,kpath,Nkpath,file,pbc,jacobi)
!   integer                              :: Nlat,Norb
!   real(8),dimension(:,:)               :: kpath
!   integer                              :: Nkpath,Nktot
!   logical                              :: pbc
!   logical,optional                     :: jacobi
!   character(len=*),optional            :: file
!   character(len=256)                   :: file_,xtics
!   integer                              :: Npts
!   integer                              :: ipts,ik,ic,unit,iorb,ilat,io,nrot
!   real(8),dimension(size(kpath,2))     :: kstart,kstop,kpoint,kdiff
!   complex(8)                           :: h(Nlat*Norb,Nlat*Norb),u(Nlat*Norb,Nlat*Norb)
!   real(8),dimension(Nlat*Norb)         :: Efoo
!   real(8),dimension(:,:,:),allocatable :: Eval
!   interface 
!      function hkr_model(kpoint,Nlat,Norb,pbc)
!        real(8),dimension(:)                      :: kpoint
!        integer                                   :: Nlat,Norb
!        logical                                   :: pbc
!        complex(8),dimension(Nlat*Norb,Nlat*Norb) :: hkr_model
!      end function hkr_model
!   end interface
!   file_="Eigenbands.tb";if(present(file))file_=file
!   Npts=size(kpath,1)
!   Nktot=(Npts-1)*Nkpath
!   allocate(Eval(Nlat,Norb,Nktot))
!   ic=0
!   do ipts=1,Npts-1
!      kstart = kpath(ipts,:)
!      kstop  = kpath(ipts+1,:)
!      kdiff  = (kstop-kstart)/Nkpath
!      do ik=1,Nkpath
!         ic=ic+1
!         kpoint = kstart + (ik-1)*kdiff
!         h = hkr_model(kpoint,Nlat,Norb,pbc)
!         if(present(jacobi).AND.jacobi)then
!            call eigh_jacobi(h,Efoo,u,nrot);h=u
!         else
!            call eigh(h,Efoo)
!         endif
!         do ilat=1,Nlat
!            do iorb=1,Norb
!               io=iorb + (ilat-1)*Norb
!               Eval(ilat,iorb,ic)=Efoo(io)
!            enddo
!         enddo
!      enddo
!   enddo
!   if(ic/=Nktot)stop "solve_HkR_along_BZpath error: bad counting of the k-points along the path"
!   do iorb=1,Norb
!      open(free_unit(unit),file="l_"//reg(txtfy(iorb))//"_"//reg(file_))
!      do ilat=1,Nlat
!         do ik=1,Nktot
!            write(unit,"(I4,4F18.9)")ik,Eval(ilat,iorb,ik)
!         enddo
!         write(unit,*)""
!      enddo
!      close(unit)
!   enddo
! end subroutine solve_HkR_along_BZpath
