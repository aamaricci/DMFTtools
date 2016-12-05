subroutine solve_Hk_along_BZpath(hk_model,Norb,kpath,Nk,colors_name,points_name,file)
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)                 :: kpoint
       complex(8),dimension(N,N)            :: hk_model
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
  file_="Eigenbands.tb";if(present(file))file_=file
  Npts=size(kpath,1)
  do iorb=1,Norb
     corb(iorb) = colors_name(iorb)
  enddo
  unit=free_unit()
  open(unit,file=reg(file_))
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
