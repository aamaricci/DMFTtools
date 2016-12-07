subroutine write_hk_w90_func(hk_model,file,No,Nd,Np,Nineq,Nkvec)
  character(len=*)                              :: file
  integer                                       :: No,Nd,Np,Nineq
  integer                                       :: Nkvec(:)
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
  integer                                       :: Nktot,unit
  integer                                       :: i,ik,iorb,jorb
  complex(8),dimension(No,No)                   :: hk
  interface 
     function hk_model(kpoint,N)
       real(8),dimension(:)      :: kpoint
       integer                   :: N
       complex(8),dimension(N,N) :: hk_model
     end function hk_model
  end interface
  !
  kgrid = TB_build_kgrid(Nkvec,.true.)
  !
  open(free_unit(unit),file=reg(file))
  write(unit,'(1A1,1x,3(A2,1x))')"#",reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
  write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,size(Nkvec))
  Nktot  = product(Nkvec)
  do ik=1,Nktot
     Hk(:,:) = hk_model(kgrid(ik,:),No)
     write(unit,"(3(F15.9,1x))")(kgrid(ik,i),i=1,size(Nkvec))
     do iorb=1,No
        write(unit,"(20(2F15.9,1x))")(Hk(iorb,jorb),jorb=1,No)
     enddo
  enddo
  close(unit)
  !
end subroutine write_hk_w90_func



subroutine write_hk_w90_array(Hk,file,No,Nd,Np,Nineq,Nkvec)
  character(len=*)                              :: file
  integer                                       :: No,Nd,Np,Nineq
  integer                                       :: Nkvec(:)
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
  integer                                       :: Nktot,unit
  integer                                       :: i,ik,iorb,jorb
  complex(8),dimension(No,No,product(Nkvec))    :: Hk
  !
  !
  kgrid = TB_build_kgrid(Nkvec,.true.)
  !
  Nktot  = product(Nkvec)!==size(Hk,3)
  !
  open(free_unit(unit),file=reg(file))
  write(unit,'(1A1,1x,3(A2,1x))')"#",reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
  write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,size(Nkvec))
  do ik=1,Nktot
     write(unit,"(3(F15.9,1x))")(kgrid(ik,i),i=1,size(Nkvec))
     do iorb=1,No
        write(unit,"(20(2F15.9,1x))")(Hk(iorb,jorb,ik),jorb=1,No)
     enddo
  enddo
  close(unit)
  !
end subroutine write_hk_w90_array





! subroutine write_hk_w90_func_1d(hk_model,file,No,Nd,Np,Nineq,Nkvec)
!   character(len=*)                              :: file
!   integer                                       :: No,Nd,Np,Nineq
!   integer                                       :: Nkvec(:)
!   real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
!   integer                                       :: Nktot,unit
!   integer                                       :: ik,iorb,jorb
!   complex(8)                                    :: hk
!   interface 
!      function hk_model(kpoint)
!        real(8),dimension(:)      :: kpoint
!        complex(8)                :: hk_model
!      end function hk_model
!   end interface
!   !
!   Nk=1
!   do ik=1,size(Nkvec)
!      Nk(ik)=Nkvec(ik)
!   enddo
!   Nktot  = product(Nk)
!   if(Nktot/=product(Nkvec))stop "write_hk_w90_ ERROR: product(Nkvec) != Nktot"
!   !
!   open(free_unit(unit),file=reg(file))
!   write(unit,'(1A1,1x,3(A2,1x))')"#",reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
!   write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nk(ik))),ik=1,3)
!   ik=0
!   do iz=1,Nk(3)
!      kz = dble(iz-1)/Nk(3)
!      do iy=1,Nk(2)
!         ky = dble(iy-1)/Nk(2)
!         do ix=1,Nk(1)
!            kx = dble(ix-1)/Nk(1)
!            kvec = kx*bk_x + ky*bk_y + kz*bk_z
!            ik = ik+1
!            Hk = hk_model(kvec)
!            write(unit,"(3(F15.9,1x))")kvec(1),kvec(2),kvec(3)
!            write(unit,"(20(2F15.9,1x))")Hk
!         enddo
!      enddo
!   enddo
!   close(unit)
!   !
! end subroutine write_hk_w90_func_1d

! subroutine write_hk_w90_array_1d(hk,file,No,Nd,Np,Nineq,Nkvec)
!   character(len=*)                     :: file
!   integer                              :: No,Nd,Np,Nineq
!   integer                              :: Nkvec(:)
!   integer                              :: Nktot,unit,Nk(3)
!   integer                              :: ik,ix,iy,iz,iorb,jorb
!   real(8)                              :: kx,ky,kz,kvec(3)
!   complex(8),dimension(product(Nkvec)) :: Hk    
!   !
!   Nk=1
!   do ik=1,size(Nkvec)
!      Nk(ik)=Nkvec(ik)
!   enddo
!   Nktot  = product(Nk)
!   if(Nktot/=product(Nkvec))stop "write_hk_w90_ ERROR: product(Nkvec) != Nktot"
!   !
!   open(free_unit(unit),file=reg(file))
!   write(unit,'(1A1,1x,3(A2,1x))')"#",reg(txtfy(Nd)),reg(txtfy(Np)),reg(txtfy(Nineq))
!   write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nk(ik))),ik=1,3)
!   ik=0
!   do iz=1,Nk(3)
!      kz = dble(iz-1)/Nk(3)
!      do iy=1,Nk(2)
!         ky = dble(iy-1)/Nk(2)
!         do ix=1,Nk(1)
!            kx = dble(ix-1)/Nk(1)
!            kvec = kx*bk_x + ky*bk_y + kz*bk_z
!            ik = ik+1
!            write(unit,"(3(F15.9,1x))")kx,ky,kz
!            write(unit,"((2F15.9,1x))")Hk(ik)
!         enddo
!      enddo
!   enddo
!   close(unit)
! end subroutine write_hk_w90_array_1d







!##################################################################
!##################################################################
!##################################################################
!##################################################################







subroutine read_hk_w90_array(hk,file,No,Nd,Np,Nineq,Nkvec,kgrid)
  character(len=*)                           :: file
  integer                                    :: No,Nd,Np,Nineq
  integer                                    :: Nkvec(:)
  real(8),dimension(3,product(Nkvec))        :: kgrid
  integer                                    :: Nktot,unit,Nk(3),Nk_
  integer                                    :: ik,ix,iy,iz,iorb,jorb
  real(8)                                    :: kx,ky,kz,kvec(3)
  complex(8),dimension(No,No,product(Nkvec)) :: Hk
  logical                                    :: ioexist
  character(len=1)                           :: achar
  inquire(file=reg(file),exist=ioexist)
  if(.not.ioexist)then
     write(*,*)"can not find file:"//reg(file)
     stop
  endif
  !
  open(free_unit(unit),file=reg(file))
  read(unit,'(1A1,1x,3(I2,1x))')achar,Nd,Np,Nineq
  read(unit,'(1A1,3(I12,1x))')achar,( Nk(ik),ik=1,3 )
  Nktot  = product(Nk)
  if(Nktot/=product(Nkvec))stop "read_hk_w90_ ERROR: product(Nkvec) != Nktot"
  !
  ik=0
  do iz=1,Nk(3)
     do iy=1,Nk(2)
        do ix=1,Nk(1)
           ik = ik+1
           read(unit,"(3(F15.9,1x))")kx,ky,kz
           kgrid(:,ik) = [kx,ky,kz]
           do iorb=1,No
              read(unit,"(20(2F15.9,1x))")(Hk(iorb,jorb,ik),jorb=1,No)
           enddo
        enddo
     enddo
  enddo
  close(unit)
  !
end subroutine read_hk_w90_array



! subroutine read_hk_w90_array_1d(hk,file,No,Nd,Np,Nineq,Nkvec,kgrid)
!   character(len=*)                    :: file
!   integer                             :: No,Nd,Np,Nineq
!   integer                             :: Nkvec(:)
!   real(8),dimension(3,product(Nkvec)) :: kgrid
!   integer                             :: Nktot,unit,Nk(3),Nk_
!   integer                             :: ik,ix,iy,iz,iorb,jorb
!   real(8)                             :: kx,ky,kz,kvec(3)
!   logical                             :: ioexist
!   complex(8),dimension(:)             :: Hk
!   character(len=1)                    :: achar
!   !
!   inquire(file=reg(file),exist=ioexist)
!   if(.not.ioexist)then
!      write(*,*)"can not find file:"//reg(file)
!      stop
!   endif
!   !
!   open(free_unit(unit),file=reg(file))
!   read(unit,'(1A1,1x,3(I2,1x))')achar,Nd,Np,Nineq
!   read(unit,'(1A1,3(I12,1x))')achar,( Nk(ik),ik=1,3 )
!   Nktot  = product(Nk)
!   if(Nktot/=product(Nkvec))stop "read_hk_w90_ ERROR: product(Nkvec) != Nktot"
!   !
!   ik=0
!   do iz=1,Nk(3)
!      do iy=1,Nk(2)
!         do ix=1,Nk(1)
!            ik = ik+1
!            read(unit,"(3(F15.9,1x))")kx,ky,kz
!            kgrid(:,ik) = [kx,ky,kz]
!            read(unit,"(20(2F15.9,1x))")Hk(ik)
!         enddo
!      enddo
!   enddo
!   close(unit)
!   !
! end subroutine read_hk_w90_array_1d
