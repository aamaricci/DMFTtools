subroutine write_hk_w90_func(hk_model,file,No,Nd,Np,Nineq,Nkvec)
  character(len=*)                              :: file
  integer                                       :: No,Nd,Np,Nineq
  integer                                       :: Nkvec(:)
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
  real(8),dimension(3)                          :: kvec
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
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)
  !
  open(free_unit(unit),file=reg(file))
  write(unit,'(10(A10,1x))')str(Nktot),str(No),str(Nd),str(Np),str(Nineq)
  write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,size(Nkvec))
  do ik=1,Nktot
     Hk(:,:) = hk_model(kgrid(ik,:),No)
     !this is to print all 3(x,y,z) components also for lower-dimensional systems
     !to be consistent with Wannier90 format:
     kvec=0d0
     kvec(:size(Nkvec)) = kgrid(ik,:size(Nkvec))
     write(unit,"(3(F15.9,1x))")(kvec(i),i=1,3) 
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
  real(8),dimension(3)                          :: kvec
  integer                                       :: Nktot,unit
  integer                                       :: i,ik,iorb,jorb
  complex(8),dimension(No,No,product(Nkvec))    :: Hk
  !
  !
  call TB_build_kgrid(Nkvec,kgrid,.true.)
  !
  Nktot  = product(Nkvec)!==size(Hk,3)
  !
  open(free_unit(unit),file=reg(file))
  write(unit,'(10(A10,1x))')str(Nktot),str(No),str(Nd),str(Np),str(Nineq)
  write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,size(Nkvec))
  do ik=1,Nktot
     !this is to print all 3(x,y,z) components also for lower-dimensional systems
     !to be consistent with Wannier90 format:
     kvec=0d0
     kvec(:size(Nkvec)) = kgrid(ik,:size(Nkvec))
     write(unit,"(3(F15.9,1x))")(kvec(i),i=1,3) 
     do iorb=1,No
        write(unit,"(1000(2F15.9,1x))")(Hk(iorb,jorb,ik),jorb=1,No)
     enddo
  enddo
  close(unit)
  !
end subroutine write_hk_w90_array











!##################################################################
!##################################################################
!##################################################################
!##################################################################







subroutine read_hk_w90_array(hk,file,No,Nd,Np,Nineq,Nkvec,kgrid)
  character(len=*)                              :: file
  integer                                       :: No,Nd,Np,Nineq
  integer                                       :: Nkvec(:)
  real(8),dimension(product(Nkvec),size(Nkvec)) :: kgrid ![Nk][Ndim]
  integer                                       :: Nktot,unit,Nk(3),Nk_
  integer                                       :: ik,ix,iy,iz,iorb,jorb
  real(8)                                       :: kx,ky,kz,kvec(3)
  complex(8),dimension(No,No,product(Nkvec))    :: Hk
  logical                                       :: ioexist
  character(len=1)                              :: achar
  inquire(file=reg(file),exist=ioexist)
  if(.not.ioexist)then
     write(*,*)"can not find file:"//reg(file)
     stop
  endif
  !
  open(free_unit(unit),file=reg(file))
  read(unit,'(10(A10,1x))')!Nktot,No,Nd,Np,Nineq
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
           kgrid(ik,:) = [kx,ky,kz]
           do iorb=1,No
              read(unit,"(1000(2F15.9,1x))")(Hk(iorb,jorb,ik),jorb=1,No)
           enddo
        enddo
     enddo
  enddo
  close(unit)
  !
end subroutine read_hk_w90_array

