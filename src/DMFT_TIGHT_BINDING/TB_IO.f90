module TB_IO
  USE TB_COMMON
  USE TB_BASIS
  implicit none



contains




  !< read/write the Hamiltonian matrix H(k) and its local part 
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
    call build_kgrid(Nkvec,kgrid,.true.)
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
          write(unit,"(1000(2F15.9,1x))")(Hk(iorb,jorb),jorb=1,No)
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
    call build_kgrid(Nkvec,kgrid,.true.)
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



  subroutine write_hk_w90_path(Hk,file,No,Nd,Np,Nineq,Nkpath,kpath)
    real(8),dimension(:,:)                                         :: kpath ![Npts][Ndim]
    character(len=*)                                               :: file
    integer                                                        :: No,Nd,Np,Nineq,Nkpath,Npts
    real(8),dimension((size(kpath,1)-1)*Nkpath,size(kpath,2))      :: kgrid ![(Npts-1)*Nkpath][Ndim]
    integer                                                        :: Nktot,unit
    integer                                                        :: i,ik,iorb,jorb
    complex(8),dimension(No,No,(size(kpath,1)-1)*Nkpath)           :: Hk
    !
    !
    call kgrid_from_path_grid(kpath,Nkpath,kgrid)
    !
    Npts   = size(kpath,1)
    Nktot  = (Npts-1)*Nkpath
    !
    open(free_unit(unit),file=reg(file))
    write(unit,'(2I12)')Nkpath,Npts
    write(unit,'(1A1,100(3F12.7,1x))')"#",(kpath(i,:),i=1,Npts)
    do ik=1,Nktot
       !this is to print all 3(x,y,z) components also for lower-dimensional systems
       !to be consistent with Wannier90 format:
       write(unit,'(3(F15.9,1x))')(kgrid(ik,i),i=1,3) 
       do iorb=1,No
          write(unit,'(1000(2F15.9,1x))')(Hk(iorb,jorb,ik),jorb=1,No)
       enddo
    enddo
    close(unit)
    !
  end subroutine write_hk_w90_path











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




  subroutine read_hk_w90_path(hk,file,No,Nkpath,kpath,kgrid)
    real(8)   ,intent(inout),dimension(:,:)                        :: kpath ![Npts][Ndim]
    real(8)   ,intent(inout),dimension(:,:)                        :: kgrid ![(Npts-1)*Nkpath][Ndim]
    complex(8),intent(inout),dimension(:,:,:)                      :: Hk    ![(Npts-1)*Nkpath][Ndim]
    character(len=*)                                               :: file
    integer                                                        :: No,Nkpath,Npts
    real(8)                                                        :: kx,ky,kz
    integer                                                        :: Nktot,unit
    integer                                                        :: i,ik,iorb,jorb
    logical                                                        :: ioexist
    character(len=1)                                               :: achar
    inquire(file=reg(file),exist=ioexist)
    if(.not.ioexist)then
       write(*,*)"can not find file:"//reg(file)
       stop
    endif
    !
    open(free_unit(unit),file=reg(file))
    read(unit,'(2I12)')Nkpath,Npts
    read(unit,'(1A1,10(3F12.7,1x))')achar,(kpath(i,:),i=1,Npts)
    Nktot  = (Npts-1)*Nkpath
    !
    do ik=1,(Npts-1)*Nkpath
       read(unit,"(3(F15.9,1x))")kx,ky,kz
       kgrid(ik,:) = [kx,ky,kz]
       do iorb=1,No
          read(unit,"(1000(2F15.9,1x))")(Hk(iorb,jorb,ik),jorb=1,No)
       enddo
    enddo
    close(unit)
    !
  end subroutine read_hk_w90_path



  !< read/write the local part of the Hamiltonian to a file
  subroutine write_hloc_1(hloc,file) ![Nlso][Nlso]
    complex(8),dimension(:,:) :: Hloc
    character(len=*),optional :: file
    integer                   :: iorb,jorb,Ni,Nj,unit
    unit=6;
    if(present(file))then
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Ni=size(Hloc,1)
    Nj=size(Hloc,2)
    if(present(file))then
       do iorb=1,Ni
          write(unit,"(9000F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       do iorb=1,Ni
          write(unit,"(9000F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
       enddo
       write(unit,*)""
       close(unit)
    else
       do iorb=1,Ni
          write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
               ('(',dreal(Hloc(iorb,jorb)),',',dimag(Hloc(iorb,jorb)),')',jorb =1,Nj)
       enddo
       write(unit,*)""
    endif
  end subroutine write_hloc_1

  subroutine write_hloc_2(hloc,file) ![Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:,:,:) :: Hloc
    character(len=*),optional     :: file
    integer                       :: iorb,jorb,ispin,jspin
    integer                       :: Norb,Nspin,unit
    unit=6;
    if(present(file))then
       unit=free_unit()
       open(unit,file=reg(file))
    endif
    Nspin=size(Hloc,1)
    Norb=size(Hloc,3)
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"Write_Hloc_2","Hloc")
    if(present(file))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(9000F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(9000F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       close(unit)
    else
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
                  (&
                  (&
                  '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
                  jorb =1,Norb),&
                  jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    endif
  end subroutine write_hloc_2




  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################




  subroutine read_hloc_1(hloc,file)
    complex(8),dimension(:,:)                    :: Hloc
    character(len=*)                             :: file
    integer                                      :: iorb,jorb,Ni,Nj,unit
    real(8),dimension(size(Hloc,1),size(Hloc,2)) :: reHloc,imHloc
    unit=free_unit()   
    open(unit,file=reg(file))
    Ni=size(Hloc,1)
    Nj=size(Hloc,2)
    do iorb=1,Ni
       read(unit,"(9000F12.6)")(reHloc(iorb,jorb),jorb=1,Nj)
    enddo
    read(unit,*)
    do iorb=1,Ni
       read(unit,"(9000F12.6)")(imHloc(iorb,jorb),jorb=1,Nj)
    enddo
    close(unit)
    Hloc = dcmplx(reHloc,imHloc)
  end subroutine read_hloc_1

  subroutine read_hloc_2(hloc,file)
    complex(8),dimension(:,:,:,:)                                          :: Hloc
    character(len=*)                                                       :: file
    integer                                                                :: iorb,jorb,ispin,jspin,unit
    integer                                                                :: Nspin,Norb
    real(8),dimension(size(Hloc,1),size(Hloc,2),size(Hloc,3),size(Hloc,4)) :: reHloc,imHloc
    unit=free_unit()   
    open(unit,file=reg(file))
    Nspin=size(Hloc,1)
    if(size(Hloc,2)/=Nspin)stop "read_hloc error: size[Hloc,2] != size[Hloc,1] == Nspin"
    Norb=size(Hloc,3)
    if(size(Hloc,4)/=Nspin)stop "read_hloc error: size[Hloc,4] != size[Hloc,3] == Norb"
    do ispin=1,Nspin
       do iorb=1,Norb
          read(unit,"(9000F12.6)")((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    read(unit,*)
    do ispin=1,Nspin
       do iorb=1,Norb
          read(unit,"(9000F12.6)")((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    close(unit)
    Hloc = dcmplx(reHloc,imHloc)
  end subroutine read_hloc_2








END MODULE TB_IO











