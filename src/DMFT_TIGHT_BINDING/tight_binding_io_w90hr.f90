  subroutine hk_from_w90_hr(Rx,Ry,Rz,ham_k,w90_file,Nspin,Norb,Nlat,Nkvec,Porder,bandwidth,kpt_latt,Hkfile,Kpointfile)
   implicit none
   real(8)               ,intent(in)            ::   Rx(:),Ry(:),Rz(:)
   complex(8),allocatable,intent(out)           ::   ham_k(:,:,:)         !(num_wann*nspin,num_wann*nspin,num_kpts)
   character(len=*)      ,intent(in)            ::   w90_file             !"seedname_hr.dat"
   integer               ,intent(in)            ::   Nspin,Norb,Nlat
   integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
   integer               ,intent(in) ,optional  ::   Porder(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
   complex(8)            ,intent(out),optional  ::   bandwidth(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
   real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)        ![ik,3]
   character(len=*)      ,intent(in) ,optional  ::   Hkfile
   character(len=*)      ,intent(in) ,optional  ::   Kpointfile
   !
   integer                                      ::   i,j,ndx1,ndx2
   integer                                      ::   ispin,jspin,iorb,jorb,ilat,iktot
   integer                                      ::   Nkx,Nky,Nkz,num_kpts
   integer                                      ::   inrpts
   integer                                      ::   P(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
   real(8)                                      ::   bk_x(size(Rx)),bk_y(size(Ry)),bk_z(size(Rz))
   real(8),dimension(product(Nkvec),size(Nkvec))::   kpt_x,kpt_y,kpt_z
   !
   !---- W90 specific ----
   !
   !Number of Wannier orbitals
   integer                                      ::   num_wann       !=Norb*Nlat
   !Wigner-Seitz grid points
   integer                                      ::   nrpts          !=147
   !Degeneracy of the Wigner-Seitz grid points
   integer(4),allocatable                       ::   ndegen(:)      !(nrpts)
   !real space vector
   real(8)   ,allocatable                       ::   irvec(:,:)     !(3,nrpts)
   !real-space Hamiltonian
   complex(8),allocatable                       ::   ham_r(:,:,:)   !(num_wann*nspin,num_wann*nspin,nrpts)
   complex(8),allocatable                       ::   ham_aux(:,:,:)
   !local Hamiltonian
   complex(8),allocatable                       ::   Hloc(:,:)      !(num_wann*nspin,num_wann*nspin)
   !dummy vars
   real(8)                                      ::   a,b,factor_hr,rdotk
   integer                                      ::   rst,qst
   factor_hr=6.28318530717959
   !
   !
   Nkx=Nkvec(1)
   Nky=Nkvec(2)
   Nkz=Nkvec(3)
   num_kpts=Nkx*Nky*Nkz
   !
   open(unit=106,file=w90_file,status="unknown",action="read")
   read(106,*)
   read(106,*) num_wann
   read(106,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
   write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
   write(*,'(A,I6,A,I6)') "  rows:",qst,"  last row elements:",rst
   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
   !
   if(allocated(kpt_latt))deallocate(kpt_latt);allocate(kpt_latt(num_kpts,3))                           ;kpt_latt=0d0
   if(allocated(ndegen))  deallocate(ndegen)  ;allocate(ndegen(nrpts))                                  ;ndegen=0
   if(allocated(irvec))   deallocate(irvec)   ;allocate(irvec(nrpts,3))                                 ;irvec=0
   if(allocated(ham_r))   deallocate(ham_r)   ;allocate(ham_r(num_wann*Nspin,num_wann*Nspin,nrpts))     ;ham_r=zero
   if(allocated(ham_k))   deallocate(ham_k)   ;allocate(ham_k(num_wann*Nspin,num_wann*Nspin,num_kpts))  ;ham_k =zero
   if(allocated(ham_aux)) deallocate(ham_aux) ;allocate(ham_aux(num_wann*Nspin,num_wann*Nspin,num_kpts));ham_aux =zero
   if(allocated(Hloc))    deallocate(Hloc)    ;allocate(Hloc(num_wann*Nspin,num_wann*Nspin))            ;Hloc =zero
   !
   !1) k-points mesh
   call TB_set_ei(Rx,Ry,Rz)
   call TB_get_bk(bk_x,bk_y,bk_z)
   call TB_set_bk(bk_x,bk_y,bk_z)
   call TB_build_kgrid(Nkvec,kpt_x,kpt_y,kpt_z,.true.)
   if(present(kpt_latt))kpt_latt=kpt_x+kpt_y+kpt_z
   !
   !2) read WS degeneracies
   do i=1,qst
      read(106,*)(ndegen(j+(i-1)*15),j=1,15)
   enddo
   read(106,*)(ndegen(j+qst*15),j=1,rst)
   write(*,'(A)')"  degen readed"
   !
   !3) read real-space Hamiltonian (no spinup-spindw hybridizations assumed)
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            read(106,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1,ndx2,a,b
            !spin up
            ham_r(ndx1,ndx2,inrpts)=dcmplx(a,b)
            !spin dw
            if(Nspin==2)ham_r(ndx1+num_wann,ndx2+num_wann,inrpts)=dcmplx(a,b)
         enddo
      enddo
   enddo
   close(106)
   write(*,'(2A)')"  H(R) readed from: ",w90_file
   !
   !4) Fourier Transform
   do iktot=1,num_kpts 
      do i=1,num_wann*nspin
         do j=1,num_wann*nspin
            do inrpts=1,nrpts
               rdotk=0.d0
               rdotk= ( irvec(inrpts,1)*dot_product(kpt_x(iktot,:),Rx) +  &
                        irvec(inrpts,2)*dot_product(kpt_y(iktot,:),Ry) +  &
                        irvec(inrpts,3)*dot_product(kpt_z(iktot,:),Rz) )
               !
               ham_k(i,j,iktot)=ham_k(i,j,iktot)+ham_r(i,j,inrpts)*dcmplx(cos(rdotk),-sin(rdotk))/ndegen(inrpts)
               !
            enddo
         enddo
      enddo
   enddo
   if (present(Porder))then
      P=Porder
   else
      P=eye(Nspin*Norb*Nlat)
   endif
   ham_aux=zero;ham_aux=ham_k;ham_k=zero
   do iktot=1,num_kpts
      ham_k(:,:,iktot)=matmul(transpose(dble(P)),matmul(ham_aux(:,:,iktot),dble(P)))
   enddo
   !
   write(*,'(A)')"  H(k) produced"
   if(present(Hkfile))then
      call TB_write_hk(ham_k,Hkfile,Nspin*Norb*Nlat,1,1,Nlat,[Nkx,Nky,Nkz])
      write(*,'(2A)')"  H(k) written on: ",Hkfile
   endif
   !
   if(present(Kpointfile))then
      open(unit=107,file=Kpointfile,status="unknown",action="write",position="rewind")
      do iktot=1,num_kpts
         write(107,'(3F15.7)') (kpt_latt(iktot,i),i=1,3)
      enddo
      close(107)
      write(*,'(2A)')"  Kpoints used written on: ",Kpointfile
   endif
   !
   deallocate(ndegen)
   deallocate(irvec)
   deallocate(ham_r)
   deallocate(ham_aux)
   deallocate(Hloc)
   if(.not.present(kpt_latt))deallocate(kpt_latt)
   !
end subroutine hk_from_w90_hr



subroutine read_Hr_w90_solve_Hk_along_BZpath(       w90_file           &      !output file of w90
                                                ,   Nspin,Norb,Nlat    &      !dimensions
                                                ,   kpath,Nk           &      !name of Kpoint on path and mesh between
                                                ,   colors_name        &      !Colors of the Norb orbitals
                                                ,   points_name        &      !Name of the point on the path
                                                ,   Porder             &      !Ordering operator, Identity if w90 is in the correct order  (optional)
                                                ,   file_eigenband     &      !Name of the file where to save the bands                    (optional)
                                                ,   Hkpathfile         &      !Name of the file where to save the H(k) on the path         (optional)
                                                ,   Kpointpathfile     &      !Name of the file where to save the Kpoints on the path      (optional)
                                                ,   S_correction_      &      !Re{impSigma(w=0)}                                           (optional)
                                                ,   ham_k              &      !k-space hamiltonian on path                                 (optional)
                                                ,   kpt_latt           )      !Kpoint path                                                 (optional)

  implicit none
  character(len=*)      ,intent(in)            ::   w90_file
  integer               ,intent(in)            ::   Nspin,Norb,Nlat
  real(8)   ,allocatable,intent(in)            ::   kpath(:,:)
  integer               ,intent(in)            ::   Nk
  type(rgb_color)       ,intent(in)            ::   colors_name(Norb*Nlat)
  character(len=*),dimension(size(kpath,1))    ::   points_name
  integer               ,intent(in) ,optional  ::   Porder(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
  character(len=*)      ,intent(in) ,optional  ::   file_eigenband
  character(len=*)      ,intent(in) ,optional  ::   Hkpathfile
  character(len=*)      ,intent(in) ,optional  ::   Kpointpathfile
  real(8)   ,allocatable,intent(in) ,optional  ::   S_correction_(:,:) ![Nspin*Norb*Nlat,Nspin*Norb*Nlat]
  complex(8),allocatable,intent(out),optional  ::   ham_k(:,:,:)   !(num_wann*nspin,num_wann*nspin,num_kpts)
  real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)
  integer                                      ::   i,j,ndx1,ndx2
  integer                                      ::   ipts,ik,ic,u1,u2
  integer                                      ::   ispin,jspin,iorb,jorb,ilat,iktot
  integer                                      ::   inrpts,num_kpts
  character(len=256)                           ::   file_,file_s1,file_s2,file_corr_s1,file_corr_s2,xtics
  integer                                      ::   Npts
  integer                                      ::   unit
  real(8),dimension(size(kpath,2))             ::   kstart,kstop,kdiff,kpoint
  integer                                      ::   P(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
  real(8)                                      ::   U(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
  real(8)                                      ::   Udag(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
  type(rgb_color)                              ::   corb(Norb*Nlat),c(Norb*Nlat)
  character(len=10)                            ::   chpoint
  character(len=32) :: fmt
  !
  !---- W90 specific ----
  !
  !Number of Wannier orbitals
  integer                                      ::   num_wann       !=Norb*Nlat
  !Wigner-Seitz grid points
  integer                                      ::   nrpts          !=147
  !Degeneracy of the Wigner-Seitz grid points
  integer(4),allocatable                       ::   ndegen(:)      !(nrpts)
  !real space vector
  real(8)   ,allocatable                       ::   irvec(:,:)     !(3,nrpts)
  !real-space Hamiltonian
  complex(8),allocatable                       ::   ham_r(:,:,:)   !(num_wann*nspin,num_wann*nspin,nrpts)
  !k-space Hamiltonian
  complex(8),allocatable                       ::   ham_aux(:,:,:) !(num_wann*nspin,num_wann*nspin,num_kpts)
  !dummy vars
  real(8)                                      ::   a,b,factor_hr,rdotk
  integer                                      ::   rst,qst
  !
  Npts=size(kpath,1)
  num_kpts=Nk*(Npts-1)
  !
  write(*,'(A)')"  Solving model along the path:"
  write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
  do ipts=1,Npts
     write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
  enddo
  !
  open(unit=106,file=w90_file,status="unknown",action="read")
  read(106,*)
  read(106,*) num_wann
  read(106,*) nrpts
  rst=mod(nrpts,15)
  qst=int(nrpts/15)
  write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
  write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
  write(*,'(A,I6,A,I6)') "  rows:",qst,"  last row elements:",rst
  if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
  !
  if(allocated(kpt_latt))deallocate(kpt_latt);allocate(kpt_latt(num_kpts,3))                           ;kpt_latt=0d0
  if(allocated(ndegen))  deallocate(ndegen)  ;allocate(ndegen(nrpts))                                  ;ndegen=0
  if(allocated(irvec))   deallocate(irvec)   ;allocate(irvec(nrpts,3))                                 ;irvec=0
  if(allocated(ham_r))   deallocate(ham_r)   ;allocate(ham_r(num_wann*Nspin,num_wann*Nspin,nrpts))     ;ham_r=zero
  if(allocated(ham_k))   deallocate(ham_k)   ;allocate(ham_k(num_wann*Nspin,num_wann*Nspin,num_kpts))  ;ham_k =zero
  if(allocated(ham_aux)) deallocate(ham_aux) ;allocate(ham_aux(num_wann*Nspin,num_wann*Nspin,num_kpts));ham_aux =zero
  !
  !1) k-points mesh on path
  ic=0
  do ipts=1,Npts-1
     kstart = kpath(ipts,:)
     kstop  = kpath(ipts+1,:)
     kdiff  = (kstop-kstart)/dble(Nk)
     do ik=1,Nk
        ic=ic+1
        kpoint = kstart + (ik-1)*kdiff
        kpt_latt(ic,:)=kpoint
     enddo
  enddo
  !
  !2) read WS degeneracies
  do i=1,qst
     read(106,*)(ndegen(j+(i-1)*15),j=1,15)
  enddo
  read(106,*)(ndegen(j+qst*15),j=1,rst)
  write(*,'(A)')"  degen readed"
  !
  !3) read real-space Hamiltonian (no spinup-spindw hybridizations assumed change sub otherwise)
  do inrpts=1,nrpts
     do i=1,num_wann
        do j=1,num_wann
           read(106,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1,ndx2,a,b
           !spin up
           ham_r(ndx1,ndx2,inrpts)=dcmplx(a,b)
           !spin dw
           if(Nspin==2)ham_r(ndx1+num_wann,ndx2+num_wann,inrpts)=dcmplx(a,b)
        enddo
     enddo
  enddo
  close(106)
  write(*,'(2A)')"  H(R) readed from: ",w90_file
  !
  !4) Fourier Transform on path
  do iktot=1,num_kpts 
     do i=1,num_wann*nspin
        do j=1,num_wann*nspin
           do inrpts=1,nrpts
              rdotk=0.d0
              rdotk= ( kpt_latt(iktot,1)*irvec(inrpts,1) + &
                       kpt_latt(iktot,2)*irvec(inrpts,2) + &
                       kpt_latt(iktot,3)*irvec(inrpts,3) )
              !
              ham_k(i,j,iktot)=ham_k(i,j,iktot)+ham_r(i,j,inrpts)*dcmplx(cos(rdotk),-sin(rdotk))/ndegen(inrpts)
              !
           enddo
        enddo
     enddo
  enddo
  write(*,'(A)')"  H(k) on path produced"
  !
  !5) re-ordering from w90 user defined to the code standard [[[Norb],Nspin],Nlat]
  if (present(Porder))then
     P=Porder
     U=dble(P);Udag=transpose(U)
  else
     P=eye(Nspin*Norb*Nlat)
     U=dble(P);Udag=transpose(U)
  endif
  ham_aux=zero;ham_aux=ham_k;ham_k=zero
  do iktot=1,num_kpts
     ham_k(:,:,iktot)=matmul(Udag,matmul(ham_aux(:,:,iktot),U))
  enddo
  !
  if(present(Hkpathfile))then
     call TB_write_hk(ham_k,Hkpathfile,Nspin*Norb*Nlat,1,1,Nlat,Nk,kpath)
     write(*,'(2A)')"  H(k) on path written on: ",Hkpathfile
  endif
  !
  if(present(Kpointpathfile))then
     open(unit=107,file=Kpointpathfile,status="unknown",action="write",position="rewind")
     do iktot=1,num_kpts
        write(107,'(3F15.7)') (kpt_latt(iktot,i),i=1,3)
     enddo
     close(107)
     write(*,'(2A)')"  Kpoints on path used written on: ",Kpointpathfile
  endif
  !
  !6) Coloured Eigenbands on path for the two different spins
  do i=1,num_wann
     corb(i) = colors_name(i)
  enddo
  !
  file_="Eigenbands"
  if(present(file_eigenband))file_=file_eigenband
  !
  file_s1=reg(file_)//"_s1"
  call solve_bands(reg(file_s1),1)
  !
  if(Nspin.gt.1)then
     !
     file_s2=reg(file_)//"_s2"
     call solve_bands(reg(file_s2),2)
     !
     if(present(S_correction_))then
        !
        file_corr_s1=reg(file_)//"_Sreal_s1"
        call solve_bands(reg(file_corr_s1),1,S_correction_)
        !
        file_corr_s2=reg(file_)//"_Sreal_s2"
        call solve_bands(reg(file_corr_s2),2,S_correction_)
        !
     endif
  else
     if(present(S_correction_))then
        !
        file_corr_s1=reg(file_)//"_Sreal_s1"
        call solve_bands(reg(file_corr_s1),1,S_correction_)
        !
     endif
  endif



  contains

  subroutine solve_bands(file_print_,spindx,S)
    implicit none
    character(len=*),intent(in)                  ::   file_print_
    integer         ,intent(in)                  ::   spindx
    real(8)   ,allocatable,intent(in) ,optional  ::   S(:,:)
    character(len=256)                           ::   file_print
    integer                                      ::   i
    real(8)                                      ::   S_correction(Norb*Nspin*Nlat,Norb*Nspin*Nlat)
    real(8)                                      ::   eval(Norb*Nlat),coeff(Norb*Nlat)
    complex(8)                                   ::   h(Norb*Nspin*Nlat,Norb*Nspin*Nlat)
    complex(8) ,dimension(Norb*Nlat,Norb*Nlat)   ::   h1,h2,hdiag
    !
    file_print=file_print_//".dat"
    write(*,'(2A)')"  Printing eigenbands on: ",reg(file_print)
    unit=free_unit()
    open(unit,file=reg(file_print))
    !
    S_correction=0d0
    if(present(S))S_correction=S
    !
    do iktot=1,num_kpts
       h=zero
       h=ham_k(:,:,iktot)+S_correction
       !reshape with spin external block (usual for w90) so as to diagonalize a block matrix
       h=matmul(U,matmul(h,Udag))
       !identify different spins
       h1(:,:)=h(1:num_wann,1:num_wann)
       h2(:,:)=h(1+num_wann:Nspin*num_wann,1+num_wann:Nspin*num_wann)
       if(spindx==1)hdiag=h1
       if(spindx==2)hdiag=h2
       call eigh(hdiag,Eval)
       do i=1,num_wann
          coeff(:)=hdiag(:,i)*conjg(hdiag(:,i))
          c(i) = coeff.dot.corb
       enddo
       write(unit,'(I12,100(F18.12,I18))')iktot,(Eval(i),rgb(c(i)),i=1,num_wann)
    enddo
    close(unit)
    xtics="'"//reg(points_name(1))//"'1,"
    do ipts=2,Npts-1
       xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//reg(txtfy((ipts-1)*Nk+1))//","
    enddo
    xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//reg(txtfy((Npts-1)*Nk))//""
       open(unit,file=reg(file_print)//".gp")
       write(unit,*)"set term wxt"
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_print)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_print)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_print)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid noytics xtics"
       !
       write(unit,*)"plot '"//reg(file_print)//"' u 1:2:3 w l lw 3 lc rgb variable,\"
       do i=2,num_wann-1
          u1=2+(i-1)*2
          u2=3+(i-1)*2
          write(unit,*)"'"//reg(file_print)//"' u 1:"//reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable,\"
       enddo
       u1=2+(num_wann-1)*2
       u2=3+(num_wann-1)*2
       write(unit,*)"'"//reg(file_print)//"' u 1:"//reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable"
       !
       close(unit)
    call system("chmod +x "//reg(file_print)//".gp")
  end subroutine solve_bands
  !
  !
end subroutine read_Hr_w90_solve_Hk_along_BZpath











