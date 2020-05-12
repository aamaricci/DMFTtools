  subroutine hk_from_w90_hr(R1,R2,R3,ham_k,ham_loc,w90_file,Nspin,Norb,Nlat,Nkvec,kpt_latt,Hkfile,Kpointfile)
   implicit none
   real(8)               ,intent(in)            ::   R1(:),R2(:),R3(:)
   complex(8),allocatable,intent(inout)         ::   ham_k(:,:,:)
   complex(8),allocatable,intent(inout)         ::   ham_loc(:,:)
   character(len=*)      ,intent(in)            ::   w90_file
   integer               ,intent(in)            ::   Nspin,Norb,Nlat
   integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
   real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)        ![ik,3]
   character(len=*)      ,intent(in) ,optional  ::   Hkfile
   character(len=*)      ,intent(in) ,optional  ::   Kpointfile
   !
   logical                                      ::   IOfile
   integer                                      ::   unitIO
   integer                                      ::   i,j,ndx1,ndx2
   integer                                      ::   Nkx,Nky,Nkz
   integer                                      ::   iktot,num_kpts
   integer                                      ::   inrpts
   real(8)                                      ::   bk1(size(R1)),bk2(size(R2)),bk3(size(R3))
   real(8),dimension(product(Nkvec),size(Nkvec))::   kpt1,kpt2,kpt3
   real(8)                                      ::   a,b,rdotk
   integer                                      ::   rst,qst
   !---- W90 specific ----
   integer                                      ::   num_wann       !=Norb*Nlat
   integer                                      ::   nrpts
   integer(4),allocatable                       ::   ndegen(:)      !(nrpts)
   integer   ,allocatable                       ::   irvec(:,:)     !(3,nrpts)
   complex(8),allocatable                       ::   ham_r(:,:,:)
   complex(8),allocatable                       ::   ham_aux(:,:,:)
   !
   !
   Nkx=Nkvec(1)
   Nky=Nkvec(2)
   Nkz=Nkvec(3)
   num_kpts=Nkx*Nky*Nkz
   !
   unitIO=free_unit()
   open(unit=unitIO,file=w90_file,status="old",action="read")
   read(unitIO,*)
   read(unitIO,*) num_wann
   read(unitIO,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   write(*,*)
   write(*,'(1A)')         "-------------- H_LDA --------------"
   write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
   write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
   write(*,'(A,I6,A,I6)') "  rows:",qst,"  last row elements:",rst
   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
   !
   if(allocated(kpt_latt))deallocate(kpt_latt);allocate(kpt_latt(num_kpts,3))                           ;kpt_latt=0d0
   if(allocated(ndegen))  deallocate(ndegen)  ;allocate(ndegen(nrpts))                                  ;ndegen=0
   if(allocated(irvec))   deallocate(irvec)   ;allocate(irvec(nrpts,3))                                 ;irvec=0
   if(allocated(ham_r))   deallocate(ham_r)   ;allocate(ham_r(num_wann*Nspin,num_wann*Nspin,nrpts))     ;ham_r=zero
   if(allocated(ham_aux)) deallocate(ham_aux) ;allocate(ham_aux(num_wann*Nspin,num_wann*Nspin,num_kpts));ham_aux=zero
   !
   !1) k-points mesh
   call TB_set_ei(R1,R2,R3)
   call TB_get_bk(bk1,bk2,bk3)
   call TB_set_bk(bk1,bk2,bk3)
   call build_kgrid_generic(Nkvec,kpt1,kpt2,kpt3,.true.)
   if(present(kpt_latt))kpt_latt=kpt1+kpt2+kpt3
   !
   !2) read WS degeneracies
   do i=1,qst
      read(unitIO,*)(ndegen(j+(i-1)*15),j=1,15)
   enddo
   if(rst.ne.0)read(unitIO,*)(ndegen(j+qst*15),j=1,rst)
   write(*,'(1A)')"  degen readed"
   !
   !3) read real-space Hamiltonian (no spinup-spindw hybridizations assumed)
   ham_loc=zero
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            read(unitIO,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1,ndx2,a,b
            !spin up
            ham_r(ndx1,ndx2,inrpts)=dcmplx(a,b)
            !spin dw
            if(Nspin==2)ham_r(ndx1+num_wann,ndx2+num_wann,inrpts)=dcmplx(a,b)
            !Hloc
            if(irvec(inrpts,1)==0.and.irvec(inrpts,2)==0.and.irvec(inrpts,3)==0)then
               ham_loc(ndx1,ndx2)=dcmplx(a,b)
               if(Nspin==2)ham_loc(ndx1+num_wann,ndx2+num_wann)=dcmplx(a,b)
            endif
         enddo
      enddo
   enddo
   close(unitIO)
   write(*,'(2A)')"  H(R) readed from: ",w90_file
   !
   !4) Fourier Transform
   do iktot = 1,num_kpts
      do inrpts=1,nrpts
         rdotk=0.d0
         rdotk= ( irvec(inrpts,1)*dot_product(kpt1(iktot,:),R1) +  &
                  irvec(inrpts,2)*dot_product(kpt2(iktot,:),R2) +  &
                  irvec(inrpts,3)*dot_product(kpt3(iktot,:),R3) )
         do i=1,num_wann*nspin
            do j=1,num_wann*nspin
               !
               ham_aux(i,j,iktot)=ham_aux(i,j,iktot)+ham_r(i,j,inrpts)*dcmplx(cos(rdotk),-sin(rdotk))/ndegen(inrpts)
               !
            enddo
         enddo
      enddo
   enddo
   ham_k=ham_aux
   !
   !5) Reordering & hermicity check
   if(Nspin==2)then
      ham_aux=zero;ham_aux=ham_k;ham_k=zero
      do iktot=1,num_kpts
         ham_k(:,:,iktot)=slo2lso(ham_aux(:,:,iktot),Nlat,Nspin,Norb)
      enddo
      ham_loc=slo2lso(ham_loc,Nlat,Nspin,Norb)
   endif
   do iktot=1,num_kpts
      call herm_check(ham_k(:,:,iktot))
   enddo
   call herm_check(ham_loc)
   !
   write(*,'(1A)')"  H(k) produced"
   if(present(Hkfile))then
      call write_hk_w90_array(ham_k,Hkfile,Nlat,Nspin,Norb,[Nkx,Nky,Nkz])
      write(*,'(2A)')"  H(k) written on: ",Hkfile
   endif
   !
   if(present(Kpointfile))then
      unitIO=free_unit()
      open(unit=unitIO,file=Kpointfile,status="unknown",action="write",position="rewind")
      do iktot=1,num_kpts
         write(unitIO,'(3F15.7)') (kpt_latt(iktot,i),i=1,3)
      enddo
      close(unitIO)
      write(*,'(2A)')"  Kpoints used written on: ",Kpointfile
   endif
   !
   deallocate(ndegen)
   deallocate(irvec)
   deallocate(ham_r)
   deallocate(ham_aux)
   if(.not.present(kpt_latt))deallocate(kpt_latt)
   !
  end subroutine hk_from_w90_hr






  subroutine hkt_from_w90_hr(field,gauge,R1,R2,R3,Ruc,ham_kt,w90_file,dipole_file,Nspin,Norb,Nlat,Nt,Nkvec)
   implicit none
   real(8)               ,intent(in)            ::   field(:,:,:) ![Nt,dim,2] 1=Efield 2=Afield
   character(len=*)      ,intent(in)            ::   gauge
   real(8)               ,intent(in)            ::   R1(:),R2(:),R3(:)
   real(8)               ,intent(in)            ::   Ruc(:,:)
   complex(8),allocatable,intent(inout)         ::   ham_kt(:,:,:,:)
   character(len=*)      ,intent(in)            ::   w90_file
   character(len=*)      ,intent(in)            ::   dipole_file
   integer               ,intent(in)            ::   Nspin,Norb,Nlat,Nt
   integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
   logical                                      ::   IOfile
   integer                                      ::   unitIO1,unitIO2
   integer                                      ::   i,j,iorb,jorb,io,jo,it
   integer                                      ::   ndx1_H,ndx2_H,ndx1_D,ndx2_D
   integer                                      ::   Nkx,Nky,Nkz
   integer                                      ::   iktot,num_kpts
   integer                                      ::   inrpts
   real(8)                                      ::   bk1(size(R1)),bk2(size(R2)),bk3(size(R3))
   real(8),dimension(product(Nkvec),size(Nkvec))::   kpt1,kpt2,kpt3
   real(8)                                      ::   a,b,Dx,Dy,Dz,rdotk,exparg
   integer                                      ::   rst,qst,limit,kvec_ndx
   integer                                      ::   auxndx,dumR1,dumR2,dumR3
   !---- light matter ----
   integer   ,allocatable,dimension(:)          ::   Kvec
   integer   ,allocatable,dimension(:,:)        ::   site_ndx
   integer   ,allocatable,dimension(:,:,:)      ::   veclist
   complex(8),allocatable,dimension(:,:,:)      ::   ham_r,StructFact
   complex(8),allocatable,dimension(:,:,:,:)    ::   lightmat_r,dip_r,ham_auxt,ham_rt
   !---- W90 specific ----
   integer                                      ::   num_wann       !=Norb*Nlat
   integer                                      ::   nrpts
   integer(4),allocatable                       ::   ndegen(:)      !(nrpts)
   integer   ,allocatable                       ::   irvec(:,:)     !(3,nrpts)
   !
   !
   Nkx=Nkvec(1)
   Nky=Nkvec(2)
   Nkz=Nkvec(3)
   num_kpts=Nkx*Nky*Nkz
   !
   unitIO1=free_unit()
   open(unit=unitIO1,file=w90_file,status="old",action="read")
   unitIO2=free_unit()
   open(unit=unitIO2,file=dipole_file,status="old",action="read")
   read(unitIO1,*)
   read(unitIO1,*) num_wann
   read(unitIO1,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   write(*,*)
   write(*,'(1A)')         "-------------- Ht_LDA --------------"
   write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
   write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
   write(*,'(A,I6,1A,I6)') "  rows:",qst,"  last row elements:",rst
   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
   !
   if(allocated(ndegen))    deallocate(ndegen)    ;allocate(ndegen(nrpts))                                       ;ndegen=0
   if(allocated(irvec))     deallocate(irvec)     ;allocate(irvec(nrpts,3))                                      ;irvec=0
   if(allocated(Kvec))      deallocate(Kvec)      ;allocate(Kvec(3))                                             ;Kvec=0
   if(allocated(veclist))   deallocate(veclist)   ;allocate(veclist(-10:10,-10:10,-10:10))                       ;veclist=0
   if(allocated(site_ndx))  deallocate(site_ndx)  ;allocate(site_ndx(nrpts,2))                                   ;site_ndx=0
   if(allocated(ham_r))     deallocate(ham_r)     ;allocate(ham_r(num_wann,num_wann,nrpts))                      ;ham_r=zero
   if(allocated(dip_r))     deallocate(dip_r)     ;allocate(dip_r(num_wann,num_wann,nrpts,3))                    ;dip_r=zero
   if(allocated(StructFact))deallocate(StructFact);allocate(StructFact(num_wann,num_wann,Nt))                    ;StructFact=zero
   if(allocated(lightmat_r))deallocate(lightmat_r);allocate(lightmat_r(num_wann,num_wann,nrpts,3))               ;lightmat_r=zero
   !
   if(allocated(ham_rt))    deallocate(ham_rt)    ;allocate(ham_rt(num_wann,num_wann,nrpts,Nt))                  ;ham_rt=zero
   if(allocated(ham_auxt))  deallocate(ham_auxt)  ;allocate(ham_auxt(num_wann*Nspin,num_wann*Nspin,num_kpts,Nt)) ;ham_auxt=zero
   !
   !1) k-points mesh
   call TB_set_ei(R1,R2,R3)
   call TB_get_bk(bk1,bk2,bk3)
   call TB_set_bk(bk1,bk2,bk3)
   call build_kgrid_generic(Nkvec,kpt1,kpt2,kpt3,.true.)
   !
   !2) read WS degeneracies
   do i=1,qst
      read(unitIO1,*)(ndegen(j+(i-1)*15),j=1,15)
   enddo
   if(rst.ne.0)read(unitIO1,*)(ndegen(j+qst*15),j=1,rst)
   write(*,'(1A)')"  degen readed"
   !
   !3) read real-space quantities
   limit=0
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            !
            !read H(R) & D(R)
            read(unitIO1,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H,a,b
            read(unitIO2,*)           dumR1,         dumR2,          dumR3,ndx1_D,ndx2_D,Dx,Dy,Dz
            !
            !consistency check
            auxndx = sum([irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H]-[dumR1,dumR2,dumR3,ndx1_D,ndx2_D])
            if(auxndx.ne.0)then
               write(*,'(10A)') "  Something is wrong between ",w90_file," and ",dipole_file," indexing"
               write(*,'(10I5)')irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H
               write(*,'(10I5)')dumR1,dumR2,dumR3,ndx1_D,ndx2_D
               stop
            endif
            !
            if(abs(dumR1).gt.limit)limit=abs(dumR1)
            veclist(irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3))=inrpts
            !
            site_ndx(inrpts,1)=floor((ndx1_H-0.01)/Norb)+1
            site_ndx(inrpts,2)=floor((ndx2_H-0.01)/Norb)+1
            !
            ham_r(ndx1_H,ndx2_H,inrpts)=dcmplx(a,b)
            dip_r(ndx1_D,ndx2_D,inrpts,1)=dcmplx(Dx,0.d0)
            dip_r(ndx1_D,ndx2_D,inrpts,2)=dcmplx(Dy,0.d0)
            dip_r(ndx1_D,ndx2_D,inrpts,3)=dcmplx(Dz,0.d0)
            !
         enddo
      enddo
   enddo
   close(unitIO1)
   close(unitIO2)
   write(*,'(1A)')"  H(R) and D(R) readed"
   !
   !4) build light-matter interaction
   if(gauge=="A")then
      !
      do i=1,nrpts
         do j=1,nrpts
            !
            Kvec(1)=irvec(i,1)-irvec(j,1)
            Kvec(2)=irvec(i,2)-irvec(j,2)
            Kvec(3)=irvec(i,3)-irvec(j,3)
            !
            if((abs(Kvec(1)).gt.limit) .or. &
               (abs(Kvec(2)).gt.limit) .or. &
               (abs(Kvec(3)).gt.limit)      )cycle
            !
            Kvec_ndx=veclist(Kvec(1),Kvec(2),Kvec(3))
            !
            lightmat_r(:,:,i,1) = lightmat_r(:,:,i,1) + matmul(dip_r(:,:,Kvec_ndx,1),ham_r(:,:,j))
            lightmat_r(:,:,i,2) = lightmat_r(:,:,i,2) + matmul(dip_r(:,:,Kvec_ndx,2),ham_r(:,:,j))
            lightmat_r(:,:,i,3) = lightmat_r(:,:,i,3) + matmul(dip_r(:,:,Kvec_ndx,3),ham_r(:,:,j))
            !
         enddo
      enddo
      !
      do i=1,nrpts
         do j=1,nrpts
            !
            Kvec(1)=irvec(i,1)-irvec(j,1)
            Kvec(2)=irvec(i,2)-irvec(j,2)
            Kvec(3)=irvec(i,3)-irvec(j,3)
            !
            if((abs(Kvec(1)).gt.limit) .or. &
               (abs(Kvec(2)).gt.limit) .or. &
               (abs(Kvec(3)).gt.limit)      )cycle
            !
            Kvec_ndx=veclist(Kvec(1),Kvec(2),Kvec(3))
            !
            lightmat_r(:,:,i,1) = lightmat_r(:,:,i,1) - matmul(ham_r(:,:,Kvec_ndx),dip_r(:,:,j,1))
            lightmat_r(:,:,i,2) = lightmat_r(:,:,i,2) - matmul(ham_r(:,:,Kvec_ndx),dip_r(:,:,j,2))
            lightmat_r(:,:,i,3) = lightmat_r(:,:,i,3) - matmul(ham_r(:,:,Kvec_ndx),dip_r(:,:,j,3))
            !
         enddo
         call herm_check(lightmat_r(:,:,i,1))
         call herm_check(lightmat_r(:,:,i,2))
         call herm_check(lightmat_r(:,:,i,3))
      enddo
      !
   elseif(gauge=="E")then
      !
      do i=1,nrpts
         !
         lightmat_r(:,:,i,1) = dip_r(:,:,i,1) ; call herm_check(lightmat_r(:,:,i,1))
         lightmat_r(:,:,i,2) = dip_r(:,:,i,2) ; call herm_check(lightmat_r(:,:,i,2))
         lightmat_r(:,:,i,3) = dip_r(:,:,i,3) ; call herm_check(lightmat_r(:,:,i,3))
         !
      enddo
      !
   endif
   deallocate(veclist,Kvec,dip_r)
   write(*,'(2A)')"  light-matter interaction built in gauge: ",gauge
   !
   !5) build prefactor
   if(gauge=="A")then
      !
      StructFact=dcmplx(1.d0,0.d0)
      !
   elseif(gauge=="E")then
      !
      StructFact=dcmplx(1.d0,0.d0)
      do it=1,Nt
         do i=1,Nlat
            do j=1,Nlat
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     io = iorb + (i-1)*Norb
                     jo = jorb + (j-1)*Norb
                     !
                     exparg = ( -field(it,1,2) * ( Ruc(i,1) - Ruc(j,1) ) &
                                -field(it,2,2) * ( Ruc(i,2) - Ruc(j,2) ) &
                                -field(it,3,2) * ( Ruc(i,3) - Ruc(j,3) ) )
                     !
                     StructFact(io,jo,it) = dcmplx(cos(exparg),sin(exparg))
                     !
                  enddo
               enddo
            enddo
         enddo
         call herm_check(StructFact(:,:,it))
      enddo
      !
   endif
   write(*,'(1A)')"  prefactor built"
   !
   !6) build interacting hamilt in real space
   if(gauge=="A")then
      !
      do inrpts=1,nrpts
         do it=1,Nt
            !
            ham_rt(:,:,inrpts,it) = ( ham_r(:,:,inrpts) + field(it,1,2) * Xi * lightmat_r(:,:,inrpts,1) &
                                                        + field(it,2,2) * Xi * lightmat_r(:,:,inrpts,2) &
                                                        + field(it,3,2) * Xi * lightmat_r(:,:,inrpts,3) )
            !
         enddo
      enddo
      !
   elseif(gauge=="E")then
      !
      do it=1,Nt
         do inrpts=1,nrpts
            !
            exparg = ( -field(it,1,2) * ( irvec(inrpts,1)*R1(1) + irvec(inrpts,2)*R2(1) + irvec(inrpts,3)*R3(1) ) &
                       -field(it,2,2) * ( irvec(inrpts,1)*R1(2) + irvec(inrpts,2)*R2(2) + irvec(inrpts,3)*R3(2) ) &
                       -field(it,3,2) * ( irvec(inrpts,1)*R1(3) + irvec(inrpts,2)*R2(3) + irvec(inrpts,3)*R3(3) ) )
            !
            ham_rt(:,:,inrpts,it) = StructFact(:,:,it)*dcmplx(cos(exparg),sin(exparg)) *          &
                                  ( ham_r(:,:,inrpts)+ field(it,1,1) * lightmat_r(:,:,inrpts,1)   &
                                                     + field(it,2,1) * lightmat_r(:,:,inrpts,2)   &
                                                     + field(it,3,1) * lightmat_r(:,:,inrpts,3)   )
            !
            !
         enddo
      enddo
      !
   endif
   deallocate(StructFact,site_ndx,ham_r,lightmat_r)
   write(*,'(1A)')"  real-space H(R,t) built"
   !
   !7) Fourier Transform
   ham_kt=zero
   do iktot=1,num_kpts
      do it=1,Nt
         do inrpts=1,nrpts
            rdotk=0.d0
            rdotk= (   irvec(inrpts,1)*dot_product(kpt1(iktot,:),R1) &
                     + irvec(inrpts,2)*dot_product(kpt2(iktot,:),R2) &
                     + irvec(inrpts,3)*dot_product(kpt3(iktot,:),R3) )
            do i=1,num_wann
               do j=1,num_wann
                  !
                  ham_auxt(i,j,iktot,it) = ham_auxt(i,j,iktot,it) + &
                                           ham_rt(i,j,inrpts,it)  * dcmplx(cos(rdotk),-sin(rdotk)) / ndegen(inrpts)
                  !
                  if(Nspin==2)ham_auxt(i+num_wann,j+num_wann,iktot,it)=ham_auxt(i,j,iktot,it)
                  !
               enddo
            enddo
         enddo
      enddo
   enddo
   !
   deallocate(ndegen,irvec,ham_rt)
   ham_kt=ham_auxt
   write(*,'(1A)')"  K-space H(K,t) built"
   !
   !5) Reordering & hermicity check
   if(Nspin==2)then
      ham_auxt=zero;ham_auxt=ham_kt;ham_kt=zero
      do iktot=1,num_kpts
         do it=1,Nt
            !
            ham_kt(:,:,iktot,it)=slo2lso(ham_auxt(:,:,inrpts,it),Nlat,Nspin,Norb)
            !
         enddo
      enddo
   endif
   do iktot=1,num_kpts
      do it=1,Nt
         call herm_check(ham_kt(:,:,iktot,it))
      enddo
   enddo
   deallocate(ham_auxt)
   write(*,'(1A)')"  H(k,t) written"
   !
  end subroutine hkt_from_w90_hr



  subroutine hloct_from_w90_hr(field,gauge,R1,R2,R3,Ruc,                      &
                               put_dipole,put_local_dipole,absorbdiagonal,    &
                               Einleads,leadlimit,Hloct,w90_file,dipole_file, &
                               Nspin,Norb,Nlat,Nt)
   implicit none
   real(8)               ,intent(in)            ::   field(:,:,:) ![Nt,dim,2] 1=Efield 2=Afield
   character(len=*)      ,intent(in)            ::   gauge
   real(8)               ,intent(in)            ::   R1(:),R2(:),R3(:)
   real(8)               ,intent(in)            ::   Ruc(:,:)
   logical               ,intent(in)            ::   put_dipole
   logical               ,intent(in)            ::   put_local_dipole
   logical               ,intent(in)            ::   absorbdiagonal
   logical               ,intent(in)            ::   Einleads
   integer               ,intent(in)            ::   leadlimit
   complex(8),allocatable,intent(inout)         ::   Hloct(:,:,:,:)
   character(len=*)      ,intent(in)            ::   w90_file
   character(len=*)      ,intent(in)            ::   dipole_file
   integer               ,intent(in)            ::   Nspin,Norb,Nlat,Nt
   logical                                      ::   IOfile
   integer                                      ::   unitIO1,unitIO2
   integer                                      ::   i,j,k,iorb,jorb,io,jo,it,ilat,jlat
   integer                                      ::   ndx1_H,ndx2_H,ndx1_D,ndx2_D
   integer                                      ::   inrpts
   real(8)                                      ::   a,b,exparg,expargR,expargD,absorbswitch
   real(8)                                      ::   REDx,REDy,REDz,IMDx,IMDy,IMDz
   integer                                      ::   rst,qst,limit,kvec_ndx
   integer                                      ::   auxndx,dumR1,dumR2,dumR3
   integer                                      ::   nRvec
   !---- light matter ----
   integer   ,allocatable,dimension(:)          ::   Kvec,ilist
   integer   ,allocatable,dimension(:,:)        ::   site_ndx,orb_ndx
   integer   ,allocatable,dimension(:,:,:)      ::   veclist
   real(8)   ,allocatable,dimension(:)          ::   locswitch
   complex(8),allocatable,dimension(:,:,:)      ::   Cvec,ham_r,StructFact
   complex(8),allocatable,dimension(:,:,:,:)    ::   dip_r,ham_rt,lightmat_r
   !---- W90 specific ----
   integer                                      ::   num_wann
   integer                                      ::   nrpts
   integer(4),allocatable                       ::   ndegen(:)      !(nrpts)
   integer   ,allocatable                       ::   irvec(:,:)     !(3,nrpts)
   !
   !
   unitIO1=free_unit();open(unit=unitIO1,file=w90_file,status="old",action="read")
   unitIO2=free_unit();open(unit=unitIO2,file=dipole_file,status="old",action="read")
   read(unitIO1,*)
   read(unitIO1,*) num_wann
   read(unitIO1,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   write(*,*)
   write(*,'(1A)')         "-------------- Ht_LDA --------------"
   write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
   write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
   write(*,'(A,I6,1A,I6)') "  rows:",qst,"  last row elements:",rst
   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
   !
   nRvec=1+2+2+2+2
   !
   if(allocated(ndegen))         deallocate(ndegen)         ;allocate(ndegen(nrpts))                         ;ndegen=0
   if(allocated(irvec))          deallocate(irvec)          ;allocate(irvec(nrpts,3))                        ;irvec=0
   if(allocated(Kvec))           deallocate(Kvec)           ;allocate(Kvec(3))                               ;Kvec=0
   if(allocated(veclist))        deallocate(veclist)        ;allocate(veclist(-10:10,-10:10,-10:10))         ;veclist=0
   if(allocated(ilist))          deallocate(ilist)          ;allocate(ilist(nRvec))                          ;ilist=0
   if(allocated(site_ndx))       deallocate(site_ndx)       ;allocate(site_ndx(num_wann,2))                  ;site_ndx=0
   if(allocated(orb_ndx))        deallocate(orb_ndx)        ;allocate(orb_ndx(num_wann,2))                   ;orb_ndx=0
   if(allocated(ham_r))          deallocate(ham_r)          ;allocate(ham_r(num_wann,num_wann,nrpts))        ;ham_r=zero
   if(allocated(Cvec))           deallocate(Cvec)           ;allocate(Cvec(Nlat,Norb,3))                     ;Cvec=zero
   if(allocated(dip_r))          deallocate(dip_r)          ;allocate(dip_r(num_wann,num_wann,nrpts,3))      ;dip_r=zero
   if(allocated(StructFact))     deallocate(StructFact)     ;allocate(StructFact(Nlat,Nlat,Nt))              ;StructFact=zero
   if(allocated(lightmat_r))     deallocate(lightmat_r)     ;allocate(lightmat_r(num_wann,num_wann,nRvec,3)) ;lightmat_r=zero
   !
   if(allocated(ham_rt))    deallocate(ham_rt)    ;allocate(ham_rt(num_wann,num_wann,nRvec,Nt))              ;ham_rt=zero
   !
   !1) read WS degeneracies
   do i=1,qst
      read(unitIO1,*)(ndegen(j+(i-1)*15),j=1,15)
   enddo
   if(rst.ne.0)read(unitIO1,*)(ndegen(j+qst*15),j=1,rst)
   write(*,'(1A)')"  degen readed"
   !
   !2) read real-space quantities
   limit=0
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            !
            !read H(R) & D(R)
            read(unitIO1,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H,a,b
            read(unitIO2,*)dumR1,dumR2,dumR3,ndx1_D,ndx2_D,REDx,IMDx,REDy,IMDy,REDz,IMDz
            !
            !consistency check
            auxndx = sum([irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H]-[dumR1,dumR2,dumR3,ndx1_D,ndx2_D])
            if(auxndx.ne.0)then
               write(*,'(10A)') "  Something is wrong between ",w90_file," and ",dipole_file," indexing"
               write(*,'(10I5)')irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H
               write(*,'(10I5)')dumR1,dumR2,dumR3,ndx1_D,ndx2_D
               stop
            endif
            !
            if(abs(dumR1).gt.limit)limit=abs(dumR1)
            veclist(irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3))=inrpts
            !
            site_ndx(ndx1_H,1)=floor((ndx1_H-0.01)/Norb)+1
            site_ndx(ndx2_H,2)=floor((ndx2_H-0.01)/Norb)+1
            !
            orb_ndx(ndx1_H,1)=ndx1_H-3*(site_ndx(ndx1_H,1)-1)
            orb_ndx(ndx2_H,2)=ndx2_H-3*(site_ndx(ndx2_H,1)-1)
            !
            ham_r(ndx1_H,ndx2_H,inrpts)=dcmplx(a,b)
            dip_r(ndx1_D,ndx2_D,inrpts,1)=dcmplx(REDx,IMDx)
            dip_r(ndx1_D,ndx2_D,inrpts,2)=dcmplx(REDy,IMDy)
            dip_r(ndx1_D,ndx2_D,inrpts,3)=dcmplx(REDz,IMDz)
            !
            if((ndx1_D.eq.ndx2_D).and.(irvec(inrpts,1).eq.0).and.(irvec(inrpts,2).eq.0).and.(irvec(inrpts,3).eq.0))then
               Cvec(site_ndx(ndx1_H,1),orb_ndx(ndx1_H,1),1)=dcmplx(REDx,IMDx)
               Cvec(site_ndx(ndx1_H,1),orb_ndx(ndx1_H,1),2)=dcmplx(REDy,IMDy)
               Cvec(site_ndx(ndx1_H,1),orb_ndx(ndx1_H,1),3)=dcmplx(REDz,IMDz)
            endif
            !
         enddo
      enddo
   enddo
   close(unitIO1)
   close(unitIO2)
   !
   ilist(1)=veclist( 0, 0, 0)
   ilist(2)=veclist(+1, 0, 0)
   ilist(3)=veclist(-1, 0, 0)
   ilist(4)=veclist( 0,+1, 0)
   ilist(5)=veclist( 0,-1, 0)
   ilist(6)=veclist( 0, 0,+1)
   ilist(7)=veclist( 0, 0,-1)
   ilist(8)=veclist(-1,+1, 0)
   ilist(9)=veclist(+1,-1, 0)
   if(allocated(locswitch))deallocate(locswitch);allocate(locswitch(nRvec));locswitch=0d0
   locswitch(1)=1.0d0
   !
   write(*,'(1A)')"  H(R) and D(R) readed"
   !
   !3) cleanup dipole components
   write(*,*)"  put_dipole:",put_dipole
   write(*,*)"  put_local_dipole:",put_local_dipole
   write(*,*)"  absorbdiagonal:",absorbdiagonal
   if(.not.put_dipole)then
      dip_r=zero
      Cvec=zero
      write(*,'(1A)')"  D(R) cancelled"
   endif
   if(.not.put_local_dipole)then
      do ilat=1,Nlat
         do iorb=1,Norb
            do jorb=1,Norb
               !
               io = iorb + (ilat-1)*Norb
               jo = jorb + (ilat-1)*Norb
               dip_r(io,jo,ilist(1),:)=zero
               !
            enddo
         enddo
      enddo
      Cvec=zero
      write(*,'(1A)')"  D(a) cancelled"
   endif
   if(absorbdiagonal)then
      do i=1,num_wann
         dip_r(i,i,ilist(1),:)=zero
      enddo
      write(*,'(1A)')"  diagonal D(0) cancelled"
   endif
   !
   !4) reduce dipole to few long range hoppings
   if(gauge=="A")then
      !
      do k=1,nRvec
         i=ilist(k)
         do j=1,nrpts
            !
            Kvec(1)=irvec(i,1)-irvec(j,1)
            Kvec(2)=irvec(i,2)-irvec(j,2)
            Kvec(3)=irvec(i,3)-irvec(j,3)
            !
            if((abs(Kvec(1)).gt.limit) .or. &
               (abs(Kvec(2)).gt.limit) .or. &
               (abs(Kvec(3)).gt.limit)      )cycle
            !
            Kvec_ndx=veclist(Kvec(1),Kvec(2),Kvec(3))
            !
            lightmat_r(:,:,k,1) = lightmat_r(:,:,k,1) + matmul(dip_r(:,:,Kvec_ndx,1),ham_r(:,:,j))
            lightmat_r(:,:,k,2) = lightmat_r(:,:,k,2) + matmul(dip_r(:,:,Kvec_ndx,2),ham_r(:,:,j))
            lightmat_r(:,:,k,3) = lightmat_r(:,:,k,3) + matmul(dip_r(:,:,Kvec_ndx,3),ham_r(:,:,j))
            !
         enddo
      enddo
      !
      do k=1,nRvec
         i=ilist(k)
         do j=1,nrpts
            !
            Kvec(1)=irvec(i,1)-irvec(j,1)
            Kvec(2)=irvec(i,2)-irvec(j,2)
            Kvec(3)=irvec(i,3)-irvec(j,3)
            !
            if((abs(Kvec(1)).gt.limit) .or. &
               (abs(Kvec(2)).gt.limit) .or. &
               (abs(Kvec(3)).gt.limit)      )cycle
            !
            Kvec_ndx=veclist(Kvec(1),Kvec(2),Kvec(3))
            !
            lightmat_r(:,:,k,1) = lightmat_r(:,:,k,1) - matmul(ham_r(:,:,Kvec_ndx),dip_r(:,:,j,1))
            lightmat_r(:,:,k,2) = lightmat_r(:,:,k,2) - matmul(ham_r(:,:,Kvec_ndx),dip_r(:,:,j,2))
            lightmat_r(:,:,k,3) = lightmat_r(:,:,k,3) - matmul(ham_r(:,:,Kvec_ndx),dip_r(:,:,j,3))
            !
         enddo
         call herm_check(lightmat_r(:,:,1,1))
         call herm_check(lightmat_r(:,:,1,2))
         call herm_check(lightmat_r(:,:,1,3))
      enddo
      !
   elseif(gauge=="E")then
      !
      do k=1,nRvec
         i=ilist(k)
         !
         lightmat_r(:,:,k,1) = dip_r(:,:,i,1) ; call herm_check(lightmat_r(:,:,1,1))
         lightmat_r(:,:,k,2) = dip_r(:,:,i,2) ; call herm_check(lightmat_r(:,:,1,2))
         lightmat_r(:,:,k,3) = dip_r(:,:,i,3) ; call herm_check(lightmat_r(:,:,1,3))
         !
      enddo
      !
   endif
   deallocate(veclist,Kvec,dip_r)
   write(*,'(2A)')"  light-matter interaction built in gauge: ",gauge
   !
   !5) build prefactor
   if(gauge=="A")then
      !
      StructFact=dcmplx(1.d0,0.d0)
      !
   elseif(gauge=="E")then
      !
      StructFact=dcmplx(1.d0,0.d0)
      do it=1,Nt
         do ilat=1,Nlat
            do jlat=1,Nlat
               exparg = ( -field(it,1,2) * ( Ruc(ilat,1) - Ruc(jlat,1) ) &
                          -field(it,2,2) * ( Ruc(ilat,2) - Ruc(jlat,2) ) &
                          -field(it,3,2) * ( Ruc(ilat,3) - Ruc(jlat,3) ) )
               !
               StructFact(ilat,jlat,it) = dcmplx(cos(exparg),sin(exparg))
               !
            enddo
         enddo
         call herm_check(StructFact(:,:,it))
      enddo
      !
   endif
   write(*,'(1A)')"  StructFact built"
   !
   !6) build interacting hamilt in real space
   if(gauge=="A")then
      !
      do it=1,Nt
         do k=1,nRvec
            inrpts=ilist(k)
            !
            ham_rt(:,:,k,it) = ( ham_r(:,:,inrpts) + field(it,1,2) * Xi * lightmat_r(:,:,k,1) &
                                                   + field(it,2,2) * Xi * lightmat_r(:,:,k,2) &
                                                   + field(it,3,2) * Xi * lightmat_r(:,:,k,3) )
            !
         enddo
      enddo
      !
   elseif(gauge=="E")then
      !
      absorbswitch=0d0
      if(absorbdiagonal)then
         absorbswitch=1.0d0
         write(*,'(1A)')"  Da(0) absorbed"
      endif
      !
      do ilat=1,Nlat
         do jlat=1,Nlat
            !
            if((.not.Einleads).and.(ilat.ge.leadlimit).and.(jlat.ge.leadlimit))then
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     io = iorb + (ilat-1)*Norb
                     jo = jorb + (jlat-1)*Norb
                     !
                     do it=1,Nt
                        do k=1,nRvec
                           inrpts=ilist(k)
                           !
                           ham_rt(io,jo,k,it)=ham_r(io,jo,inrpts)
                           !
                        enddo
                     enddo
                  enddo
               enddo
            else
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     io = iorb + (ilat-1)*Norb
                     jo = jorb + (jlat-1)*Norb
                     !
                     do it=1,Nt
                        do k=1,nRvec
                           inrpts=ilist(k)
                           !
                           expargR = ( &
                           +field(it,1,2) * ( irvec(inrpts,1)*R1(1) + irvec(inrpts,2)*R2(1) + irvec(inrpts,3)*R3(1) ) &
                           +field(it,2,2) * ( irvec(inrpts,1)*R1(2) + irvec(inrpts,2)*R2(2) + irvec(inrpts,3)*R3(2) ) &
                           +field(it,3,2) * ( irvec(inrpts,1)*R1(3) + irvec(inrpts,2)*R2(3) + irvec(inrpts,3)*R3(3) ) )
                           !
                           expargD = ( &
                           -field(it,1,2) * ( Cvec(ilat,iorb,1) - Cvec(jlat,Jorb,1)*locswitch(i) ) &
                           -field(it,2,2) * ( Cvec(ilat,iorb,2) - Cvec(jlat,Jorb,2)*locswitch(i) ) &
                           -field(it,3,2) * ( Cvec(ilat,iorb,3) - Cvec(jlat,Jorb,3)*locswitch(i) ) ) * absorbswitch
                           !
                           ham_rt(io,jo,k,it) =                                                                            &
                           dcmplx(cos(expargR),sin(expargR))*dcmplx(cos(expargD),sin(expargD))*StructFact(ilat,jlat,it) *  &
                          ( ham_r(io,jo,inrpts) + field(it,1,1) * lightmat_r(io,jo,k,1)                                    &
                                                + field(it,2,1) * lightmat_r(io,jo,k,2)                                    &
                                                + field(it,3,1) * lightmat_r(io,jo,k,3) )
                           !
                        enddo
                     enddo
                  enddo
               enddo
            endif
            !
         enddo
      enddo
      !
   endif
   deallocate(StructFact,site_ndx,ham_r,lightmat_r,locswitch)
   write(*,'(1A)')"  real-space H(R,t) built"
   !
   !6) Reordering & hermicity check
   Hloct=zero
   do k=1,nRvec
      do it=1,Nt
         Hloct(:,:,k,it)=ham_rt(:,:,k,it)
         if(Nspin==2)then
            Hloct(1+num_wann:2*num_wann,1+num_wann:2*num_wann,k,it)=ham_rt(:,:,k,it)
            Hloct(:,:,k,it)=slo2lso(Hloct(:,:,k,it),Nlat,Nspin,Norb)
         endif
      enddo
   enddo
   deallocate(ham_rt)
   write(*,'(1A)')"  Hloc(K,t) passed to main"
   !
  end subroutine hloct_from_w90_hr



! subroutine read_Hr_w90_solve_Hk_along_BZpath(       w90_file           &      !output file of w90
!                                                 ,   Nspin,Norb,Nlat    &      !dimensions
!                                                 ,   kpath,Nk           &      !name of Kpoint on path and mesh between
!                                                 ,   colors_name        &      !Colors of the Norb orbitals
!                                                 ,   points_name        &      !Name of the point on the path
!                                                 ,   file_eigenband     &      !Name of the file where to save the bands                    (optional)
!                                                 ,   Hkpathfile         &      !Name of the file where to save the H(k) on the path         (optional)
!                                                 ,   Kpointpathfile     &      !Name of the file where to save the Kpoints on the path      (optional)
!                                                 ,   S_correction_      &      !Re{impSigma(w=0)}                                           (optional)
!                                                 ,   ham_k              &      !k-space hamiltonian on path                                 (optional)
!                                                 ,   kpt_latt           )      !Kpoint path                                                 (optional)

!   implicit none
!   character(len=*)      ,intent(in)            ::   w90_file
!   integer               ,intent(in)            ::   Nspin,Norb,Nlat
!   real(8)   ,allocatable,intent(in)            ::   kpath(:,:)
!   integer               ,intent(in)            ::   Nk
!   type(rgb_color)       ,intent(in)            ::   colors_name(Norb*Nlat)
!   character(len=*),dimension(size(kpath,1))    ::   points_name
!   character(len=*)      ,intent(in) ,optional  ::   file_eigenband
!   character(len=*)      ,intent(in) ,optional  ::   Hkpathfile
!   character(len=*)      ,intent(in) ,optional  ::   Kpointpathfile
!   real(8)   ,allocatable,intent(in) ,optional  ::   S_correction_(:,:) ![Nspin*Norb*Nlat,Nspin*Norb*Nlat]
!   complex(8),allocatable,intent(out),optional  ::   ham_k(:,:,:)   !(num_wann*nspin,num_wann*nspin,num_kpts)
!   real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)
!   integer                                      ::   i,j,ndx1,ndx2
!   integer                                      ::   ipts,ik,ic,u1,u2
!   integer                                      ::   iktot,num_kpts
!   integer                                      ::   inrpts
!   character(len=256)                           ::   file_,file_s1,file_s2,file_corr_s1,file_corr_s2,xtics
!   integer                                      ::   Npts
!   integer                                      ::   unit
!   real(8),dimension(size(kpath,2))             ::   kstart,kstop,kdiff,kpoint
!   real(8)                                      ::   U(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
!   real(8)                                      ::   Udag(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
!   type(rgb_color)                              ::   corb(Norb*Nlat),c(Norb*Nlat)
!   character(len=10)                            ::   chpoint
!   character(len=32) :: fmt
!   real(8)                                      ::   a,b,rdotk
!   integer                                      ::   rst,qst
!   !---- W90 specific ----
!   integer                                      ::   num_wann       !=Norb*Nlat
!   integer                                      ::   nrpts
!   integer(4),allocatable                       ::   ndegen(:)      !(nrpts)
!   integer   ,allocatable                       ::   irvec(:,:)     !(3,nrpts)
!   complex(8),allocatable                       ::   ham_r(:,:,:)
!   complex(8),allocatable                       ::   ham_aux(:,:,:)
!   !
!   Npts=size(kpath,1)
!   num_kpts=Nk*(Npts-1)
!   !
!   write(*,'(1A)')"  Solving model along the path:"
!   write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
!   do ipts=1,Npts
!      write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
!   enddo
!   !
!   open(unit=106,file=w90_file,status="unknown",action="read")
!   read(106,*)
!   read(106,*) num_wann
!   read(106,*) nrpts
!   rst=mod(nrpts,15)
!   qst=int(nrpts/15)
!   write(*,'(1A)')         "-------------- H_LDA --------------"
!   write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
!   write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
!   write(*,'(A,I6,A,I6)') "  rows:",qst,"  last row elements:",rst
!   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
!   !
!   if(allocated(kpt_latt))deallocate(kpt_latt);allocate(kpt_latt(num_kpts,3))                           ;kpt_latt=0d0
!   if(allocated(ndegen))  deallocate(ndegen)  ;allocate(ndegen(nrpts))                                  ;ndegen=0
!   if(allocated(irvec))   deallocate(irvec)   ;allocate(irvec(nrpts,3))                                 ;irvec=0
!   if(allocated(ham_r))   deallocate(ham_r)   ;allocate(ham_r(num_wann*Nspin,num_wann*Nspin,nrpts))     ;ham_r=zero
!   if(allocated(ham_k))   deallocate(ham_k)   ;allocate(ham_k(num_wann*Nspin,num_wann*Nspin,num_kpts))  ;ham_k =zero
!   if(allocated(ham_aux)) deallocate(ham_aux) ;allocate(ham_aux(num_wann*Nspin,num_wann*Nspin,num_kpts));ham_aux =zero
!   !
!   !1) k-points mesh on path
!   ic=0
!   do ipts=1,Npts-1
!      kstart = kpath(ipts,:)
!      kstop  = kpath(ipts+1,:)
!      kdiff  = (kstop-kstart)/dble(Nk)
!      do ik=1,Nk
!         ic=ic+1
!         kpoint = kstart + (ik-1)*kdiff
!         kpt_latt(ic,:)=kpoint
!      enddo
!   enddo
!   !
!   !2) read WS degeneracies
!   do i=1,qst
!      read(106,*)(ndegen(j+(i-1)*15),j=1,15)
!   enddo
!   if(rst.ne.0)read(106,*)(ndegen(j+qst*15),j=1,rst)
!   write(*,'(1A)')"  degen readed"
!   !
!   !3) read real-space Hamiltonian (no spinup-spindw hybridizations assumed change sub otherwise)
!   do inrpts=1,nrpts
!      do i=1,num_wann
!         do j=1,num_wann
!            read(106,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1,ndx2,a,b
!            !spin up
!            ham_r(ndx1,ndx2,inrpts)=dcmplx(a,b)
!            !spin dw
!            if(Nspin==2)ham_r(ndx1+num_wann,ndx2+num_wann,inrpts)=dcmplx(a,b)
!         enddo
!      enddo
!   enddo
!   close(106)
!   write(*,'(2A)')"  H(R) readed from: ",w90_file
!   !
!   !4) Fourier Transform on path
!   do iktot=1,num_kpts
!      do i=1,num_wann*nspin
!         do j=1,num_wann*nspin
!            do inrpts=1,nrpts
!               rdotk=0.d0
!               rdotk= ( kpt_latt(iktot,1)*irvec(inrpts,1) + &
!                        kpt_latt(iktot,2)*irvec(inrpts,2) + &
!                        kpt_latt(iktot,3)*irvec(inrpts,3) )
!               !
!               ham_k(i,j,iktot)=ham_k(i,j,iktot)+ham_r(i,j,inrpts)*dcmplx(cos(rdotk),-sin(rdotk))/ndegen(inrpts)
!               !
!            enddo
!         enddo
!      enddo
!   enddo
!   write(*,'(1A)')"  H(k) on path produced"
!   !
!   !5) re-ordering from w90 user defined to the code standard [[[Norb],Nspin],Nlat]
!   if (Nspin==2)then
!   ham_aux=zero;ham_aux=ham_k;ham_k=zero
!      do iktot=1,num_kpts
!         ham_k(:,:,iktot)=slo2lso(ham_aux(:,:,iktot),Nlat,Nspin,Norb)
!      enddo
!   endif
!   !
!   if(present(Hkpathfile))then
!      call TB_write_hk(ham_k,Hkpathfile,Nspin*Norb*Nlat,1,1,Nlat,Nk,kpath)
!      write(*,'(2A)')"  H(k) on path written on: ",Hkpathfile
!   endif
!   !
!   if(present(Kpointpathfile))then
!      open(unit=107,file=Kpointpathfile,status="unknown",action="write",position="rewind")
!      do iktot=1,num_kpts
!         write(107,'(3F15.7)') (kpt_latt(iktot,i),i=1,3)
!      enddo
!      close(107)
!      write(*,'(2A)')"  Kpoints on path used written on: ",Kpointpathfile
!   endif
!   !
!   !6) Coloured Eigenbands on path for the two different spins
!   do i=1,num_wann
!      corb(i) = colors_name(i)
!   enddo
!   !
!   file_="Eigenbands"
!   if(present(file_eigenband))file_=file_eigenband
!   !
!   file_s1=reg(file_)//"_s1"
!   call solve_bands(reg(file_s1),1)
!   !
!   if(Nspin.gt.1)then
!      !
!      file_s2=reg(file_)//"_s2"
!      call solve_bands(reg(file_s2),2)
!      !
!      if(present(S_correction_))then
!         !
!         file_corr_s1=reg(file_)//"_Sreal_s1"
!         call solve_bands(reg(file_corr_s1),1,S_correction_)
!         !
!         file_corr_s2=reg(file_)//"_Sreal_s2"
!         call solve_bands(reg(file_corr_s2),2,S_correction_)
!         !
!      endif
!   else
!      if(present(S_correction_))then
!         !
!         file_corr_s1=reg(file_)//"_Sreal_s1"
!         call solve_bands(reg(file_corr_s1),1,S_correction_)
!         !
!      endif
!   endif



!   contains

!   subroutine solve_bands(file_print_,spindx,S)
!     implicit none
!     character(len=*),intent(in)                  ::   file_print_
!     integer         ,intent(in)                  ::   spindx
!     real(8)   ,allocatable,intent(in) ,optional  ::   S(:,:)
!     character(len=256)                           ::   file_print
!     integer                                      ::   i
!     real(8)                                      ::   S_correction(Norb*Nspin*Nlat,Norb*Nspin*Nlat)
!     real(8)                                      ::   eval(Norb*Nlat),coeff(Norb*Nlat)
!     complex(8)                                   ::   h(Norb*Nspin*Nlat,Norb*Nspin*Nlat)
!     complex(8) ,dimension(Norb*Nlat,Norb*Nlat)   ::   h1,h2,hdiag
!     !
!     file_print=file_print_//".dat"
!     write(*,'(2A)')"  Printing eigenbands on: ",reg(file_print)
!     unit=free_unit()
!     open(unit,file=reg(file_print))
!     !
!     S_correction=0d0
!     if(present(S))S_correction=S
!     !
!     do iktot=1,num_kpts
!        h=zero
!        h=ham_k(:,:,iktot)+S_correction
!        !reshape with spin external block (usual for w90) so as to diagonalize a block matrix
!        h=matmul(U,matmul(h,Udag))
!        !identify different spins
!        h1(:,:)=h(1:num_wann,1:num_wann)
!        h2(:,:)=h(1+num_wann:Nspin*num_wann,1+num_wann:Nspin*num_wann)
!        if(spindx==1)hdiag=h1
!        if(spindx==2)hdiag=h2
!        call eigh(hdiag,Eval)
!        do i=1,num_wann
!           coeff(:)=hdiag(:,i)*conjg(hdiag(:,i))
!           c(i) = coeff.dot.corb
!        enddo
!        write(unit,'(I12,100(F18.12,I18))')iktot,(Eval(i),rgb(c(i)),i=1,num_wann)
!     enddo
!     close(unit)
!     xtics="'"//reg(points_name(1))//"'1,"
!     do ipts=2,Npts-1
!        xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//reg(txtfy((ipts-1)*Nk+1))//","
!     enddo
!     xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//reg(txtfy((Npts-1)*Nk))//""
!        open(unit,file=reg(file_print)//".gp")
!        write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
!        write(unit,*)"#set out '"//reg(file_print)//".png'"
!        write(unit,*)""
!        write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
!        write(unit,*)"#set out '"//reg(file_print)//".svg'"
!        write(unit,*)""
!        write(unit,*)"#set term postscript eps enhanced color 'Times'"
!        write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_print)//".pdf'"
!        write(unit,*)"unset key"
!        write(unit,*)"set xtics ("//reg(xtics)//")"
!        write(unit,*)"set grid noytics xtics"
!        !
!        write(unit,*)"plot '"//reg(file_print)//"' u 1:2:3 w l lw 3 lc rgb variable,\"
!        do i=2,num_wann-1
!           u1=2+(i-1)*2
!           u2=3+(i-1)*2
!           write(unit,*)"'"//reg(file_print)//"' u 1:"//reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable,\"
!        enddo
!        u1=2+(num_wann-1)*2
!        u2=3+(num_wann-1)*2
!        write(unit,*)"'"//reg(file_print)//"' u 1:"//reg(txtfy(u1))//":"//reg(txtfy(u2))//" w l lw 3 lc rgb variable"
!        !
!        close(unit)
!     call system("chmod +x "//reg(file_print)//".gp")
!   end subroutine solve_bands
!   !
!   !
! end subroutine read_Hr_w90_solve_Hk_along_BZpath
