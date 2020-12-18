subroutine hk_from_w90_hr_mpi(MpiComm,R1,R2,R3,ham_k,ham_loc,w90_file,Nspin,Norb,Nlat,Nkvec,kpt_latt,Hkfile,Kpointfile)
  implicit none
  integer               ,intent(in)            ::   MpiComm
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
  integer                                      ::   mpi_ierr
  integer                                      ::   mpi_rank
  integer                                      ::   mpi_size
  logical                                      ::   mpi_master
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
  !MPI setup:
  mpi_size  = MPI_Get_size(MpiComm)
  mpi_rank =  MPI_Get_rank(MpiComm)
  mpi_master= MPI_Get_master(MpiComm)
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
  if(mpi_master)then
     write(*,*)
     write(*,'(1A)')         "-------------- H_LDA --------------"
     write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
     write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
     write(*,'(A,I6,A,I6)') "  rows:",qst,"  last row elements:",rst
  endif
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
  call build_kgrid_generic(Nkvec,kpt1,kpt2,kpt3)
  if(present(kpt_latt))kpt_latt=kpt1+kpt2+kpt3
  !
  !2) read WS degeneracies
  do i=1,qst
     read(unitIO,*)(ndegen(j+(i-1)*15),j=1,15)
  enddo
  if(rst.ne.0)read(unitIO,*)(ndegen(j+qst*15),j=1,rst)
  if(mpi_master)write(*,'(1A)')"  degen readed"
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
  if(mpi_master)write(*,'(2A)')"  H(R) readed from: ",w90_file
  !
  !4) Fourier Transform
  do iktot = 1+mpi_rank, num_kpts, mpi_size
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
  call Mpi_AllReduce(ham_aux,ham_k, size(ham_k), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
  call MPI_Barrier(MpiComm,mpi_ierr)
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
  if(mpi_master)then
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
  endif
  !
  deallocate(ndegen)
  deallocate(irvec)
  deallocate(ham_r)
  deallocate(ham_aux)
  if(.not.present(kpt_latt))deallocate(kpt_latt)
  call MPI_Barrier(MpiComm,mpi_ierr)
  !
end subroutine hk_from_w90_hr_mpi






  subroutine hkt_from_w90_hr_mpi(MpiComm,field,gauge,R1,R2,R3,Ruc,dipole_flag,ham_kt,w90_file,dipole_file,Nspin,Norb,Nlat,Nt,Nkvec)
   implicit none
   integer               ,intent(in)            ::   MpiComm
   real(8)               ,intent(in)            ::   field(:,:,:) ![Nt,dim,2] 1=Efield 2=Afield
   character(len=*)      ,intent(in)            ::   gauge
   real(8)               ,intent(in)            ::   R1(:),R2(:),R3(:)
   real(8)               ,intent(in)            ::   Ruc(:,:)
   logical               ,intent(in)            ::   dipole_flag
   complex(8),allocatable,intent(inout)         ::   ham_kt(:,:,:,:)
   character(len=*)      ,intent(in)            ::   w90_file
   character(len=*)      ,intent(in)            ::   dipole_file
   integer               ,intent(in)            ::   Nspin,Norb,Nlat,Nt
   integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
   !
   integer                                      ::   mpi_ierr
   integer                                      ::   mpi_rank
   integer                                      ::   mpi_size
   logical                                      ::   mpi_master
   !
   logical                                      ::   IOfile
   integer                                      ::   unitIO1,unitIO2
   integer                                      ::   i,j,iorb,jorb,io,jo,it
   integer                                      ::   ndx1_H,ndx2_H,ndx1_D,ndx2_D
   integer                                      ::   Nkx,Nky,Nkz
   integer                                      ::   iktot,num_kpts
   integer                                      ::   inrpts
   real(8)                                      ::   bk1(size(R1)),bk2(size(R2)),bk3(size(R3))
   real(8),dimension(product(Nkvec),size(Nkvec))::   kpt1,kpt2,kpt3
   real(8)                                      ::   a,b,exparg,rdotk
   real(8)                                      ::   REDx,REDy,REDz,IMDx,IMDy,IMDz
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
   !MPI setup:
   mpi_size  = MPI_Get_size(MpiComm)
   mpi_rank =  MPI_Get_rank(MpiComm)
   mpi_master= MPI_Get_master(MpiComm)
   !
   !
   Nkx=Nkvec(1)
   Nky=Nkvec(2)
   Nkz=Nkvec(3)
   num_kpts=Nkx*Nky*Nkz
   !
   unitIO1=free_unit()
   open(unit=unitIO1,file=w90_file,status="old",action="read")
   if(dipole_flag)then
      unitIO2=free_unit()
      open(unit=unitIO2,file=dipole_file,status="old",action="read")
   endif
   read(unitIO1,*)
   read(unitIO1,*) num_wann
   read(unitIO1,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   if(mpi_master)then
      write(*,*)
      write(*,'(1A)')         "-------------- Ht_LDA --------------"
      write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
      write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
      write(*,'(A,I6,1A,I6)') "  rows:",qst,"  last row elements:",rst
   endif
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
   call build_kgrid_generic(Nkvec,kpt1,kpt2,kpt3)
   !
   !2) read WS degeneracies
   do i=1,qst
      read(unitIO1,*)(ndegen(j+(i-1)*15),j=1,15)
   enddo
   if(rst.ne.0)read(unitIO1,*)(ndegen(j+qst*15),j=1,rst)
   if(mpi_master)write(*,'(1A)')"  degen readed"
   !
   !3) read real-space quantities
   limit=0
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            !
            !read H(R) & D(R)
            read(unitIO1,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H,a,b
            if(dipole_flag)read(unitIO2,*)dumR1,dumR2,dumR3,ndx1_D,ndx2_D,REDx,IMDx,REDy,IMDy,REDz,IMDz
            !
            !consistency check
            if(dipole_flag)then
               auxndx = sum([irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H]-[dumR1,dumR2,dumR3,ndx1_D,ndx2_D])
               if(auxndx.ne.0)then
                  write(*,'(10A)') "  Something is wrong between ",w90_file," and ",dipole_file," indexing"
                  write(*,'(10I5)')irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1_H,ndx2_H
                  write(*,'(10I5)')dumR1,dumR2,dumR3,ndx1_D,ndx2_D
                  stop
               endif
            endif
            !
            if(abs(dumR1).gt.limit)limit=abs(dumR1)
            veclist(irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3))=inrpts
            !
            site_ndx(inrpts,1)=floor((ndx1_H-0.01)/Norb)+1
            site_ndx(inrpts,2)=floor((ndx2_H-0.01)/Norb)+1
            !
            ham_r(ndx1_H,ndx2_H,inrpts)=dcmplx(a,b)
            if(dipole_flag)then
               dip_r(ndx1_D,ndx2_D,inrpts,1)=dcmplx(REDx,IMDx)
               dip_r(ndx1_D,ndx2_D,inrpts,2)=dcmplx(REDy,IMDy)
               dip_r(ndx1_D,ndx2_D,inrpts,3)=dcmplx(REDz,IMDz)
            endif
            !
         enddo
      enddo
   enddo
   close(unitIO1)
   if(dipole_flag)close(unitIO2)
   if(mpi_master)write(*,'(1A)')"  H(R) and D(R) readed"
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
   if(mpi_master)write(*,'(2A)')"  light-matter interaction built in gauge: ",gauge
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
   if(mpi_master)write(*,'(1A)')"  prefactor built"
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
   if(mpi_master)write(*,'(1A)')"  real-space H(R,t) built"
   !
   !7) Fourier Transform
   ham_kt=zero
   do iktot=1+mpi_rank,num_kpts,mpi_size
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
   call Mpi_AllReduce(ham_auxt,ham_kt,size(ham_kt),MPI_Double_Complex,MPI_Sum,MpiComm,mpi_ierr)
   call MPI_Barrier(MpiComm,mpi_ierr)
   if(mpi_master)write(*,'(1A)')"  K-space H(K,t) built"
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
   if(mpi_master)write(*,'(1A)')"  H(k,t) written"
   call MPI_Barrier(MpiComm,mpi_ierr)
   !
  end subroutine hkt_from_w90_hr_mpi
