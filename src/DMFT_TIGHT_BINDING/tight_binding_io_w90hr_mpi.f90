  subroutine hk_from_w90_hr_mpi(MpiComm,Rx,Ry,Rz,ham_k,w90_file,Nspin,Norb,Nlat,Nkvec,Porder,kpt_latt,Hkfile,Kpointfile)
   implicit none
   integer               ,intent(in)            ::   MpiComm
   real(8)               ,intent(in)            ::   Rx(:),Ry(:),Rz(:)
   complex(8),allocatable,intent(out)           ::   ham_k(:,:,:)         !(num_wann*nspin,num_wann*nspin,num_kpts)
   character(len=*)      ,intent(in)            ::   w90_file             !"seedname_hr.dat"
   integer               ,intent(in)            ::   Nspin,Norb,Nlat
   integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
   integer               ,intent(in) ,optional  ::   Porder(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
   real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)        ![ik,3]
   character(len=*)      ,intent(in) ,optional  ::   Hkfile
   character(len=*)      ,intent(in) ,optional  ::   Kpointfile
   !
   integer                                      ::   mpi_ierr
   integer                                      ::   mpi_rank
   integer                                      ::   mpi_size
   logical                                      ::   mpi_master
   integer                                      ::   i,j,ndx1,ndx2
   integer                                      ::   ispin,jspin,iorb,jorb,ilat,iktot
   integer                                      ::   Nkx,Nky,Nkz,num_kpts
   integer                                      ::   inrpts
   integer                                      ::   P(Nspin*Norb*Nlat,Nspin*Norb*Nlat)
   real(8)                                      ::   bk_x(size(Rx)),bk_y(size(Ry)),bk_z(size(Rz))
   real(8),dimension(product(Nkvec),size(Nkvec))::   kpt_x,kpt_y,kpt_z
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
   real(8)                                      ::   a,b,rdotk
   integer                                      ::   rst,qst
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
   open(unit=106,file=w90_file,status="unknown",action="read")
   read(106,*)
   read(106,*) num_wann
   read(106,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   if(mpi_master)then
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
   if(mpi_master)write(*,'(A)')"  degen readed"
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
   if(mpi_master)write(*,'(2A)')"  H(R) readed from: ",w90_file
   !
   !4) Fourier Transform
   do iktot = 1+mpi_rank, num_kpts, mpi_size
      do inrpts=1,nrpts
         rdotk=0.d0
         rdotk= ( irvec(inrpts,1)*dot_product(kpt_x(iktot,:),Rx) +  &
                  irvec(inrpts,2)*dot_product(kpt_y(iktot,:),Ry) +  &
                  irvec(inrpts,3)*dot_product(kpt_z(iktot,:),Rz) )
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
   !5) Reordering
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
   if(mpi_master)then
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
   endif
   !
   deallocate(ndegen)
   deallocate(irvec)
   deallocate(ham_r)
   deallocate(ham_aux)
   deallocate(Hloc)
   if(.not.present(kpt_latt))deallocate(kpt_latt)
   !
  end subroutine hk_from_w90_hr_mpi

