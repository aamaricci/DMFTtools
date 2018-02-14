  subroutine hk_from_w90_hr(ham_k,w90_file,Nspin,Norb,Nlat,Nkvec,kpt_latt,Hkfile,Kpointfile)
   implicit none
   complex(8),allocatable,intent(out)           ::   ham_k(:,:,:)         !(num_wann*nspin,num_wann*nspin,num_kpts)
   character(len=*)      ,intent(in)            ::   w90_file             !"LVO_hr.dat"
   integer            ,intent(in)               ::   Nspin,Norb,Nlat
   integer(4),allocatable,intent(in)            ::   Nkvec(:)             ![Nkx,Nky,Nkz]
   real(8)   ,allocatable,intent(out),optional  ::   kpt_latt(:,:)        ![ik,3]
   character(len=*)      ,intent(in) ,optional  ::   Hkfile
   character(len=*)      ,intent(in) ,optional  ::   Kpointfile
   !
   integer                                      ::   i,j,ndx1,ndx2
   integer                                      ::   ispin,jspin,iorb,jorb,ilat
   integer                                      ::   ikx,iky,ikz,iktot
   integer                                      ::   Nkx,Nky,Nkz,num_kpts
   integer                                      ::   inrpts,ndim
   real(8)   ,allocatable                       ::   bkx(:),bky(:),bkz(:)
   complex(8),parameter                         ::   zero=dcmplx(0.d0,0.d0)
   complex(8),parameter                         ::   xi  =dcmplx(0.d0,1.d0)
   real(8)   ,parameter                         ::   pi  =3.14159265358979323846
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
   !local Hamiltonian
   complex(8),allocatable                       ::   Hloc(:,:)      !(num_wann*nspin,num_wann*nspin)
   !dummy vars
   real(8)                                      ::   a,b,factor_hr,rdotk
   integer                                      ::   rst,qst
   factor_hr=6.28318530717959
   !
   !
   ndim=3!shape(Nkvec)
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
   allocate(kpt_latt(num_kpts,3))                          ;kpt_latt=0d0
   allocate(ndegen(nrpts))                                 ;ndegen=0
   allocate(irvec(nrpts,3))                                ;irvec=0
   allocate(ham_r(num_wann*Nspin,num_wann*Nspin,nrpts))    ;ham_r=zero
   allocate(ham_k(num_wann*Nspin,num_wann*Nspin,num_kpts)) ;ham_k =zero
   allocate(Hloc(num_wann*Nspin,num_wann*Nspin))           ;Hloc =zero
   allocate(bkx(3))                                        ;bkx=0d0
   allocate(bky(3))                                        ;bky=0d0
   allocate(bkz(3))                                        ;bkz=0d0
   !
   !1) k-points mesh
   call TB_build_kgrid(Nkvec,kpt_latt,.true.)
   call TB_get_bk(bkx,bky,bkz)
   !
   !2) read WS degeneracies
   do i=1,qst
      !write(*,*) "read line:", i
      read(106,*)(ndegen(j+(i-1)*15),j=1,15)
   enddo
   read(106,*)(ndegen(j+qst*15),j=1,rst)
   write(*,*)"  degen readed"
   !
   !3) read real-space Hamiltonian
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            read(106,*)irvec(inrpts,1),irvec(inrpts,2),irvec(inrpts,3),ndx1,ndx2,a,b
            !spin up
            ham_r(ndx1,ndx2,inrpts)=dcmplx(a,b)
            !spin dw
            ham_r(ndx1+num_wann,ndx2+num_wann,inrpts)=dcmplx(a,b)
         enddo
      enddo
   enddo
   close(106)
   write(*,*)"  H(R) readed from: ",w90_file
   !
   !4) Fourier Transform
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
   write(*,*)"  H(k) produced"
   if(present(Hkfile))then
      call TB_write_hk(ham_k,Hkfile,Nspin*Norb*Nlat,1,1,Nlat,[Nkx,Nky,Nkz])
      write(*,*)"  H(k) written on: ",Hkfile
   endif
   !
   if(present(Kpointfile))then
      open(unit=107,file=Kpointfile,status="unknown",action="write",position="rewind")
      do iktot=1,num_kpts
         write(107,'(3F15.7)') (kpt_latt(iktot,i),i=1,3)
      enddo
      close(107)
      write(*,*)"  Kpoints used written on: ",Hkfile
   endif
   !
!   Hloc=sum(ham_k,dim=3)/num_kpts
!   write(*,*)"  Hloc produced"
!   call TB_write_Hloc(Hloc,"Hloc.dat")
!   write(*,*)"  Hloc written on","Hlocfile.dat"
   !
end subroutine hk_from_w90_hr
