  subroutine dipole_t2g_LDA(Rx,Ry,Rz,Ruc,w90_file,Norb,Nlat,outputfile,Nr,varfact,t_thresh_)
   implicit none
   real(8)               ,intent(in)            ::   Rx(:),Ry(:),Rz(:)
   real(8)               ,intent(in)            ::   Ruc(:,:)
   character(len=*)      ,intent(in)            ::   w90_file
   integer               ,intent(in)            ::   Norb,Nlat
   character(len=*)      ,intent(in)            ::   outputfile
   integer               ,intent(in)            ::   Nr
   real(8)               ,intent(in)            ::   varfact
   real(8)               ,intent(in),optional   ::   t_thresh_
   !
   !
   integer                                      ::   i,j,rst,qst
   integer                                      ::   ix,iy,iz
   integer                                      ::   ndx,vecndx
   integer                                      ::   totdim,totIntegrals
   integer                                      ::   bra_orb_ndx,ket_orb_ndx
   integer                                      ::   Ruc_bra_ndx,Ruc_ket_ndx
   real(8)                                      ::   x,y,z,dr
   real(8)                                      ::   R_thresh,R_single
   real(8)                                      ::   dist,t_thresh
   real(8)                                      ::   wfc_bra,wfc_ket
   real(8)                                      ::   Dx,Dy,Dz
   real(8),allocatable                          ::   Rbra(:,:),Rket(:,:)
   real(8),allocatable                          ::   dipole(:,:)
   integer,allocatable                          ::   Integral_ndx(:)
   integer,allocatable                          ::   orbital_ndx(:,:)
   integer,allocatable                          ::   site_ndx(:,:)
   !---- W90 specific ----
   integer                                      ::   num_wann
   integer                                      ::   nrpts,inrpts
   real(8)                                      ::   a,b
   integer                                      ::   n1,n2,n3
   integer                                      ::   ndx1,ndx2
   !
   !
   t_thresh=1e-3
   if(present(t_thresh_))t_thresh=t_thresh_
   !
   !
   open(unit=106,file=w90_file,status="unknown",action="read")
   read(106,*)
   read(106,*) num_wann
   read(106,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   totdim = nrpts*num_wann*num_wann
   !
   write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
   write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
   write(*,'(A,I6,A,I6)') "  rows:",qst,"  last row elements:",rst
   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
   if(size(Ruc,dim=1).ne.Nlat)stop "mismatch between number of lattice sites and their position"
   !
   !
   if(allocated(Rbra))         deallocate(Rbra)        ;allocate(Rbra(totdim,3))        ;Rbra =0d0
   if(allocated(Rket))         deallocate(Rket)        ;allocate(Rket(totdim,3))        ;Rket =0d0
   if(allocated(dipole))       deallocate(dipole)      ;allocate(dipole(totdim,3))      ;dipole=zero
   if(allocated(Integral_ndx)) deallocate(Integral_ndx);allocate(Integral_ndx(totdim))  ;Integral_ndx=0
   if(allocated(site_ndx))     deallocate(site_ndx)    ;allocate(site_ndx(totdim,2))    ;site_ndx=0
   if(allocated(orbital_ndx))  deallocate(orbital_ndx) ;allocate(orbital_ndx(totdim,2)) ;orbital_ndx=0
   !
   !
   !1) read WS degeneracies spaces
   do i=1,qst
      read(106,*)
   enddo
   read(106,*)
   write(*,'(A)')"  degen skipped"
   !
   !
   R_single = norm2(Rx+Ry+Rz)
   R_thresh = R_single
   !
   !
   !2) read real-space vectors
   vecndx=0
   totIntegrals=0
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            vecndx=vecndx+1
            !
            read(106,*) n1,n2,n3,ndx1,ndx2,a,b
            !
            Ruc_bra_ndx=floor((ndx1-0.01)/Norb)+1
            Ruc_ket_ndx=floor((ndx2-0.01)/Norb)+1
            !
            Rbra(vecndx,:) = Ruc(Ruc_bra_ndx,:) + n1*Rx + n2*Ry + n3*Rz
            Rket(vecndx,:) = Ruc(Ruc_ket_ndx,:)
            !
            !write(997,*)Rbra(vecndx,1),Rbra(vecndx,2),Rbra(vecndx,3)
            !
            if(abs(dcmplx(a,b)).gt.t_thresh .and. (n1+n2+n3).ne.0)then
               !
               dist = norm2(Rbra(vecndx,:))
               if(dist.ge.R_thresh)R_thresh=dist
               !
               totIntegrals=totIntegrals+1
               !
               Integral_ndx(totIntegrals)=vecndx
               !
               site_ndx(totIntegrals,1)=Ruc_bra_ndx
               site_ndx(totIntegrals,2)=Ruc_ket_ndx
               !
               orbital_ndx(totIntegrals,1)=ndx1-floor((ndx1-0.01)/Norb)*Norb
               orbital_ndx(totIntegrals,2)=ndx2-floor((ndx2-0.01)/Norb)*Norb
               !
            endif
            !
         enddo
      enddo
   enddo
   close(106)
   write(*,'(2A)')      "  positions readed from: ",w90_file
   write(*,'(A,I8,A)')  "  computation of: ",totIntegrals," integrals"
   write(*,'(A,F12.6)') "  integration radius (Bohr): ",R_thresh
   write(*,'(A,F12.6)') "  ratio with u.c. radius: ",R_thresh/R_single
   !
   !
   !3) compute integrals
   dr = 2.d0*R_thresh/Nr
   do ndx = 1, totIntegrals
      !
      Dx = 0d0
      Dy = 0d0
      Dz = 0d0
      !
      vecndx = Integral_ndx(ndx)
      !
      !write(998,*)Rbra(vecndx,1),Rbra(vecndx,2),Rbra(vecndx,3)
      !
      Ruc_bra_ndx=site_ndx(ndx,1)
      Ruc_ket_ndx=site_ndx(ndx,2)
      !
      bra_orb_ndx=orbital_ndx(ndx,1)
      ket_orb_ndx=orbital_ndx(ndx,2)
      !
      do iz=1,Nr+1
         z = -R_thresh+dble(iz-1)*dr
         do iy=1,Nr+1
            y = -R_thresh+dble(iy-1)*dr
            do ix=1,Nr+1
               x = -R_thresh+dble(ix-1)*dr
               !
               wfc_bra = atomic_wfc(bra_orb_ndx,(x-Rbra(vecndx,1)),(y-Rbra(vecndx,2)),(z-Rbra(vecndx,3)),2.d0*varfact*R_single)
               wfc_ket = atomic_wfc(ket_orb_ndx,(x-Rket(vecndx,1)),(y-Rket(vecndx,2)),(z-Rket(vecndx,3)),2.d0*varfact*R_single)
               !
               Dx = Dx +  wfc_bra * ( x - Ruc(Ruc_bra_ndx,1) ) * wfc_ket * dr**3
               Dy = Dy +  wfc_bra * ( y - Ruc(Ruc_bra_ndx,2) ) * wfc_ket * dr**3
               Dz = Dz +  wfc_bra * ( z - Ruc(Ruc_bra_ndx,3) ) * wfc_ket * dr**3
               !
               !if(ndx==1)write(999,*)x,y,z
               !
            enddo
         enddo
      enddo
      !
      dipole(vecndx,1) = Dx
      dipole(vecndx,2) = Dy
      dipole(vecndx,3) = Dz
      !
   enddo
   !
   !
   !4) write to file with same indexing of W90
   open(unit=106,file=w90_file  ,status="unknown",action="read")
   open(unit=107,file=outputfile,status="unknown",action="write")
   do i=1,qst+4
      read(106,*)
   enddo
   do vecndx=1,totdim
      read (106,*) n1,n2,n3,ndx1,ndx2
      write(107,'(5I5,3F12.6)') n1,n2,n3,ndx1,ndx2,dipole(vecndx,:)
   enddo
   !
   deallocate(Rbra,Rket,dipole,Integral_ndx,site_ndx,orbital_ndx)
   !
   contains


   pure function atomic_wfc(orb,x,y,z,var) result(phi)
   integer,intent(in)     :: orb
   real(8),intent(in)     :: x,y,z,var
   real(8)                :: phi
   real(8)                :: r,rsq,radial
   !
   rsq = x**2 + y**2 + z**2
   r   = sqrt(rsq)
   radial = exp((-(r/var)**2)/3.)/(2*pi*3*var) 
   !
   ! xz
   if(orb==1) phi = sqrt(15/(4*pi))*(x*z/rsq)*radial
   ! yz
   if(orb==2) phi = sqrt(15/(4*pi))*(y*z/rsq)*radial
   ! xy
   if(orb==3) phi = sqrt(15/(4*pi))*(x*y/rsq)*radial
   !
   end function atomic_wfc

end subroutine dipole_t2g_LDA



!   pure function orb_index(ndx,ilat,Norb) result(orb)
!   integer,intent(in)     :: ndx,ilat,Norb
!   integer                :: orb
!   integer                :: app
!   !
!   app = ndx-floor((ndx-0.01)/Norb)*Norb
!   !
!   orb = app
!   if(app.le.2)  orb = app + (abs(app-Norb)-app)*mod(ilat+1,2)
!   !
!   end function orb_index



