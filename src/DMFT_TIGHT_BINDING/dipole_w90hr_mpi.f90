  subroutine dipole_t2g_LDA_mpi(MpiComm,R1,R2,R3,Ruc,w90_file,Norb,Nlat,outputfile,Nr,var,optimize,t_thresh_,cg_niter_,cg_Ftol_)
   implicit none
   integer               ,intent(in)            ::   MpiComm
   real(8)               ,intent(in)            ::   R1(:),R2(:),R3(:)
   real(8)               ,intent(in)            ::   Ruc(:,:)
   character(len=*)      ,intent(in)            ::   w90_file
   integer               ,intent(in)            ::   Norb,Nlat
   character(len=*)      ,intent(in)            ::   outputfile
   integer               ,intent(in)            ::   Nr
   real(8),allocatable   ,intent(inout)         ::   var(:)
   logical               ,intent(in)            ::   optimize
   real(8)               ,intent(in),optional   ::   t_thresh_
   integer               ,intent(in),optional   ::   cg_niter_
   real(8)               ,intent(in),optional   ::   cg_Ftol_
   !
   integer                                      ::   mpi_ierr
   integer                                      ::   mpi_rank
   integer                                      ::   mpi_size
   logical                                      ::   mpi_master
   !
   integer                                      ::   i,j,rst,qst
   integer                                      ::   ix,iy,iz
   integer                                      ::   ndx,vecndx
   integer                                      ::   totdim,totIntegrals,totIntegrals_opt
   integer                                      ::   bra_orb_ndx,ket_orb_ndx
   integer                                      ::   Ruc_bra_ndx,Ruc_ket_ndx
   real(8)                                      ::   x,y,z,dr
   real(8)                                      ::   R_thresh,R_single
   real(8)                                      ::   dist,t_thresh
   real(8)                                      ::   wfc_bra,wfc_ket
   real(8)                                      ::   Dx,Dy,Dz
   real(8),allocatable                          ::   Rbra(:,:),Rket(:,:)
   complex(8),allocatable                       ::   dipole(:,:),dipole_aux(:,:)
   integer,allocatable                          ::   Integral_ndx(:)
   integer,allocatable                          ::   orbital_ndx(:,:)
   integer,allocatable                          ::   site_ndx(:,:)
   !---- optimization ----
   logical                                      ::   IOfile
   integer                                      ::   unitIO1,unitIO2
   integer                                      ::   dim_opt,iter,cg_niter
   real(8)                                      ::   chi,cg_Ftol
   integer   ,allocatable                       ::   Integral_ndx_opt(:)
   integer   ,allocatable                       ::   orbital_ndx_opt(:,:)
   integer   ,allocatable                       ::   site_ndx_opt(:,:)
   complex(8),allocatable                       ::   tab_R(:)
   complex(8),allocatable                       ::   tab_o(:,:)
   !---- W90 specific ----
   integer                                      ::   num_wann
   integer                                      ::   nrpts,inrpts
   real(8)                                      ::   a,b
   integer                                      ::   n1,n2,n3
   integer                                      ::   ndx1,ndx2
   !
   !
   !MPI setup:
   mpi_size  = MPI_Get_size(MpiComm)
   mpi_rank =  MPI_Get_rank(MpiComm)
   mpi_master= MPI_Get_master(MpiComm)
   !
   !
   t_thresh=1e-4
   if(present(t_thresh_))t_thresh=t_thresh_
   cg_niter=800
   if(present(cg_niter_))cg_niter=cg_niter_
   cg_Ftol=1e-5
   if(present(cg_Ftol_))cg_Ftol=cg_Ftol_
   !
   !
   if(mpi_master)then
      write(*,*)
      write(*,'(1A)')         "-------------- DIPOLE --------------"
   endif
   !
   inquire(file="variance.dat",exist=IOfile)
   if(IOfile)then
      unitIO1=free_unit()
      open(unit=unitIO1,file="variance.dat",status="unknown",action="read")
      read(unitIO1,'(3F12.8)')var
      close(unitIO1)
      if(mpi_master)write(*,'(1A)') "  Replacing variance"
   endif
   !
   !
   unitIO1=free_unit()
   open(unit=unitIO1,file=w90_file,status="old",action="read")
   read(unitIO1,*)
   read(unitIO1,*) num_wann
   read(unitIO1,*) nrpts
   rst=mod(nrpts,15)
   qst=int(nrpts/15)
   totdim = nrpts*num_wann*num_wann
   dim_opt = 9*num_wann*num_wann
   !dim_opt = 2*num_wann*num_wann
   !
   if(mpi_master)then
      write(*,'(A,I6)')      "  number of Wannier functions:   ",num_wann
      write(*,'(A,I6)')      "  number of Wigner-Seitz vectors:",nrpts
      write(*,'(A,I6,1A,I6)')"  rows:",qst,"  last row elements:",rst
   endif
   if(num_wann.ne.Nlat*Norb)stop "hk_from_w90_hr. Something is wrong"
   if(size(Ruc,dim=1).ne.Nlat)stop "mismatch between number of lattice sites and their position"
   !
   !
   if(allocated(Rbra))             deallocate(Rbra)            ;allocate(Rbra(totdim,3))            ;Rbra =0d0
   if(allocated(Rket))             deallocate(Rket)            ;allocate(Rket(totdim,3))            ;Rket =0d0
   if(allocated(dipole))           deallocate(dipole)          ;allocate(dipole(totdim,3))          ;dipole=zero
   if(allocated(dipole_aux))       deallocate(dipole_aux)      ;allocate(dipole_aux(totdim,3))      ;dipole_aux=zero
   if(allocated(Integral_ndx))     deallocate(Integral_ndx)    ;allocate(Integral_ndx(totdim))      ;Integral_ndx=0
   if(allocated(site_ndx))         deallocate(site_ndx)        ;allocate(site_ndx(totdim,2))        ;site_ndx=0
   if(allocated(orbital_ndx))      deallocate(orbital_ndx)     ;allocate(orbital_ndx(totdim,2))     ;orbital_ndx=0
   !
   if(optimize)then
      if(allocated(Integral_ndx_opt)) deallocate(Integral_ndx_opt);allocate(Integral_ndx_opt(dim_opt)) ;Integral_ndx_opt=0
      if(allocated(site_ndx_opt))     deallocate(site_ndx_opt)    ;allocate(site_ndx_opt(dim_opt,2))   ;site_ndx_opt=0
      if(allocated(orbital_ndx_opt))  deallocate(orbital_ndx_opt) ;allocate(orbital_ndx_opt(dim_opt,2));orbital_ndx_opt=0
      if(allocated(tab_R))            deallocate(tab_R)           ;allocate(tab_R(dim_opt))            ;tab_R=zero
      if(allocated(tab_o))            deallocate(tab_o)           ;allocate(tab_o(num_wann,num_wann))  ;tab_o=zero
   endif
   !
   !
   !1) read WS degeneracies spaces
   do i=1,qst
      read(unitIO1,*)
   enddo
   if(rst.ne.0)read(unitIO1,*)
   if(mpi_master)write(*,'(1A)')"  degen skipped"
   !
   !
   R_single = norm2(R1+R2+R3)
   R_thresh = R_single
   !
   !
   !2) read real-space vectors
   vecndx=0
   totIntegrals=0
   totIntegrals_opt=0
   do inrpts=1,nrpts
      do i=1,num_wann
         do j=1,num_wann
            vecndx=vecndx+1
            !
            read(unitIO1,'(5I5,2F12.6)') n1,n2,n3,ndx1,ndx2,a,b
            !
            Ruc_bra_ndx=floor((ndx1-0.01)/Norb)+1
            Ruc_ket_ndx=floor((ndx2-0.01)/Norb)+1
            !
            Rbra(vecndx,:) = Ruc(Ruc_bra_ndx,:) + n1*R1 + n2*R2 + n3*R3
            Rket(vecndx,:) = Ruc(Ruc_ket_ndx,:)
            !
            dist = norm2(Rbra(vecndx,:))
            !
            if(abs(dcmplx(a,b)).gt.t_thresh)then
               if(.not.(dot_product(abs([n1,n2,n3]),[1,1,1])==0 .and. Ruc_bra_ndx==Ruc_ket_ndx))then
               !
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
            endif
            !
            if(optimize)then
               if(dot_product(abs([n1,n2,n3]),[1,1,1])==1) then
                  !
                  totIntegrals_opt=totIntegrals_opt+1
                  !
                  Integral_ndx_opt(totIntegrals_opt)=vecndx
                  !
                  site_ndx_opt(totIntegrals_opt,1)=Ruc_bra_ndx
                  site_ndx_opt(totIntegrals_opt,2)=Ruc_ket_ndx
                  !
                  orbital_ndx_opt(totIntegrals_opt,1)=ndx1-floor((ndx1-0.01)/Norb)*Norb
                  orbital_ndx_opt(totIntegrals_opt,2)=ndx2-floor((ndx2-0.01)/Norb)*Norb
                  !
                  tab_R(totIntegrals_opt)=dcmplx(a,b)
                  !
               elseif(dot_product(abs([n1,n2,n3]),[1,1,1])==0) then
                  tab_o(ndx1,ndx2)=dcmplx(a,b)
               endif
            endif
            !
         enddo
      enddo
   enddo
   close(unitIO1)
   !
   !
   if(mpi_master)then
      write(*,'(A,30F12.7)')  "  starting variance vec: ",var
      !
      if(optimize)then
         !
         write(*,'(A,I8,1A)') "  optimization of: ",totIntegrals_opt," integrals"
         !
         call fmin_cgminimize(var,chi2_hopping_ratio,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
         !
         write(*,'(A,30F12.7)')  "  final variance vec: ",var
         write(*,'(A,30F12.7)')  "  final chi: ",chi
         !
         unitIO1=free_unit()
         open(unit=unitIO1,file="variance.dat",status="unknown",action="write",position="rewind")
         write(unitIO1,'(4F12.8,1I6)')var,chi,iter
         close(unitIO1)
         !
      endif
   endif
   if(optimize)call MPI_Bcast(var,size(var),MPI_DOUBLE_PRECISION,mpi_rank,MpiComm,mpi_ierr)
   !
   if(mpi_master)then
      write(*,'(2A)')      "  positions readed from: ",w90_file
      write(*,'(A,I8,1A)') "  computation of: ",totIntegrals," integrals"
      write(*,'(A,F12.6)') "  integration radius (Bohr): ",R_thresh
      write(*,'(A,F12.6)') "  ratio with u.c. radius: ",R_thresh/R_single
   endif
   !
   !
   !3) compute integrals
   dr = 2.d0*R_thresh/Nr
   do ndx = 1+mpi_rank, totIntegrals, mpi_size
      !
      Dx = 0d0
      Dy = 0d0
      Dz = 0d0
      !
      vecndx = Integral_ndx(ndx)
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
               wfc_bra = atomic_wfc(bra_orb_ndx,(x-Rbra(vecndx,1)),(y-Rbra(vecndx,2)),(z-Rbra(vecndx,3)),var(bra_orb_ndx)*R_single)
               wfc_ket = atomic_wfc(ket_orb_ndx,(x-Rket(vecndx,1)),(y-Rket(vecndx,2)),(z-Rket(vecndx,3)),var(ket_orb_ndx)*R_single)
               !
               Dx = Dx +  wfc_bra * ( x - (Rbra(vecndx,1)+Ruc(Ruc_ket_ndx,1))/2.d0 ) * wfc_ket * dr**3
               Dy = Dy +  wfc_bra * ( y - (Rbra(vecndx,2)+Ruc(Ruc_ket_ndx,2))/2.d0 ) * wfc_ket * dr**3
               Dz = Dz +  wfc_bra * ( z - (Rbra(vecndx,3)+Ruc(Ruc_ket_ndx,3))/2.d0 ) * wfc_ket * dr**3
               !
            enddo
         enddo
      enddo
      !
      dipole_aux(vecndx,1) = dcmplx(Dx,0.d0)
      dipole_aux(vecndx,2) = dcmplx(Dy,0.d0)
      dipole_aux(vecndx,3) = dcmplx(Dz,0.d0)
      !
   enddo
   call Mpi_AllReduce(dipole_aux,dipole, size(dipole), MPI_Double_Complex, MPI_Sum, MpiComm, mpi_ierr)
   call MPI_Barrier(MpiComm,mpi_ierr)
   !
   !
   !4) write to file with same indexing of W90
   if(mpi_master)then
      unitIO1=free_unit()
      open(unit=unitIO1,file=w90_file  ,status="unknown",action="read")
      unitIO2=free_unit()
      open(unit=unitIO2,file=outputfile,status="unknown",action="write")
      do i=1,qst+3
         read(unitIO1,*)
      enddo
     if(rst.ne.0)read(unitIO1,*)
      do vecndx=1,totdim
         read (unitIO1,'(5I5)') n1,n2,n3,ndx1,ndx2
         write(unitIO2,'(5I5,6F15.9)') n1,n2,n3,ndx1,ndx2,real(dipole(vecndx,1)),aimag(dipole(vecndx,1)) &
                                                         ,real(dipole(vecndx,2)),aimag(dipole(vecndx,2)) &
                                                         ,real(dipole(vecndx,3)),aimag(dipole(vecndx,3))
      enddo
      close(unitIO1)
      close(unitIO2)
   endif
   !
   deallocate(Rbra,Rket,dipole,dipole_aux,Integral_ndx,site_ndx,orbital_ndx)
   if(optimize)deallocate(Integral_ndx_opt,site_ndx_opt,orbital_ndx_opt,tab_R,tab_o)
   !
   contains



   pure function atomic_wfc(orb,x,y,z,var) result(phi)
   implicit none
   integer,intent(in)     :: orb
   real(8),intent(in)     :: x,y,z,var
   real(8)                :: phi
   real(8)                :: r,rsq,radial
   !
   !this is to avoid the 0/0 condition below
   rsq = x**2 + y**2 + z**2 + 0.0000001
   r   = sqrt(rsq)
   radial = exp((-(r/var)**2)/3.)/(2*pi*3*abs(var)) 
   !
   ! xz
   if(orb==1) phi = sqrt(15/(4*pi))*(x*z/rsq)*radial
   ! yz
   if(orb==2) phi = sqrt(15/(4*pi))*(y*z/rsq)*radial
   ! xy
   if(orb==3) phi = sqrt(15/(4*pi))*(x*y/rsq)*radial
   !
   end function atomic_wfc



   function chi2_hopping_ratio(variance) result(chi2)
     implicit none
     real(8),dimension(:)         ::  variance
     integer                      ::  Rucbra,Rucket
     integer                      ::  orbbra,orbket
     integer                      ::  Nrsmall
     real(8)                      ::  rabsq_o(num_wann,num_wann)
     real(8)                      ::  rsq,rsq_o,t_R,t_o
     real(8)                      ::  chi2
     integer,save                 ::  intern_iter=0
     !
     Nrsmall=Nr-6
     !
     dr = 2.d0*R_single/Nrsmall
     do ndx1=1,num_wann
        do ndx2=1,num_wann
           !
           Rucbra = ndx1-floor((ndx1-0.01)/Norb)*Norb
           Rucket = ndx2-floor((ndx2-0.01)/Norb)*Norb
           !
           orbbra=floor((ndx1-0.01)/Norb)+1
           orbket=floor((ndx2-0.01)/Norb)+1
           !
           rsq_o=0d0
           do iz=1,Nrsmall+1
              z = -R_single+dble(iz-1)*dr
              do iy=1,Nrsmall+1
                 y = -R_single+dble(iy-1)*dr
                 do ix=1,Nrsmall+1
                    x = -R_single+dble(ix-1)*dr
                    !
                    wfc_bra = atomic_wfc(orbbra,(x-Ruc(Rucbra,1)),(y-Ruc(Rucbra,2)),(z-Ruc(Rucbra,3)),variance(orbbra)*R_single)
                    wfc_ket = atomic_wfc(orbket,(x-Ruc(Rucket,1)),(y-Ruc(Rucket,2)),(z-Ruc(Rucket,3)),variance(orbket)*R_single)
                    !
                    rsq_o = rsq_o +  wfc_bra * ( x - (Ruc(Rucbra,1)+Ruc(Rucket,1))/2.d0 )**2 * wfc_ket * dr**3  &
                                  +  wfc_bra * ( y - (Ruc(Rucbra,2)+Ruc(Rucket,2))/2.d0 )**2 * wfc_ket * dr**3  &
                                  +  wfc_bra * ( z - (Ruc(Rucbra,3)+Ruc(Rucket,3))/2.d0 )**2 * wfc_ket * dr**3
                    !
                 enddo
              enddo
           enddo
           rabsq_o(ndx1,ndx2) = rsq_o 
        enddo
     enddo
     !
     chi2=zero
     dr = 2.d0*(1.5d0*R_single)/Nrsmall
     do ndx = 1, totIntegrals_opt
        !
        vecndx = Integral_ndx_opt(ndx)
        !
        Rucbra=site_ndx_opt(ndx,1)
        Rucket=site_ndx_opt(ndx,2)
        !
        orbbra=orbital_ndx_opt(ndx,1)
        orbket=orbital_ndx_opt(ndx,2)
        !
        t_R   = real( tab_R( ndx ))
        t_o   = real( tab_o( orbbra + (Rucbra-1)*Norb , orbket + (Rucket-1)*Norb ) )
        rsq_o =     rabsq_o( orbbra + (Rucbra-1)*Norb , orbket + (Rucket-1)*Norb )
        !
        if((t_o.le.t_thresh).and.(rsq_o.le.t_thresh))then
           cycle
        endif
        !
        rsq=0d0
        do iz=1,Nrsmall+1
           z = -(1.5d0*R_single)+dble(iz-1)*dr
           do iy=1,Nrsmall+1
              y = -(1.5d0*R_single)+dble(iy-1)*dr
              do ix=1,Nrsmall+1
                 x = -(1.5d0*R_single)+dble(ix-1)*dr
                 !
                 wfc_bra = atomic_wfc(orbbra,(x-Rbra(vecndx,1)),(y-Rbra(vecndx,2)),(z-Rbra(vecndx,3)),variance(orbbra)*R_single)
                 wfc_ket = atomic_wfc(orbket,(x-Rket(vecndx,1)),(y-Rket(vecndx,2)),(z-Rket(vecndx,3)),variance(orbket)*R_single)
                 !
                 rsq = rsq +  wfc_bra * ( x - (Rbra(vecndx,1)+Rket(vecndx,1))/2.d0 )**2 * wfc_ket * dr**3  &
                           +  wfc_bra * ( y - (Rbra(vecndx,2)+Rket(vecndx,2))/2.d0 )**2 * wfc_ket * dr**3  &
                           +  wfc_bra * ( z - (Rbra(vecndx,3)+Rket(vecndx,3))/2.d0 )**2 * wfc_ket * dr**3
                 !
              enddo
           enddo
        enddo
        !
        if((t_o.gt.t_thresh).and.(rsq_o.gt.t_thresh))then
           chi2  = chi2 + abs( rsq / rsq_o - t_R / t_o ) / (totIntegrals_opt*Nlat*Norb)
           !write(*,'(2I8,10F10.6)') ndx,vecndx, rsq,rsq_o,t_R,t_o
        endif
        !
     enddo
     !
     intern_iter=intern_iter+1
     if(intern_iter==1) write(*,'(A,1F12.6)') "  initail chi2: ",chi2
     write(*,'(A,1I6,30F12.6)') "  minimize ",intern_iter,chi2,variance
     !
   end function chi2_hopping_ratio

end subroutine dipole_t2g_LDA_mpi

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



