subroutine dmft_get_weiss_superc_main(Gloc,Floc,Smats,SAmats,Weiss,Hloc)
  complex(8),dimension(:,:,:,:,:),intent(in)      :: Gloc         ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(in)      :: Floc         ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(in)      :: Smats        !
  complex(8),dimension(:,:,:,:,:),intent(in)      :: SAmats        !
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        !
  complex(8),dimension(:,:,:,:),intent(in)        :: Hloc         ! [Nspin][Nspin][Norb][Norb]
  !aux
  complex(8),dimension(:,:,:),allocatable         :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable         :: Smats_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable         :: invGloc_site ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable         :: calG0_site   ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
  integer                                         :: Nspin,Norb,Nso,Nso2,Lmats
  integer                                         :: i,iorb,jorb,ispin,jspin,io,jo
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  !Testing part:
  Nspin = size(Gloc,1)
  Norb  = size(Gloc,3)
  Lmats = size(Gloc,5)
  Nso   = Nspin*Norb
  Nso2  = 2*Nso
  call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_main","Gloc")
  call assert_shape(Floc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_main","Floc")
  call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_main","Smats")
  call assert_shape(SAmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_main","SAmats")
  call assert_shape(Weiss,[2,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_main","Weiss")
  call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_weiss_superc_main","Hloc")
  !
  if(allocated(wm))deallocate(wm)
  allocate(wm(Lmats))
  allocate(zeta_site(Nso2,Nso2,Lmats))
  allocate(Smats_site(Nso2,Nso2,Lmats))
  allocate(invGloc_site(Nso2,Nso2,Lmats))
  allocate(calG0_site(Nso2,Nso2,Lmats))
  !
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
  zeta_site=zero
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        zeta_site(io,io,:)         = xi*wm(:) + xmu 
        zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu 
     enddo
  enddo
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              zeta_site(io,jo,:)           = zeta_site(io,jo,:)         &
                   - Hloc(ispin,jspin,iorb,jorb) - Smats(ispin,jspin,iorb,jorb,:)
              zeta_site(io,jo+Nso,:)       =     - SAmats(ispin,jspin,iorb,jorb,:)
              zeta_site(io+Nso,jo,:)       =     - SAmats(ispin,jspin,iorb,jorb,:)
              zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + &
                   Hloc(ispin,jspin,iorb,jorb) + conjg(Smats(ispin,jspin,iorb,jorb,:))
              !
              invGloc_site(io,jo,:)        = Gloc(ispin,jspin,iorb,jorb,:)
              invGloc_site(io,jo+Nso,:)    = Floc(ispin,jspin,iorb,jorb,:)
              invGloc_site(io+Nso,jo,:)    = Floc(ispin,jspin,iorb,jorb,:)
              invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(ispin,jspin,iorb,jorb,:))
              !
              Smats_site(io,jo,:)          = Smats(ispin,jspin,iorb,jorb,:)
              Smats_site(io,jo+Nso,:)      = SAmats(ispin,jspin,iorb,jorb,:)
              Smats_site(io+Nso,jo,:)      = SAmats(ispin,jspin,iorb,jorb,:)
              Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(ispin,jspin,iorb,jorb,:))
           enddo
        enddo
     enddo
  enddo
  !
  !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
  do i=1,Lmats
     call inv(invGloc_site(:,:,i))
  enddo
  !
  ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
  calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
  do i=1,Lmats
     call inv(calG0_site(:,:,i))
  enddo
  !
  !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
  !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              io = iorb + (ispin-1)*Norb
              jo = jorb + (jspin-1)*Norb
              Weiss(1,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
              Weiss(2,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
           enddo
        enddo
     enddo
  enddo
  !
end subroutine dmft_get_weiss_superc_main




subroutine dmft_get_weiss_superc_ineq(Gloc,Floc,Smats,SAmats,Weiss,Hloc)
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Floc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: Smats        !
  complex(8),dimension(:,:,:,:,:,:),intent(in)      :: SAmats        ! 
  complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Weiss        ! 
  complex(8),dimension(:,:,:,:,:),intent(in)        :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
  !aux
  complex(8),dimension(:,:,:),allocatable           :: zeta_site    ![2*Nspin*Norb][2*Nspin*Norb][Lmats]
  complex(8),dimension(:,:,:),allocatable           :: Smats_site   !
  complex(8),dimension(:,:,:),allocatable           :: invGloc_site !
  complex(8),dimension(:,:,:),allocatable           :: calG0_site   !
  integer                                           :: Nlat,Nspin,Norb,Nso,Nso2,Nlso,Lmats
  integer                                           :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js,inambu,jnambu
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  !Testing part:
  Nlat  = size(Gloc,1)
  Nspin = size(Gloc,2)
  Norb  = size(Gloc,4)
  Lmats = size(Gloc,6)
  Nso   = Nspin*Norb
  Nso2  = 2*Nso
  Nlso  = Nlat*Nspin*Norb
  call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_ineq","Gloc")
  call assert_shape(Floc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_ineq","Floc")
  call assert_shape(Smats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_ineq","Smats")
  call assert_shape(SAmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_ineq","SAmats")
  call assert_shape(Weiss,[2,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_superc_ineq","Weiss")
  call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_weiss_superc_ineq","Hloc")
  !
  if(allocated(wm))deallocate(wm)
  allocate(wm(Lmats))
  allocate(zeta_site(Nso2,Nso2,Lmats))
  allocate(Smats_site(Nso2,Nso2,Lmats))
  allocate(invGloc_site(Nso2,Nso2,Lmats))
  allocate(calG0_site(Nso2,Nso2,Lmats))
  !
  wm = pi/beta*(2*arange(1,Lmats)-1)
  Weiss       = zero
  do ilat=1,Nlat
     !Dump the Gloc and the Smats for the ilat-th site into a [Norb*Nspin]^2 matrix and create the zeta_site
     zeta_site=zero
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           zeta_site(io,io,:)         = xi*wm(:) + xmu - Hloc(ilat,ispin,ispin,iorb,iorb)
           zeta_site(io+Nso,io+Nso,:) = xi*wm(:) - xmu + Hloc(ilat,ispin,ispin,iorb,iorb)
        enddo
     enddo
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 !
                 zeta_site(io,jo,:)           = zeta_site(io,jo,:)         - Smats(ilat,ispin,jspin,iorb,jorb,:)
                 zeta_site(io,jo+Nso,:)       =-SAmats(ilat,ispin,jspin,iorb,jorb,:)
                 zeta_site(io+Nso,jo,:)       =-SAmats(ilat,ispin,jspin,iorb,jorb,:)
                 zeta_site(io+Nso,jo+Nso,:)   = zeta_site(io+Nso,jo+Nso,:) + conjg(Smats(ilat,ispin,jspin,iorb,jorb,:))
                 !
                 invGloc_site(io,jo,:)        = Gloc(ilat,ispin,jspin,iorb,jorb,:)
                 invGloc_site(io,jo+Nso,:)    = Floc(ilat,ispin,jspin,iorb,jorb,:)
                 invGloc_site(io+Nso,jo,:)    = Floc(ilat,ispin,jspin,iorb,jorb,:)
                 invGloc_site(io+Nso,jo+Nso,:)=-conjg(Gloc(ilat,ispin,jspin,iorb,jorb,:))
                 !
                 Smats_site(io,jo,:)          = Smats(ilat,ispin,jspin,iorb,jorb,:)
                 Smats_site(io,jo+Nso,:)      = SAmats(ilat,ispin,jspin,iorb,jorb,:)
                 Smats_site(io+Nso,jo,:)      = SAmats(ilat,ispin,jspin,iorb,jorb,:)
                 Smats_site(io+Nso,jo+Nso,:)  =-conjg(Smats(ilat,ispin,jspin,iorb,jorb,:))
              enddo
           enddo
        enddo
     enddo
     !Invert the ilat-th site [Norb*Nspin]**2 Gloc block matrix 
     do i=1,Lmats
        call inv(invGloc_site(:,:,i))
     enddo
     !
     ![calG0]_ilat = [ [Gloc]_ilat^-1 + [Smats]_ilat ]^-1
     calG0_site(:,:,1:Lmats) = invGloc_site(:,:,1:Lmats) + Smats_site(:,:,1:Lmats)
     do i=1,Lmats
        call inv(calG0_site(:,:,i))
     enddo
     !
     !Dump back the [Norb*Nspin]**2 block of the ilat-th site into the 
     !output structure of [Nlat,Nspsin,Nspin,Norb,Norb] matrix
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 Weiss(1,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo,:)
                 Weiss(2,ilat,ispin,jspin,iorb,jorb,:) = calG0_site(io,jo+Nso,:)
              enddo
           enddo
        enddo
     enddo
  end do
  !
end subroutine dmft_get_weiss_superc_ineq
