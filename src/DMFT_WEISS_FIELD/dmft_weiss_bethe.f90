subroutine dmft_get_weiss_normal_bethe(Gloc,Weiss,Hloc,Wbands)
  complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
  real(8),dimension(:),intent(in)               :: Wbands ![Nspin*Norb]
  !aux
  complex(8),dimension(size(Gloc,5))            :: invWeiss ![Lmats]
  integer                                       :: Nspin,Norb,Nso,Lmats
  integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
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
  call assert_shape(Gloc,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_bethe","Gloc")
  call assert_shape(Weiss,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_bethe","Weiss")
  call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_bethe","Hloc")
  call assert_shape(Wbands,[Nspin*Norb],"dmft_get_weiss_normal_bethe","Wbands")
  !
  if(allocated(wm))deallocate(wm)
  allocate(wm(Lmats))
  !
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  !\calG0^{-1}_aa = iw + mu - H_0 - d**2/4*Gmats
  Weiss=zero
  do ispin=1,Nspin
     do iorb=1,Norb
        invWeiss = (xi*wm(:)+xmu) - Hloc(ispin,ispin,iorb,iorb) - 0.25d0*Wbands(iorb+(ispin-1)*Norb)*Gloc(ispin,ispin,iorb,iorb,:)
        Weiss(ispin,ispin,iorb,iorb,:) = one/invWeiss
     enddo
  enddo
  !
  deallocate(wm)
end subroutine dmft_get_weiss_normal_bethe



subroutine dmft_get_weiss_normal_bethe_ineq(Gloc,Weiss,Hloc,Wbands)
  complex(8),dimension(:,:,:,:,:,:),intent(in)    :: Gloc         ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Weiss        ! [Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),dimension(:,:,:,:,:),intent(in)      :: Hloc         ! [Nlat][Nspin][Nspin][Norb][Norb]
  real(8),dimension(:),intent(in)                 :: Wbands       ! [Nlat*Nspin*Norb]
  !aux
  complex(8),dimension(size(Gloc,6))              :: invWeiss     ![Lmats]
  integer                                         :: Nlat,Nspin,Norb,Nso,Nlso,Lmats
  integer                                         :: i,j,iorb,jorb,ispin,jspin,ilat,jlat,io,jo,js
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
  Nlso  = Nlat*Nspin*Norb
  call assert_shape(Gloc,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_bethe_ineq","Gloc")
  call assert_shape(Weiss,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_get_weiss_normal_bethe_ineq","Weiss")
  call assert_shape(Hloc,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_get_weiss_normal_ineq","Hloc")
  call assert_shape(Wbands,[Nlat*Nspin*Norb],"dmft_get_weiss_normal_ineq","Wbands")
  !
  if(allocated(wm))deallocate(wm)
  allocate(wm(Lmats))
  !
  wm = pi/beta*(2*arange(1,Lmats)-1)
  Weiss     = zero
  do ilat=1,Nlat
     !\calG0^{-1}_aa_i = iw + mu - H_0 - d**2/4*Gmats_aa_i
     Weiss=zero
     do ispin=1,Nspin
        do iorb=1,Norb
           io       = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
           invWeiss = (xi*wm(:)+xmu) - Hloc(ilat,ispin,ispin,iorb,iorb) - 0.25d0*Wbands(io)*Gloc(ilat,ispin,ispin,iorb,iorb,:)
           Weiss(ilat,ispin,ispin,iorb,iorb,:) = one/invWeiss
        enddo
     enddo
  end do
  !
  deallocate(wm)
  !
end subroutine dmft_get_weiss_normal_bethe_ineq
