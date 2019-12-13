subroutine TB_reset_ei
  ei_x=[1d0,0d0,0d0]
  ei_y=[0d0,1d0,0d0]
  ei_z=[0d0,0d0,1d0]
  set_eivec=.false.
end subroutine TB_reset_ei


subroutine TB_reset_bk
  bk_x=[1d0,0d0,0d0]*pi2
  bk_y=[0d0,1d0,0d0]*pi2
  bk_z=[0d0,0d0,1d0]*pi2
  set_bkvec=.false.
end subroutine TB_reset_bk


subroutine TB_set_ei(eix,eiy,eiz)
  real(8),dimension(:),intent(in)                  :: eix
  real(8),dimension(size(eix)),intent(in),optional :: eiy
  real(8),dimension(size(eix)),intent(in),optional :: eiz
  ei_x = 0d0
  ei_x(1:size(eix)) = eix
  !
  if(present(eiy))then
     ei_y = 0d0
     ei_y(1:size(eiy)) = eiy
  endif
  !
  if(present(eiz))then
     ei_z = 0d0
     ei_z(1:size(eiz)) = eiz
  endif
  !
  set_eivec=.true.
end subroutine TB_set_ei


subroutine TB_set_bk(bkx,bky,bkz)
  real(8),dimension(:),intent(in)                  :: bkx
  real(8),dimension(size(bkx)),intent(in),optional :: bky
  real(8),dimension(size(bkx)),intent(in),optional :: bkz
  bk_x = 0d0
  bk_x(1:size(bkx)) = bkx
  !
  if(present(bky))then
     bk_y = 0d0
     bk_y(1:size(bky)) = bky
  endif
  !
  if(present(bkz))then
     bk_z = 0d0
     bk_z(1:size(bkz)) = bkz
  endif
  !
  set_bkvec=.true.
end subroutine TB_set_bk





subroutine TB_get_bk(bkx,bky,bkz)
  real(8),dimension(:),intent(inout)                  :: bkx
  real(8),dimension(size(bkx)),intent(inout),optional :: bky
  real(8),dimension(size(bkx)),intent(inout),optional :: bkz
  real(8),dimension(3)                                :: b1,b2,b3
  !
  call TB_reciprocal_basis(a1=ei_x,a2=ei_y,a3=ei_z,b1=b1,b2=b2,b3=b3)
  !
  bkx = b1(1:size(bkx))
  if(present(bky))bky = b2(1:size(bky))
  if(present(bkz))bkz = b3(1:size(bkz))
  !
end subroutine TB_get_bk

subroutine TB_get_ei(eix,eiy,eiz)
  real(8),dimension(:),intent(inout)                  :: eix
  real(8),dimension(size(eix)),intent(inout),optional :: eiy
  real(8),dimension(size(eix)),intent(inout),optional :: eiz
  real(8),dimension(3)                                :: a1,a2,a3
  !
  call TB_reciprocal_basis(a1=bk_x,a2=bk_y,a3=bk_z,b1=a1,b2=a2,b3=a3)
  !
  eix = a1(1:size(eix))
  if(present(eiy))eiy = a2(1:size(eiy))
  if(present(eiz))eiz = a3(1:size(eiz))
  !
end subroutine TB_get_ei





subroutine TB_build_bk(verbose)
  logical,optional :: verbose
  logical          :: verbose_
  if(.not.set_eivec)stop "TB_build_bk ERROR: Direct basis not set, set_eivec=F"
  verbose_=.false.;if(present(verbose))verbose_=verbose
  call TB_reciprocal_basis(&
       a1=ei_x,a2=ei_y,a3=ei_z,&
       b1=bk_x,b2=bk_y,b3=bk_z)
  set_bkvec=.true.
  if(verbose_)call print_bk
end subroutine TB_build_bk

subroutine TB_build_ei(verbose)
  logical,optional :: verbose
  logical          :: verbose_
  if(.not.set_bkvec)stop "TB_build_ei ERROR: Reciprocal basis not set, set_bkvec=F"
  verbose_=.false.;if(present(verbose))verbose_=verbose
  !note that we exchange the direct basis with the reciprocal
  call TB_reciprocal_basis(&
       b1=ei_x,b2=ei_y,b3=ei_z,&
       a1=bk_x,a2=bk_y,a3=bk_z)
  set_eivec=.true.
  if(verbose_)call print_ei
end subroutine TB_build_ei





subroutine print_ei(pfile)
  character(len=*),optional :: pfile
  integer                   :: unit,i
  if(io_eivec)return
  unit=6
  if(present(pfile))open(free_unit(unit),file=reg(pfile))
  write(unit,"(A)")"Using Direct Lattice vectors:"
  write(unit,"(A,3F8.4,A1)")"ei_x = [",(ei_x(i),i=1,3),"]"
  write(unit,"(A,3F8.4,A1)")"ei_y = [",(ei_y(i),i=1,3),"]"
  write(unit,"(A,3F8.4,A1)")"ei_z = [",(ei_z(i),i=1,3),"]"
  if(present(pfile))close(unit)
  io_eivec=.true.
end subroutine print_ei

subroutine print_bk(pfile)
  character(len=*),optional :: pfile
  integer                   :: unit,i
  if(io_bkvec)return
  unit=6
  if(present(pfile))open(free_unit(unit),file=reg(pfile))
  write(unit,"(A)")"Using Reciprocal Lattice vectors:"
  write(unit,"(A,3F8.4,A1)")"bk_x = [",(bk_x(i),i=1,3),"]"
  write(unit,"(A,3F8.4,A1)")"bk_y = [",(bk_y(i),i=1,3),"]"
  write(unit,"(A,3F8.4,A1)")"bk_z = [",(bk_z(i),i=1,3),"]"
  if(present(pfile))close(unit)
  io_bkvec=.true.
end subroutine print_bk





!-------------------------------------------------------------------------------------------
!PURPOSE:  Build the reciprocal space {b1[,b2,b3]} given the direct basis {a1[,a2,a3]}
!  in dimensions up to 3
!-------------------------------------------------------------------------------------------
subroutine TB_reciprocal_basis(a1,a2,a3, b1,b2,b3)
  !This routine generates the reciprocal lattice vectors b1,b2,b3
  !given the space vectors a1,a2,a3.
  !the vectors are in units of the lattice constant (a=1)
  real(8),dimension(:),intent(in)                 :: a1
  real(8),dimension(size(a1)),intent(in),optional :: a2
  real(8),dimension(size(a1)),intent(in),optional :: a3
  real(8),dimension(size(a1))                     :: b1
  real(8),dimension(size(a1)),optional            :: b2
  real(8),dimension(size(a1)),optional            :: b3
  !
  real(8),dimension(3)                            :: ar1,ar2,ar3
  real(8),dimension(3)                            :: bk1,bk2,bk3
  real(8)                                         :: den, s
  integer                                         :: iperm, i, j, k, l, ipol
  integer                                         :: N
  real(8),dimension(3,3)                          :: Mat
  !
  N = size(a1)
  !
  ar1=[1d0,0d0,0d0]
  ar2=[0d0,1d0,0d0]
  ar3=[0d0,0d0,1d0]
  !
  ar1(:N)=a1
  if(present(a2))ar2(:N)=a2
  if(present(a3))ar3(:N)=a3
  !
  den = det(transpose(reshape([ar1,ar2,ar3],shape=[3,3])))
  !
  !    here we compute the reciprocal vectors
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     bk1(ipol) = (ar2(j)*ar3(k) - ar2(k)*ar3 (j) )/den*pi2
     bk2(ipol) = (ar3(j)*ar1(k) - ar3(k)*ar1 (j) )/den*pi2
     bk3(ipol) = (ar1(j)*ar2(k) - ar1(k)*ar2 (j) )/den*pi2
     l = i
     i = j
     j = k
     k = l
  enddo
  b1=bk1(:N)
  if(present(b2))b2=bk2(:N)
  if(present(b3))b3=bk3(:N)
  return
end subroutine TB_reciprocal_basis
