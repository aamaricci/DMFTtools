subroutine dmft_interaction_setKanamori_serial(Utensor,Uloc,Ust,Jh,Jx,Jp)
  implicit none
  real(8),dimension(:,:,:,:,:),intent(inout)    :: Utensor  ! [Norb][Norb][Norb][Norb][2]
  real(8),dimension(:),intent(in)               :: Uloc
  real(8),intent(in)                            :: Ust
  real(8),intent(in)                            :: Jh,Jx,Jp
  !
  integer                                       :: Norb
  integer                                       :: iorb,jorb,korb,lorb
  integer                                       :: io,jo
  logical                                       :: test=.true.
  !
  !Retrieve parameters:
  call get_ctrl_var(Norb,"NORB")
  !
  !
  !Testing part:
  call assert_shape(Utensor,[Norb,Norb,Norb,Norb,2],"dmft_interaction_setKanamori_mpi","Utensor")
  !
  Utensor=0d0
  !
  !density-density interaction: same orbital, opposite spins
  do iorb=1,Norb
     Utensor(iorb,iorb,iorb,iorb,1) = Uloc(iorb)
  enddo
  !
  !density-density interaction: different orbitals, opposite spins
  do iorb=1,Norb
     do jorb=1+iorb,Norb
        Utensor(iorb,iorb,jorb,jorb,1) = Ust
        Utensor(jorb,jorb,iorb,iorb,1) = Ust
     enddo
  enddo
  !
  !density-density interaction: different orbitals, parallel spins
  do iorb=1,Norb
     do jorb=1+iorb,Norb
        Utensor(iorb,iorb,jorb,jorb,2) = (Ust-Jh)/2.d0
        Utensor(jorb,jorb,iorb,iorb,2) = (Ust-Jh)/2.d0
     enddo
  enddo
  !
  !spin-exchange and pair-hopping
  do iorb=1,Norb
     do jorb=1+iorb,Norb
        Utensor(iorb,jorb,jorb,iorb,1) = Jx
        Utensor(jorb,iorb,iorb,jorb,1) = Jx
        Utensor(iorb,jorb,iorb,jorb,1) = Jp
        Utensor(jorb,iorb,jorb,iorb,1) = Jp
     enddo
  enddo
  !
  !These elements are unessential in the impurity problem (tested) due to non existence of connecting states
  !but mandatory if the rotational invariance of the interaction has to be preserved
  if(test)then
     do iorb=1,Norb
        do jorb=1+iorb,Norb
           Utensor(iorb,jorb,jorb,iorb,2) = -(Ust-Jh)/2.d0
           Utensor(jorb,iorb,iorb,jorb,2) = -(Ust-Jh)/2.d0
        enddo
     enddo
  endif
  !
end subroutine dmft_interaction_setKanamori_serial



subroutine dmft_interaction_rotate_serial(Utensor,rot)
  implicit none
  real(8),dimension(:,:,:,:,:),intent(inout)    :: Utensor  ! [Norb][Norb][Norb][Norb][2]
  complex(8),dimension(:,:),intent(in)          :: rot
  !
  real(8),dimension(:,:,:),allocatable          :: Umatrix
  complex(8),dimension(:,:,:),allocatable       :: Nab
  complex(8),dimension(:,:),allocatable         :: rotdag
  integer                                       :: Norb
  integer                                       :: iorb,jorb,korb,lorb
  integer                                       :: io,jo,count
  !
  !Retrieve parameters:
  call get_ctrl_var(Norb,"NORB")
  !
  !Testing part:
  call assert_shape(Utensor,[Norb,Norb,Norb,Norb,2],"dmft_interaction_rotate_mpi","Utensor")
  if(size(rot,1).ne.Norb) stop "Rotation currently implemented only in the orbital space"
  !
  allocate(Umatrix(Norb*Norb,Norb*Norb,2));Umatrix=0d0
  allocate(Nab(Norb,Norb,Norb*Norb));Nab=zero
  !
  allocate(rotdag(Norb,Norb));rotdag=zero
  rotdag=transpose(conjg(rot))
  !
  !N=c^+,c = {col,row}
  count=0
  do iorb=1,Norb
     do jorb=1,Norb
        count=count+1
        Nab(jorb,iorb,count) = cmplx(1.d0,0d0)
     enddo
  enddo
  !
  do io=1,Norb*Norb
     Nab(:,:,io) = matmul(rotdag,matmul(Nab(:,:,io),rot))
  enddo
  !
  do io=1,Norb*Norb
     do jo=1,Norb*Norb
        call getMKstride(io,iorb,jorb,Norb)
        call getMKstride(jo,korb,lorb,Norb)
        Umatrix(:,:,1) = Umatrix(:,:,1) + kronecker_product(Nab(:,:,io),Nab(:,:,jo))*Utensor(iorb,jorb,korb,lorb,1)
        Umatrix(:,:,2) = Umatrix(:,:,2) + kronecker_product(Nab(:,:,io),Nab(:,:,jo))*Utensor(iorb,jorb,korb,lorb,2)
     enddo
  enddo
  !
  forall(io=1:Norb*Norb,jo=1:Norb*Norb, abs(Umatrix(io,jo,1)).lt.1e-4) Umatrix(io,jo,1)=0d0
  forall(io=1:Norb*Norb,jo=1:Norb*Norb, abs(Umatrix(io,jo,2)).lt.1e-4) Umatrix(io,jo,2)=0d0
  !
  do iorb=1,Norb !row of external submatrix
     do jorb=1,Norb !col of external submatrix
        do korb=1,Norb !row of internal submatrix
           do lorb=1,Norb !col of internal submatrix
              !
              io = korb + (iorb-1)*Norb
              jo = lorb + (jorb-1)*Norb
              !
              Utensor(iorb,jorb,korb,lorb,1) = Umatrix(io,jo,1)
              Utensor(iorb,jorb,korb,lorb,2) = Umatrix(io,jo,2)
           enddo
        enddo
     enddo
  enddo
  !
contains

   subroutine getUKstride(ndx,iorb,jorb,Norb)
     implicit none
     integer,intent(in)                          :: iorb,jorb,Norb
     integer,intent(out)                         :: ndx
     integer,dimension(Norb,Norb)                :: stride
     integer                                     :: count
     integer                                     :: io,jo
     !
     count=0
     do io=1,Norb
        do jo=1,Norb
           count=count+1
           stride(io,jo)=count
        enddo
     enddo
     !
     ndx = stride(iorb,jorb)
     !
   end subroutine getUKstride

   subroutine getMKstride(ndx,iorb,jorb,Norb)
     implicit none
     integer,intent(in)                          :: ndx,Norb
     integer,intent(out)                         :: iorb,jorb
     integer,dimension(Norb,Norb)                :: stride
     integer                                     :: count
     integer                                     :: io,jo
     !
     count=0
     do io=1,Norb
        do jo=1,Norb
           count=count+1
           stride(io,jo)=count
        enddo
     enddo
     !
     strideloop:do iorb=1,Norb
        do jorb=1,Norb
           if(stride(iorb,jorb)==ndx) exit strideloop
        enddo
     enddo strideloop
     !
   end subroutine getMKstride

end subroutine dmft_interaction_rotate_serial
