!> MATSUBARA: single site
subroutine dmft_gf_print_matsubara_main(Gmats,fname,iprint)
  complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats
  character(len=*),intent(in)                :: fname
  integer,intent(in)                         :: iprint
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  !    
  Nspin = size(Gmats,1)
  Norb  = size(Gmats,3)
  Lmats = size(Gmats,5)
  call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_print_matsubara",reg(fname)//"_mats")
  !
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  select case(iprint)
  case(1)                  !print only diagonal elements
     write(*,"(A,1x,A)") reg(fname),"matsubara: write spin-orbital diagonal elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_iw"//str(gf_suffix)
           call splot(reg(suffix),wm,Gmats(ispin,ispin,iorb,iorb,:))
        enddo
     enddo
     !
  case(2)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)") reg(fname),"matsubara: write spin diagonal and all orbitals elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_iw"//reg(gf_suffix)
              call splot(reg(suffix),wm,Gmats(ispin,ispin,iorb,jorb,:))
           enddo
        enddo
     enddo
     !
  case default                  !print all off-diagonals
     write(*,"(A,1x,A)") reg(fname),"matsubara: write all elements."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_iw"//reg(gf_suffix)
                 call splot(reg(suffix),wm,Gmats(ispin,jspin,iorb,jorb,:))
              enddo
           enddo
        enddo
     enddo
     !
  end select
end subroutine dmft_gf_print_matsubara_main

!> MATSUBARA: ineq sites
subroutine dmft_gf_print_matsubara_ineq(Gmats,fname,iprint,ineq_index,ineq_pad,itar)
  complex(8),dimension(:,:,:,:,:,:),intent(in) :: Gmats
  character(len=*),intent(in)                  :: fname
  integer,intent(in)                           :: iprint
  character(len=*),optional                    :: ineq_index
  integer,optional                             :: ineq_pad
  logical,optional                             :: itar
  character(len=:),allocatable                 :: index
  integer                                      :: pad
  !
  if(present(ineq_index))then
     index=trim(adjustl(trim(ineq_index)))
  else
     index='indx'
  endif
  pad=6
  if(present(ineq_pad))pad=ineq_pad

  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  !    
  !
  Nlat  = size(Gmats,1)
  Nspin = size(Gmats,2)
  Norb  = size(Gmats,4)
  Lmats = size(Gmats,6)
  call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gf_print_matsubara_ineq",reg(fname)//"_mats")
  !
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  select case(iprint)
  case(1)                  !print only diagonal elements
     write(*,"(A,1x,A)")reg(fname),"matsubara: write spin-orbital diagonal elements. No Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_iw"//reg(gf_suffix)
           call splot(reg(suffix),wm,Gmats(:,ispin,ispin,iorb,iorb,:))
           call file_gzip(reg(suffix))
        enddo
     enddo
     !
  case(2)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. No Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_iw"//reg(gf_suffix)
              call splot(reg(suffix),wm,Gmats(:,ispin,ispin,iorb,jorb,:))
              call file_gzip(reg(suffix))
           enddo
        enddo
     enddo
     !
  case(3)                  !print all off-diagonals
     write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. No Split."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_iw"//reg(gf_suffix)
                 call splot(reg(suffix),wm,Gmats(:,ispin,jspin,iorb,jorb,:))
                 call file_gzip(reg(suffix))
              enddo
           enddo
        enddo
     enddo
     !
  case(4)                  !print only diagonal elements
     write(*,"(A,1x,A)")reg(fname),"matsubara: write spin-orbital diagonal elements. Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           do ilat=1,Nlat
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(iorb)//&
                   "_s"//str(ispin)//&
                   "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
              call splot(reg(suffix),wm,Gmats(ilat,ispin,ispin,iorb,iorb,:))
           enddo
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_iw"
           call file_targz(tarball=reg(suffix),&
                pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
        enddo
     enddo
     !
  case(5)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements. Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              do ilat=1,Nlat
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//&
                      "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                 call splot(reg(suffix),wm,Gmats(ilat,ispin,ispin,iorb,jorb,:))
              enddo
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_iw"
              call file_targz(tarball=reg(suffix),&
                   pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
           enddo
        enddo
     enddo
     !
  case default
     write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements. Split."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 do ilat=1,Nlat
                    suffix=reg(fname)//&
                         "_l"//str(iorb)//str(jorb)//&
                         "_s"//str(ispin)//str(jspin)//&
                         "_iw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                    call splot(reg(suffix),wm,Gmats(ilat,ispin,jspin,iorb,jorb,:))
                 enddo
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_iw"
                 call file_targz(tarball=reg(suffix),&
                      pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
              enddo
           enddo
        enddo
     enddo
     !
  end select
end subroutine dmft_gf_print_matsubara_ineq



!> MATSUBARA: full Gij
subroutine dmft_gij_print_matsubara(Gmats,fname,iprint)
  complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Gmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats/Lreal]
  real(8),dimension(size(Gmats,7))               :: w
  character(len=*),intent(in)                    :: fname
  integer,intent(in)                             :: iprint
  !
  !Retrieve parameters:
  call get_ctrl_var(beta,"BETA")
  !
  Nlat  = size(Gmats,1)
  Nspin = size(Gmats,3)
  Norb  = size(Gmats,5)
  Lmats = size(Gmats,7)
  call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats],"dmft_gij_print_matsubara",reg(fname)//"_mats")
  !
  if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  !
  select case(iprint)
  case (0)
     write(*,"(A,1x,A)")reg(fname),"matsubara: not written on file."
  case(1)                  !print only diagonal elements     
     write(*,"(A,1x,A)")reg(fname),"matsubara: spin-orbital diagonal elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_iw"//reg(gf_suffix)
           call splot(reg(suffix),wm,Gmats(:,:,ispin,ispin,iorb,iorb,:))
           call file_gzip(reg(suffix))
        enddo
     enddo
  case(2)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)")reg(fname),"matsubara: write spin diagonal and all orbitals elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_iw"//reg(gf_suffix)
              call splot(reg(suffix),wm,Gmats(:,:,ispin,ispin,iorb,jorb,:))
              call file_gzip(reg(suffix))
           enddo
        enddo
     enddo
  case default              !print all off-diagonals
     write(*,"(A,1x,A)")reg(fname),"matsubara: write all elements."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_iw"//reg(gf_suffix)
                 call splot(reg(suffix),wm,Gmats(:,:,ispin,jspin,iorb,jorb,:))
                 call file_gzip(reg(suffix))
              enddo
           enddo
        enddo
     enddo
  end select
end subroutine dmft_gij_print_matsubara





!##################################################################
!##################################################################
!##################################################################







!> REALAXIS: single site
subroutine dmft_gf_print_realaxis_main(Greal,fname,iprint)
  complex(8),dimension(:,:,:,:,:),intent(in) :: Greal
  character(len=*),intent(in)                :: fname
  integer,intent(in)                         :: iprint
  !
  !Retrieve parameters:
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  !
  Nspin = size(Greal,1)
  Norb  = size(Greal,3)
  Lreal = size(Greal,5)
  call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_print_realaxis",reg(fname)//"_real")
  !
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  select case(iprint)
  case(1)                  !print only diagonal elements
     write(*,"(A,1x,A)") reg(fname),"real: write spin-orbital diagonal elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_realw"//reg(gf_suffix)
           call splot(reg(suffix),wr,Greal(ispin,ispin,iorb,iorb,:))
        enddo
     enddo
     !
  case(2)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)") reg(fname),"real: write spin diagonal and all orbitals elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_realw"//reg(gf_suffix)
              call splot(reg(suffix),wr,Greal(ispin,ispin,iorb,jorb,:))
           enddo
        enddo
     enddo
     !
  case default                  !print all off-diagonals
     write(*,"(A,1x,A)") reg(fname),"real: write all elements."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_realw"//reg(gf_suffix)
                 call splot(reg(suffix),wr,Greal(ispin,jspin,iorb,jorb,:))
              enddo
           enddo
        enddo
     enddo
     !
  end select
end subroutine dmft_gf_print_realaxis_main



!> REALAXIS: ineq sites
subroutine dmft_gf_print_realaxis_ineq(Greal,fname,iprint,ineq_index,ineq_pad)
  complex(8),dimension(:,:,:,:,:,:),intent(in) :: Greal
  character(len=*),intent(in)                  :: fname
  integer,intent(in)                           :: iprint
  character(len=*),optional                    :: ineq_index
  integer,optional                             :: ineq_pad
  character(len=:),allocatable                 :: index
  integer                                      :: pad
  !
  if(present(ineq_index))then
     index=trim(adjustl(trim(ineq_index)))
  else
     index='indx'
  endif
  pad=6
  if(present(ineq_pad))pad=ineq_pad

  !
  !Retrieve parameters:
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  !
  Nlat  = size(Greal,1)
  Nspin = size(Greal,2)
  Norb  = size(Greal,4)
  Lreal = size(Greal,6)
  call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gf_print_realaxis_ineq",reg(fname)//"_real")
  !
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  select case(iprint)
  case(1)                  !print only diagonal elements
     write(*,"(A,1x,A)")reg(fname),"real: write spin-orbital diagonal elements. No Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_realw"//reg(gf_suffix)
           call splot(reg(suffix),wr,Greal(:,ispin,ispin,iorb,iorb,:))
           call file_gzip(reg(suffix))
        enddo
     enddo
     !
  case(2)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. No Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_realw"//reg(gf_suffix)
              call splot(reg(suffix),wr,Greal(:,ispin,ispin,iorb,jorb,:))
              call file_gzip(reg(suffix))
           enddo
        enddo
     enddo
     !
  case(3)
     write(*,"(A,1x,A)")reg(fname),"real: write all elements. No Split."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_realw"//reg(gf_suffix)
                 call splot(reg(suffix),wr,Greal(:,ispin,jspin,iorb,jorb,:))
                 call file_gzip(reg(suffix))
              enddo
           enddo
        enddo
     enddo
     !
  case(4)                  !print only diagonal elements
     write(*,"(A,1x,A)")reg(fname),"real: write spin-orbital diagonal elements. Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           do ilat=1,Nlat
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(iorb)//&
                   "_s"//str(ispin)//&
                   "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
              call splot(reg(suffix),wr,Greal(ilat,ispin,ispin,iorb,iorb,:))
           enddo
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_realw"
           call file_targz(tarball=reg(suffix),&
                pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
        enddo
     enddo
     !
  case(5)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements. Split."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              do ilat=1,Nlat
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//&
                      "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                 call splot(reg(suffix),wr,Greal(ilat,ispin,ispin,iorb,jorb,:))
              enddo
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_realw"
              call file_targz(tarball=reg(suffix),&
                   pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
           enddo
        enddo
     enddo
     !
  case default
     write(*,"(A,1x,A)")reg(fname),"real: write all elements. Split."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 do ilat=1,Nlat
                    suffix=reg(fname)//&
                         "_l"//str(iorb)//str(jorb)//&
                         "_s"//str(ispin)//str(jspin)//&
                         "_realw_"//reg(index)//str(ilat,pad)//reg(gf_suffix)
                    call splot(reg(suffix),wr,Greal(ilat,ispin,jspin,iorb,jorb,:))
                 enddo
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_realw"
                 call file_targz(tarball=reg(suffix),&
                      pattern=reg(suffix)//"_"//reg(index)//"*"//reg(gf_suffix))
              enddo
           enddo
        enddo
     enddo
     !
  end select
end subroutine dmft_gf_print_realaxis_ineq


!> REALAXIS: full GF
subroutine dmft_gij_print_realaxis(Greal,fname,iprint)
  complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Greal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal/Lreal]
  real(8),dimension(size(Greal,7))               :: w
  character(len=*),intent(in)                    :: fname
  integer,intent(in)                             :: iprint
  integer                                        :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb,ilat,jlat
  !
  !Retrieve parameters:
  call get_ctrl_var(wini,"WINI")
  call get_ctrl_var(wfin,"WFIN")
  !
  Nlat  = size(Greal,1)
  Nspin = size(Greal,3)
  Norb  = size(Greal,5)
  Lreal = size(Greal,7)
  call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal],"dmft_gij_print_realaxis",reg(fname)//"_real")
  !
  if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  !
  select case(iprint)
  case (0)
     write(*,"(A,1x,A)")reg(fname),"real: not written on file."
  case(1)                  !print only diagonal elements     
     write(*,"(A,1x,A)")reg(fname),"real: spin-orbital diagonal elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           suffix=reg(fname)//&
                "_l"//str(iorb)//str(iorb)//&
                "_s"//str(ispin)//&
                "_realw"//reg(gf_suffix)
           call splot(reg(suffix),wr,Greal(:,:,ispin,ispin,iorb,iorb,:))
           call file_gzip(reg(suffix))
        enddo
     enddo
  case(2)                  !print spin-diagonal, all orbitals 
     write(*,"(A,1x,A)")reg(fname),"real: write spin diagonal and all orbitals elements."
     do ispin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              suffix=reg(fname)//&
                   "_l"//str(iorb)//str(jorb)//&
                   "_s"//str(ispin)//&
                   "_realw"//reg(gf_suffix)
              call splot(reg(suffix),wr,Greal(:,:,ispin,ispin,iorb,jorb,:))
              call file_gzip(reg(suffix))
           enddo
        enddo
     enddo
  case default              !print all off-diagonals
     write(*,"(A,1x,A)")reg(fname),"real: write all elements."
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 suffix=reg(fname)//&
                      "_l"//str(iorb)//str(jorb)//&
                      "_s"//str(ispin)//str(jspin)//&
                      "_realw"//reg(gf_suffix)
                 call splot(reg(suffix),wr,Greal(:,:,ispin,jspin,iorb,jorb,:))
                 call file_gzip(reg(suffix))
              enddo
           enddo
        enddo
     enddo
  end select
end subroutine dmft_gij_print_realaxis


