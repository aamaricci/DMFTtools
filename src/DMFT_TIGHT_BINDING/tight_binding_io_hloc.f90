subroutine write_hloc_1(hloc,file) ![Nlso][Nlso]
  complex(8),dimension(:,:) :: Hloc
  character(len=*),optional :: file
  integer                   :: iorb,jorb,Ni,Nj,unit
  unit=6;
  if(present(file))then
     unit=free_unit()
     open(unit,file=reg(file))
  endif
  Ni=size(Hloc,1)
  Nj=size(Hloc,2)
  if(present(file))then
     do iorb=1,Ni
        write(unit,"(9000F12.6)")(dreal(Hloc(iorb,jorb)),jorb=1,Nj)
     enddo
     write(unit,*)""
     do iorb=1,Ni
        write(unit,"(9000F12.6)")(dimag(Hloc(iorb,jorb)),jorb=1,Nj)
     enddo
     write(unit,*)""
     close(unit)
  else
     do iorb=1,Ni
        write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
             ('(',dreal(Hloc(iorb,jorb)),',',dimag(Hloc(iorb,jorb)),')',jorb =1,Nj)
     enddo
     write(unit,*)""
  endif
end subroutine write_hloc_1

subroutine write_hloc_2(hloc,file) ![Nspin][Nspin][Norb][Norb]
  complex(8),dimension(:,:,:,:) :: Hloc
  character(len=*),optional     :: file
  integer                       :: iorb,jorb,ispin,jspin
  integer                       :: Norb,Nspin,unit
  unit=6;
  if(present(file))then
     unit=free_unit()
     open(unit,file=reg(file))
  endif
  Nspin=size(Hloc,1)
  Norb=size(Hloc,3)
  call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"Write_Hloc_2","Hloc")
  if(present(file))then
     do ispin=1,Nspin
        do iorb=1,Norb
           write(unit,"(9000F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
        enddo
     enddo
     write(unit,*)""
     do ispin=1,Nspin
        do iorb=1,Norb
           write(unit,"(9000F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
        enddo
     enddo
     write(unit,*)""
     close(unit)
  else
     do ispin=1,Nspin
        do iorb=1,Norb
           write(unit,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
                (&
                (&
                '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
                jorb =1,Norb),&
                jspin=1,Nspin)
        enddo
     enddo
     write(unit,*)""
  endif
end subroutine write_hloc_2




!##################################################################
!##################################################################
!##################################################################
!##################################################################




subroutine read_hloc_1(hloc,file)
  complex(8),dimension(:,:)                    :: Hloc
  character(len=*)                             :: file
  integer                                      :: iorb,jorb,Ni,Nj,unit
  real(8),dimension(size(Hloc,1),size(Hloc,2)) :: reHloc,imHloc
  unit=free_unit()   
  open(unit,file=reg(file))
  Ni=size(Hloc,1)
  Nj=size(Hloc,2)
  do iorb=1,Ni
     read(unit,"(9000F12.6)")(reHloc(iorb,jorb),jorb=1,Nj)
  enddo
  read(unit,*)
  do iorb=1,Ni
     read(unit,"(9000F12.6)")(imHloc(iorb,jorb),jorb=1,Nj)
  enddo
  close(unit)
  Hloc = dcmplx(reHloc,imHloc)
end subroutine read_hloc_1

subroutine read_hloc_2(hloc,file)
  complex(8),dimension(:,:,:,:)                                          :: Hloc
  character(len=*)                                                       :: file
  integer                                                                :: iorb,jorb,ispin,jspin,unit
  integer                                                                :: Nspin,Norb
  real(8),dimension(size(Hloc,1),size(Hloc,2),size(Hloc,3),size(Hloc,4)) :: reHloc,imHloc
  unit=free_unit()   
  open(unit,file=reg(file))
  Nspin=size(Hloc,1)
  if(size(Hloc,2)/=Nspin)stop "read_hloc error: size[Hloc,2] != size[Hloc,1] == Nspin"
  Norb=size(Hloc,3)
  if(size(Hloc,4)/=Nspin)stop "read_hloc error: size[Hloc,4] != size[Hloc,3] == Norb"
  do ispin=1,Nspin
     do iorb=1,Norb
        read(unit,"(9000F12.6)")((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
     enddo
  enddo
  read(unit,*)
  do ispin=1,Nspin
     do iorb=1,Norb
        read(unit,"(9000F12.6)")((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
     enddo
  enddo
  close(unit)
  Hloc = dcmplx(reHloc,imHloc)
end subroutine read_hloc_2




