module DMFT_GFIO
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS,   only: reg,txtfy,splot,file_gzip
  USE SF_MISC,      only: assert_shape
  implicit none
  private


  interface dmft_print_gf_matsubara
     module procedure :: dmft_gf_print_matsubara_main
     module procedure :: dmft_gf_print_matsubara_ineq
  end interface dmft_print_gf_matsubara


  interface dmft_print_gf_realaxis
     module procedure :: dmft_gf_print_realaxis_main
     module procedure :: dmft_gf_print_realaxis_ineq
  end interface dmft_print_gf_realaxis


  public :: set_gf_suffix
  !
  public :: dmft_print_gf_matsubara
  public :: dmft_print_gf_realaxis
  !
  public :: dmft_print_gij_matsubara
  public :: dmft_print_gij_realaxis


  character(len=128)                        :: suffix
  character(len=128)                        :: gf_suffix='.dat'
  integer                                   :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats
  integer                                   :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is,js
  !
  real(8)                                   :: beta
  real(8)                                   :: xmu,eps
  real(8)                                   :: wini,wfin 




contains



  subroutine set_gf_suffix(string)
    character(len=*) :: string
    gf_suffix=reg(string)
  end subroutine set_gf_suffix




  !PRINT GLOC SINGLE SITE
  subroutine dmft_gf_print_matsubara_main(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats
    real(8),dimension(size(Gmats,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    Nspin = size(Gmats,1)
    Norb  = size(Gmats,3)
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,size(Gmats,5)],"dmft_gloc_print_matsubara",reg(fname)//"_mats")
    !
    select case(iprint)
    case(1)                  !print only diagonal elements
       write(*,"(A)") "Gloc matsubara: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gf_suffix)
             call splot(reg(suffix),w,Gmats(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,"(A)") "Gloc matsubara: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gf_suffix)
                call splot(reg(suffix),w,Gmats(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default                  !print all off-diagonals
       write(*,"(A)") "Gloc matsubara: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gf_suffix)
                   call splot(reg(suffix),w,Gmats(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gf_print_matsubara_main


  subroutine dmft_gf_print_matsubara_ineq(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Gmats
    real(8),dimension(size(Gmats,6))             :: w        
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,2)
    Norb  = size(Gmats,4)
    call assert_shape(Gmats,[Nlat,Nspin,Nspin,Norb,Norb,size(Gmats,6)],"dmft_gloc_print_matsubara_ineq",reg(fname)//"_mats")
    !
    select case(iprint)
    case(1)                  !print only diagonal elements
       write(*,*)"Gloc matsubara: write spin-orbital diagonal elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gf_suffix)
             call splot(reg(suffix),w,Gmats(:,ispin,ispin,iorb,iorb,:))
             call file_gzip(reg(suffix))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc matsubara: write spin diagonal and all orbitals elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gf_suffix)
                call splot(reg(suffix),w,Gmats(:,ispin,ispin,iorb,jorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       enddo
       !
    case(3)                  !print all off-diagonals
       write(*,*)"Gloc matsubara: write all elements. No Split."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gf_suffix)
                   call splot(reg(suffix),w,Gmats(:,ispin,jspin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       enddo
       !
    case(4)                  !print only diagonal elements
       write(*,*)"Gloc matsubara: write spin-orbital diagonal elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw_indx"//reg(txtfy(ilat,6))//reg(gf_suffix)
                call splot(reg(suffix),w,Gmats(ilat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
       !
    case(5)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc matsubara: write spin diagonal and all orbitals elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//&
                        "_iw_indx"//reg(txtfy(ilat,6))//reg(gf_suffix)
                   call splot(reg(suffix),w,Gmats(ilat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    case default
       write(*,*)"Gloc matsubara: write all elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                           "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                           "_iw_indx"//reg(txtfy(ilat,6))//reg(gf_suffix)
                      call splot(reg(suffix),w,Gmats(ilat,ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gf_print_matsubara_ineq







  subroutine dmft_gf_print_realaxis_main(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Greal
    real(8),dimension(size(Greal,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    Nspin = size(Greal,1)
    Norb  = size(Greal,3)
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,size(Greal,5)],"dmft_gloc_print_realaxis",reg(fname)//"_real")
    !
    select case(iprint)
    case(1)                  !print only diagonal elements
       write(*,"(A)") "Gloc real: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gf_suffix)
             call splot(reg(suffix),w,Greal(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,"(A)") "Gloc real: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gf_suffix)
                call splot(reg(suffix),w,Greal(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default                  !print all off-diagonals
       write(*,"(A)") "Gloc real: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gf_suffix)
                   call splot(reg(suffix),w,Greal(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gf_print_realaxis_main

  subroutine dmft_gf_print_realaxis_ineq(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:),intent(in) :: Greal
    real(8),dimension(size(Greal,6))             :: w
    character(len=*),intent(in)                  :: fname
    integer,intent(in)                           :: iprint
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,2)
    Norb  = size(Greal,4)
    call assert_shape(Greal,[Nlat,Nspin,Nspin,Norb,Norb,size(Greal,6)],"dmft_gloc_print_realaxis_ineq",reg(fname)//"_real")
    !
    select case(iprint)
    case(1)                  !print only diagonal elements
       write(*,*)"Gloc real: write spin-orbital diagonal elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gf_suffix)
             call splot(reg(suffix),w,Greal(:,ispin,ispin,iorb,iorb,:))
             call file_gzip(reg(suffix))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc real: write spin diagonal and all orbitals elements. No Split."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gf_suffix)
                call splot(reg(suffix),w,Greal(:,ispin,ispin,iorb,jorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       enddo
       !
    case(3)
       write(*,*)"Gloc real: write all elements. No Split."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gf_suffix)
                   call splot(reg(suffix),w,Greal(:,ispin,jspin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       enddo
       !
    case(4)                  !print only diagonal elements
       write(*,*)"Gloc real: write spin-orbital diagonal elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw_indx"//reg(txtfy(ilat,6))//reg(gf_suffix)
                call splot(reg(suffix),w,Greal(ilat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
       !
    case(5)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc real: write spin diagonal and all orbitals elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//&
                        "_realw_indx"//reg(txtfy(ilat,6))//reg(gf_suffix)
                   call splot(reg(suffix),w,Greal(ilat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    case default
       write(*,*)"Gloc real: write all elements. Split."
       do ilat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      suffix=reg(fname)//&
                           "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                           "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                           "_realw_indx"//reg(txtfy(ilat,6))//reg(gf_suffix)
                      call splot(reg(suffix),w,Greal(ilat,ispin,jspin,iorb,jorb,:))
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gf_print_realaxis_ineq







  !PRINT GLOC FULL LATTICE
  subroutine dmft_print_gij_matsubara(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Gmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats/Lreal]
    real(8),dimension(size(Gmats,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    !
    Nlat  = size(Gmats,1)
    Nspin = size(Gmats,3)
    Norb  = size(Gmats,5)
    call assert_shape(Gmats,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Gmats,7)],"dmft_gloc_print_matsubara_gij",reg(fname)//"_mats")
    !
    select case(iprint)
    case (0)
       write(*,*)"Gloc matsubara: not written on file."
    case(1)                  !print only diagonal elements     
       write(*,*)"Gloc matsubara: spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gf_suffix)
             call splot(reg(suffix),w,Gmats(:,:,ispin,ispin,iorb,iorb,:))
             call file_gzip(reg(suffix))
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc matsubara: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gf_suffix)
                call splot(reg(suffix),w,Gmats(:,:,ispin,ispin,iorb,jorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       enddo
    case default              !print all off-diagonals
       write(*,*)"Gloc matsubara: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gf_suffix)
                   call splot(reg(suffix),w,Gmats(:,:,ispin,jspin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine dmft_print_gij_matsubara

  subroutine dmft_print_gij_realaxis(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:,:,:),intent(in) :: Greal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal/Lreal]
    real(8),dimension(size(Greal,7))               :: w
    character(len=*),intent(in)                    :: fname
    integer,intent(in)                             :: iprint
    integer                                        :: Nlat,Nspin,Norb,ispin,jspin,iorb,jorb,ilat,jlat
    !
    Nlat  = size(Greal,1)
    Nspin = size(Greal,3)
    Norb  = size(Greal,5)
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Greal,7)],"dmft_gloc_print_realaxis_gij",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(*,*)"Gloc real: not written on file."
    case(1)                  !print only diagonal elements     
       write(*,*)"Gloc real: spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gf_suffix)
             call splot(reg(suffix),w,Greal(:,:,ispin,ispin,iorb,iorb,:))
             call file_gzip(reg(suffix))
          enddo
       enddo
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,*)"Gloc real: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gf_suffix)
                call splot(reg(suffix),w,Greal(:,:,ispin,ispin,iorb,jorb,:))
                call file_gzip(reg(suffix))
             enddo
          enddo
       enddo
    case default              !print all off-diagonals
       write(*,*)"Gloc real: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gf_suffix)
                   call splot(reg(suffix),w,Greal(:,:,ispin,jspin,iorb,jorb,:))
                   call file_gzip(reg(suffix))
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine dmft_print_gij_realaxis







end module DMFT_GFIO
