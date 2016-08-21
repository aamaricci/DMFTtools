module DMFT_MISC
  USE SF_TIMER
  USE SF_CONSTANTS, only: pi,zero
  implicit none
  private


  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block


  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape


  !ARRAY RESHAPE for DMFT:
  public :: blocks_to_matrix
  public :: matrix_to_blocks
  public :: select_block
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape


  !LOOP:
  public :: start_loop
  public :: end_loop

  !SORT 2D:
  public :: find2Dmesh

  !OTHER:
  ! public :: get_local_density
  public :: order_of_magnitude
  public :: get_density_from_matsubara_gf
  public :: get_matsubara_gf_from_dos
  public :: get_free_dos
  public :: sum_overk_zeta
  public :: finalize_run


contains


  subroutine start_loop(loop,max,name,unit,id)
    integer                   :: loop
    integer,optional          :: max,unit,id
    character(len=*),optional :: name
    character(len=16)         :: loop_name
    integer                   :: unit_,id_
    loop_name="main-loop";if(present(name))loop_name=name
    unit_    =6          ;if(present(unit))unit_=unit
    id_      =0          ;if(present(id))id_=id
    write(unit_,*)
    if(.not.present(max))then
       write(unit_,"(A,I5)")"-----"//trim(adjustl(trim(loop_name))),loop,"-----"
    else
       write(unit_,"(A,I5,A,I5,A)")"-----"//trim(adjustl(trim(loop_name))),loop,&
            " (max:",max,")-----"
    endif
    call start_timer
  end subroutine start_loop


  subroutine end_loop(unit,id)
    integer,optional :: unit,id
    integer          :: unit_,id_
    unit_=6 ; if(present(unit))unit_=unit
    id_  =0 ; if(present(id))id_=id
    write(unit_,"(A)")"====================================="
    call stop_timer
    write(unit_,*)
    write(unit_,*)
  end subroutine end_loop


  subroutine finalize_run(iter,error,filename)
    integer,intent(in)                   :: iter
    real(8),intent(in)                   :: error
    character(len=*),intent(in),optional :: filename
    character(len=100)                   :: filename_
    integer                              :: unit
    filename_="job_done.out";if(present(filename))filename_=filename
    unit=897
    open(unit,file=trim(filename),status="new")
    write(unit,*)iter,error
    close(unit)
  end subroutine finalize_run



  function order_of_magnitude(x) result(norder)
    integer :: norder
    real(8) :: x
    norder = floor(log10(abs(x)))
  end function order_of_magnitude



  ! function get_local_density(giw,beta) result(n)
  !   complex(8),dimension(:) :: giw
  !   real(8)                 :: gtau(0:size(giw))
  !   real(8)                 :: beta,n
  !   call fftgf_iw2tau(giw,gtau,beta)
  !   n = -2.d0*gtau(size(giw))
  ! end function get_local_density



  function get_density_from_matsubara_gf(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    complex(8)              :: tail
    real(8)                 :: sum
    real(8)                 :: n,wmax,beta,mues,At,w
    integer                 :: i,Liw
    Liw=size(giw)
    wmax = pi/beta*real(2*Liw-1,8)
    mues =-dreal(giw(Liw))*wmax**2
    sum=0.d0
    do i=1,Liw
       w=pi/beta*real(2*i-1,8)
       tail=-cmplx(mues,w,8)/(mues**2+w**2)
       sum=sum + dreal(giw(i)-tail)
    enddo
    At = -1.d0/(1.d0 + exp(-beta*mues))
    if((mues*beta) >  30.d0)At = -1.d0
    if((mues*beta) < -30.d0)At = -exp(mues*beta)
    n=sum*2.d0/beta+At+1.d0
  end function get_density_from_matsubara_gf


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate sum over k-points (very simple method)
  !+-------------------------------------------------------------------+
  function sum_overk_zeta(zeta,ek,wk) result(fg)
    complex(8)                    :: zeta,fg
    real(8),dimension(:)          :: ek
    real(8),dimension(:),optional :: wk
    real(8),dimension(size(ek))   :: wk_
    wk_ = 1.d0/size(ek);if(present(wk))wk_=wk
    fg=sum(wk(:)/(zeta-ek(:)))
  end function sum_overk_zeta


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine get_free_dos(ek,wk,dos,file,store,wmin,wmax,eps)
    real(8),dimension(:)          :: ek,wk
    real(8),dimension(:),optional :: dos
    character(len=*),optional     :: file
    logical,optional              :: store
    real(8),optional              :: wmin,wmax,eps
    character(len=32)             :: file_     
    logical                       :: store_
    real(8)                       :: wini,wfin,eta
    integer,parameter             :: M=2048
    integer                       :: i,ik,Lk,L
    real(8)                       :: w,dew
    complex(8)                    :: gf,iw
    file_="DOSfree.lattice" ; if(present(file))file_=file
    L=M           ; if(present(dos))L=size(dos)
    store_=.true. ; if(present(store))store_=store
    wini = -10.d0 ; if(present(wmin))wini=wmin
    wfin =  10.d0 ; if(present(wmax))wfin=wmax
    eta  = 0.01d0 ; if(present(eps))eta=eps
    Lk =size(ek)  ; dew=abs(wfin-wini)/real(L,8)
    open(70,file=trim(adjustl(trim(file_))))
    do i=1,L
       w  = wini + dble(i-1)*dew;iw=cmplx(w,eta)
       gf = sum(wk/(iw-ek))
       if(present(dos))dos(i)=-aimag(gf)/pi
       write(70,*)w,-aimag(gf)/pi
    enddo
    close(70)
  end subroutine get_free_dos




  !+-----------------------------------------------------------------+
  !PURPOSE  : calculate the Matsubara GF given the spectral density  
  !G^M(\iw)= int_\RR d\e A(\e)/iw+\e=-\iw\int_\RR d\e A(\e)/w^2+\e^2 
  !+-----------------------------------------------------------------+
  subroutine get_matsubara_gf_from_dos(win,gin,gout,beta)
    implicit none
    integer :: i,j,Nin,Nout
    real(8) :: w,wm,beta,A,wini,dw
    real(8) :: gmats_re,gmats_im
    real(8),dimension(:)    :: win
    complex(8),dimension(:) :: gin
    complex(8),dimension(:) :: gout
    complex(8),dimension(:),allocatable :: dummy_out
    wini=minval(win)
    dw=abs(win(2)-win(1)) !Assume constant step in x-grid
    Nin=size(gin)
    Nout=size(gout)
    allocate(dummy_out(4*Nout))

    do j=1,Nout
       wm=pi/beta*real(2*j-1,8)
       gmats_re=zero;gmats_im=zero
       do i=1,Nin
          w=wini+dble(i)*dw
          A=aimag(gin(i))/pi !DOS
          gmats_im=gmats_im + dw*A*wm/(wm**2+w**2)
          gmats_re=gmats_re + dw*A*w/(wm**2+w**2)
       enddo
       dummy_out(j)=cmplx(gmats_re,gmats_im)
    enddo
    gout = dummy_out(1:Nout)
  end subroutine get_matsubara_gf_from_dos







  ! SORTING 2D:
  !###################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine find2Dmesh(gridX,gridY,xin,iout)
    implicit none
    integer :: i,j,k
    real(8) :: dist
    integer,dimension(2) :: iout
    real(8),dimension(2) :: xin
    real(8),dimension(:) :: gridX,gridY
    !Find the index of the nearest point along X-axis
    iout(1)=fastsearchreal(xin(1),gridX(:))
    !Find the index of the nearest point along Y-axis
    iout(2)=fastsearchreal(xin(2),gridY(:))

    !Further checks on X:
    k=min(iout(1)+1,size(gridX))
    dist=abs(xin(1)-gridX(iout(1)) )
    if(abs(xin(1)-gridX(k)) < dist ) iout(1)=k
    k=max(iout(1)-1,1)
    if(abs(xin(1)-gridX(k)) < dist ) iout(1)=k

    !Further checks on Y:
    k=min(iout(2)+1,size(gridY))
    dist=abs(xin(2)-gridY(iout(2)))
    if(abs(xin(2)-gridY(k)) < dist ) iout(2)=k
    k=max(iout(2)-1,1)
    if(abs(xin(2)-gridY(1)) < dist ) iout(2)=k
    return
  end subroutine find2Dmesh
  !---------------------------------------------!
  function fastsearchreal(xx,tab)
    integer :: i1,i2,is,siz
    real(8) :: xx,tab(:)
    integer :: fastsearchreal
    siz=size(tab)
    is=siz/2
    i1=1
    i2=siz
    fastsearchreal=0
    if(tab(1)>xx)then
       fastsearchreal=siz
       return
    endif
    if(tab(siz)<=xx)then
       fastsearchreal=siz
       return
    endif
    do
       if(tab(is)<=xx) then
          i1=is
          is=i1+max(1,(i2-i1)/2)
          goto 28
       endif
       if(tab(is)>xx)then
          i2=is
          is=i1+(i2-i1)/2
          goto 28
       endif
28     continue
       if(is==siz.and.tab(is)<=xx)then
          fastsearchreal=is
          return
       endif
       if(is+1<=siz)then
          if(tab(is)<=xx.and.tab(is+1)>xx)then
             fastsearchreal=is
             return
          endif
       endif
    enddo
  end function fastsearchreal





  !####################################################################
  !                    computational ROUTINES
  !####################################################################
  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix,Nlat,Nso) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks






  !--------------------------------------------------------------------!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !--------------------------------------------------------------------!
  function select_block_Nlso(ip,Matrix,Nlat,Nso) result(Vblock)
    integer                                 :: ip
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nso,Nso)           :: Vblock
    integer                                 :: i,j
    Vblock=zero
    i = 1+(ip-1)*Nso
    j =       ip*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ip,Matrix,Nlat,Nspin,Norb) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block_nnn








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso




end module DMFT_MISC
