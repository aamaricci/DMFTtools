module TB_EFERMI
  USE TB_COMMON
  USE DMFT_CTRL_VARS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none



contains



  subroutine TB_FermiLevel(Hk,filling,Ef,Nspin,verbose)
    complex(8),dimension(:,:,:)                 :: Hk
    real(8)                                     :: filling
    real(8)                                     :: Ef,Dmin,Dmax
    integer,optional                            :: Nspin
    logical,optional                            :: verbose
    complex(8),dimension(:,:),allocatable       :: Uvec
    real(8),dimension(:),allocatable            :: Eval
    complex(8),dimension(:,:,:),allocatable     :: Rk,Rk_tmp
    real(8),dimension(:,:),allocatable          :: Ek,Ek_tmp
    integer                                     :: Nk,Nso,ik,io,info,Nspin_
    logical                                     :: verbose_
    !
    verbose_ = .true. ;if(present(verbose))verbose_=verbose
    Nspin_   = 1      ;if(present(Nspin))Nspin_=Nspin
    !
    Nso = size(Hk,1)
    Nk = size(Hk,3)
    call assert_shape(Hk,[Nso,Nso,Nk],"TB_FermiLevel","Hk")
    !
    !MPI setup:
#ifdef _MPI 
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    allocate(Uvec(Nso,Nso))
    allocate(Eval(Nso))
    allocate(Rk_tmp(Nso,Nso,Nk),Rk(Nso,Nso,Nk))
    allocate(Ek_tmp(Nso,Nk),Ek(Nso,Nk))
    Rk_tmp = zero
    Ek_tmp = 0d0
    do ik=1+mpi_rank,Nk,mpi_size
       Uvec = Hk(:,:,ik)
       call eigh(Uvec,Eval)
       Rk_tmp(:,:,ik) = Uvec
       Ek_tmp(:,ik)   = Eval
    enddo
    !
    !
#ifdef _MPI
    if(check_MPI())then
       Rk = zero ; Ek = 0d0
       call AllReduce_MPI(MPI_COMM_WORLD,Ek_tmp,Ek)
       call AllReduce_MPI(MPI_COMM_WORLD,Rk_tmp,Rk)
       deallocate(Ek_tmp,Rk_tmp)
    else
       Rk = Rk_tmp
       Ek = Ek_tmp
    endif
#else
    Rk = Rk_tmp
    Ek = Ek_tmp
#endif
    !
    !
    Dmin = minval(Ek)
    Dmax = maxval(Ek)
    Ef   = Dmin
    if(mpi_master.AND.verbose_)write(*,*)"D_min,D_max:",Dmin,Dmax
    call fzero(get_dens,Ef,Dmax,info,rguess=Dmin+0.5d0*(Dmax-Dmin))
    if(info/=1)then
       if(mpi_master)write(*,*)"ERROR TB_get_Fermi: fzero returned info>1 ",info
       stop
    endif
    !
    if(mpi_master)write(*,*)"Fermi Level: ",Ef
    !
  contains
    !
    function get_dens(ef) result(ndens)
      real(8),intent(in)      :: ef
      real(8)                 :: dens,dens_tmp
      real(8)                 :: ndens
      real(8),dimension(Nso)  :: rtmp,wt
      dens_tmp = 0d0
      do ik=1+mpi_rank,Nk,mpi_size
         do io=1,Nso
            Rtmp(io) = sum( abs(Rk(io,:,ik))**2*fermi(Ek(:,ik)-ef,1000d0))
         enddo
         dens_tmp= dens_tmp + sum(Rtmp)/Nk
      enddo
#ifdef _MPI
      if(check_MPI())then
         dens = 0d0
         call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens)
      else
         dens = dens_tmp
      endif
#else
      dens = dens_tmp
#endif
      dens  = (3-Nspin_)*dens
      ndens = dens-filling
      if(mpi_master.AND.verbose_)write(*,"(A9,3G18.9)")"Ef,N0   =",ef,dens,filling
    end function get_dens
    !
  end subroutine TB_FermiLevel


END MODULE TB_EFERMI








  !   subroutine TB_FermiLevel(Hk,filling,Ef,nspin,verbose)
  !     complex(8),dimension(:,:,:)           :: Hk
  !     real(8)                               :: filling
  !     real(8)                               :: Ef
  !     integer,optional                            :: Nspin
  !     logical,optional                            :: verbose
  !     complex(8),dimension(:,:),allocatable :: Uk
  !     real(8),dimension(:),allocatable      :: Ek
  !     real(8),dimension(:),allocatable      :: Ek_all,Ek_tmp
  !     integer,dimension(:),allocatable      :: Ek_indx
  !     integer                               :: stride,Nk,Nso,ik,indx
  !     Nk = size(Hk,3)
  !     Nso = size(Hk,1)
  !     call assert_shape(Hk,[Nso,Nso,Nk],"TB_FermiLevel","Hk")
  !     allocate(Uk(Nso,Nso),Ek(Nso),Ek_all(Nso*Nk),Ek_indx(Nk*Nso),Ek_tmp(Nk*Nso))
  !     !MPI setup:
  ! #ifdef _MPI    
  !     if(check_MPI())then
  !        mpi_size  = get_size_MPI()
  !        mpi_rank =  get_rank_MPI()
  !        mpi_master= get_master_MPI()
  !     else
  !        mpi_size=1
  !        mpi_rank=0
  !        mpi_master=.true.
  !     endif
  ! #else
  !     mpi_size=1
  !     mpi_rank=0
  !     mpi_master=.true.
  ! #endif
  !     Ek_tmp = 0d0
  !     stride = 0
  !     do ik=1+mpi_rank,Nk,mpi_size
  !        stride = (ik-1)*Nso
  !        Uk = Hk(:,:,ik)
  !        call eigh(Uk,Ek)
  !        Ek_tmp(stride+1:stride+Nso) = Ek
  !     enddo
  ! #ifdef _MPI
  !     if(check_MPI())then
  !        Ek_all = 0d0
  !        call AllReduce_MPI(MPI_COMM_WORLD,Ek_tmp,Ek_all)
  !     else
  !        Ek_all = Ek_tmp
  !     endif
  ! #else
  !     Ek_all = Ek_tmp
  ! #endif
  !     call sort_array(Ek_all,Ek_indx)
  !     if(filling==0d0)filling=Nso
  !     indx  = ceiling(filling*Nk/2d0)
  !     Ef    = Ek_all(indx)
  !   end subroutine TB_FermiLevel











