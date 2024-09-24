program test
  USE SCIFOR
  USE DMFT_CTRL_VARS
  USE DMFT_GF
  USE ASSERTING
  implicit none
  integer,parameter                                  :: Nspin = 2
  integer,parameter                                  :: Norb  = 2
  integer,parameter                                  :: Nx    = 6, Nlat=Nx*Nx
  integer,parameter                                  :: Nkx   = 6, Nk=Nkx*Nkx
  integer,parameter                                  :: Nso   = Nspin*Norb
  integer,parameter                                  :: Nlso  = Nlat*Nso
  integer,parameter                                  :: L     = 512
  integer,parameter                                  :: Le    = 5000
  real(8)                                            :: xmu   = 0d0
  real(8)                                            :: beta  = 100d0
  real(8)                                            :: wini  = -5d0
  real(8)                                            :: wfin  = 5d0
  real(8)                                            :: eps   = 0.02d0  
  real(8)                                            :: Mh    = 1d0
  real(8)                                            :: lambda= 0.3d0
  real(8)                                            :: delta = 0.3d0
  complex(8),dimension(Nso,Nso,L)                    :: S_rank2=zero
  complex(8),dimension(Nlat,Nso,Nso,L)               :: S_rank3=zero
  complex(8),dimension(Nspin,Nspin,Norb,Norb,L)      :: S_rank4=zero
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L) :: S_rank5=zero
  complex(8),dimension(Nso,Nso,L)                    :: mG0_rank2,mG_rank2
  complex(8),dimension(Nlat,Nso,Nso,L)               :: mG0_rank3,mG_rank3
  complex(8),dimension(Nspin,Nspin,Norb,Norb,L)      :: mG0_rank4,mG_rank4
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L) :: mG0_rank5,mG_rank5
  complex(8),dimension(Nso,Nso,L)                    :: rG0_rank2,rG_rank2
  complex(8),dimension(Nlat,Nso,Nso,L)               :: rG0_rank3,rG_rank3
  complex(8),dimension(Nspin,Nspin,Norb,Norb,L)      :: rG0_rank4,rG_rank4
  complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L) :: rG0_rank5,rG_rank5
  !
  complex(8),dimension(2,Norb,Norb,L)                :: scS_rank2=zero
  complex(8),dimension(2,Nlat,Norb,Norb,L)           :: scS_rank3=zero
  complex(8),dimension(2,1,1,Norb,Norb,L)            :: scS_rank4=zero
  complex(8),dimension(2,Nlat,1,1,Norb,Norb,L)       :: scS_rank5=zero
  complex(8),dimension(2,Norb,Norb,L)                :: scmG0_rank2,scmG_rank2
  complex(8),dimension(2,Nlat,Norb,Norb,L)           :: scmG0_rank3,scmG_rank3
  complex(8),dimension(2,1,1,Norb,Norb,L)            :: scmG0_rank4,scmG_rank4
  complex(8),dimension(2,Nlat,1,1,Norb,Norb,L)       :: scmG0_rank5,scmG_rank5
  complex(8),dimension(2,Norb,Norb,L)                :: scrG0_rank2,scrG_rank2
  complex(8),dimension(2,Nlat,Norb,Norb,L)           :: scrG0_rank3,scrG_rank3
  complex(8),dimension(2,1,1,Norb,Norb,L)            :: scrG0_rank4,scrG_rank4
  complex(8),dimension(2,Nlat,1,1,Norb,Norb,L)       :: scrG0_rank5,scrG_rank5
  integer,dimension(:,:),allocatable                 :: Links
  complex(8),dimension(:,:,:,:),allocatable          :: Hlat
  complex(8),dimension(:,:,:),allocatable            :: Hk,Hij
  complex(8),dimension(:,:,:,:),allocatable          :: scHk,scHij
  complex(8),dimension(:,:),allocatable              :: Hloc
  real(8),dimension(Nso,Le)                          :: Ebands,Dbands
  real(8),dimension(Nso)                             :: Wband
  real(8),dimension(Nso)                             :: De,H0
  complex(8),dimension(Nso,Nso)                      :: Gamma1,Gamma2,Gamma5
  real(8),dimension(3)                               :: ei_1=[1d0,0d0,0d0],ei_2=[0d0,1d0,0d0],ei_3=[0d0,0d0,1d0]
  real(8),dimension(3)                               :: bk_1=[pi2,0d0,0d0],bk_2=[0d0,pi2,0d0],bk_3=[0d0,0d0,pi2]
  integer                                            :: i,j,iso,jso,ii,jj
  integer                                            :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo
  complex(8)                                         :: zeta
  real(8),dimension(L)                               :: wm,wr



  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")


  gamma1=kron( pauli_tau_z, pauli_sigma_x)
  gamma2=kron( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron( pauli_tau_0, pauli_sigma_z)



  !######################################################################
  !NORMAL
  !######################################################################
  !Setup BHZ model 2d,
  ! H(k)
  allocate(Hk(Nso,Nso,Nk))
  allocate(Hloc(Nso,Nso))
  call TB_build_Hk(Hk,hk_bhz,Nso,[Nkx,Nkx])
  Hloc = sum(Hk,3)/Nk;where(abs(Hloc)<1d-4)Hloc=zero
  call print_matrix(Hloc)
  ! H(ij)
  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  allocate(Links(4,2))
  Links(1,:) = [1 ,0]
  Links(2,:) = [0 ,1]
  Links(3,:) = [-1,0]
  Links(4,:) = [0,-1]
  call TB_build_Hlat(Hlat,ts_bhz,Nso,[Nx,Nx],Links,pbc=.true.)
  allocate(Hij(Nlso,Nlso,1))
  do concurrent(i=1:Nlat,j=1:Nlat,iso=1:Nso,jso=1:Nso)
     ii = iso + (i-1)*Nso
     jj = jso + (j-1)*Nso
     Hij(ii,jj,1) = Hlat(iso,jso,i,j)
  enddo
  !
  !Get exact Green's functions:
  call get_gf_hk(Hk,mG0_rank2,'mats')
  call get_gf_hk(Hk,rG0_rank2,'real')
  call get_gf_hlat(Hij(:,:,1),mG0_rank3,'mats')
  call get_gf_hlat(Hij(:,:,1),rG0_rank3,'real')
  do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
     io = iorb+(ispin-1)*Norb
     jo = jorb+(jspin-1)*Norb
     mG0_rank4(ispin,jspin,iorb,jorb,:)   = mG0_rank2(io,jo,:)
     mG0_rank5(:,ispin,jspin,iorb,jorb,:) = mG0_rank3(:,io,jo,:)
     rG0_rank4(ispin,jspin,iorb,jorb,:)   = rG0_rank2(io,jo,:)
     rG0_rank5(:,ispin,jspin,iorb,jorb,:) = rG0_rank3(:,io,jo,:)
  enddo
  !
  !
  !> CHECK: Build all ranks  GF, Mats and Real
  call get_gloc(Hk, mG_rank2,S_rank2,axis='m')
  call get_gloc(Hk, mG_rank4,S_rank4,axis='m')
  call get_gloc(Hij,mG_rank3,S_rank3,axis='m')
  call get_gloc(Hij,mG_rank5,S_rank5,axis='m')
  !
  call get_gloc(Hk, rG_rank2,S_rank2,axis='r')
  call get_gloc(Hk, rG_rank4,S_rank4,axis='r')
  call get_gloc(Hij,rG_rank3,S_rank3,axis='r')
  call get_gloc(Hij,rG_rank5,S_rank5,axis='r')
  !
  !> ASSERT EQUALITY WITH EXACT ONES:
  call assert(mG_rank2,mG0_rank2,"Matsubara H(k): get_gloc_rank2",tol=1d-12)
  call assert(mG_rank3,mG0_rank3,"Matsubara H(k): get_gloc_rank3",tol=1d-12)
  call assert(mG_rank4,mG0_rank4,"Matsubara H(k): get_gloc_rank4",tol=1d-12)
  call assert(mG_rank5,mG0_rank5,"Matsubara H(k): get_gloc_rank5",tol=1d-12)
  call assert(rG_rank2,rG0_rank2,"Realaxis  H(k): get_gloc_rank2",tol=1d-12)
  call assert(rG_rank3,rG0_rank3,"Realaxis  H(k): get_gloc_rank3",tol=1d-12)
  call assert(rG_rank4,rG0_rank4,"Realaxis  H(k): get_gloc_rank4",tol=1d-12)
  call assert(rG_rank5,rG0_rank5,"Realaxis  H(k): get_gloc_rank5",tol=1d-12)
  call wait(500)
  

  ! DOS:
  Wband = [1d0,0.5d0,1d0,0.5d0]
  H0    = [0d0,0d0,0d0,0d0]
  do io=1,Nso
     Ebands(io,:) = linspace(-Wband(io),Wband(io),Le,mesh=de(io))
     Dbands(io,:) = dens_bethe(Ebands(io,:),Wband(io))*de(io)
  enddo
  !
  !> Get exact Green's functions 
  wm = pi/beta*(2*arange(1,L)-1)
  wr = linspace(wini,wfin,L)
  mG0_rank2=zero
  S_rank2 =zero
  do io=1,Nso
     do i=1,L
        zeta = xi*wm(i) - S_rank2(io,io,i) - H0(io)
        mG0_rank2(io,io,i) = gfbethe(wm(i),zeta,Wband(io))
        zeta = dcmplx(wr(i),eps) - S_rank2(io,io,i) - H0(io)
        rG0_rank2(io,io,i) = gfbether(wr(i),zeta,Wband(io))
     enddo
  enddo
  do concurrent(ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
     io = iorb+(ispin-1)*Norb
     jo = jorb+(jspin-1)*Norb
     mG0_rank4(ispin,jspin,iorb,jorb,:)   = mG0_rank2(io,jo,:)
     rG0_rank4(ispin,jspin,iorb,jorb,:)   = rG0_rank2(io,jo,:)
  enddo
  !
  !> CHECK: Build all ranks  GF, Mats and Real 
  call get_gloc(Ebands,Dbands,H0,mG_rank2,S_rank2,axis='m')
  call get_gloc(Ebands,Dbands,H0,mG_rank4,S_rank4,axis='m')
  call get_gloc(Ebands,Dbands,H0,rG_rank2,S_rank2,axis='r')
  call get_gloc(Ebands,Dbands,H0,rG_rank4,S_rank4,axis='r')
  !> ASSERT EQUALITY WITH EXACT ONES:
  call assert(mG_rank2,mG0_rank2,"Matsubara DOS: get_gloc_rank2",tol=1d-4)
  call assert(mG_rank4,mG0_rank4,"Matsubara DOS: get_gloc_rank4",tol=1d-4)
  call assert(rG_rank2,rG0_rank2,"Realaxis  DOS: get_gloc_rank2",tol=1d-2)
  call assert(rG_rank4,rG0_rank4,"Realaxis  DOS: get_gloc_rank4",tol=1d-2)
  call wait(500)
  !######################################################################
  !######################################################################



  deallocate(Hk,Hloc,Hlat,Links,Hij)


  !######################################################################
  !SUPERC
  !######################################################################
  !Setup 2bHM model 2d,
  !H(k)
  allocate(Hk(Norb,Norb,Nk))
  allocate(Hloc(Norb,Norb))
  call TB_build_Hk(Hk,hk_hm2b,Norb,[Nkx,Nkx])
  Hloc = sum(Hk,3)/Nk;where(abs(Hloc)<1d-4)Hloc=zero
  call print_matrix(Hloc)
  allocate(scHk(2,Norb,Norb,Nk));scHk=zero
  scHk(1,:,:,:) = Hk 
  scHk(2,:,:,:) =-conjg(Hk)
  !
  !H(ij)
  allocate(Hlat(Norb,Norb,Nlat,Nlat))
  allocate(Links(4,2))
  Links(1,:) = [1 ,0]
  Links(2,:) = [0 ,1]
  Links(3,:) = [-1,0]
  Links(4,:) = [0,-1]
  call TB_build_Hlat(Hlat,ts_hm2b,Norb,[Nx,Nx],Links,pbc=.true.)
  allocate(Hij(Nlat*Norb,Nlat*Norb,1))
  do concurrent(i=1:Nlat,j=1:Nlat,iso=1:Norb,jso=1:Norb)
     ii = iso + (i-1)*Norb
     jj = jso + (j-1)*Norb
     Hij(ii,jj,1) = Hlat(iso,jso,i,j)
  enddo
  allocate(scHij(2,Nlat*Norb,Nlat*Norb,1));scHij=zero
  scHij(1,:,:,1) = Hij(:,:,1)
  scHij(2,:,:,1) =-conjg(Hij(:,:,1))
  !
  !Set SC Sigma (BCS-like)
  scS_rank2 = zero; forall(io=1:Norb)scS_rank2(2,io,io,:) = -delta
  scS_rank3 = zero; forall(io=1:Norb)scS_rank3(2,:,io,io,:) = -delta
  scS_rank4 = zero; forall(io=1:Norb)scS_rank4(2,1,1,io,io,:) = -delta
  scS_rank5 = zero; forall(io=1:Norb)scS_rank5(2,:,1,1,io,io,:) = -delta
  !
  !Get exact BCS Green's functions:
  call get_gf_schk(Hk,scS_rank2,scmG0_rank2,'mats')
  call get_gf_schk(Hk,scS_rank2,scrG0_rank2,'real')
  call get_gf_schlat(Hij(:,:,1),scS_rank3,scmG0_rank3,'mats')
  call get_gf_schlat(Hij(:,:,1),scS_rank3,scrG0_rank3,'real')
  do concurrent(iorb=1:Norb,jorb=1:Norb)
     scmG0_rank4(:,1,1,iorb,jorb,:)   = scmG0_rank2(:,iorb,jorb,:)
     scrG0_rank4(:,1,1,iorb,jorb,:)   = scrG0_rank2(:,iorb,jorb,:)
     do ilat=1,Nlat
        scmG0_rank3(:,ilat,iorb,jorb,:)       = scmG0_rank3(:,ilat,iorb,jorb,:)
        scrG0_rank3(:,ilat,iorb,jorb,:)       = scrG0_rank3(:,ilat,iorb,jorb,:)
        scmG0_rank5(:,ilat,1,1,iorb,jorb,:)   = scmG0_rank3(:,ilat,iorb,jorb,:)
        scrG0_rank5(:,ilat,1,1,iorb,jorb,:)   = scrG0_rank3(:,ilat,iorb,jorb,:)
     enddo
  enddo
  !
  !> CHECK: Build all ranks Nambu GF, Mats and Real
  call get_gloc(scHk, scmG_rank2,scS_rank2,axis='m')
  call get_gloc(scHk, scmG_rank4,scS_rank4,axis='m')
  call get_gloc(scHij,scmG_rank3,scS_rank3,axis='m')
  call get_gloc(scHij,scmG_rank5,scS_rank5,axis='m')
  call get_gloc(scHk, scrG_rank2,scS_rank2,axis='r')
  call get_gloc(scHk, scrG_rank4,scS_rank4,axis='r')
  call get_gloc(scHij,scrG_rank3,scS_rank3,axis='r')
  call get_gloc(scHij,scrG_rank5,scS_rank5,axis='r')
  !
  !> ASSERT EQUALITY WITH EXACT ONES:
  call assert(scmG_rank2,scmG0_rank2,"SC Matsubara H(k): get_gloc_rank2",tol=1d-12)
  call assert(scmG_rank3,scmG0_rank3,"SC Matsubara H(k): get_gloc_rank3",tol=1d-12)
  call assert(scmG_rank4,scmG0_rank4,"SC Matsubara H(k): get_gloc_rank4",tol=1d-12)
  call assert(scmG_rank5,scmG0_rank5,"SC Matsubara H(k): get_gloc_rank5",tol=1d-12)
  call assert(scrG_rank2,scrG0_rank2,"SC Realaxis  H(k): get_gloc_rank2",tol=1d-12)
  call assert(scrG_rank3,scrG0_rank3,"SC Realaxis  H(k): get_gloc_rank3",tol=1d-12)
  call assert(scrG_rank4,scrG0_rank4,"SC Realaxis  H(k): get_gloc_rank4",tol=1d-12)
  call assert(scrG_rank5,scrG0_rank5,"SC Realaxis  H(k): get_gloc_rank5",tol=1d-12)
  call wait(500)


contains


  function hk_bhz(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
  end function hk_bhz

  function hk_hm2b(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=2)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = ek*pauli_tau_z  !+ Mh*pauli_tau_z
  end function hk_hm2b



  function ts_bhz(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*Gamma5
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    case (2) !UP HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    case (3) !LEFT HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    case (4) !DOWN HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    case default 
       stop "ts_model ERROR: link != {0,...,4}"
    end select
  end function ts_bhz

  function ts_hm2b(link,N) result(Hts)
    integer                       :: link
    integer                       :: N
    complex(8),dimension(N,N) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  0d0!Mh*pauli_tau_z
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*pauli_tau_z !+ xi*0.5d0*lambda*Gamma1
    case (2) !UP HOPPING
       Hts = -0.5d0*pauli_tau_z !+ xi*0.5d0*lambda*Gamma2
    case (3) !LEFT HOPPING
       Hts = -0.5d0*pauli_tau_z !- xi*0.5d0*lambda*Gamma1
    case (4) !DOWN HOPPING
       Hts = -0.5d0*pauli_tau_z !- xi*0.5d0*lambda*Gamma2
    case default 
       stop "ts_model ERROR: link != {0,...,4}"
    end select
  end function ts_hm2b


  subroutine get_gf_hk(Hk,Gloc,axis)
    complex(8),dimension(Nso,Nso,Nk) :: Hk
    complex(8),dimension(Nso,Nso,L)  :: Gloc
    character(len=*)                 :: axis
    complex(8),dimension(L)          :: wfreq
    complex(8),dimension(Nso,L)      :: csi
    complex(8),dimension(Nso,Nso)    :: U
    real(8),dimension(Nso)           :: E
    integer                          :: i,ilat,io,jo,ik
    !
    write(*,"(A)")"Get local Green's function Hk, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    Gloc = zero
    do ik=1,Nk
       U = Hk(:,:,ik)
       call eigh(U,E)
       forall(i=1:L)csi(:,i)=one/(wfreq(i)+xmu-E(:))
       do concurrent(io=1:Nso,jo=1:Nso,i=1:L)
          Gloc(io,jo,i) = Gloc(io,jo,i) + sum(U(io,:)*conjg(U(jo,:))*csi(:,i))/Nk!can use matmul
       enddo
    enddo
  end subroutine get_gf_hk

  subroutine get_gf_schk(Hk,S,Gloc,axis)
    complex(8),dimension(Norb,Norb,Nk)  :: Hk
    complex(8),dimension(2,Norb,Norb,L) :: S
    complex(8),dimension(2,Norb,Norb,L) :: Gloc
    character(len=*)                    :: axis
    complex(8),dimension(L)             :: wfreq
    complex(8),dimension(2*Norb,L)      :: csi
    complex(8),dimension(2*Norb,2*Norb) :: U
    real(8),dimension(2*Norb)           :: E
    integer                             :: i,j,io,jo,ik,is,js,N
    !
    write(*,"(A)")"Get local Green's function Hk, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    N = Norb
    Gloc = zero
    do ik=1,Nk
       U(1:N,1:N)         = Hk(:,:,ik)+S(1,:,:,1)
       U(1:N,N+1:2*N)     = S(2,:,:,1)
       U(N+1:2*N,1:N)     = S(2,:,:,1)
       U(N+1:2*N,N+1:2*N) =-conjg(Hk(:,:,ik)+S(1,:,:,L))
       call eigh(U,E)
       forall(i=1:L)csi(:,i)=one/(wfreq(i)+xmu-E(:))
       do concurrent(io=1:N,jo=1:N,i=1:L)
          Gloc(1,io,jo,i) = Gloc(1,io,jo,i) + sum(U(io,:)*conjg(U(jo,:))*csi(:,i))/Nk!can use matmul
          Gloc(2,io,jo,i) = Gloc(2,io,jo,i) + sum(U(io,:)*conjg(U(N+jo,:))*csi(:,i))/Nk!can use matmul
       enddo
    enddo
  end subroutine get_gf_schk



  subroutine get_gf_hlat(Hij,Gloc,axis)
    complex(8),dimension(Nlso,Nlso)      :: Hij
    complex(8),dimension(Nlat,Nso,Nso,L) :: Gloc
    character(len=*)                     :: axis
    complex(8),dimension(L)              :: wfreq
    complex(8),dimension(Nlso,L)         :: csi
    complex(8),dimension(Nlat,Nso,Nso,L) :: Gtmp
    complex(8),dimension(Nlso,Nlso)      :: U
    real(8),dimension(Nlso)              :: E
    integer                              :: i,ilat,io,jo,is,js
    !
    write(*,"(A)")"Get local Green's function Hij, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    U = Hij
    call eigh(U,E)
    !
    !Allocate and setup the Matsubara freq.
    forall(i=1:L)csi(:,i)=one/(wfreq(i)+xmu-E(:))
    !
    Gloc = zero
    do i=1,L
       do concurrent(ilat=1:Nlat,io=1:Nso,jo=1:Nso)
          is = io + (ilat-1)*Nso
          js = jo + (ilat-1)*Nso
          Gloc(ilat,io,jo,i) = sum(U(is,:)*conjg(U(js,:))*csi(:,i))!can use matmul
       enddo
    enddo
  end subroutine get_gf_hlat

  subroutine get_gf_schlat(Hij,S,Gloc,axis)
    complex(8),dimension(Nlat*Norb,Nlat*Norb)     :: Hij
    complex(8),dimension(2,Nlat,Norb,Norb,L)      :: S
    complex(8),dimension(2,Nlat,Norb,Norb,L)      :: Gloc
    character(len=*)                              :: axis
    complex(8),dimension(L)                       :: wfreq
    complex(8),dimension(2*Nlat*Norb,L)           :: csi
    complex(8),dimension(2,Nlat,Norb,Norb,L)      :: Gtmp
    complex(8),dimension(2*Nlat*Norb,2*Nlat*Norb) :: U
    real(8),dimension(2*Nlat*Norb)                :: E
    integer                                       :: i,ilat,io,jo,is,js,N
    !
    write(*,"(A)")"Get local Green's function Hij, axis:"//str(axis)
    wfreq = build_frequency_array(axis)
    !
    N = Nlat*Norb
    U(1:N,1:N)         = Hij     + reshape_rank3_to_matrix(S(1,:,:,:,1),Nlat,Norb)
    U(1:N,N+1:2*N)     =           reshape_rank3_to_matrix(S(2,:,:,:,1),Nlat,Norb)
    U(N+1:2*N,1:N)     =           reshape_rank3_to_matrix(S(2,:,:,:,1),Nlat,Norb)
    U(N+1:2*N,N+1:2*N) =-conjg(Hij+reshape_rank3_to_matrix(S(1,:,:,:,L),Nlat,Norb))
    !
    call eigh(U,E)
    !
    !Allocate and setup the Matsubara freq.
    forall(i=1:L)csi(:,i)=one/(wfreq(i)+xmu-E(:))
    !
    Gloc = zero
    do i=1,L
       do concurrent(ilat=1:Nlat,io=1:Norb,jo=1:Norb)
          is = io + (ilat-1)*Norb
          js = jo + (ilat-1)*Norb
          Gloc(1,ilat,io,jo,i) = sum(U(is,:)*conjg(U(js,:))*csi(:,i))
          Gloc(2,ilat,io,jo,i) = sum(U(is,:)*conjg(U(N+js,:))*csi(:,i))
       enddo
    enddo
  end subroutine get_gf_schlat


  function build_frequency_array(axis) result(wfreq)
    character(len=*)                    :: axis
    complex(8),dimension(:),allocatable :: wfreq
    if(allocated(wfreq))deallocate(wfreq)
    allocate(wfreq(L))
    !
    select case(to_lower(axis))
    case default;
       stop "build_frequency_array ERROR: axis undefined. axis=[matsubara,realaxis]"
    case("matsubara","mats","m")
       wfreq = dcmplx(0d0,pi/beta*(2*arange(1,L)-1))
    case("realaxis","real","r")
       wfreq = dcmplx(linspace(wini,wfin,L),eps)
    end select
    return
  end function build_frequency_array




  function i2indices(istate,Nvec) result(ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end function i2indices


  subroutine build_kgrid(Nkvec,kgrid)
    integer,dimension(:)             :: Nkvec
    real(8),dimension(:,:)           :: kgrid ![Nk][Ndim]
    real(8),dimension(size(Nkvec))   :: kvec
    real(8),dimension(:),allocatable :: grid_x,grid_y,grid_z
    integer                          :: ik,Ivec(3),Nk(3),ndim,Nktot,i
    real(8),dimension(3)             :: ktmp,wtmp
    !
    Nktot = product(Nkvec)
    Ndim  = size(Nkvec)          !dimension of the grid to be built
    wtmp      = 1d0
    call assert_shape(kgrid,[Nktot,Ndim],"build_kgrid","kgrid")
    !
    Nk=1
    do ik=1,Ndim
       Nk(ik)=Nkvec(ik)
    enddo
    if(product(Nk)/=product(Nkvec))stop "TB_build_grid ERROR: product(Nkvec) != product(Nk)"
    !
    allocate(grid_x(Nk(1)))
    allocate(grid_y(Nk(2)))   
    allocate(grid_z(Nk(3)))
    !
    grid_x = linspace(0d0,wtmp(1),Nk(1),iend=.false.)
    grid_y = linspace(0d0,wtmp(2),Nk(2),iend=.false.)
    grid_z = linspace(0d0,wtmp(3),Nk(3),iend=.false.)
    !
    do ik=1,Nktot
       ivec = i2indices(ik,Nk)
       ktmp = [grid_x(ivec(1)), grid_y(ivec(2)), grid_z(ivec(3))]
       kgrid(ik,:) = ktmp(1)*bk_1(:ndim) + ktmp(2)*bk_2(:ndim) + ktmp(3)*bk_3(:ndim)
    end do
  end subroutine build_kgrid




  !< build the \hat{H}({\mathbf k}) or \hat{H}({\mathbf k};i,j) Hamiltonian matrix
  ! from the function user defined hk_model procedure.
  subroutine TB_build_Hk(Hk,hk_model,Norb,Nkvec)
    integer,dimension(:),intent(in)                :: Nkvec
    integer                                        :: Norb
    real(8),dimension(product(Nkvec),size(Nkvec))  :: kgrid ![Nk][Ndim]
    integer                                        :: Nktot
    integer                                        :: ik
    complex(8),dimension(Norb,Norb,product(Nkvec)) :: hk,haux
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)                      :: kpoint
         integer                                   :: N
         complex(8),dimension(N,N)                 :: hk_model
       end function hk_model
    end interface
    !
    !
    call build_kgrid(Nkvec,kgrid)
    !
    Nktot  = product(Nkvec)
    Hk   = zero
    do ik=1,Nktot
       Hk(:,:,ik) = hk_model(kgrid(ik,:),Norb)
    enddo
  end subroutine TB_build_Hk



  subroutine TB_build_Hlat(Hij,ts_model,Nso,Nrvec,Links,pbc)
    integer                                                     :: Nso
    integer,dimension(:),intent(in)                             :: Nrvec
    integer,dimension(:,:),intent(in)                           :: Links ![Nlink][dim]
    logical,optional                                            :: pbc
    logical                                                     :: pbc_
    integer,dimension(product(Nrvec),size(Nrvec))               :: RNgrid
    integer                                                     :: Nlat,Nlink
    integer                                                     :: ilat,jlat
    integer,dimension(size(Nrvec))                              :: Ri,Rj
    integer                                                     :: i,ilink
    complex(8),dimension(Nso,Nso,product(Nrvec),product(Nrvec)) :: Hij
    integer                                                     :: ir,ix,iy,iz,Nr(3)
    !
    interface
       function ts_model(link,Nso)
         integer                       :: link
         integer                       :: Nso
         complex(8),dimension(Nso,Nso) :: ts_model
       end function ts_model
    end interface
    !
    pbc_ = .true. ; if(present(pbc))pbc_=pbc
    !
    if(size(Nrvec)/=size(Links,2))stop "TB_build_Hij ERROR: size(Nvec) != size(Links,2)"
    !
    Nlat  = product(Nrvec)
    Nlink = size(Links,1)
    !
    call TB_build_CoordGrid(Nrvec,RNgrid)
    !
    Hij = zero
    !
    do_lattice: do ilat = 1,Nlat
       Hij(:,:,ilat,ilat) = Hij(:,:,ilat,ilat) + ts_model(0,Nso)
       Ri = RNgrid(ilat,:)
       do_links: do ilink=1,Nlink
          Rj = Ri + Links(ilink,:)
          if(pbc_)then
             do i=1,size(Nrvec)
                if( Rj(i)==0 )Rj(i)=Nrvec(i)
                if( Rj(i)==Nrvec(i)+1)Rj(i)=1
             enddo
          endif
          jlat = TB_find_IndxCoord(Rj,Nrvec)
          if(jlat==0)cycle do_links
          !
          !>Build Hij
          Hij(:,:,ilat,jlat) = Hij(:,:,ilat,jlat) + ts_model(ilink,Nso)
       enddo do_links
    enddo do_lattice
    return
  end subroutine TB_build_Hlat

  subroutine TB_build_CoordGrid(Nrvec,RNgrid)
    integer,dimension(:)                          :: Nrvec  ! Nr=product(Nrvec);Ndim=size(Nrvec)
    integer,dimension(product(Nrvec),size(Nrvec)) :: RNgrid ![Nr][Ndim]
    integer                                       :: ir,ix,iy,iz,Nr(3),Rvec(3)
    !
    Nr=1
    do ir=1,size(Nrvec)
       Nr(ir)=Nrvec(ir)
    enddo
    !
    ir=0
    do iz=1,Nr(3)
       do iy=1,Nr(2)
          do ix=1,Nr(1)
             ir=ir+1
             Rvec = [ix,iy,iz]
             RNgrid(ir,:) = Rvec(1:size(Nrvec))
          enddo
       enddo
    enddo
    ir=0
  end subroutine TB_build_CoordGrid

  function TB_find_IndxCoord(RNvec,Nrvec) result(ir)
    integer,dimension(:)           :: Nrvec
    integer,dimension(size(Nrvec)) :: RNvec
    integer                        :: Nr(3),Rvec(3)
    integer                        :: ir,ix,iy,iz
    logical :: bool
    !
    Nr=1
    Rvec=1
    do ir=1,size(Nrvec)
       Nr(ir)=Nrvec(ir)
       Rvec(ir)=RNvec(ir)
    enddo
    !
    ir=0
    do iz=1,Nr(3)
       do iy=1,Nr(2)
          do ix=1,Nr(1)
             ir=ir+1
             bool = Rvec(1)==ix.AND.Rvec(2)==iy.AND.Rvec(3)==iz
             if(bool)return
          enddo
       enddo
    enddo
    ir=0
  end function TB_find_IndxCoord




end program test
  
