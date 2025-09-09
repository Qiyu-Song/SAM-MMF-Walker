module module_hostmodel
  use grid, only: nsx, nzm, nz, adz, adzw, dx, dz, dx_hm, dt_hm
  use vars, only: rho, rhow, hm_step
  implicit none
  private
  public :: host_model_init, host_model_finalize, host_model_evolve, nudging_hm

 

contains

subroutine host_model_init()
  use vars
  implicit none

  if (.not. allocated(wsub_map)) then
    allocate(wsub_map(nsx, nz))
  end if

  if (.not. wsub_inited) then
    wsub_map = 0.0     ! 设定初值
    wsub_inited = .true.
    hm_step = 0
  end if
end subroutine host_model_init


subroutine host_model_finalize()  !暂时不打算调用
  use vars
  implicit none
  if (allocated(wsub_map)) deallocate(wsub_map)
  wsub_inited = .false.
end subroutine host_model_finalize


subroutine host_model_evolve( &
  u0_in, v0_in, wsub_in, t0_in, q0_in,  &
  u_hm_map, v_hm_map, w_hm_map, t_hm_map, q_hm_map)

  implicit none
  ! -------- 输入（不含 ghost） --------
  real, intent(in) :: u0_in(nsx, nzm)
  real, intent(in) :: v0_in(nsx, nzm)
  real, intent(in) :: wsub_in(nsx, nz)
  real, intent(in) :: t0_in(nsx, nzm)
  real, intent(in) :: q0_in(nsx, nzm)
  

  ! -------- 输出（不含 ghost） --------
  real, intent(out) :: u_hm_map(nsx, nzm)
  real, intent(out) :: v_hm_map(nsx, nzm)
  real, intent(out) :: w_hm_map(nsx, nz)
  real, intent(out) :: t_hm_map(nsx, nzm)
  real, intent(out) :: q_hm_map(nsx, nzm)

  ! -------- 局部 --------
  real :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)
  ! for advection of scalars
  real :: u1_hm_map(nsx, nzm), v1_hm_map(nsx, nzm), w1_hm_map(nsx, nz)
  real :: p_phys(nsx, nzm)    ! 压力势（诊断用，可不输出）
  real :: tmp(nsx, nzm)

  ! 拷贝初值
  u_hm_map = u0_in
  v_hm_map = v0_in
  w_hm_map = wsub_in
  t_hm_map = t0_in
  q_hm_map = q0_in

  ! interpolate for u,v since using Arakawa C-type grid
  ! u,v should be on the left boundary of grid box
  tmp(1,     :) = 0.5 * (u_hm_map(1,     :) + u_hm_map(nsx,     :))
  tmp(2:nsx, :) = 0.5 * (u_hm_map(2:nsx, :) + u_hm_map(1:nsx-1, :))
  u_hm_map = tmp

  tmp(1,     :) = 0.5 * (v_hm_map(1,     :) + v_hm_map(nsx,     :))
  tmp(2:nsx, :) = 0.5 * (v_hm_map(2:nsx, :) + v_hm_map(1:nsx-1, :))
  v_hm_map = tmp

  dudt_hm = 0.0; dvdt_hm = 0.0; dwdt_hm = 0.0

  ! 1) 动量平流（2D，二阶中心）
  call advect_mom_hm(u_hm_map, v_hm_map, w_hm_map, rho, rhow, adz, adzw, &
                     dudt_hm, dvdt_hm, dwdt_hm, dx_hm, dz)

  ! 2) 压力投影
  call pressure_hm(u_hm_map, w_hm_map, rho, rhow, adz, adzw, dt_hm, dx_hm, dz, &
                            dudt_hm, dwdt_hm, p_phys)

  ! 3) AB 时间推进
  call adams_hm(u_hm_map, v_hm_map, w_hm_map, dudt_hm, dvdt_hm, dwdt_hm, &
                  u1_hm_map, v1_hm_map, w1_hm_map, dt_hm, dx_hm, dz)

  ! 4) 标量平流（上风，正定）
  call advect_scalars_hm(t_hm_map, u1_hm_map, w1_hm_map, rho, rhow, adz, adzw)
  call advect_scalars_hm(q_hm_map, u1_hm_map, w1_hm_map, rho, rhow, adz, adzw)

  ! interpolate back for u,v (due to Arakawa C-type grid)
  tmp(1:nsx-1, :) = 0.5 * (u_hm_map(1:nsx-1, :) + u_hm_map(2:nsx, :))
  tmp(nsx,     :) = 0.5 * (u_hm_map(nsx,     :) + u_hm_map(1,     :))
  u_hm_map = tmp

  tmp(1:nsx-1, :) = 0.5 * (v_hm_map(1:nsx-1, :) + v_hm_map(2:nsx, :))
  tmp(nsx,     :) = 0.5 * (v_hm_map(nsx,     :) + v_hm_map(1,     :))
  v_hm_map = tmp
 
end subroutine host_model_evolve

!================== 动量平流：2D 二阶中心 ==================

subroutine advect_mom_hm(u_hm_map, v_hm_map, w_hm_map, rho, rhow, adz, adzw, dudt_hm, dvdt_hm, dwdt_hm, dx_hm, dz_hm)
  implicit none
  ! 输入 
  real, intent(in)  :: u_hm_map(nsx, nzm), v_hm_map(nsx, nzm), w_hm_map(nsx, nz)
  real, intent(in)  :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)
  real, intent(in)  :: dx_hm, dz_hm                ! dx 为 host 水平间距（列间距）
  ! 输出
  real, intent(inout) :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)

  ! 局部
  real :: dx25, dz25, irho_w, irho_k, irhow_k
  real fu(nsx,nzm), fv(nsx,nzm), fw(nsx, nzm)     ! x 向通量
  real fuz(nsx, nz), fvz(nsx, nz), fwz(nsx, nzm)  ! z 向通量（注意 f*u/f*v 在 w 层，大小 nz）
  integer :: i, ic, ib, k, kc, kb, kcu

  !---- 初始化清零 ----
  fu = 0.; fv = 0.; fw = 0.
  ! fuz / fvz 在 w 层定义：k=1…nz；边界 k=1,nz 设 0
  fuz(:,1) = 0.; fvz(:,1) = 0.; fwz(:,1) = 0.
  fuz(:,nz) = 0.; fvz(:,nz) = 0.; fwz(:,nzm) = 0.

  dx25 = 0.25 / dx_hm         
  dz25 = 1.   / (4.*dz_hm)

  !==================== x 向通量 ====================
  do k = 1, nzm
    kc  = k + 1
    kcu = min(kc, nzm)
    irho_w = 1.0 / ( rhow(kc) * adzw(kc) )
    do i = 1, nsx
      ic = i + 1
      if (ic > nsx) ic = ic - nsx
      fu(i,k) = dx25 * (u_hm_map(ic,k)+u_hm_map(i,k)) * (u_hm_map(i,k)+u_hm_map(ic,k))
      fv(i,k) = dx25 * (u_hm_map(ic,k)+u_hm_map(i,k)) * (v_hm_map(i,k)+v_hm_map(ic,k))   ! advect v by u
      fw(i,k) = dx25 * ( u_hm_map(ic,k)*rho(k)*adz(k) + u_hm_map(ic,kcu)*rho(kcu)*adz(kcu) ) * &
                        ( w_hm_map(i,kc) + w_hm_map(ic,kc) )
    end do
    do i = 1, nsx
      ib = i - 1
      if (ib < 1) ib = nsx + ib
      dudt_hm(i,k)   = dudt_hm(i,k)   - ( fu(i,k) - fu(ib,k) )
      dvdt_hm(i,k)   = dvdt_hm(i,k)   - ( fv(i,k) - fv(ib,k) )
      dwdt_hm(i,kc)  = dwdt_hm(i,kc)  - irho_w * ( fw(i,k) - fw(ib,k) )
    end do
  end do

  !==================== z 向通量 ====================
  
  do k = 2, nzm
    kb = k - 1
    do i = 1, nsx
      ib = i - 1; if (ib < 1) ib = nsx + ib
      fuz(i,k) = dz25 * rhow(k) * ( w_hm_map(i,k) + w_hm_map(ib,k) ) * ( u_hm_map(i,k) + u_hm_map(i,kb) )
      fvz(i,k) = dz25 * rhow(k) * ( w_hm_map(i,k) + w_hm_map(ib,k) ) * ( v_hm_map(i,k) + v_hm_map(i,kb) )
    end do
  end do
  ! 顶层接口 k = nz (=nzm+1) 自然保持 0

  ! 应用到 dudt_hm/dvdt_hm；fwz 给 dwdt_hm
  do k = 1, nzm
    kc = k + 1
    irho_k = 1.0 / ( rho(k) * adz(k) )
    do i = 1, nsx
      dudt_hm(i,k) = dudt_hm(i,k) - ( fuz(i,kc) - fuz(i,k) ) * irho_k
      dvdt_hm(i,k) = dvdt_hm(i,k) - ( fvz(i,kc) - fvz(i,k) ) * irho_k
      fwz(i,k)  = dz25 * ( w_hm_map(i,kc)*rhow(kc) + w_hm_map(i,k)*rhow(k) ) * ( w_hm_map(i,kc) + w_hm_map(i,k) )
    end do
  end do

  do k = 2, nzm
    kb = k - 1
    irhow_k = 1.0 / ( rhow(k) * adzw(k) )
    do i = 1, nsx
      dwdt_hm(i,k) = dwdt_hm(i,k) - ( fwz(i,k) - fwz(i,kb) ) * irhow_k
    end do
  end do
 
end subroutine advect_mom_hm 



!================== 压力投影 ==================
subroutine pressure_hm(u_hm_map, w_hm_map, rho, rhow, adz, adzw, dt_hm, dx_hm, dz_hm, &
                                dudt_hm,  dwdt_hm, p_phys)
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  real, intent(in)    :: u_hm_map(1:nsx, nzm)
  real, intent(in)    :: w_hm_map(1:nsx, nz)
  real, intent(in)    :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)
  real, intent(in)    :: dt_hm, dx_hm, dz_hm
  

  real, intent(inout) :: dudt_hm(nsx, nzm),  dwdt_hm(nsx, nz)

  real, intent(out)   :: p_phys(nsx, nzm)

  real :: rhs(1:nsx, nzm)
  integer, parameter :: nx2 = nsx + 2

  ! x 向变换的工作数组
  real(8) :: F(nx2, nzm), WORK(nx2,1), trigx(3*nsx/2+1)
  integer :: ifax(100)

  ! 竖直三对角系数与谱特征值
  real(8) :: a(nzm), c(nzm), eigx, ddx2, pii, factx, xnx
  real(8) :: alfa(nzm-1), beta(nzm-1), fline(nzm), denom

  ! AB 系数与 press_rhs 系数
  real :: atc, btc, ctc, dta, btat, ctat, rdx, rdz, rup, rdn

  integer :: i, k, kx, ip, im, id
  integer :: stage


  ! --------- AB 系数（沿用 adams_hm 的阶段，但此处还未旋转历史，所以为 hm_step+1）---------
  ! stage = min(3, hm_step + 1)
  ! select case(stage)
  ! case (1)
  !    atc = 1.0 ; btc = 0.0         ; ctc = 0.0
  ! case (2)
  !    atc = 1.5 ; btc = -0.5        ; ctc = 0.0
  ! case default
  !    atc = 23.0/12.0 ; btc = -16.0/12.0 ; ctc = 5.0/12.0
  ! end select

  ! dta  = 1.0/dt_hm/atc
  ! btat = btc/atc
  ! ctat = ctc/atc
  rdx  = 1.0/dx_hm

  ! --------- RHS (press_rhs) — 2D(x,z) 与 SAM 一致的形式 ---------
  do k = 1, nzm
    rdz = 1.0/(adz(k)*dz_hm)
    rup = rhow(k+1)/rho(k) * rdz    ! 上界面系数 (kc=k+1)
    rdn = rhow(k  )/rho(k) * rdz    ! 下界面系数 (k)

    do i = 1, nsx
      ip = i + 1; if (ip > nsx) ip = ip-nsx   

      ! 有历史几步加速度信息时
      ! rhs(i,k) = &
      !   ( rdx*(u_hm_map(ip,k) - u_hm_map(i,k)) + ( w_hm_map(i,k+1)*rup - w_hm_map(i,k)*rdn ) )*dta  &
      ! + ( rdx*(dudt_hm(ip,k) - dudt_hm(i,k)) + ( dwdt_hm(i,k+1)*rup - dwdt_hm(i,k)*rdn ) )  &
      ! + btat*( rdx*(dudt_hm_hist(ip,k,1) - dudt_hm_hist(i,k,1))  &
      !        + ( dwdt_hm_hist(i,k+1,1)*rup - dwdt_hm_hist(i,k,1)*rdn ) )              &
      ! + ctat*( rdx*(dudt_hm_hist(ip,k,2) - dudt_hm_hist(i,k,2))  &
      !        + ( dwdt_hm_hist(i,k+1,2)*rup - dwdt_hm_hist(i,k,2)*rdn ) )

      rhs(i,k) = &
        ( rdx*(u_hm_map(ip,k) - u_hm_map(i,k)) + ( w_hm_map(i,k+1)*rup - w_hm_map(i,k)*rdn ) ) /dt_hm  &
      + ( rdx*(dudt_hm(ip,k) - dudt_hm(i,k)) + ( dwdt_hm(i,k+1)*rup - dwdt_hm(i,k)*rdn ) )
      
    end do
  end do
  ! --------- end RHS ---------   

  call fftfax_crm(nsx, ifax, trigx)   

  ! --------- x 正变换 ---------
  do k = 1, nzm  ! 认为nzm = nzslab
    F(1:nsx,      k) = rhs(1:nsx,k)
    call fft991_crm(F(1,k), WORK, trigx, ifax, 1, nx2, nsx, 1, -1) 
  end do

  ! --------- 构造 z 向三对角系数 ---------  
  ! assuming dowallx = .false., dowally = .false.
  do k = 1, nzm
    a(k) = rhow(k  ) /( rho(k)*adz(k)*adzw(k  ) * dz_hm*dz_hm )
    c(k) = rhow(k+1) /( rho(k)*adz(k)*adzw(k+1) * dz_hm*dz_hm )
  end do

  ddx2 = 1._8/(dx_hm*dx_hm)
  pii  = acos(-1._8)
  xnx=pii/nsx
  factx= 2.d0

  ! --------- 对每个 kx 解竖直三对角 ---------
  do i = 1, nsx+1
    id = (i-0.1)/2.
    eigx = (2._8*cos(factx*xnx*id) - 2._8)*ddx2

    fline(1:nzm) = F(i, 1:nzm)

    if (id.eq.0) then
      beta(1) = fline(1)/(eigx - a(1) - c(1))
      alfa(1) = -c(1)   /(eigx - a(1) - c(1))
    else
      beta(1) = fline(1)/(eigx - c(1))
      alfa(1) = -c(1)   /(eigx - c(1))
    end if

    do k = 2, nzm-1
      denom   = eigx - a(k) - c(k) + a(k)*alfa(k-1)
      alfa(k) = -c(k) / denom
      beta(k) = (fline(k) - a(k)*beta(k-1)) / denom
    end do

    fline(nzm) = (fline(nzm) - a(nzm)*beta(nzm-1)) / (eigx - a(nzm) + a(nzm)*alfa(nzm-1))
    do k = nzm-1, 1, -1
      fline(k) = alfa(k)*fline(k+1) + beta(k)
    end do

    F(i,1:nzm) = fline(1:nzm)
  end do

  ! --------- x 逆变换 → 物理空间 φ(x,z) ---------  
  do k = 1, nzm 
    call fft991_crm(F(1,k), WORK, trigx, ifax, 1, nx2, nsx, 1, +1)
  end do

  do k = 1, nzm
    do i = 1, nsx
      p_phys(i,k) = F(i,k)
    end do
  end do

  ! --------- 压力梯度修正 RHS（2D：仅 u,w）---------
  do k = 1, nzm
    do i = 1, nsx
      im = i - 1; if (im < 1)   im = nsx + im
      dudt_hm(i,k) = dudt_hm(i,k) - (p_phys(i,k) - p_phys(im,k))/dx_hm  !dvdt_hm不做更新
      dwdt_hm(i,k) = dwdt_hm(i,k) - (p_phys(i,k)-p_phys(i,max(1,k-1)))/(dz_hm*adzw(k))
    end do
  end do

  do k = 1, nzm
    do i = 1, nsx
      p_phys(i,k) = p_phys(i,k) * rho(k)   ! convert p'/rho to p'
    end do
  end do

end subroutine pressure_hm


!================== Adams–Bashforth 时间推进 ==================
subroutine adams_hm(u, v, w, dudt_hm, dvdt_hm, dwdt_hm, u1, v1, w1, dt_hm, dx_hm, dz_hm)
  implicit none
  real, intent(inout) :: u(nsx, nzm), v(nsx, nzm), w(nsx, nz)
  real, intent(in)    :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)
  real, intent(out)   :: u1(nsx, nzm), v1(nsx, nzm), w1(nsx, nz)
  real, intent(in)    :: dt_hm, dx_hm, dz_hm
  real :: at, bt, ct
  real :: dtdx, dtdz, rhox, rhoy, rhoz , a1, a2
  integer :: i,k

  hm_step = hm_step + 1

  ! 有历史几步加速度信息时
  ! if (hm_step == 1) then
  !   at=1.0; bt=0.0; ct=0.0
  ! else if (hm_step == 2) then
  !   at=1.5; bt=-0.5; ct=0.0
  ! else
  !   at=23.0/12.0; bt=-16.0/12.0; ct=5.0/12.0
  ! end if

  ! do k=1,nzm; do i=1,nsx
  !   dudt_hm_hist(i,k,3) = dudt_hm_hist(i,k,2)
  !   dudt_hm_hist(i,k,2) = dudt_hm_hist(i,k,1)
  !   dudt_hm_hist(i,k,1) = dudt_hm(i,k)
  !   dvdt_hm_hist(i,k,3) = dvdt_hm_hist(i,k,2)
  !   dvdt_hm_hist(i,k,2) = dvdt_hm_hist(i,k,1)
  !   dvdt_hm_hist(i,k,1) = dvdt_hm(i,k)
  ! end do; end do
  ! do k=1,nz; do i=1,nsx
  !   dwdt_hm_hist(i,k,3) = dwdt_hm_hist(i,k,2)
  !   dwdt_hm_hist(i,k,2) = dwdt_hm_hist(i,k,1)
  !   dwdt_hm_hist(i,k,1) = dwdt_hm(i,k)
  ! end do; end do

  ! do k=1,nzm; do i=1,nsx
  !   u(i,k) = u(i,k) + dt_hm*( at*dudt_hm_hist(i,k,1) + bt*dudt_hm_hist(i,k,2) + ct*dudt_hm_hist(i,k,3) )
  !   v(i,k) = v(i,k) + dt_hm*( at*dvdt_hm_hist(i,k,1) + bt*dvdt_hm_hist(i,k,2) + ct*dvdt_hm_hist(i,k,3) )
  !   w(i,k) = w(i,k) + dt_hm*( at*dwdt_hm_hist(i,k,1) + bt*dwdt_hm_hist(i,k,2) + ct*dwdt_hm_hist(i,k,3) )  !原代码只更新到nzm
  ! end do; end do

  u1(:,:) = u(:,:)
  v1(:,:) = v(:,:)
  w1(:,:) = w(:,:)

  do k=1,nzm; do i=1,nsx
    u(i,k) = u1(i,k) + dt_hm*dudt_hm(i,k) 
    v(i,k) = v1(i,k) + dt_hm*dvdt_hm(i,k)
    w(i,k) = w1(i,k) + dt_hm*dwdt_hm(i,k)
  end do; end do
  
  ! compute time averaged velocties for second-order advection of scalars:
  dtdx = dt_hm/dx_hm
  dtdz = dt_hm/dz_hm
  a1 = 0.5
  a2 = 0.5
  if(hm_step.eq.1) then
    a1 = 1.
    a2 = 0.
  end if

  do k=1,nzm
    rhox = rho(k)*dtdx
    rhoz = rhow(k)*dtdz
    do i=1,nsx
      u1(i,k) = (a1*u(i,k)+a2*u1(i,k))*rhox
      v1(i,k) = (a1*v(i,k)+a2*v1(i,k))*rhox ! assume dx_hm = dy_hm
      w1(i,k) = (a1*w(i,k)+a2*w1(i,k))*rhoz
    end do
  end do

end subroutine adams_hm


!================== 标量平流：质量通量上风 ==================
! using MPDATA method
subroutine advect_scalars_hm(f, u_hm_map, w_hm_map, rho, rhow, adz, adzw)
  implicit none
  
  real, intent(inout) :: f(1:nsx, nzm)
  real, intent(in)    :: u_hm_map(1:nsx, nzm), w_hm_map(1:nsx, nz)
  real, intent(in)    :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)

  ! -------- 局部变量 --------
  real :: mx (1:nsx, nzm)
  real :: mn (1:nsx, nzm)
  real :: uuu(1:nsx, nzm)         ! x向面通量
  real :: www(1:nsx, nz)          ! z向界面通量
  real :: irho(nzm), iadz(nzm), irhow(nz)
  real :: dd, eps
  integer :: i, ib, ic, k, kb, kc
  logical nonos

  real x1, x2, a, b, a1, a2, y
  real andiff,across,pp,pn
  andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
  across(x1,a1,a2)=0.03125*a1*a2*x1
  pp(y)= max(0.,y)
  pn(y)=-min(0.,y)


!========================================================
  nonos = .true.
  eps = 1.e-10


  www(:, nz) = 0.

  if(nonos) then

    do k = 1, nzm
      kc = min(nzm, k+1)
      kb = max(1, k-1)
      do i = 1, nsx
        ib = i-1
        if (ib < 1) ib = nsx + ib
        ic = i + 1
        if (ic > nsx) ic = ic - nsx
        mx(i, k) = max(f(ib, k), f(ic, k), f(i, kb), f(i, kc), f(i, k))
        mn(i, k) = min(f(ib, k), f(ic, k), f(i, kb), f(i, kc), f(i, k))
      end do
    end do

  end if  ! nonos

  !========================
  ! 第 1 步：低阶（迎风）通量
  !========================
  do k = 1, nzm
    kb = max(1, k-1)
    do i = 1, nsx
      ib = i - 1
      if (ib < 1) ib = nsx + ib
      uuu(i, k) = max(0., u_hm_map(i, k))*f(ib, k) + &
                  min(0., u_hm_map(i, k))*f(i, k)
    end do

    do i = 1, nsx
      www(i, k) = max(0., w_hm_map(i, k))*f(i, kb) + &
                  min(0., w_hm_map(i, k))*f(i, k )
    end do
  end do

  
  do k = 1, nzm
    irho(k) = 1. / rho(k)
    iadz(k) = 1. / adz(k)
    do i = 1, nsx
      ic = i + 1
      if (ic > nsx) ic = ic - nsx
      f(i, k) = f(i, k) - ( uuu(ic, k) - uuu(i, k)   &
                   + ( www(i, k+1) - www(i, k) )*iadz(k) ) * irho(k)
    end do
  end do

  !========================
  ! 第 2 步：反扩散通量
  !========================
  do k = 1, nzm
    kc = min(nzm, k+1)
    kb = max(1,   k-1)
    dd = 2. / (kc - kb) / adz(k)
    irhow(k) = 1. / ( rhow(k) * adz(k) )

    
    do i = 1, nsx
      ib = i-1
      if (ib < 1) ib = nsx + ib
      uuu(i, k) = andiff( f(ib,k), f(i,k), u_hm_map(i,k), irho(k) )  &
        - across( dd*( f(ib,kc)+f(i,kc) - f(ib,kb)-f(i,kb) ), &
                  u_hm_map(i,k), w_hm_map(ib,k)+w_hm_map(ib,kc)+w_hm_map(i,k)+w_hm_map(i,kc) ) * irho(k)
    end do

    do i = 1, nsx
      ib = i - 1
      if (ib < 1) ib = nsx + ib
      ic = i + 1
      if (ic > nsx) ic = ic - nsx
      www(i, k) = andiff( f(i,kb), f(i,k), w_hm_map(i,k), irhow(k) )  &
        - across( f(ic,kb)+f(ic,k) - f(ib,kb)-f(ib,k), &
                  w_hm_map(i,k), ( u_hm_map(i,kb)+u_hm_map(i,k)+u_hm_map(ic,k)+u_hm_map(ic,kb) ) ) * irhow(k)
    end do
  end do

  ! 底部边界
  www(:, 1) = 0.

  !---------- non-osscilatory option ---------------

  if (nonos) then

    do k = 1, nzm
      kc = min(nzm, k+1)
      kb = max(1, k-1)
      do i = 1, nsx
        ib = i - 1
        if (ib < 1) ib = nsx + ib
        ic = i + 1
        if (ic > nsx) ic = ic - nsx
        mx(i, k) = max(f(ib, k), f(ic, k), f(i, kb), f(i, kc), f(i, k), mx(i, k))
        mn(i, k) = min(f(ib, k), f(ic, k), f(i, kb), f(i, kc), f(i, k), mn(i, k))
      end do
    end do

    do k = 1, nzm
      kc = min(nzm, k+1)
      do i = 1, nsx
        ic = i + 1
        if (ic > nsx) ic = ic - nsx
        mx(i, k) = rho(k) * (mx(i,k)-f(i,k)) / (pn(uuu(ic,k))+pp(uuu(i,k)) + &
                  iadz(k)*(pn(www(i,kc))+pp(www(i,k))) + eps)	
        mn(i, k) = rho(k) * (f(i,k)-mn(i,k)) / (pp(uuu(ic,k))+pn(uuu(i,k)) + &
                  iadz(k)*(pp(www(i,kc))+pn(www(i,k))) + eps)	
      end do
    end do

    do k = 1, nzm
      kb = max(1, k-1)
      do i = 1, nsx
        ib = i - 1
        if (ib < 1) ib = nsx + ib
        uuu(i,k) = pp(uuu(i,k)) * min(1.,mx(i,k), mn(ib,k)) &
                 - pn(uuu(i,k)) * min(1.,mx(ib,k),mn(i,k))
      end do
      do i=1,nsx
        www(i,k) = pp(www(i,k)) * min(1.,mx(i,k), mn(i,kb)) &
                 - pn(www(i,k)) * min(1.,mx(i,kb),mn(i,k))
      end do
    end do

  endif ! nonos

  do k = 1, nzm 
    do i = 1, nsx
      ic = i+1
      if (ic > nsx) ic = ic - nsx
      f(i,k)= max(0., f(i,k) - (uuu(ic,k)-uuu(i,k) &
                      +(www(i,k+1)-www(i,k))*iadz(k))*irho(k))
    end do
 end do 
  
end subroutine advect_scalars_hm


!==========================================================================================
subroutine nudging_hm()
! keep nudging terms constant in each host model time step
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
implicit none

real coef, coef1
integer i,j,k
	
call t_startf ('nudging_hm')

tnudge = 0.
qnudge = 0.
unudge = 0.
vnudge = 0.

coef = 1./dt_hm

if(donudging_uv) then
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0_local_hm(k)-ug0_hm(k))*coef
        vnudge(k)=vnudge(k) - (v0_local_hm(k)-vg0_hm(k))*coef
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0_local_hm(k)-ug0_hm(k))*coef
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0_local_hm(k)-vg0_hm(k))*coef
          end do
        end do
      end if
    end do
endif

! no minus gamaz here since both t0_local_hm and tg0_hm include gamaz
if(donudging_tq.or.donudging_t) then
    coef1 = dtn / dt_hm
    do k=1,nzm
      if(z(k).ge.nudging_t_z1.and.z(k).le.nudging_t_z2) then
        tnudge(k)=tnudge(k) -(t0_local_hm(k)-tg0_hm(k))*coef
        do j=1,ny
          do i=1,nx
             t(i,j,k)=t(i,j,k)-(t0_local_hm(k)-tg0_hm(k))*coef1
          end do
        end do
      end if
    end do
endif

if(donudging_tq.or.donudging_q) then
    coef1 = dtn / dt_hm
    do k=1,nzm
      if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
        qnudge(k)=qnudge(k) -(q0_local_hm(k)-qg0_hm(k))*coef
        do j=1,ny
          do i=1,nx
             micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)-(q0_local_hm(k)-qg0_hm(k))*coef1
          end do
        end do
      end if
    end do
endif

call t_stopf('nudging_hm')

end subroutine nudging_hm




end module module_hostmodel
