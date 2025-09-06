module module_hostmodel
  use vars
  implicit none
  private
  public :: host_model_evolve, nudging_hm

  integer :: hm_step = 0

contains



subroutine host_model_evolve( &
  u0_map, v0_map, wsub_map, t0_map, q0_map,  &
  rho, rhow, adz, adzw,                    &
  u_hm_map, v_hm_map, w_hm_map, t_hm_map, q_hm_map)

  implicit none
  ! -------- 输入（含 ghost） --------
  real, intent(in) :: u0_map(1:nsx, nzm)
  real, intent(in) :: v0_map(1:nsx, nzm)
  real, intent(in) :: wsub_map(1:nsx, nz)
  real, intent(in) :: t0_map(1:nsx, nzm)
  real, intent(in) :: q0_map(1:nsx, nzm)
  real, intent(in) :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)
  

  ! -------- 输出（含 ghost） --------
  real, intent(out) :: u_hm_map(1:nsx, nzm)
  real, intent(out) :: v_hm_map(1:nsx, nzm)
  real, intent(out) :: w_hm_map(1:nsx, nz)
  real, intent(out) :: t_hm_map(1:nsx, nzm)
  real, intent(out) :: q_hm_map(1:nsx, nzm)

  ! -------- 局部 --------
  real :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)
  real :: p_phys(1:nsx, nzm)    ! 压力势（诊断用，可不输出）

  ! 拷贝初值
  u_hm_map = u0_map
  v_hm_map = v0_map
  w_hm_map = wsub_map
  t_hm_map = t0_map
  q_hm_map = q0_map

  ! 1) 动量平流（2D，二阶中心）
  call advect_mom_hm(u_hm_map, v_hm_map, w_hm_map, rho, rhow, adz, adzw, &
                     dudt_hm, dvdt_hm, dwdt_hm, dx_hm, dz_hm)

  ! 2) 压力投影（谱法：x-FFT/DCT + z-三对角）
  call pressure_hm_spectral(u_hm_map, w_hm_map, rho, rhow, adz, adzw, dt_hm, dx_hm, dz_hm, &
                            dudt_hm, dwdt_hm, p_phys)

  ! 3) AB 时间推进
  call adams_hm(u_hm_map, v_hm_map, w_hm_map, dudt_hm, dvdt_hm, dwdt_hm)

  ! 4) 标量平流（上风，正定）
  call advect_scalars_hm(t_hm_map, u_hm_map, w_hm_map, rho, rhow, adz, adzw, dx_hm, dz_hm, dt_hm, .false.)
  call advect_scalars_hm(q_hm_map, u_hm_map, w_hm_map, rho, rhow, adz, adzw, dx_hm, dz_hm, dt_hm, .true.)

 
end subroutine host_model_evolve

!================== 动量平流：2D 二阶中心 ==================

subroutine advect_mom_hm(u_hm_map, v_hm_map, w_hm_map, rho, rhow, adz, adzw, dudt_hm, dvdt_hm, dwdt_hm, dx_hm, dz_hm)
  implicit none
  ! 输入 
  real, intent(in)  :: u_hm_map(1:nsx, nzm), v_hm_map(1:nsx, nzm), w_hm_map(1:nsx, nz)
  real, intent(in)  :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)
  real, intent(in)  :: dx_hm, dz_hm                ! dx 为 host 水平间距（列间距）
  ! 输出
  real, intent(out) :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)

  ! 局部
  real :: dx25, dz25, irho_w, irho_k, irhow_k
  real, allocatable :: fu(:,:), fv(:,:), fw(:,:)     ! x 向通量
  real, allocatable :: fuz(:,:), fvz(:,:), fwz(:,:)  ! z 向通量（注意 f*u/f*v 在 w 层，大小 nz）
  integer :: i, ic, k, kc, kb, kcu

  !---- 分配并清零
  allocate(fu(0:nsx, nzm), fv(0:nsx, nzm), fw(0:nsx, nzm))
  allocate(fuz(1:nsx, nz ), fvz(1:nsx, nz ), fwz(1:nsx, nzm))
  fu = 0.0; fv = 0.0; fw = 0.0
  fuz = 0.0; fvz = 0.0; fwz = 0.0

  dudt_hm = 0.0; dvdt_hm = 0.0; dwdt_hm = 0.0

  dx25 = 0.25 / dx_hm         
  dz25 = 1.0  / (4.0*dz_hm)

  !==================== x 向通量 ====================
  do k = 1, nzm
    kc  = k + 1
    kcu = min(kc, nzm)
    irho_w = 1.0 / ( rhow(kc) * adzw(kc) )
    do i = 0, nsx
      if (i < 1) i = nsx + i
      ic = i + 1
      if (ic > nsx) ic = ic - nsx
      fu(i,k) = dx25 * (u_hm_map(ic,k)+u_hm_map(i,k)) * (u_hm_map(i,k)+u_hm_map(ic,k))
      fv(i,k) = dx25 * (u_hm_map(ic,k)+u_hm_map(i,k)) * (v_hm_map(i,k)+v_hm_map(ic,k))   ! advect v by u
      fw(i,k) = dx25 * ( u_hm_map(ic,k)*rho(k)*adz(k) + u_hm_map(ic,kcu)*rho(kcu)*adz(kcu) ) * &
                        ( w_hm_map(i,kc) + w_hm_map(ic,kc) )
    end do
    do i = 1, nsx
      ic = i - 1
      if (ic < 1) ic = nsx + ic
      dudt_hm(i,k)   = dudt_hm(i,k)   - ( fu(i,k) - fu(ic,k) )
      dvdt_hm(i,k)   = dvdt_hm(i,k)   - ( fv(i,k) - fv(ic,k) )
      dwdt_hm(i,kc)  = dwdt_hm(i,kc)  - irho_w * ( fw(i,k) - fw(ic,k) )
    end do
  end do

  !==================== z 向通量 ====================
  ! fuz / fvz 在 w 层定义：k=1…nz；边界 k=1,nz 设 0
  fuz(:,1) = 0.0; fuz(:,nz) = 0.0
  fvz(:,1) = 0.0; fvz(:,nz) = 0.0
  fwz(:,1) = 0.0; fwz(:,nzm) = 0.0
  do k = 2, nzm
    kb = k - 1
    do i = 1, nsx
      fuz(i,k) = dz25 * rhow(k) * ( w_hm_map(i,k) + w_hm_map(i-1,k) ) * ( u_hm_map(i,k) + u_hm_map(i,kb) )
      fvz(i,k) = dz25 * rhow(k) * ( w_hm_map(i,k) + w_hm_map(i-1,k) ) * ( v_hm_map(i,k) + v_hm_map(i,kb) )
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
 

  ! 释放
  deallocate(fu, fv, fw, fuz, fvz, fwz)
end subroutine advect_mom_hm 



!================== 压力投影：谱法（x 变换 + z 三对角） ==================
subroutine pressure_hm_spectral(u_hm_map, w_hm_map, rho, rhow, adz, adzw, dt_hm, dx_hm, dz_hm, &
                                dudt_hm,  dwdt_hm, p_phys)
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  real, intent(in)    :: u_hm_map(1:nsx, nzm)
  real, intent(in)    :: w_hm_map(1:nsx, nz)
  real, intent(in)    :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)
  real, intent(in)    :: dt_hm, dx_hm, dz_hm
  

  real, intent(inout) :: dudt_hm(nsx, nzm),  dwdt_hm(nsx, nz)

  real, intent(out)   :: p_phys(0:nsx, nzm)

  real :: rhs(1:nsx, nzm)
  integer, parameter :: nx2 = nsx + 2

  ! x 向变换的工作数组
  real(real64) :: F(nx2, nzm), WORK(nx2,1), trigx(3*nsx/2+1)
  integer :: ifax(100)

  ! 竖直三对角系数与谱特征值
  real(real64) :: a(nzm), c(nzm), eigx, ddx2, pii, factx
  real(real64) :: alfa(nzm-1), beta(nzm-1), fline(nzm), denom

  ! AB 系数与 press_rhs 系数
  real :: atc, btc, ctc, dta, btat, ctat, rdx, rdz, rup, rdn

  integer :: i, k, kx, ip, im
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
  ! rdx  = 1.0/dx_hm

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
        ( rdx*(u_hm_map(ip,k) - u_hm_map(i,k)) + ( w_hm_map(i,k+1)*rup - w_hm_map(i,k)*rdn ) )*dt_hm  &
      + ( rdx*(dudt_hm(ip,k) - dudt_hm(i,k)) + ( dwdt_hm(i,k+1)*rup - dwdt_hm(i,k)*rdn ) )  &
      
    end do
  end do
  ! --------- end RHS ---------   


  ! --------- 构造 z 向三对角系数 ---------  
  do k = 1, nzm
    a(k) = rhow(k  ) /( rho(k)*adz(k)*adzw(k  ) * dz_hm*dz_hm )
    c(k) = rhow(k+1) /( rho(k)*adz(k)*adzw(k+1) * dz_hm*dz_hm )
  end do

  ddx2 = 1._8/(dx_hm*dx_hm)
  pii  = acos(-1._8)
  factx= 2.d0

  call fftfax_crm(nsx, ifax, trigx)   

  ! --------- x 正变换 ---------
  do k = 1, nzm  ! 认为nzm = nzslab
    F(1:nsx,      k) = rhs(1:nsx,k)
    call fft991_crm(F(1,k), WORK, trigx, ifax, 1, nx2, nsx, 1, -1) 
  end do

  ! --------- 对每个 kx 解竖直三对角 ---------
  do kx = 0, nsx-1                      !这里和SAM 不太一样
    eigx = (2._8*cos(factx*pii*kx/nsx) - 2._8)*ddx2

    do k = 1, nzm
      fline(k) = F(kx+1, k)
    end do

    if (kx == 0) then
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

    do k = 1, nzm
      F(kx+1, k) = fline(k)
    end do
  end do

  ! --------- x 逆变换 → 物理空间 φ(x,z) ---------  
  do k = 1, nzm 
    call fft991_crm(F(1,k), WORK, trigx, ifax, 1, nx2, nsx, 1, +1)
  end do

  do k = 1, nzm
    do i = 1, nsx
      p_phys(i,k) = real(F(i,k))
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



end subroutine pressure_hm_spectral


!================== Adams–Bashforth 时间推进 ==================
subroutine adams_hm(u, v, w, dudt_hm, dvdt_hm, dwdt_hm)
  implicit none
  real, intent(inout) :: u(1:nsx, nzm), v(1:nsx, nzm), w(1:nsx, nz)
  real, intent(in)    :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)
  real :: at, bt, ct
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

  do k=1,nzm; do i=1,nsx
    u(i,k) = u(i,k) + dt_hm*dudt_hm(i,k) 
    v(i,k) = v(i,k) + dt_hm*dvdt_hm(i,k)
    w(i,k) = w(i,k) + dt_hm*dwdt_hm(i,k)
  end do; end do


end subroutine adams_hm


!================== 标量平流：质量通量上风 ==================
subroutine advect_scalars_hm(f, u_hm_map, w_hm_map, rho, rhow, adz, adzw, dx_hm, dz_hm, dtloc, clip_nonneg)   !!!原代码里索引有3个ghost
  implicit none
  ! -------- 接口保持不变 --------
  real, intent(inout) :: f(1:nsx, nzm)
  real, intent(in)    :: u_hm_map(1:nsx, nzm), w_hm_map(1:nsx, nz)
  real, intent(in)    :: rho(nzm), rhow(nz), adz(nzm), adzw(nz)
  real, intent(in)    :: dx_hm, dz_hm, dtloc
  logical, intent(in) :: clip_nonneg

  ! -------- 局部变量 --------
  real :: uuu(0:nsx+1, nzm)         ! x向面通量
  real :: www(1:nsx, nz)          ! z向界面通量
  real :: irho(nzm), iadz(nzm), irhow(nz)
  real :: dd, eps
  integer :: i, ib, ic, k, kb, kc

  real x1, x2, a, b, a1, a2
  real andiff,across,pp,pn
  andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
  across(x1,a1,a2)=0.03125*a1*a2*x1



!========================================================
  eps = 1.e-10

  ! 顶部 w 通量为 0（与 SAM 一致）
  do i = 1, nsx
    www(i, nz) = 0.
  end do

  !========================
  ! 第 1 步：低阶（迎风）通量
  !========================
  do k = 1, nzm
    kb = max(1, k-1)
    do i = 1, nsx+1
      uuu(i, k) = max(0., u_hm_map(i, k))*f(i-1, k) + &
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
    irhow(k) = 1. / ( rhow(k) * adzw(k) )

    
    do i = 1, nsx
      ic = i-1
      if (ic < 1) ic = nsx + ic
      uuu(i, k) = andiff( f(ic,k), f(i,k), u_hm_map(i,k), 1./rho(k) )  &
        - across( dd*( f(ic,kc)+f(i,kc) - f(ic,kb)-f(i,kb) ), &
                  u_hm_map(i,k), w_hm_map(ic,k)+w_hm_map(ic,kc)+w_hm_map(i,k)+w_hm_map(i,kc) ) * (1./rho(k))
    end do

    do i = 1, nsx
      ic = i-1
      if (ic < 1) ic = nsx + ic
      ib = i+1
      if (ib > nsx) ib = ib - nsx
      www(i, k) = andiff( f(i,kb), f(i,k), w_hm_map(i,k), irhow(k) )  &
        - across( f(ib,kb)+f(ib,k) - f(ic,kb)-f(ic,k), &
                  w_hm_map(i,k), ( u_hm_map(i,kb)+u_hm_map(i,k)+u_hm_map(ib,k)+u_hm_map(ib,kb) ) ) * irhow(k)
    end do
  end do

  ! 底部边界（与 SAM 同步）
  do i = 1, nsx
    www(i, 1) = 0.0
  end do

  do k=1,nzm 
    do i=1,nsx
      ib = i+1
      if (ib > nsx) ib = ib - nsx
      f(i,k)= max(0., f(i,k) - (uuu(ib,k)-uuu(i,k) &
                      +(www(i,k+1)-www(i,k))*iadz(k))*irho(k))
    end do
 end do 
  
end subroutine advect_scalars_hm


!==========================================================================================
subroutine nudging_hm()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
implicit none

real coef, coef1
integer i,j,k
	
call t_startf ('nudging')

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

coef = 1./tautqls

if(donudging_tq.or.donudging_t) then
    coef1 = dtn / tautqls
    do k=1,nzm
      if(z(k).ge.nudging_t_z1.and.z(k).le.nudging_t_z2) then
        tnudge(k)=tnudge(k) -(t0_local_hm(k)-tg0_hm(k)-gamaz(k))*coef
        do j=1,ny
          do i=1,nx
             t(i,j,k)=t(i,j,k)-(t0_local_hm(k)-tg0_hm(k)-gamaz(k))*coef1
          end do
        end do
      end if
    end do
endif

if(donudging_tq.or.donudging_q) then
    coef1 = dtn / tautqls
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

call t_stopf('nudging')

end subroutine nudging_hm




end module module_hostmodel
