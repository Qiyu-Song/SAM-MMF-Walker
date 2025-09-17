module module_hostmodel
  use grid, only: nsx, nzm, nz, adz, adzw, dz, dx_hm, dt_hm
  use vars, only: rho, rhow, hm_step
  use domain, only: nsubdomains_x
  implicit none
  private
  public :: host_model_init, host_model_finalize, host_model_evolve, nudging_hm

 

contains

subroutine host_model_init()
  use vars
  implicit none

  if (.not. allocated(wsub_map)) then
    allocate(wsub_map(nsubdomains_x, nz))
  end if

  if (.not. wsub_inited) then
    wsub_map = 0.0     ! 设定初值
    wsub_inited = .true.
    hm_step = 0
  end if
  ! 平均subdomain里的sst
  sstxy(1:nx, 1:ny) = sum(sstxy(:,:))/(nx*ny)
end subroutine host_model_init


subroutine host_model_finalize()  !暂时不打算调用
  use vars
  implicit none
  if (allocated(wsub_map)) deallocate(wsub_map)
  wsub_inited = .false.
end subroutine host_model_finalize


subroutine host_model_evolve( &
  u0_in, v0_in, wsub_in, t0_in, q0_in,  &
  tabs0_in, qv0_in, qn0_in, qp0_in, &
  u_out_map, v_out_map, w_out_map, t_out_map, q_out_map)

  implicit none
  ! -------- 输入（不含 ghost） --------
  real, intent(in) :: u0_in(nsubdomains_x, nzm)
  real, intent(in) :: v0_in(nsubdomains_x, nzm)
  real, intent(in) :: wsub_in(nsubdomains_x, nz)
  real, intent(in) :: t0_in(nsubdomains_x, nzm)
  real, intent(in) :: q0_in(nsubdomains_x, nzm)
  real, intent(in) :: tabs0_in(nsubdomains_x, nzm)
  real, intent(in) :: qv0_in(nsubdomains_x, nzm)
  real, intent(in) :: qn0_in(nsubdomains_x, nzm)
  real, intent(in) :: qp0_in(nsubdomains_x, nzm)
  

  ! -------- 输出（不含 ghost） --------
  real, intent(out) :: u_out_map(nsubdomains_x, nzm)
  real, intent(out) :: v_out_map(nsubdomains_x, nzm)
  real, intent(out) :: w_out_map(nsubdomains_x, nz)
  real, intent(out) :: t_out_map(nsubdomains_x, nzm)
  real, intent(out) :: q_out_map(nsubdomains_x, nzm)

  ! -------- 局部 --------
  real :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)
  ! for advection of scalars
  real :: u1_hm_map(nsx, nzm), v1_hm_map(nsx, nzm), w1_hm_map(nsx, nz)
  real :: p_phys(nsx, nzm)    ! 压力势（诊断用，可不输出）
  real :: tmp(nsx, nzm)
  ! real :: u_hm_map_edge(nsx, nzm)
  ! real :: v_hm_map_edge(nsx, nzm)
  real :: u_hm_map(nsx, nzm), v_hm_map(nsx, nzm), w_hm_map(nsx, nz), t_hm_map(nsx, nzm), q_hm_map(nsx, nzm)
  real :: tabs0_hm_map(nsx, nzm), qv0_hm_map(nsx, nzm), qn0_hm_map(nsx, nzm), qp0_hm_map(nsx, nzm)
  
  ! 拷贝初值
  ! u_hm_map = u0_in
  ! v_hm_map = v0_in
  ! w_hm_map = wsub_in
  ! t_hm_map = t0_in
  ! q_hm_map = q0_in

  ! interpolate for u,v since using Arakawa C-type grid
  ! u,v should be on the left boundary of grid box
  ! tmp(1,     :) = 0.5 * (u_hm_map(1,     :) + u_hm_map(nsx,     :))
  ! tmp(2:nsx, :) = 0.5 * (u_hm_map(2:nsx, :) + u_hm_map(1:nsx-1, :))
  ! u_hm_map = tmp

  ! tmp(1,     :) = 0.5 * (v_hm_map(1,     :) + v_hm_map(nsx,     :))
  ! tmp(2:nsx, :) = 0.5 * (v_hm_map(2:nsx, :) + v_hm_map(1:nsx-1, :))
  ! v_hm_map = tmp

  ! 平均变量到粗网格
  call pair_avg_U(u0_in, u_hm_map)
  call pair_avg_U(v0_in, v_hm_map)
  call pair_avg_T(t0_in, t_hm_map)
  call pair_avg_T(q0_in, q_hm_map)
  call pair_avg_w(wsub_in, w_hm_map)
  call pair_avg_T(tabs0_in, tabs0_hm_map)
  call pair_avg_T(qv0_in, qv0_hm_map)
  call pair_avg_T(qn0_in, qn0_hm_map)
  call pair_avg_T(qp0_in, qp0_hm_map)

  v_hm_map = 0.

  call output_host_model_single_variable(u_hm_map, 'U00', 'U_after_1st_interpolation' , 'm/s')
  call output_host_model_single_variable(v_hm_map, 'V00', 'V_after_1st_interpolation' , 'm/s')
  call output_host_model_single_variable(t_hm_map, 'T00', 'T_after_1st_interpolation' , 'K')
  call output_host_model_single_variable(q_hm_map, 'Q00', 'Q_after_1st_interpolation' , 'kg/kg')
  tmp(:, :) = w_hm_map(:,1:nzm)
  call output_host_model_single_variable(tmp, 'W00', 'W_after_buoyancy' , 'm/s2')

  dudt_hm = 0.0; dvdt_hm = 0.0; dwdt_hm = 0.0

  ! 1) 浮力
  call buoyancy_hm(tabs0_hm_map, qv0_hm_map, qn0_hm_map, qp0_hm_map, dwdt_hm)

  tmp(:, :) = dwdt_hm(:,1:nzm)
  call output_host_model_single_variable(tmp, 'dwdt1', 'dwdt_after_buoyancy' , 'm/s2')
  call output_host_model_single_variable(qv0_hm_map, 'qv01', 'qv0_hm_map_for_buoyancy' , 'kg/kg')
  call output_host_model_single_variable(qp0_hm_map, 'qp01', 'qp0_hm_map_for_buoyancy' , 'kg/kg')
  call output_host_model_single_variable(qn0_hm_map, 'qn01', 'qn0_hm_map_for_buoyancy' , 'kg/kg')
  call output_host_model_single_variable(tabs0_hm_map, 'tabs01', 'tabs0_hm_map_for_buoyancy' , 'K')

  ! 1-1) damping
  call damping_hm(u_hm_map, v_hm_map, w_hm_map, dudt_hm, dvdt_hm, dwdt_hm)

  call output_host_model_single_variable(dudt_hm, 'dudt_dam', 'dudt_after_damping' , 'm/s2')
  call output_host_model_single_variable(dvdt_hm, 'dvdt_dam', 'dvdt_after_damping' , 'm/s2')
  tmp(:, :) = dwdt_hm(:,1:nzm)
  call output_host_model_single_variable(tmp, 'dwdt_dam', 'dwdt_after_damping' , 'm/s2')

  ! 2) 动量平流（2D，二阶中心）
  call advect_mom_hm(u_hm_map, v_hm_map, w_hm_map, &
                     dudt_hm, dvdt_hm, dwdt_hm)

  call output_host_model_single_variable(dudt_hm, 'dudt2', 'dudt_after_advect_mom' , 'm/s2')
  call output_host_model_single_variable(dvdt_hm, 'dvdt2', 'dvdt_after_advect_mom' , 'm/s2')
  tmp(:, :) = dwdt_hm(:,1:nzm)
  call output_host_model_single_variable(tmp, 'dwdt2', 'dwdt_after_advect_mom' , 'm/s2')

  ! 3) 压力投影
  call pressure_hm(u_hm_map, w_hm_map, &
                            dudt_hm, dwdt_hm, p_phys)

  call output_host_model_single_variable(dudt_hm, 'dudt3', 'dudt_after_pressure' , 'm/s2')
  tmp(:, :) = dwdt_hm(:,1:nzm)
  call output_host_model_single_variable(tmp, 'dwdt3', 'dwdt_after_pressure' , 'm/s2')
  call output_host_model_single_variable(p_phys, 'p_phys3', 'Pressure_Perturbation' , 'Pa')

  ! 4) AB 时间推进
  call adams_hm(u_hm_map, v_hm_map, w_hm_map, dudt_hm, dvdt_hm, dwdt_hm, &
                  u1_hm_map, v1_hm_map, w1_hm_map)

  call output_host_model_single_variable(u_hm_map, 'U4', 'U_after_adams' , 'm/s')
  call output_host_model_single_variable(v_hm_map, 'V4', 'V_after_adams' , 'm/s')
  tmp(:, :) = w_hm_map(:,1:nzm)
  call output_host_model_single_variable(tmp, 'W4', 'W_after_adams' , 'm/s')
  call output_host_model_single_variable(u1_hm_map, 'U41', 'U1_after_adams' , 'm/s')
  call output_host_model_single_variable(v1_hm_map, 'V41', 'V1_after_adams' , 'm/s')
  tmp(:, :) = w1_hm_map(:,1:nzm)
  call output_host_model_single_variable(tmp, 'W41', 'W1_after_adams' , 'm/s')

  ! 5) 标量平流（上风，正定）
  call advect_scalars_hm(t_hm_map, u1_hm_map, w1_hm_map)
  call advect_scalars_hm(q_hm_map, u1_hm_map, w1_hm_map)

  ! u_hm_map_edge = u_hm_map
  ! v_hm_map_edge = v_hm_map

  ! ! interpolate back for u,v (due to Arakawa C-type grid)
  ! tmp(1:nsx-1, :) = 0.5 * (u_hm_map(1:nsx-1, :) + u_hm_map(2:nsx, :))
  ! tmp(nsx,     :) = 0.5 * (u_hm_map(nsx,     :) + u_hm_map(1,     :))
  ! u_hm_map = tmp

  ! tmp(1:nsx-1, :) = 0.5 * (v_hm_map(1:nsx-1, :) + v_hm_map(2:nsx, :))
  ! tmp(nsx,     :) = 0.5 * (v_hm_map(nsx,     :) + v_hm_map(1,     :))
  ! v_hm_map = tmp

  ! 插值回去
  call copy_interp_U(u_hm_map, u_out_map)
  call copy_interp_U(v_hm_map, v_out_map)
  call copy_interp_T(t_hm_map, t_out_map)
  call copy_interp_T(q_hm_map, q_out_map)
  call copy_interp_w(w_hm_map, w_out_map)
  ! 输出in, out (nsubdomains_x,)
  call output_host_model(u0_in, v0_in, t0_in, q0_in,  &
                          u_out_map, v_out_map, t_out_map, q_out_map, w_out_map,&
                          tabs0_in, qv0_in, qn0_in, qp0_in)

end subroutine host_model_evolve



subroutine buoyancy_hm(tabs0_in, qv0_in, qn0_in, qp0_in, dwdt_hm)
  use vars
  use params
  implicit none

  real, intent(in)  :: tabs0_in(nsx, nzm), qv0_in(nsx, nzm), qn0_in(nsx, nzm), qp0_in(nsx, nzm)
  real, intent(inout) :: dwdt_hm(nsx, nz)
  real :: tabs0_entire_domain(nzm), qv0_entire_domain(nzm), qn0_entire_domain(nzm), qp0_entire_domain(nzm)

  integer i,k,kb
  real betu, betd
	
  do k = 1, nzm
      tabs0_entire_domain(k) = sum( tabs0_in(:,k) ) / nsx
      qv0_entire_domain(k) = sum( qv0_in(:,k) ) / nsx
      qn0_entire_domain(k) = sum( qn0_in(:,k) ) / nsx
      qp0_entire_domain(k) = sum( qp0_in(:,k) ) / nsx  
  end do

  do k=2,nzm	
    kb=k-1
    betu=adz(kb)/(adz(k)+adz(kb))
    betd=adz(k)/(adz(k)+adz(kb))
    
    do i=1,nsx
      dwdt_hm(i,k)=dwdt_hm(i,k) +  &
          bet(k)*betu* &
        ( tabs0_entire_domain(k)*(epsv*(qv0_in(i,k)-qv0_entire_domain(k))-(qn0_in(i,k)-qn0_entire_domain(k)+qp0_in(i,k)-qp0_entire_domain(k))) &
          +(tabs0_in(i,k)-tabs0_entire_domain(k))*(1.+epsv*qv0_entire_domain(k)-qn0_entire_domain(k)-qp0_entire_domain(k)) ) &
        + bet(kb)*betd* &
        ( tabs0_entire_domain(kb)*(epsv*(qv0_in(i,kb)-qv0_entire_domain(kb))-(qn0_in(i,kb)-qn0_entire_domain(kb)+qp0_in(i,kb)-qp0_entire_domain(kb))) &
          +(tabs0_in(i,kb)-tabs0_entire_domain(kb))*(1.+epsv*qv0_entire_domain(kb)-qn0_entire_domain(kb)-qp0_entire_domain(kb)) )  

    end do ! i
  end do ! k

end subroutine buoyancy_hm



subroutine damping_hm(u_hm_map, v_hm_map, w_hm_map, dudt_hm, dvdt_hm, dwdt_hm)
    use vars
    implicit none

    real, intent(in)  :: u_hm_map(nsx, nzm), v_hm_map(nsx, nzm), w_hm_map(nsx, nz)
    real, intent(inout) :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)

    real :: u0_entire_domain(nzm), v0_entire_domain(nzm)

    real tau_min	! minimum damping time-scale (at the top)
    real tau_max    ! maxim damping time-scale (base of damping layer)
    real damp_depth ! damping depth as a fraction of the domain height
    parameter(tau_min=1800., tau_max=3600., damp_depth=0.3)
    real tau(nzm)   
    integer i, k, n_damp

   

    do k=nzm,1,-1
        if(z(nzm)-z(k).lt.damp_depth*z(nzm)) then 
            n_damp=nzm-k+1
        endif
    end do
    ! print*, 'n_damp', n_damp
    do k=nzm,nzm-n_damp,-1
        tau(k) = tau_min *(tau_max/tau_min)**((z(nzm)-z(k))/(z(nzm)-z(nzm-n_damp)))
        tau(k)=1./tau(k)
    end do
    ! print*, 'finish tau'
    do k = 1, nzm
        u0_entire_domain(k) = sum( u_hm_map(:,k) ) / nsx
        v0_entire_domain(k) = sum( v_hm_map(:,k) ) / nsx
    end do
    ! print*, 'horizontal mean'
    do k = nzm, nzm-n_damp, -1
        do i=1,nsx
            dudt_hm(i,k)= dudt_hm(i,k)-(u_hm_map(i,k)-u0_entire_domain(k)) * tau(k)
            dvdt_hm(i,k)= dvdt_hm(i,k)-(v_hm_map(i,k)-v0_entire_domain(k)) * tau(k)
            dwdt_hm(i,k)= dwdt_hm(i,k)-w_hm_map(i,k) * tau(k)
        end do
    end do 

    
end subroutine damping_hm



!================== 动量平流：2D 二阶中心 ==================

subroutine advect_mom_hm(u_hm_map, v_hm_map, w_hm_map, dudt_hm, dvdt_hm, dwdt_hm)
  implicit none
  ! 输入 
  real, intent(in)  :: u_hm_map(nsx, nzm), v_hm_map(nsx, nzm), w_hm_map(nsx, nz)

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
  dz25 = 1.   / (4.*dz)

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
subroutine pressure_hm(u_hm_map, w_hm_map,  &
                                dudt_hm,  dwdt_hm, p_phys)
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  real, intent(in)    :: u_hm_map(1:nsx, nzm)
  real, intent(in)    :: w_hm_map(1:nsx, nz)

  

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
    rdz = 1.0/(adz(k)*dz)
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
    a(k) = rhow(k  ) /( rho(k)*adz(k)*adzw(k  ) * dz*dz )
    c(k) = rhow(k+1) /( rho(k)*adz(k)*adzw(k+1) * dz*dz )
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
      dwdt_hm(i,k) = dwdt_hm(i,k) - (p_phys(i,k)-p_phys(i,max(1,k-1)))/(dz*adzw(k))
    end do
  end do

  do k = 1, nzm
    do i = 1, nsx
      p_phys(i,k) = p_phys(i,k) * rho(k)   ! convert p'/rho to p'
    end do
  end do

end subroutine pressure_hm


!================== Adams–Bashforth 时间推进 ==================
subroutine adams_hm(u, v, w, dudt_hm, dvdt_hm, dwdt_hm, u1, v1, w1)
  implicit none
  real, intent(inout) :: u(nsx, nzm), v(nsx, nzm), w(nsx, nz)
  real, intent(in)    :: dudt_hm(nsx, nzm), dvdt_hm(nsx, nzm), dwdt_hm(nsx, nz)
  real, intent(out)   :: u1(nsx, nzm), v1(nsx, nzm), w1(nsx, nz)

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
  dtdz = dt_hm/dz
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
subroutine advect_scalars_hm(f, u_hm_map, w_hm_map)
  implicit none
  
  real, intent(inout) :: f(1:nsx, nzm)
  real, intent(in)    :: u_hm_map(1:nsx, nzm), w_hm_map(1:nsx, nz)

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



     
subroutine output_host_model(u0_in, v0_in, t0_in, q0_in,  &
  u_out_map, v_out_map, t_out_map, q_out_map, w_out_map,  &
  tabs0_in, qv0_in, qn0_in, qp0_in)         ! 输出形状是（nsubdomains_x，nzm ）的变量
	
    use vars

    implicit none
    real, intent(in) :: u0_in(nsubdomains_x, nzm), v0_in(nsubdomains_x, nzm), t0_in(nsubdomains_x, nzm), q0_in(nsubdomains_x, nzm)
    real, intent(in) :: u_out_map(nsubdomains_x, nzm), v_out_map(nsubdomains_x, nzm), t_out_map(nsubdomains_x, nzm), q_out_map(nsubdomains_x, nzm), w_out_map(nsubdomains_x, nz)
    real, intent(in) :: tabs0_in(nsubdomains_x, nzm), qv0_in(nsubdomains_x, nzm), qn0_in(nsubdomains_x, nzm), qp0_in(nsubdomains_x, nzm)
    character *120 filename
    character *80 long_name
    character *8 name
    character *10 timechar
    character *4 rankchar
    character *5 sepchar
    character *6 filetype
    character *10 units

    integer i,k,nfields_hm,nfields1_hm
    real(4) tmp(nsubdomains_x,1,nzm)
    integer, external :: lenstr


    nfields_hm=13 ! number of 3D fields to save
    nfields1_hm=0

    sepchar=""

    write(rankchar,'(i4)') 1 !nsubdomains
    write(timechar,'(i10)') nstep
    do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
    end do

    print*, 'Rank=', rank, '*************begin output_host_model***************'

    filetype = '.bin2D'
    filename='./OUT_3D/host_model_output_new_damping_2/'//trim(case)//'_'//trim(caseid)//&
    filetype//sepchar
    if(nrestart.eq.0.and.notopened3D) then
        open(46,file=filename,status='unknown',form='unformatted')	
    else
        open(46,file=filename,status='unknown', &
                            form='unformatted', position='append')
    end if
    notopened3D=.false.

    write(46) nsubdomains_x,1,nzm,1,1,1,nfields_hm
 
    do k=1,nzm
        write(46) real(z(k),4)
    end do
 
    do k=1,nzm
        write(46) real(pres(k),4)
    end do
   
    write(46) real(dx_hm,4)
   
    write(46) real(dy,4)
    write(46) real(float(nstep)*dt/(3600.*24.)+day0,4)

  
    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=u0_in(i,k)
        end do
    end do
    name='U0_In'
    long_name='Input X Wind Component For Host Model'
    units='m/s'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

   
    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=v0_in(i,k)
        end do
    end do
    name='V0_In'
    long_name='Input Y Wind Component For Host Model'
    units='m/s'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=t0_in(i,k)
        end do
    end do
    name='T0_In'
    long_name='Input Liquid/Ice Water Static Energy For Host Model'
    units='K'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=q0_in(i,k)
        end do
    end do
    name='Q0_In'
    long_name='Input Water Vapor For Host Model'
    units='kg/kg'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=tabs0_in(i,k)
        end do
    end do
    name='Tabs0_In'
    long_name='Input Tabs For Host Model'
    units='K'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=qv0_in(i,k)
        end do
    end do
    name='Qv0_In'
    long_name='Input Qv0 For Host Model'
    units='kg/kg'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)
    
    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=qn0_in(i,k)
        end do
    end do
    name='Qn0_In'
    long_name='Input Qn0 For Host Model'
    units='kg/kg'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=qp0_in(i,k)
        end do
    end do
    name='Qp0_In'
    long_name='Input Qp0 For Host Model'
    units='kg/kg'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=u_out_map(i,k)
        end do
    end do
    name='U0_Out'
    long_name='Output X Wind Component For Host Model'
    units='m/s'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=v_out_map(i,k)
        end do
    end do
    name='V0_Out'
    long_name='Output Y Wind Component For Host Model'
    units='m/s'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=t_out_map(i,k)
        end do
    end do
    name='T0_Out'
    long_name='Output Liquid/Ice Water Static Energy For Host Model'
    units='K'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=q_out_map(i,k)
        end do
    end do
    name='Q0_Out'
    long_name='Output Water Vapor For Host Model'
    units='kg/kg'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsubdomains_x
            tmp(i,1,k)=w_out_map(i,k)
        end do
    end do
    name='W_SUB_OUT'
    long_name='Output Large Scale Z Wind For Host Model'
    units='m/s'
    call compress3D_hm(tmp,nsubdomains_x,1,nzm,name,long_name,units)


    if(nfields_hm.ne.nfields1_hm) then
        print*,'host model write_fields3D error: nfields_hm=',nfields_hm,'  nfields1_hm=',nfields1_hm
        call task_abort()
    end if

    close (46)

    print*, 'Appending 3D data. file:'//filename

end subroutine output_host_model



subroutine output_host_model_single_variable(u0_in, v_name,v_longname,v_unit) ! 输出形状是（nsx， ）的变量
	
    use vars

    implicit none
    real, intent(in) :: u0_in(nsx, nzm)
    character(*), intent(in) :: v_name, v_longname, v_unit
    character *120 filename
    character *80 long_name
    character *8 name
    character *10 timechar
    character *4 rankchar
    character *5 sepchar
    character *6 filetype
    character *10 units
   

    integer i,k,nfields_hm,nfields1_hm
    real(4) tmp(nsx,1,nzm)
    integer, external :: lenstr
    name = v_name
    long_name = v_longname
    units = v_unit
    nfields_hm=1
    nfields1_hm=0

    sepchar=""

    write(rankchar,'(i4)') 1 !nsubdomains
    write(timechar,'(i10)') nstep
    do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
    end do

    print*, 'Rank=', rank, '*************begin output_hm_single_variable***************'

    filetype = '.bin2D'
    filename='./OUT_3D/host_model_output_new_damping_2/'//trim(case)//'_'//trim(caseid)//&
    '_'//trim(name)//filetype//sepchar
    if(nrestart.eq.0.and.notopened3D) then
        open(46,file=filename,status='unknown',form='unformatted')	
    else
        open(46,file=filename,status='unknown', &
                            form='unformatted', position='append')
    end if
    notopened3D=.false.
    write(46) nsx,1,nzm,1,1,1,nfields_hm
    do k=1,nzm
        write(46) real(z(k),4)
    end do
    do k=1,nzm
        write(46) real(pres(k),4)
    end do
    write(46) real(dx_hm,4)  
    write(46) real(dy,4)
    write(46) real(float(nstep)*dt/(3600.*24.)+day0,4)

  
    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=u0_in(i,k)
        end do
    end do
    call compress3D_hm(tmp,nsx,1,nzm, name,long_name,units)

    if(nfields_hm.ne.nfields1_hm) then
        print*,'host model write_fields3D error: nfields_hm=',nfields_hm,'  nfields1_hm=',nfields1_hm
        call task_abort()
    end if

    close (46)

    print*, 'Appending 3D data. file:'//filename

end subroutine output_host_model_single_variable



subroutine compress3D_hm (f,nx,ny,nz,name, long_name, units)
    implicit none

    integer nx,ny,nz
    real(4) f(nx,1,nz)
    character*(*) name,long_name,units

    integer(2), allocatable :: byte(:)
    real(4), allocatable :: byte4(:)
    integer size,count

    character(10) value_min(nz), value_max(nz)
    character(7) form
    integer int_fac, integer_max, integer_min
    parameter (int_fac=2,integer_min=-32000, integer_max=32000)
    !	parameter (int_fac=1,integer_min=-127, integer_max=127)
    real(4) f_max(1),f_min(1), f_max1(1), f_min1(1), scale
    integer i,j,k,req

    size=nx*ny*nz 
    allocate (byte4(size))
    count = 0

    do k=1,nz
        do i=1,nx
            count = count+1
            byte4(count) = f(i,1,k)
        end do
    end do

    write(46) name,' ',long_name,' ',units
    write(46) (byte4(k),k=1,count)

    deallocate(byte4)

end subroutine compress3D_hm
	

! 如果要输出txt
! subroutine write_host_diag(u0_in, v0_in, t0_in, q0_in, &
!                            u_hm_map, v_hm_map, t_hm_map, q_hm_map, w_hm_map, hm_step)

!   use grid, only: nsx, nzm
!   implicit none
!   ! ===== 输入参数 =====
!   real, intent(in) :: u0_in(nsx, nzm), v0_in(nsx, nzm), t0_in(nsx, nzm), q0_in(nsx, nzm)
!   real, intent(in) :: u_hm_map(nsx, nzm), v_hm_map(nsx, nzm), t_hm_map(nsx, nzm), q_hm_map(nsx, nzm)
!   real, intent(in) :: w_hm_map(nsx, nz)
!   integer, intent(in) :: hm_step   ! 当前步数

!   ! ===== 局部变量 =====
!   integer :: iunit, i, k
!   character(len=120) :: diagfile

!   diagfile = './OUT_3D/host_model_output/host_model_diag.txt'

!   ! 打开文件：存在则追加，不存在则新建
!   open(newunit=iunit, file=diagfile, status='unknown', position='append', &
!        action='write', form='formatted')

!   write(iunit,*) '===== Host Model Diagnostic Output ====='
!   write(iunit,*) 'Step:', hm_step

!   ! -------- 写 9 个场 --------
!   ! 输入
!   call write_field(iunit, 'U0_In', 'Input X Wind Component', 'm/s', u0_in)
!   call write_field(iunit, 'V0_In', 'Input Y Wind Component', 'm/s', v0_in)
!   call write_field(iunit, 'T0_In', 'Input Liquid/Ice Water Static Energy', 'K', t0_in)
!   call write_field(iunit, 'Q0_In', 'Input Water Vapor', 'g/kg', q0_in)

!   ! 输出
!   call write_field(iunit, 'U0_Out', 'Output X Wind Component', 'm/s', u_hm_map)
!   call write_field(iunit, 'V0_Out', 'Output Y Wind Component', 'm/s', v_hm_map)
!   call write_field(iunit, 'T0_Out', 'Output Liquid/Ice Water Static Energy', 'K', t_hm_map)
!   call write_field(iunit, 'Q0_Out', 'Output Water Vapor', 'g/kg', q_hm_map)

!   ! W 特殊：nz 层
!   write(iunit,*) '--- W_SUB_OUT (Output Large Scale Z Wind, m/s) ---'
!   do k = 1, nz
!      write(iunit,'(100f12.5)') (w_hm_map(i,k), i=1,nsx)
!   end do

!   close(iunit)

! end subroutine write_host_diag
! subroutine write_field(iunit, name, long_name, units, f)
!   use grid, only: nsx, nzm
!   implicit none
!   integer, intent(in) :: iunit
!   character(*), intent(in) :: name, long_name, units
!   real, intent(in) :: f(nsx, nzm)
!   integer :: i, k

!   write(iunit,*) '--- ', trim(name), ' (', trim(long_name), ', ', trim(units), ') ---'
!   do k = 1, nzm
!      write(iunit,'(100f12.5)') (f(i,k), i=1,nsx)
!   end do
! end subroutine write_field

subroutine pair_avg_U(u_hi_res, u_lo_res)
    use grid, only: nsx, nzm
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: u_hi_res(nsubdomains_x,nzm)
    real, intent(out)  :: u_lo_res(nsx,nzm)

    u_lo_res(1,:) = 0.5*(u_hi_res(nsubdomains_x,:) + u_hi_res(1,:))
    u_lo_res(2:nsx,:) = 0.5*( u_hi_res(2:nsubdomains_x-2:2,:) + u_hi_res(3:nsubdomains_x-1:2,:) )
end subroutine pair_avg_U


subroutine pair_avg_T(T_hi_res, T_lo_res)
    use grid, only: nsx, nzm
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: T_hi_res(nsubdomains_x,nzm)
    real, intent(out)  :: T_lo_res(nsx,nzm)

    T_lo_res(1:nsx,:) = 0.5*( T_hi_res(1:nsubdomains_x-1:2,:) + T_hi_res(2:nsubdomains_x:2,:) )
    
end subroutine pair_avg_T

subroutine pair_avg_w(T_hi_res, T_lo_res)
    use grid, only: nsx, nz
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: T_hi_res(nsubdomains_x,nz)
    real, intent(out)  :: T_lo_res(nsx,nz)

    T_lo_res(1:nsx,:) = 0.5*( T_hi_res(1:nsubdomains_x-1:2,:) + T_hi_res(2:nsubdomains_x:2,:) )
    
end subroutine pair_avg_w


subroutine inv_pair_interp_U(u_lo_res, u_hi_res)
    use grid, only: nsx, nzm
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: u_lo_res(nsx,nzm) 
    real, intent(out)  :: u_hi_res(nsubdomains_x,nzm)
    integer p, pp

    do p = 1, nsx-1   
      pp = p + 2
      if (pp > nsx) pp = pp - nsx       
      u_hi_res(2*p, :) = 0.25*u_lo_res(p, :) + 0.75*u_lo_res(p + 1, :)
      u_hi_res(2*p + 1, :) = 0.25*u_lo_res(p + 2, :) + 0.75*u_lo_res(p + 1, :)
    end do
    u_hi_res(1, :) = 0.25*u_lo_res(2, :) + 0.75*u_lo_res(1, :)
    u_hi_res(nsubdomains_x, :) = 0.25*u_lo_res(nsx, :) + 0.75*u_lo_res(1, :)
end subroutine inv_pair_interp_U


subroutine inv_pair_interp_T(T_lo_res, T_hi_res)
    use grid, only: nsx, nzm
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: T_lo_res(nsx,nzm)
    real, intent(out)  :: T_hi_res(nsubdomains_x,nzm)
    integer p

    do p = 1, nsx-1            
      T_hi_res(2*p, :) = 0.75*T_lo_res(p, :) + 0.25*T_lo_res(p + 1, :)
      T_hi_res(2*p + 1, :) = 0.25*T_lo_res(p, :) + 0.75*T_lo_res(p + 1, :)
    end do
    T_hi_res(1, :) = 0.25*T_lo_res(nsx, :) + 0.75*T_lo_res(1, :)
    T_hi_res(nsubdomains_x, :) = 0.75*T_lo_res(nsx, :) + 0.25*T_lo_res(1, :)
end subroutine inv_pair_interp_T

subroutine inv_pair_interp_w(T_lo_res, T_hi_res)
    use grid, only: nsx, nz
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: T_lo_res(nsx,nz)
    real, intent(out)  :: T_hi_res(nsubdomains_x,nz)
    integer p

    do p = 1, nsx-1            
      T_hi_res(2*p, :) = 0.75*T_lo_res(p, :) + 0.25*T_lo_res(p + 1, :)
      T_hi_res(2*p + 1, :) = 0.25*T_lo_res(p, :) + 0.75*T_lo_res(p + 1, :)
    end do
    T_hi_res(1, :) = 0.25*T_lo_res(nsx, :) + 0.75*T_lo_res(1, :)
    T_hi_res(nsubdomains_x, :) = 0.75*T_lo_res(nsx, :) + 0.25*T_lo_res(1, :)
end subroutine inv_pair_interp_w


subroutine copy_interp_U(u_lo_res, u_hi_res)
    use grid, only: nsx, nzm
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: u_lo_res(nsx,nzm) 
    real, intent(out)  :: u_hi_res(nsubdomains_x,nzm)
    integer p

    do p = 1, nsx-1         
      u_hi_res(2*p, :) = u_lo_res(p + 1, :)
      u_hi_res(2*p + 1, :) = u_lo_res(p + 1, :)
    end do
    u_hi_res(1, :) = u_lo_res(1, :)
    u_hi_res(nsubdomains_x, :) = u_lo_res(1, :)
end subroutine copy_interp_U

subroutine copy_interp_T(T_lo_res, T_hi_res)
    use grid, only: nsx, nzm
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: T_lo_res(nsx,nzm)
    real, intent(out)  :: T_hi_res(nsubdomains_x,nzm)
    integer p

    do p = 1, nsx-1            
      T_hi_res(2*p, :) = T_lo_res(p, :)
      T_hi_res(2*p + 1, :) = T_lo_res(p + 1, :)
    end do
    T_hi_res(1, :) = T_lo_res(1, :)
    T_hi_res(nsubdomains_x, :) = T_lo_res(nsx, :) 
end subroutine copy_interp_T

subroutine copy_interp_w(T_lo_res, T_hi_res)
    use grid, only: nsx, nz
    use domain, only: nsubdomains_x
    implicit none
    real, intent(in)   :: T_lo_res(nsx,nz)
    real, intent(out)  :: T_hi_res(nsubdomains_x,nz)
    integer p

    do p = 1, nsx-1            
      T_hi_res(2*p, :) = T_lo_res(p, :) 
      T_hi_res(2*p + 1, :) = T_lo_res(p + 1, :)
    end do
    T_hi_res(1, :) = T_lo_res(1, :)
    T_hi_res(nsubdomains_x, :) = T_lo_res(nsx, :)
end subroutine copy_interp_w

end module module_hostmodel
