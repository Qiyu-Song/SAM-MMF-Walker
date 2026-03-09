module module_hostmodel
  use grid, only: nsx, nzm, nz, adz, adzw, dz, dx_hm, dt_hm, dt_hm_subcycle
  use vars, only: rho, rhow, hm_step
  implicit none
  private
  public :: host_model_init, host_model_finalize, host_model_evolve, nudging_hm, nudging_hm_nouv, modify_U_for_subdomain, set_sin_x_sst
 

contains

subroutine host_model_init()
  use vars
  implicit none
  integer i, k
  real :: dudt_tmp(nsx, nzm), dwdt_tmp(nsx, nz), tmp(nsx, nzm), tmp2(nsx, nz)

  if (.not. allocated(wsub_map)) then
    allocate(wsub_map(nsx, nz))
  end if

  if (.not. wsub_inited) then
    wsub_map = 0.0     
    wsub_inited = .true.
    could_hm_nudging = .false.
    hm_step = 0
    u_hm_map_save = 0.
    u_sub_map_save = 0.
    u_hm_updated_map_save = 0.
    dudt_hm_hist = 0.
    dwdt_hm_hist = 0.
    ! ----------------添加条带的初始场-----------------
    ! do k = 1, 3
    !   do i = 1, nsx
    !     u_hm_map_save(i,k) = 1.0*(-1.0)**(i+1)
    !   end do
    ! end do
    ! tmp = 0.
    ! call pressure_hm(u_hm_map_save, wsub_map, dudt_tmp, dwdt_tmp, tmp)
    ! call adams_hm(u_hm_map_save,  wsub_map, dudt_tmp, dwdt_tmp, tmp, tmp2)
    ! u_sub_map_save = 0.
    ! u_hm_updated_map_save = u_hm_map_save
    ! ------------------------------------------------
    do i = 1, nsx
      t_hm_map_save(i,:) = t0(:)
      t_hm_updated_map_save(i,:) = t0(:)
      t_sub_map_save(i,:) = t0(:)
      q_hm_map_save(i,:) = q0(:)
      q_hm_updated_map_save(i,:) = q0(:)
      q_sub_map_save(i,:) = q0(:)
      
    end do
    
    ! if (hm_only) then
    !   call set_sin_x_sst_for_hm(t_hm_map_save)
    ! end if
  end if

end subroutine host_model_init

! subroutine set_constant_sst_hm()
!   use vars, only: sstxy,t00
!   use params, only: tabs_s, delta_sst, ocean_type
!   use grid
!   sstxy = tabs_s - t00
! end subroutine set_constant_sst_hm


! subroutine set_sin_x_sst_for_hm(t_map)
!   use vars
!   implicit none
!   real, intent(inout) :: t_map(nsx, nzm)
!   real(8) pii
!   integer i,k
! 
!   pii = atan2(0.d0,-1.d0)
!   do i = 1, nsx
!     do k = 1, 6
!       t_map(i,k) = t_map(i,k) - 2.5*cos(2.*pii*i/nsx)         ! - 1.0*((-1.)**(i))
!     end do
!   end do
! end subroutine set_sin_x_sst_for_hm

subroutine set_sin_x_sst()
  use vars, only: sstxy,t00
  use params, only: tabs_s, delta_sst, ocean_type
  use grid
  implicit none
  real(8) tmpx(nx), pii, lx
  integer i,j
  sstxy = tabs_s - t00
  lx = float(nx_gl)*dx
  do i = 1,nx
    tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
  end do
  pii = atan2(0.d0,-1.d0)
  do j=1,ny
    do i=1,nx
      sstxy(i,j) = tabs_s-delta_sst*cos(2.*pii*tmpx(i)/lx) - t00
    end do
  end do
  sstxy(1:nx, 1:ny) = sum(sstxy(:,:))/(nx*ny)
end subroutine set_sin_x_sst

! subroutine set_sin_x_sst_stripe()
!   use vars, only: sstxy,t00
!   use params, only: tabs_s, delta_sst, ocean_type
!   use grid
!   implicit none
!   real(8) tmpx(nx), pii, lx
!   integer i,j
!   sstxy = tabs_s - t00
!   lx = float(nx_gl)*dx
!   do i = 1,nx
!     tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
!   end do
!   pii = atan2(0.d0,-1.d0)
!   do j=1,ny
!     do i=1,nx
!       sstxy(i,j) = tabs_s-delta_sst*cos(2.*pii*tmpx(i)/lx) - t00 + 1.0*((-1.0)**(rank+1))
!     end do
!   end do
!   sstxy(1:nx, 1:ny) = sum(sstxy(:,:))/(nx*ny)
! end subroutine set_sin_x_sst_stripe

subroutine host_model_finalize()  !暂时不打算调用
  use vars
  implicit none
  if (allocated(wsub_map)) deallocate(wsub_map)
  wsub_inited = .false.
end subroutine host_model_finalize


subroutine host_model_evolve( &
   u0_in, wsub_in, t0_in, q0_in,  &
  tabs0_in, qn0_in, qp0_in, &
  qni0_in, qnl0_in, qpi0_in, qpl0_in, prec_flx_map, &
  u_out_map,  w_out_map, t_out_map, q_out_map, u_press_modify)
  use vars
  use params, only: fac_cond, fac_fus, fac_sub
  implicit none
  ! -------- 输入（不含 ghost） --------
  real, intent(in) :: u0_in(nsx, nzm)
  real, intent(in) :: wsub_in(nsx, nz)
  real, intent(in) :: t0_in(nsx, nzm)
  real, intent(in) :: q0_in(nsx, nzm)
  real, intent(in) :: tabs0_in(nsx, nzm)

  real, intent(in) :: qn0_in(nsx, nzm)
  real, intent(in) :: qp0_in(nsx, nzm)
  real, intent(in) :: qni0_in(nsx, nzm)
  real, intent(in) :: qnl0_in(nsx, nzm)
  real, intent(in) :: qpi0_in(nsx, nzm)
  real, intent(in) :: qpl0_in(nsx, nzm)
  real, intent(in) :: prec_flx_map(nsx, nzm)
  
  

  ! -------- 输出（不含 ghost） --------
  real, intent(out) :: u_out_map(nsx, nzm)
  real, intent(out) :: w_out_map(nsx, nz)
  real, intent(out) :: t_out_map(nsx, nzm)
  real, intent(out) :: q_out_map(nsx, nzm)
  real, intent(out) :: u_press_modify(nsx, nzm)


  ! -------- 局部 --------
  real :: dudt_hm(nsx, nzm),  dwdt_hm(nsx, nz),  dwdt_hm_tmp(nsx, nz)
  ! for advection of scalars
  real :: u1_hm_map(nsx, nzm),  w1_hm_map(nsx, nz)
  real :: p_phys(nsx, nzm)    ! 压力势（诊断用，可不输出）
  real :: tmp(nsx, nzm), tmp1(nsx, nzm), tmp2(nsx, nzm), tmp_U(nsx, nzm), tmp_dudt(nsx, nzm)
  real :: u_hm_map(nsx, nzm),  w_hm_map(nsx, nz), t_hm_map(nsx, nzm), q_hm_map(nsx, nzm) ! , qni_hm_map(nsx, nzm) , qnl_hm_map(nsx, nzm), qpi_hm_map(nsx, nzm), qpl_hm_map(nsx, nzm)
  integer :: i, k
  real :: tabs_map_hm(nsx, nzm)
  logical :: do_3step_adams_tmp
  integer :: icyc

  call output_host_model_single_variable(prec_flx_map, 'PrecFlux', '1Prec_Rate_2Sensible_heat_flux_3Latent_heat_flux' , 'mm/day_W/m2', 0)
  u_out_map = 0.
  w_out_map = 0.
  t_out_map = 0.
  q_out_map = 0.
  u_press_modify = 0.

  
  w_hm_map = wsub_in
  
  if (hm_only) then
    u_hm_map = u_hm_updated_map_save         ! u_hm_updated_map_save 上一次hm_step更新之后的u
    t_hm_map = t_hm_updated_map_save
    q_hm_map = q_hm_updated_map_save

    u_hm_map_save = u_hm_map
    t_hm_map_save = t_hm_map
    q_hm_map_save = q_hm_map

  
  else
    if (nouvchatting) then
      u_hm_map = u_hm_updated_map_save
    else !全都通信的情况
      call face2center_U((u_hm_updated_map_save-u_hm_map_save),tmp1)  !u_hm_updated_map_save 上一次hm_step更新之后的u； u_hm_map_save 上一次hm_step更新之前的u

      call output_host_model_single_variable(u_hm_updated_map_save-u_hm_map_save, 'deltaU', 'u_hm_updated-u_hm' , 'm/s', 0)

      call center2face_U((u0_in-u_sub_map_save-tmp1), tmp2)

      call output_host_model_single_variable(tmp2, 'delta2U', 'modification_to_Uhm' , 'm/s', 0)

      u_hm_map = u_hm_updated_map_save + tmp2
      tmp_U = u_hm_map

      ! --------------------------------------------修正u,w------------------------------------------------------------
      dudt_hm = 0.0
      dwdt_hm = 0.0
      
      ! the reason to do this modification is that u/w fields in the host model should hold continuity equation,
      ! but CRM updates does not guarantee this.
      ! still need to temporarily enforce do_3step_adams = .false. in pressure_hm
      do_3step_adams_tmp = do_3step_adams
      do_3step_adams = .false.
      call pressure_hm(u_hm_map, w_hm_map, &
                                dudt_hm, dwdt_hm, p_phys)
      do_3step_adams = do_3step_adams_tmp

      call output_host_model_single_variable(dudt_hm, 'dudt_R', 'dudt_after_pressure_R' , 'm/s2', 0)
      tmp(:, :) = dwdt_hm(:,1:nzm)
      call output_host_model_single_variable(tmp, 'dwdt_R', 'dwdt_after_pressure_R' , 'm/s2', 0)
      call output_host_model_single_variable(p_phys, 'p_phys_R', 'Pressure_Perturbation_R' , 'Pa', 0)

      ! add modification terms to u/w fields
      ! this will be the final state of the previous step / the initial state for the next step for the host model. 
      ! dt_hm_subcycle is used because this is the time step assumed in pressure_hm computation of tendencies.
      u_hm_map = u_hm_map + dudt_hm * dt_hm_subcycle
      w_hm_map = w_hm_map + dwdt_hm * dt_hm_subcycle

      call output_host_model_single_variable(u_hm_map, 'U_R', 'U_after_adams_R' , 'm/s', 0)
      tmp(:, :) = w_hm_map(:,1:nzm)
      call output_host_model_single_variable(tmp, 'W_R', 'W_after_adams_R' , 'm/s', 0)
      
      ! ------------------------------------------------------------------------------------------------------------------
      call face2center_U((u_hm_map - tmp_U), u_press_modify)
      call output_host_model_single_variable(u_press_modify, 'U_modify', 'U_back_to_subdomain' , 'm/s', 0)
      u_sub_map_save = u0_in + u_press_modify
    end if   ! if (nouvchatting) else

    t_hm_map = t_hm_map_save + t0_in - t_sub_map_save
    t_sub_map_save = t0_in

    q_hm_map = q_hm_map_save + q0_in - q_sub_map_save
    q_sub_map_save = q0_in


    u_hm_map_save = u_hm_map
    t_hm_map_save = t_hm_map
    q_hm_map_save = q_hm_map


  end if    ! if (hm_only) else


  ! call output_host_model_single_variable(qp0_in, 'qp0in', 'qp0_in_for_buoyancy' , 'kg/kg', 0)
  ! call output_host_model_single_variable(qn0_in, 'qn0in', 'qn0_in_for_buoyancy' , 'kg/kg', 0)
  ! call output_host_model_single_variable(tabs0_in, 'tabs0in', 'tabs0_in_for_buoyancy' , 'K', 0)


  !----------------------------------------------------------------------------------
  do icyc = 1,hm_subcycle

    dudt_hm = 0.0; dwdt_hm = 0.0

    ! 1) 浮力

    ! -------------------加个bubble------------------------------------------------------------------------------------
    
    ! if (hm_only) then
    !   if (hm_step.le.20) then
    !     call hot_bubble(hm_step, t_hm_map)
    !   end if
    ! end if
    ! call output_host_model_single_variable(t_hm_map, 't1', 't_after_bubble' , 'K', icyc)
    
    
    ! -------------------------------------------------------------------------------------------------------
    if (hm_only) then
      call buoyancy_only_in_hm(t_hm_map, q_hm_map, dwdt_hm)
    else
      do k = 1, nzm
          do i = 1, nsx
            
              tabs_map_hm(i,k) = t_hm_map(i,k) - gamaz(k)+ fac_cond * (qnl0_in(i,k)+qpl0_in(i,k)) +fac_sub *(qni0_in(i,k) + qpi0_in(i,k))    ! tabs(i,j,k) = t(i,j,k)-gamaz(k)+ fac_cond * (qcl(i,j,k)+qpl(i,j,k)) +fac_sub *(qci(i,j,k) + qpi(i,j,k))
            
          end do   
      end do
      
      call buoyancy_hm(tabs_map_hm, q_hm_map-qn0_in, qn0_in, qp0_in, dwdt_hm)
      
    end if

    tmp(:, :) = dwdt_hm(:,1:nzm)
    call output_host_model_single_variable(tmp, 'dwdt1', 'dwdt_after_buoyancy' , 'm/s2', icyc)
    call output_host_model_single_variable(q_hm_map-qn0_in, 'qv0', 'qv0_in_for_buoyancy' , 'kg/kg', icyc)   
    call output_host_model_single_variable(tabs_map_hm, 'tabshm', 'tabs_map_hm_for_buoyancy' , 'K', icyc)

    ! 1-1) damping
    call damping_hm(u_hm_map, w_hm_map, dudt_hm, dwdt_hm)

    call output_host_model_single_variable(dudt_hm, 'dudt_dam', 'dudt_after_damping' , 'm/s2', icyc)
    tmp(:, :) = dwdt_hm(:,1:nzm)
    call output_host_model_single_variable(tmp, 'dwdt_dam', 'dwdt_after_damping' , 'm/s2', icyc)


    tmp_dudt = dudt_hm
    call diffuse_u(u_hm_map, dudt_hm)
    call diffuse_w(w_hm_map, dwdt_hm)
    ! call diffuse_TQ(t_hm_map)
    ! call diffuse_TQ(q_hm_map)
    call output_host_model_single_variable(dudt_hm-tmp_dudt, 'dudt_dif', 'dudt_diffuse' , 'm/s2', icyc)


    call output_host_model_single_variable(u_hm_map, 'u_smooth', 'u_after_diffusion' , 'm/s', icyc)


    ! 2) 动量平流（2D，二阶中心）
    call advect_mom_hm(u_hm_map,  w_hm_map, &
                      dudt_hm, dwdt_hm)

    call output_host_model_single_variable(dudt_hm, 'dudt2', 'dudt_after_advect_mom' , 'm/s2', icyc)
    tmp(:, :) = dwdt_hm(:,1:nzm)
    call output_host_model_single_variable(tmp, 'dwdt2', 'dwdt_after_advect_mom' , 'm/s2', icyc)

    

    ! 3) 压力投影
    call pressure_hm(u_hm_map, w_hm_map, &
                              dudt_hm, dwdt_hm, p_phys)

    call output_host_model_single_variable(dudt_hm, 'dudt3', 'dudt_after_pressure' , 'm/s2', icyc)
    tmp(:, :) = dwdt_hm(:,1:nzm)
    call output_host_model_single_variable(tmp, 'dwdt3', 'dwdt_after_pressure' , 'm/s2', icyc)
    call output_host_model_single_variable(p_phys, 'p_phys3', 'Pressure_Perturbation' , 'Pa', icyc)

    ! 4) AB 时间推进
    call adams_hm(u_hm_map, w_hm_map, dudt_hm,dwdt_hm, &
                    u1_hm_map, w1_hm_map)

    call output_host_model_single_variable(u_hm_map, 'U4', 'U_after_adams' , 'm/s', icyc)

    tmp(:, :) = w_hm_map(:,1:nzm)
    call output_host_model_single_variable(tmp, 'W4', 'W_after_adams' , 'm/s', icyc)
    call output_host_model_single_variable(u1_hm_map, 'U41', 'U1_after_adams' , 'm/s', icyc)

    tmp(:, :) = w1_hm_map(:,1:nzm)
    call output_host_model_single_variable(tmp, 'W41', 'W1_after_adams' , 'm/s', icyc)
    ! 5) 标量平流
    call advect_scalars_hm(t_hm_map, u1_hm_map, w1_hm_map)
    call advect_scalars_hm(q_hm_map, u1_hm_map, w1_hm_map)
 
    call output_host_model_single_variable(t_hm_map, 't4', 't_after_advect' , 'K', icyc)
    call output_host_model_single_variable(q_hm_map, 'q4', 'q_after_advect' , 'kg/kg', icyc)
  end do

  u_hm_updated_map_save = u_hm_map
  t_hm_updated_map_save = t_hm_map
  q_hm_updated_map_save = q_hm_map 

  t_out_map = t_hm_map - t_hm_map_save
  q_out_map = q_hm_map - q_hm_map_save
 
  call face2center_U((u_hm_map-u_hm_map_save), u_out_map)
  w_out_map = w_hm_map


  call output_host_model(u0_in, t0_in, q0_in,  &
                        tabs0_in, qn0_in, qp0_in, &
                        qni0_in, qnl0_in, qpi0_in, qpl0_in,  &
                        u_out_map,  w_out_map, t_out_map, q_out_map)



 
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



subroutine damping_hm(u_hm_map, w_hm_map, dudt_hm, dwdt_hm)
    use vars
    implicit none

    real, intent(in)  :: u_hm_map(nsx, nzm), w_hm_map(nsx, nz)
    real, intent(inout) :: dudt_hm(nsx, nzm), dwdt_hm(nsx, nz)

    real :: u0_entire_domain(nzm)
    real :: w0_entire_domain(nz)

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

    do k=nzm,nzm-n_damp,-1
        tau(k) = tau_min *(tau_max/tau_min)**((z(nzm)-z(k))/(z(nzm)-z(nzm-n_damp)))
        tau(k)=1./tau(k)
    end do
   
    do k = 1, nzm
        u0_entire_domain(k) = sum( u_hm_map(:,k) ) / nsx
    end do

    do k = 1, nz
        w0_entire_domain(k) = sum( w_hm_map(:,k) ) / nsx
    end do



    !  ---------------------------------------damp to remove upper-layer waves-------------------------------------------------
    do k = nzm, nzm-n_damp, -1
        do i=1,nsx
            dudt_hm(i,k)= dudt_hm(i,k)- (u_hm_map(i,k) - u0_entire_domain(k)) * tau(k)
            dwdt_hm(i,k)= dwdt_hm(i,k)- (w_hm_map(i,k) - w0_entire_domain(k)) * tau(k)
        end do
    end do 


    ! ---------------------------------------damp to remove layer-mean "dapgRM"------------------------------------------------
    do k = 1,nzm
        do i = 1, nsx
            dudt_hm(i,k) = dudt_hm(i,k) - u0_entire_domain(k) /(20.0*24.0*3600.0)
        end do
    end do

    do k = 1,nz
        do i = 1, nsx
            dwdt_hm(i,k) = dwdt_hm(i,k) - w0_entire_domain(k)  /(20.0*24.0*3600.0)
        end do
    end do
  ! ----------------------------------------------------------------------------------------------------------
end subroutine damping_hm



!================== 动量平流：2D 二阶中心 ==================
subroutine advect_mom_hm(u_hm_map, w_hm_map, dudt_hm, dwdt_hm)
  implicit none
  ! 输入 
  real, intent(in)  :: u_hm_map(nsx, nzm), w_hm_map(nsx, nz)

  ! 输出
  real, intent(inout) :: dudt_hm(nsx, nzm), dwdt_hm(nsx, nz)

  ! 局部
  real :: dx25, dz25, irho_w, irho_k, irhow_k
  real fu(nsx,nzm), fw(nsx, nzm)     ! x 向通量
  real fuz(nsx, nz), fwz(nsx, nzm)  ! z 向通量（注意 f*u/f*v 在 w 层，大小 nz）
  integer :: i, ic, ib, k, kc, kb, kcu

  !---- 初始化清零 ----
  fu = 0.; fw = 0.
  ! fuz / fvz 在 w 层定义：k=1…nz；边界 k=1,nz 设 0
  fuz(:,1) = 0.; fwz(:,1) = 0.
  fuz(:,nz) = 0.; fwz(:,nzm) = 0.

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
      
      ! advect2
      ! if ((u_hm_map(i,k)+u_hm_map(ic,k)) >= 0.0) then
      !   fu(i,k) =  0.5*(u_hm_map(i,k)+u_hm_map(ic,k)) * u_hm_map(i,k)/dx_hm ! 从左（上游）
      ! else
      !   fu(i,k) = 0.5*(u_hm_map(i,k)+u_hm_map(ic,k)) * u_hm_map(ic,k)/dx_hm  ! 从右（上游）
      ! end if
     
      !原来的
      fu(i,k) = dx25 * (u_hm_map(ic,k)+u_hm_map(i,k)) * (u_hm_map(i,k)+u_hm_map(ic,k))

      !advect22
      ! if (( u_hm_map(ic,k)*rho(k)*adz(k) + u_hm_map(ic,kcu)*rho(kcu)*adz(kcu) ) >= 0.0) then
      !   fw(i,k) =  0.5*( u_hm_map(ic,k)*rho(k)*adz(k) + u_hm_map(ic,kcu)*rho(kcu)*adz(kcu) ) * &
      !                  w_hm_map(i,kc)/dx_hm ! 从左（上游）
      ! else
      !   fw(i,k) = 0.5*( u_hm_map(ic,k)*rho(k)*adz(k) + u_hm_map(ic,kcu)*rho(kcu)*adz(kcu) ) * &
      !                  w_hm_map(ic,kc)/dx_hm  ! 从右（上游）
      ! end if

      !原来的
      fw(i,k) = dx25 * ( u_hm_map(ic,k)*rho(k)*adz(k) + u_hm_map(ic,kcu)*rho(kcu)*adz(kcu) ) * &
                        ( w_hm_map(i,kc) + w_hm_map(ic,kc) )
    end do
    do i = 1, nsx
      ib = i - 1
      if (ib < 1) ib = nsx + ib
      dudt_hm(i,k)   = dudt_hm(i,k)   - ( fu(i,k) - fu(ib,k) )
      dwdt_hm(i,kc)  = dwdt_hm(i,kc)  - irho_w * ( fw(i,k) - fw(ib,k) )
    end do
  end do

  !==================== z 向通量 ====================
  
  do k = 2, nzm
    kb = k - 1
    do i = 1, nsx
      ib = i - 1; if (ib < 1) ib = nsx + ib
      !advect22z2
      ! if (( w_hm_map(i,k) + w_hm_map(ib,k) ) >= 0.0) then
      !   fuz(i,k) = 0.5*rhow(k) * ( w_hm_map(i,k) + w_hm_map(ib,k) ) * u_hm_map(i,kb)/dz ! 从下（上游）
      ! else
      !   fuz(i,k) = 0.5*rhow(k) * ( w_hm_map(i,k) + w_hm_map(ib,k) ) * u_hm_map(i,k)/dz  ! 从上（上游）
      ! end if

      !原来的
      fuz(i,k) = dz25 * rhow(k) * ( w_hm_map(i,k) + w_hm_map(ib,k) ) * ( u_hm_map(i,k) + u_hm_map(i,kb) )
    end do
  end do
  ! 顶层接口 k = nz (=nzm+1) 自然保持 0

  ! 应用到 dudt_hm/dvdt_hm；fwz 给 dwdt_hm
  do k = 1, nzm
    kc = k + 1
    irho_k = 1.0 / ( rho(k) * adz(k) )
    do i = 1, nsx
      dudt_hm(i,k) = dudt_hm(i,k) - ( fuz(i,kc) - fuz(i,k) ) * irho_k
      !advect22z22
      ! if (( w_hm_map(i,kc)*rhow(kc) + w_hm_map(i,k)*rhow(k) ) >= 0.0) then
      !   fwz(i,k) = 0.5*( w_hm_map(i,kc)*rhow(kc) + w_hm_map(i,k)*rhow(k) ) * w_hm_map(i,k)/dz ! 从下（上游）
      ! else
      !   fwz(i,k) = 0.5*( w_hm_map(i,kc)*rhow(kc) + w_hm_map(i,k)*rhow(k) ) * w_hm_map(i,kc)/dz  ! 从上（上游）
      ! end if
      !原来的
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
                                dudt_hm, dwdt_hm, p_phys)
  use, intrinsic :: iso_fortran_env, only: real64
  use vars, only: do_3step_adams, dudt_hm_hist, dwdt_hm_hist
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
  if (do_3step_adams) then
    stage = min(3, hm_step + 1)
    select case(stage)
    case (1)
      atc = 1.0 ; btc = 0.0         ; ctc = 0.0
    case (2)
      atc = 1.5 ; btc = -0.5        ; ctc = 0.0
    case default
      atc = 23.0/12.0 ; btc = -16.0/12.0 ; ctc = 5.0/12.0
    end select

    dta  = 1.0/dt_hm_subcycle/atc
    btat = btc/atc
    ctat = ctc/atc
  end if

  rdx  = 1.0/dx_hm

  ! --------- RHS (press_rhs) — 2D(x,z) 与 SAM 一致的形式 ---------
  do k = 1, nzm
    rdz = 1.0/(adz(k)*dz)
    rup = rhow(k+1)/rho(k) * rdz    ! 上界面系数 (kc=k+1)
    rdn = rhow(k  )/rho(k) * rdz    ! 下界面系数 (k)

    do i = 1, nsx
      ip = i + 1; if (ip > nsx) ip = ip-nsx   

      ! 有历史几步加速度信息时
      if (do_3step_adams) then
        rhs(i,k) = &
          ( rdx*(u_hm_map(ip,k) - u_hm_map(i,k)) + ( w_hm_map(i,k+1)*rup - w_hm_map(i,k)*rdn ) )*dta  &
        + ( rdx*(dudt_hm(ip,k) - dudt_hm(i,k)) + ( dwdt_hm(i,k+1)*rup - dwdt_hm(i,k)*rdn ) )  &
        + btat*( rdx*(dudt_hm_hist(ip,k,1) - dudt_hm_hist(i,k,1))  &
              + ( dwdt_hm_hist(i,k+1,1)*rup - dwdt_hm_hist(i,k,1)*rdn ) )              &
        + ctat*( rdx*(dudt_hm_hist(ip,k,2) - dudt_hm_hist(i,k,2))  &
              + ( dwdt_hm_hist(i,k+1,2)*rup - dwdt_hm_hist(i,k,2)*rdn ) )
      else
        rhs(i,k) = &
          ( rdx*(u_hm_map(ip,k) - u_hm_map(i,k)) + ( w_hm_map(i,k+1)*rup - w_hm_map(i,k)*rdn ) ) /dt_hm_subcycle  &
        + ( rdx*(dudt_hm(ip,k) - dudt_hm(i,k)) + ( dwdt_hm(i,k+1)*rup - dwdt_hm(i,k)*rdn ) )
      end if 

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
subroutine adams_hm(u, w, dudt_hm, dwdt_hm,u1, w1)
  use vars, only: do_3step_adams, dudt_hm_hist, dwdt_hm_hist
  implicit none
  real, intent(inout) :: u(nsx, nzm), w(nsx, nz)
  real, intent(in)    :: dudt_hm(nsx, nzm),  dwdt_hm(nsx, nz)
  real, intent(out)   :: u1(nsx, nzm), w1(nsx, nz)

  real :: at, bt, ct
  real :: dtdx, dtdz, rhox, rhoy, rhoz , a1, a2
  integer :: i,k

  hm_step = hm_step + 1

  u1(:,:) = u(:,:)
  w1(:,:) = w(:,:)

  if (do_3step_adams) then
    if (hm_step == 1) then
      at=1.0; bt=0.0; ct=0.0
    else if (hm_step == 2) then
      at=1.5; bt=-0.5; ct=0.0
    else
      at=23.0/12.0; bt=-16.0/12.0; ct=5.0/12.0
    end if

    do k=1,nzm; do i=1,nsx
      dudt_hm_hist(i,k,3) = dudt_hm_hist(i,k,2)
      dudt_hm_hist(i,k,2) = dudt_hm_hist(i,k,1)
      dudt_hm_hist(i,k,1) = dudt_hm(i,k)
    end do; end do
    do k=1,nz; do i=1,nsx
      dwdt_hm_hist(i,k,3) = dwdt_hm_hist(i,k,2)
      dwdt_hm_hist(i,k,2) = dwdt_hm_hist(i,k,1)
      dwdt_hm_hist(i,k,1) = dwdt_hm(i,k)
    end do; end do

    do k=1,nzm; do i=1,nsx
      u(i,k) = u1(i,k) + dt_hm_subcycle*( at*dudt_hm_hist(i,k,1) + bt*dudt_hm_hist(i,k,2) + ct*dudt_hm_hist(i,k,3) )
      w(i,k) = w1(i,k) + dt_hm_subcycle*( at*dwdt_hm_hist(i,k,1) + bt*dwdt_hm_hist(i,k,2) + ct*dwdt_hm_hist(i,k,3) )  !原代码只更新到nzm
    end do; end do

  else
    do k=1,nzm; do i=1,nsx
      u(i,k) = u1(i,k) + dt_hm_subcycle*dudt_hm(i,k) 
      w(i,k) = w1(i,k) + dt_hm_subcycle*dwdt_hm(i,k)
    end do; end do
  end if
    
  ! compute time averaged velocties for second-order advection of scalars:
  dtdx = dt_hm_subcycle/dx_hm
  dtdz = dt_hm_subcycle/dz
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
    
  ! call t_startf ('nudging_hm')

  tnudge = 0.
  qnudge = 0.
  unudge = 0.
  vnudge = 0.

  coef = 1./dt_hm

  if(donudging_uv) then
      do k=1,nzm
        if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
          unudge(k)=unudge(k) - (-ug0_hm(k))*coef
          vnudge(k)=vnudge(k) - (v0(k)-vg0(k))/tauls
          do j=1,ny
            do i=1,nx
              dudt(i,j,k,na)=dudt(i,j,k,na)-(-ug0_hm(k))*coef
              dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))/tauls
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
          tnudge(k)=tnudge(k) -(-tg0_hm(k))*coef
          do j=1,ny
            do i=1,nx
              t(i,j,k)=t(i,j,k)-(-tg0_hm(k))*coef1
            end do
          end do
        end if
      end do
  endif

  if(donudging_tq.or.donudging_q) then
      coef1 = dtn / dt_hm
      do k=1,nzm
        if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
          qnudge(k)=qnudge(k) -(-qg0_hm(k))*coef
          do j=1,ny
            do i=1,nx
                micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)-(-qg0_hm(k))*coef1
            end do
          end do
        end if
      end do
  endif


end subroutine nudging_hm




! 只nudge T q
subroutine nudging_hm_nouv()
  ! keep nudging terms constant in each host model time step
	
  use vars
  use params
  use microphysics, only: micro_field, index_water_vapor
  implicit none

  real coef, coef1
  integer i,j,k
    

  tnudge = 0.
  qnudge = 0.
  unudge = 0.
  vnudge = 0.

  coef = 1./dt_hm

  if(donudging_uv) then
      do k=1,nzm
        if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
          unudge(k)=unudge(k) - (u0(k)-ug0(k))/tauls
          vnudge(k)=vnudge(k) - (v0(k)-vg0(k))/tauls
          do j=1,ny
            do i=1,nx
              dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))/tauls
              dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))/tauls
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
          tnudge(k)=tnudge(k) -(-tg0_hm(k))*coef
          do j=1,ny
            do i=1,nx
              t(i,j,k)=t(i,j,k)-(-tg0_hm(k))*coef1
            end do
          end do
        end if
      end do
  endif

  if(donudging_tq.or.donudging_q) then
      coef1 = dtn / dt_hm
      do k=1,nzm
        if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
          qnudge(k)=qnudge(k) -(-qg0_hm(k))*coef
          do j=1,ny
            do i=1,nx
              micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)-(-qg0_hm(k))*coef1
            end do
          end do
        end if
      end do
  endif


end subroutine nudging_hm_nouv

     
subroutine output_host_model(u0_in, t0_in, q0_in,  &
                        tabs0_in, qn0_in, qp0_in, &
                        qni0_in, qnl0_in, qpi0_in, qpl0_in,  &
                        u_out_map,  w_out_map, t_out_map, q_out_map)        ! 输出形状是（nsx，nzm ）的变量
	
    use vars

    implicit none
    real, intent(in) :: u0_in(nsx, nzm), t0_in(nsx, nzm), q0_in(nsx, nzm)
    real, intent(in) :: tabs0_in(nsx, nzm), qn0_in(nsx, nzm), qp0_in(nsx, nzm)
    real, intent(in) :: qni0_in(nsx, nzm), qnl0_in(nsx, nzm), qpi0_in(nsx, nzm), qpl0_in(nsx, nzm)
    real, intent(in) :: u_out_map(nsx, nzm), w_out_map(nsx, nz), t_out_map(nsx, nzm), q_out_map(nsx, nzm)
    
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


    nfields_hm=14 ! number of 3D fields to save
    nfields1_hm=0

    sepchar=""

    write(rankchar,'(i4)') 1 !nsubdomains
    write(timechar,'(i10)') nstep
    do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
    end do

    ! print*, 'Rank=', rank, '*************begin output_host_model***************'

    filetype = '.bin2D'
    filename='./OUT_3D/hm_output_'//trim(case)//'_'//trim(caseid)//&
    filetype//sepchar
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
    name='U0_In'
    long_name='Input U from subdomain'
    units='m/s'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

   
    

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=t0_in(i,k)
        end do
    end do
    name='T0_In'
    long_name='Input Liquid/Ice Water Static Energy from subdomain'
    units='K'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=q0_in(i,k)
        end do
    end do
    name='Q0_In'
    long_name='Input Water Vapor from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=tabs0_in(i,k)
        end do
    end do
    name='Tabs0_In'
    long_name='Input Tabs from subdomain'
    units='K'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=qn0_in(i,k)
        end do
    end do
    name='Qn0_In'
    long_name='Input Qn0 from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)
    

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=qp0_in(i,k)
        end do
    end do
    name='Qp0_In'
    long_name='Input Qp0 from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=qni0_in(i,k)
        end do
    end do
    name='Qni0_In'
    long_name='Input Qni0 from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=qnl0_in(i,k)
        end do
    end do
    name='Qnl0_In'
    long_name='Input Qnl0 from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=qpi0_in(i,k)
        end do
    end do
    name='Qpi0_In'
    long_name='Input Qpi0 from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=qpl0_in(i,k)
        end do
    end do
    name='Qpl0_In'
    long_name='Input Qpl0 from subdomain'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)


    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=u_out_map(i,k)
        end do
    end do
    name='U_Out'
    long_name='Output U'
    units='m/s'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)


    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=w_out_map(i,k)
        end do
    end do
    name='W_OUT'
    long_name='Output W W_hm_map'
    units='m/s'

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=t_out_map(i,k)
        end do
    end do
    name='T_Out'
    long_name='Output Liquid/Ice Water Static Energy t'
    units='K'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    nfields1_hm=nfields1_hm+1
    do k=1,nzm
        do i=1,nsx
            tmp(i,1,k)=q_out_map(i,k)
        end do
    end do
    name='Q_Out'
    long_name='Output Water Vapor q'
    units='kg/kg'
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)

    
    call compress3D_hm(tmp,nsx,1,nzm,name,long_name,units)


    if(nfields_hm.ne.nfields1_hm) then
        ! print*,'host model write_fields3D error: nfields_hm=',nfields_hm,'  nfields1_hm=',nfields1_hm
        call task_abort()
    end if

    close (46)

    print*, 'Appending 3D data. file:'//filename

end subroutine output_host_model



subroutine output_host_model_single_variable(u0_in, v_name,v_longname,v_unit,icyc)
	
    use vars
    use params,only: nstephostmodel

    implicit none
    real, intent(in) :: u0_in(nsx, nzm)
    character(*), intent(in) :: v_name, v_longname, v_unit
    integer, intent(in) :: icyc
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
    write(timechar,'(i10)') nstep + icyc/hm_subcycle*nstephostmodel
    do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
    end do

    print*, 'Rank=', rank, '*************begin output_host_model***************'

    filetype = '.bin2D'
    filename='./OUT_3D/hm_output_'//trim(case)//'_'//trim(caseid)//&
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
    
    write(46) real(REAL(nstep + REAL(icyc)/REAL(hm_subcycle)*REAL(nstephostmodel))*dt/(3600.*24.)+day0, 4)

  
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



subroutine compress3D_hm(f,nx,ny,nz,name, long_name, units)
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

subroutine center2face_U(u_center_map, u_face_map)
    use grid, only: nsx, nzm
    implicit none
    real, intent(in)   :: u_center_map(nsx,nzm)
    real, intent(out)  :: u_face_map(nsx,nzm)
    u_face_map(1,     :) = 0.5 * (u_center_map(1,     :) + u_center_map(nsx,     :))
    u_face_map(2:nsx, :) = 0.5 * (u_center_map(2:nsx, :) + u_center_map(1:nsx-1, :))
end subroutine center2face_U

subroutine face2center_U(u_face_map, u_center_map)
    use grid, only: nsx, nzm
    implicit none
    real, intent(in)   :: u_face_map(nsx,nzm)
    real, intent(out)  :: u_center_map(nsx,nzm)
    u_center_map(1:nsx-1, :) = 0.5 * (u_face_map(1:nsx-1, :) + u_face_map(2:nsx, :))
    u_center_map(nsx,     :) = 0.5 * (u_face_map(nsx,     :) + u_face_map(1,     :))
end subroutine face2center_U

subroutine hot_bubble(hm_step, t)
    use grid, only: nsx, nzm
    implicit none
    integer, intent(in) :: hm_step       
    real, intent(inout) :: t(nsx, nzm)    
    
    real :: amp, x_wave, z_wave, x_center
    integer :: i, k
    real, parameter :: period = 20.0   
    real, parameter :: base_amp = 1.0  
    real, parameter :: z_scale = 8.0   
    real, parameter :: x_scale = 2.0   
    real, parameter :: pi = acos(-1.0)

    x_center = real(nsx) / 2.0
    amp = base_amp * sin(pi * min(1.,real(hm_step) / period))

 
    do k = 1, nzm
        z_wave = sin(pi * min(1.,real(k)/z_scale))
        do i = 1, nsx
            x_wave = exp(-((real(i)-x_center)**2) / (x_scale)**2)
            t(i,k) = t(i,k) + amp * x_wave * z_wave
        end do
    end do

end subroutine hot_bubble

subroutine buoyancy_only_in_hm(t_hm_map, q_hm_map, dwdt_hm)   !  qni_hm_map, qnl_hm_map, qpi_hm_map, qpl_hm_map,
  use vars
  use params
  implicit none
  real, intent(in)  :: t_hm_map(nsx, nzm), q_hm_map(nsx, nzm)
  real, intent(inout) :: dwdt_hm(nsx, nz)
  real :: tabs0_entire_domain(nzm), qv0_entire_domain(nzm)
  real :: tabs_map(nsx,nzm)
  integer i,k,kb
  real betu, betd
	
  do k = 1, nzm
      do i = 1, nsx   ! tabs(i,j,k) = t(i,j,k)-gamaz(k)+ fac_cond * (qcl(i,j,k)+qpl(i,j,k)) +fac_sub *(qci(i,j,k) + qpi(i,j,k))
        tabs_map(i,k) = t_hm_map(i,k) - gamaz(k)
      end do
      tabs0_entire_domain(k) = sum( tabs_map(:,k) ) / nsx
      qv0_entire_domain(k) = sum( q_hm_map(:,k) ) / nsx    
  end do
  do k=2,nzm	
    kb=k-1
    betu=adz(kb)/(adz(k)+adz(kb))
    betd=adz(k)/(adz(k)+adz(kb))
    do i=1,nsx
      dwdt_hm(i,k)=dwdt_hm(i,k) +  &
          bet(k)*betu* &
        ( tabs0_entire_domain(k)*(epsv*(q_hm_map(i,k)-qv0_entire_domain(k))) &
          +(tabs_map(i,k)-tabs0_entire_domain(k))*(1.+epsv*qv0_entire_domain(k)) ) &
        + bet(kb)*betd* &
        ( tabs0_entire_domain(kb)*(epsv*(q_hm_map(i,kb)-qv0_entire_domain(kb))) &
          +(tabs_map(i,kb)-tabs0_entire_domain(kb))*(1.+epsv*qv0_entire_domain(kb)) )  
    end do
  end do

end subroutine buoyancy_only_in_hm




subroutine diffuse_u(u_map,dudt_hm)
    use vars
    implicit none
    real, intent(in)   :: u_map(nsx,nzm)
    real, intent(inout)   :: dudt_hm(nsx,nzm)

    integer i,ic,ib,k
   
    do k = 1,nzm-2
      do i=1,nsx
        ic = i + 1
        if (ic > nsx) ic = ic - nsx
        ib = i - 1
        if (ib < 1) ib = ib + nsx
        dudt_hm(i,k) = dudt_hm(i,k) + diffuse_intensity*(u_map(ic,k) -2*u_map(i,k) + u_map(ib,k))/dt_hm_subcycle
      end do
    end do
end subroutine diffuse_u

subroutine diffuse_w(w_map,dwdt_hm)
    use vars
    implicit none
    real, intent(in)   :: w_map(nsx,nz)
    real, intent(inout)   :: dwdt_hm(nsx,nz)
    integer i,ic,ib,k
   
    do k = 1,nz
      do i=1,nsx
        ic = i + 1
        if (ic > nsx) ic = ic - nsx
        ib = i - 1
        if (ib < 1) ib = ib + nsx
        dwdt_hm(i,k) = dwdt_hm(i,k) + diffuse_intensity*(w_map(ic,k) -2*w_map(i,k) + w_map(ib,k))/dt_hm_subcycle
      end do
    end do
end subroutine diffuse_w

subroutine diffuse_TQ(t_map)
    use vars
    implicit none
    real, intent(inout)   :: t_map(nsx,nzm)

    integer i,ic,ib,k
   
    do k = 1,nzm
      do i=1,nsx
        ic = i + 1
        if (ic > nsx) ic = ic - nsx
        ib = i - 1
        if (ib < 1) ib = ib + nsx
        t_map(i,k) = t_map(i,k) + diffuse_intensity*(t_map(ic,k) -2*t_map(i,k) + t_map(ib,k))
      end do
    end do

end subroutine diffuse_TQ



subroutine modify_U_for_subdomain()
    use vars
    implicit none
    integer i,j,k

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          u(i,j,k) = u(i,j,k)+ ug0_press_modify(k)
        end do
      end do
    end do
end subroutine modify_U_for_subdomain

! ----------------------傅里叶变换消最高频------------------------------------------
subroutine damp_highest_wavenumber(u_map)
    use vars
    implicit none

    real, intent(inout) :: u_map(nsx, nzm)

    complex, allocatable :: u_fft(:)  
    real, allocatable :: temp_row(:)   
    integer j, max_wavenumber_index
  
    allocate(u_fft(nsx))
    allocate(temp_row(nsx))
    
    max_wavenumber_index = nsx / 2 + 1
   
    do j = 1, nzm
        temp_row = u_map(:, j)

        call dft_1d(temp_row, u_fft, nsx)
        
        u_fft(max_wavenumber_index) = cmplx(0.0, 0.0) 

        call idft_1d(u_fft, u_map(:, j), nsx) 
    end do

    deallocate(u_fft)
    deallocate(temp_row)

end subroutine damp_highest_wavenumber


subroutine dft_1d(x_in, x_out, N_size)
    implicit none
    integer, intent(in) :: N_size
    real, intent(in)    :: x_in(N_size)
    complex, intent(out) :: x_out(N_size)

    integer n, k
    real    :: arg
    real    :: pi
    pi = acos(-1.0)
    
    ! 外层循环：遍历输出频率 k (波数)
    do k = 1, N_size
        x_out(k) = cmplx(0.0, 0.0)
        
        ! 内层循环：遍历输入时间/空间点 n
        do n = 1, N_size
            ! 计算角度：-2 * pi * (k-1) * (n-1) / N
            arg = -2.0 * pi * real( (k-1) * (n-1) ) / real(N_size)
            x_out(k) = x_out(k) + x_in(n) * cmplx(cos(arg), sin(arg))

        end do
    end do
end subroutine dft_1d



subroutine idft_1d(x_in_complex, x_out_real, N_size)
    implicit none
    integer, intent(in) :: N_size
    complex, intent(in) :: x_in_complex(N_size)
    real, intent(out)   :: x_out_real(N_size)

    integer n, k
    real    :: arg
    complex :: sum_val
    real    :: pi
    pi = acos(-1.0)


    do n = 1, N_size
        sum_val = cmplx(0.0, 0.0)
        
        ! 内层循环：遍历频率 k
        do k = 1, N_size
            ! IDFT 的角度是正号：2 * pi * (k-1) * (n-1) / N
            arg = 2.0 * pi * real( (k-1) * (n-1) ) / real(N_size)
            
            ! IDFT 定义：Xn = (1/N) * SUM [ Xk * e^(i * arg) ]
            sum_val = sum_val + x_in_complex(k) * cmplx(cos(arg), sin(arg))

        end do
        
        ! IDFT 还需要除以 N
        x_out_real(n) = real(sum_val) / real(N_size)
    end do
end subroutine idft_1d

end module module_hostmodel
