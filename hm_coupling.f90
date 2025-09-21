subroutine hm_couple_step()
    use vars
    use microphysics, only: micro_field, index_water_vapor
    use module_hostmodel, only: host_model_evolve
    implicit none

    integer :: k, m, is_hm
    integer :: nother_hm, nexp_zm_hm
    integer :: ireq_hm, rf_hm, tag_hm
    integer :: jzm
    integer :: req_u, req_v, req_t, req_q
    logical :: done
    integer :: nleft
    integer :: cnt_zm
    real :: ucol(nzm), vcol(nzm), tcol(nzm), qcol(nzm)
    real :: sendbuf(nzm)
   
   

    real, allocatable :: u0_map(:,:), v0_map(:,:), t0_map(:,:), q0_map(:,:)
    real, allocatable :: u_out_map(:,:), v_out_map(:,:), t_out_map(:,:), q_out_map(:,:)
    real, allocatable :: w_out_map(:,:)
    real, allocatable :: dummy2d(:,:)

    real,    allocatable :: rbuf_zm_hm(:,:)      ! (nzm, nexp_zm_hm)
    integer, allocatable :: reqs_zm_hm(:)
    logical, allocatable :: done_zm_hm(:)

    real, allocatable :: tabs0_map(:,:), qv0_map(:,:), qn0_map(:,:), qp0_map(:,:)

    if (masterproc) then
        allocate(u0_map(nsx, nzm), v0_map(nsx, nzm),  &
                t0_map(nsx, nzm), q0_map(nsx, nzm),  &
                tabs0_map(nsx, nzm), qv0_map(nsx, nzm),  &
                qn0_map(nsx, nzm), qp0_map(nsx, nzm))
                
        allocate(u_out_map(nsx, nzm), v_out_map(nsx, nzm),  &
                t_out_map(nsx, nzm), q_out_map(nsx, nzm))
        allocate(w_out_map(nsx, nz))
    else
        allocate(dummy2d(nsx, nzm))  ! assign dummy to avoid alloc error
    end if

    call t_startf ('host_model')

    do k = 1, nzm
    !     u0_local_hm(k) = sum( u(1:nx,1:ny,k) ) / real(nx*ny)   !! nx是每个subdomain里的x格点数
    !     v0_local_hm(k) = sum( v(1:nx,1:ny,k) ) / real(nx*ny)
    !     t0_local_hm(k) = sum( t(1:nx,1:ny,k) ) / real(nx*ny)
        q0_local_hm(k) = sum( micro_field(1:nx,1:ny,k,index_water_vapor) ) / real(nx*ny)
    end do

    ! u0, v0, t0, q0 were calculated in diagnose.f90 for each subdomain
    ! here we gather them to the masterproc (rank=0) for host model coupling
    ! then the masterproc will distribute the results back to each subdomain

    u0_local_hm = u0
    v0_local_hm = v0
    t0_local_hm = t0
    ! q0_local_hm = q0
    tabs0_local_hm = tabs0
    qv0_local_hm = qv0
    qn0_local_hm = qn0
    qp0_local_hm = qp0

    
    ! gather u0
    if (masterproc) then
        call task_bgather_float_map(0, u0_local_hm(1), nzm, nsx, u0_map)
    else
        call task_bgather_float_map(0, u0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather v0
    if (masterproc) then
        call task_bgather_float_map(0, v0_local_hm(1), nzm, nsx, v0_map)
    else
        call task_bgather_float_map(0, v0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather t0
    if (masterproc) then
        call task_bgather_float_map(0, t0_local_hm(1), nzm, nsx, t0_map)
    else
        call task_bgather_float_map(0, t0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather q0
    if (masterproc) then
        call task_bgather_float_map(0, q0_local_hm(1), nzm, nsx, q0_map)
    else
        call task_bgather_float_map(0, q0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather tabs0
    if (masterproc) then
        call task_bgather_float_map(0, tabs0_local_hm(1), nzm, nsx, tabs0_map)
    else
        call task_bgather_float_map(0, tabs0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather qv0
    if (masterproc) then
        call task_bgather_float_map(0, qv0_local_hm(1), nzm, nsx, qv0_map)
    else
        call task_bgather_float_map(0, qv0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather qn0
    if (masterproc) then
        call task_bgather_float_map(0, qn0_local_hm(1), nzm, nsx, qn0_map)
    else
        call task_bgather_float_map(0, qn0_local_hm(1), nzm, nsx, dummy2d)
    end if
    ! gather qp0
    if (masterproc) then
        call task_bgather_float_map(0, qp0_local_hm(1), nzm, nsx, qp0_map)
    else
        call task_bgather_float_map(0, qp0_local_hm(1), nzm, nsx, dummy2d)
    end if
    
    
    ! if (.not. masterproc) then
    !     call task_bsend_float(0, u0_local_hm(1), nzm, 101)
    !     call task_bsend_float(0, v0_local_hm(1), nzm, 102)
    !     call task_bsend_float(0, t0_local_hm(1), nzm, 103)
    !     call task_bsend_float(0, q0_local_hm(1), nzm, 104)
    ! 
    ! end if

    if (masterproc) then

        ! allocate(u0_map(nsx, nzm), v0_map(nsx, nzm),  &
        !         t0_map(nsx, nzm), q0_map(nsx, nzm))
        ! allocate(u_hm_map(nsx, nzm), v_hm_map(nsx, nzm),  &
        !         t_hm_map(nsx, nzm), q_hm_map(nsx, nzm))
        ! allocate(w_hm_map(nsx, nz))
        

        ! ! master 自己这一列先放进去
        ! u0_map(1,:)   = u0_local_hm(:)
        ! v0_map(1,:)   = v0_local_hm(:)
        ! t0_map(1,:)   = t0_local_hm(:)
        ! q0_map(1,:)   = q0_local_hm(:)


        ! ! 其余需要接收的消息数
        ! nother_hm  = max(0, nsx-1)
        ! nexp_zm_hm = nother_hm * 4          ! u/v/t/q, 每列 4 条，长度 nzm


        ! if (nexp_zm_hm > 0) then
        !     allocate(rbuf_zm_hm(nzm, nexp_zm_hm), reqs_zm_hm(nexp_zm_hm), done_zm_hm(nexp_zm_hm))
        !     done_zm_hm = .false.
            

        !     jzm = 0
        !     do is_hm = 2, nsx
        !         if (nexp_zm_hm > 0) then
        !             jzm=jzm+1; call task_receive_float_tag(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm),101)  ! 101
        !             jzm=jzm+1; call task_receive_float_tag(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm),102)  ! 102
        !             jzm=jzm+1; call task_receive_float_tag(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm),103)  ! 103
        !             jzm=jzm+1; call task_receive_float_tag(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm),104)  ! 104
        !         end if  
        !     end do

        !     ! 2) 轮询直到全部完成；按 tag_hm+来源 rf_hm 分拣
            
        !     cnt_zm = nexp_zm_hm

        !     do while (cnt_zm > 0)
        !         if (nexp_zm_hm > 0) then
        !             do ireq_hm = 1, nexp_zm_hm
        !                 if (.not. done_zm_hm(ireq_hm)) then
        !                     call task_test(reqs_zm_hm(ireq_hm), done_zm_hm(ireq_hm), rf_hm, tag_hm)  ! 返回rf_hm(rank), tag_hm
        !                     if (done_zm_hm(ireq_hm)) then
        !                         is_hm = rf_hm + 1
        !                         select case (tag_hm)
        !                         case (101); u0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
        !                         case (102); v0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
        !                         case (103); t0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
        !                         case (104); q0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
        !                         case default
        !                             ! 忽略与 nzm 组不匹配的 tag_hm（理论上不会来 105）
        !                         end select
        !                         cnt_zm = cnt_zm - 1
        !                     end if
        !                 end if
        !             end do
        !         end if
        !     end do

        !     if (nexp_zm_hm > 0) deallocate(rbuf_zm_hm, reqs_zm_hm, done_zm_hm)

        ! end if  ! 有需要接收的

        ! ==== 至此，u0_map/v0_map/t0_map/q0_map/wsub_map 都齐了 ====



        !------------- 调用 host model -------------
        call host_model_evolve( u0_in=u0_map, v0_in=v0_map, wsub_in=wsub_map, &
                            t0_in=t0_map, q0_in=q0_map,                    &
                            tabs0_in = tabs0_map, qv0_in = qv0_map, qn0_in = qn0_map, qp0_in = qp0_map,   &
                            u_out_map=u_out_map, v_out_map=v_out_map,            &
                            w_out_map=w_out_map, t_out_map=t_out_map, q_out_map=q_out_map )
        
        wsub_map(:, :)         = w_out_map(:, :)

        ! ug0_hm(1:nzm)          = u_hm_map(1,1:nzm)
        ! vg0_hm(1:nzm)          = v_hm_map(1,1:nzm)
        ! tg0_hm(1:nzm)          = t_hm_map(1,1:nzm)
        ! qg0_hm(1:nzm)          = q_hm_map(1,1:nzm)

        ! !----------------send back to subdomains---------------------------
        ! do m = 1, nsubdomains-1  ! m=0是master本身
        !     sendbuf(:) = u_hm_map(m+1,1:nzm)
        !     call task_bsend_float(m, sendbuf(1), nzm, 301)
        !     sendbuf(:) = v_hm_map(m+1,1:nzm)
        !     call task_bsend_float(m, sendbuf(1), nzm, 302)
        !     sendbuf(:) = t_hm_map(m+1,1:nzm)
        !     call task_bsend_float(m, sendbuf(1), nzm, 303)
        !     sendbuf(:) = q_hm_map(m+1,1:nzm)
        !     call task_bsend_float(m, sendbuf(1), nzm, 304)
            
        ! end do

    end if  ! masterproc

    ! distribute ug0_hm
    if (masterproc) then
        call task_bscatter_float_map(0, u_out_map, nzm, nsx, ug0_hm(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, ug0_hm(1))
    end if
    ! distribute vg0_hm
    if (masterproc) then
        call task_bscatter_float_map(0, v_out_map, nzm, nsx, vg0_hm(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, vg0_hm(1))
    end if
    ! distribute tg0_hm
    if (masterproc) then
        call task_bscatter_float_map(0, t_out_map, nzm, nsx, tg0_hm(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, tg0_hm(1))
    end if
    ! distribute qg0_hm
    if (masterproc) then
        call task_bscatter_float_map(0, q_out_map, nzm, nsx, qg0_hm(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, qg0_hm(1))
    end if

    ! !------------------------------------------------------------
    ! ! 非 master
    ! !------------------------------------------------------------
    ! if (.not. masterproc) then
        

    !     ! 先挂五个非阻塞接收（必须先接收再发送以避免死锁；master 端按顺序 bsend）
    !     call task_receive_float_tag(ucol, nzm, req_u,301)  ! tag_hm=301
    !     call task_receive_float_tag(vcol, nzm, req_v,302)  ! tag_hm=302
    !     call task_receive_float_tag(tcol, nzm, req_t,303)  ! tag_hm=303
    !     call task_receive_float_tag(qcol, nzm, req_q,304)  ! tag_hm=304


    !     ! 轮询全部完成
    !     nleft = 4
    !     do while (nleft > 0)
    !         call task_test(req_u, done, rf_hm, tag_hm); if (done .and. req_u .ne. -1) then; req_u = -1; nleft = nleft - 1; end if
    !         call task_test(req_v, done, rf_hm, tag_hm); if (done .and. req_v .ne. -1) then; req_v = -1; nleft = nleft - 1; end if
    !         call task_test(req_t, done, rf_hm, tag_hm); if (done .and. req_t .ne. -1) then; req_t = -1; nleft = nleft - 1; end if
    !         call task_test(req_q, done, rf_hm, tag_hm); if (done .and. req_q .ne. -1) then; req_q = -1; nleft = nleft - 1; end if
            
    !     end do

        
    !     ug0_hm(1:nzm)           = ucol(1:nzm)
    !     vg0_hm(1:nzm)           = vcol(1:nzm)
    !     tg0_hm(1:nzm)           = tcol(1:nzm)
    !     qg0_hm(1:nzm)           = qcol(1:nzm)

    ! end if

    !------------------------------------------------------------
    ! clean up
    !------------------------------------------------------------
    if (masterproc) then
        if (allocated(u0_map))    deallocate(u0_map)
        if (allocated(v0_map))    deallocate(v0_map)
        if (allocated(t0_map))    deallocate(t0_map)
        if (allocated(q0_map))    deallocate(q0_map)
        if (allocated(tabs0_map))    deallocate(tabs0_map)
        if (allocated(qv0_map))    deallocate(qv0_map)
        if (allocated(qn0_map))    deallocate(qn0_map)
        if (allocated(qp0_map))    deallocate(qp0_map)
        if (allocated(u_out_map))  deallocate(u_out_map)  
        if (allocated(v_out_map))  deallocate(v_out_map)  
        if (allocated(t_out_map))  deallocate(t_out_map)  
        if (allocated(q_out_map))  deallocate(q_out_map)  
        if (allocated(w_out_map))  deallocate(w_out_map)  
    else
        if (allocated(dummy2d))  deallocate(dummy2d)
    end if

    call t_stopf('host_model')

end subroutine hm_couple_step
