subroutine hm_couple_step()
    use vars
    use microphysics, only: micro_field, index_water_vapor
    use module_hostmodel, only: host_model_evolve
    implicit none

    integer :: k, m, is_hm
    integer :: nother_hm, nexp_zm_hm, nexp_w_hm
    integer :: ireq_hm, rf_hm, tag_hm
    integer :: jzm, jw
    integer :: req_u, req_v, req_t, req_q, req_w
    logical :: done
    integer :: nleft
    integer :: cnt_zm, cnt_w
   

    real, allocatable :: u0_map(:,:), v0_map(:,:), t0_map(:,:), q0_map(:,:)
    real, allocatable :: wsub_map(:,:)
    real, allocatable :: u_hm_map(:,:), v_hm_map(:,:), t_hm_map(:,:), q_hm_map(:,:)
    real, allocatable :: w_hm_map(:,:)

    real,    allocatable :: rbuf_zm_hm(:,:)      ! (nzm, nexp_zm_hm)
    integer, allocatable :: reqs_zm_hm(:)
    logical, allocatable :: done_zm_hm(:)

    real,    allocatable :: rbuf_w_hm(:,:)       ! (nz, nexp_w_hm)
    integer, allocatable :: reqs_w_hm(:)
    logical, allocatable :: done_w_hm(:)


    do k = 1, nzm
        u0_local_hm(k) = sum( u(1:nx,1:ny,k) ) / real(nx*ny)   !! nx是每个subdomain里的x格点数
        v0_local_hm(k) = sum( v(1:nx,1:ny,k) ) / real(nx*ny)
        t0_local_hm(k) = sum( t(1:nx,1:ny,k) ) / real(nx*ny)
        q0_local_hm(k) = sum( micro_field(1:nx,1:ny,k,index_water_vapor) ) / real(nx*ny)
    end do

    
    if (.not. masterproc) then
        call task_bsend_float(0, u0_local_hm, nzm, 101)
        call task_bsend_float(0, v0_local_hm, nzm, 102)
        call task_bsend_float(0, t0_local_hm, nzm, 103)
        call task_bsend_float(0, q0_local_hm, nzm, 104)
        call task_bsend_float(0, wsub_subdomain, nz, 105)

    end if


    if (masterproc) then
        
        
        allocate(u0_map(nsx, nzm), v0_map(nsx, nzm),  &
                t0_map(nsx, nzm), q0_map(nsx, nzm))
        allocate(wsub_map(nsx, nz))
        allocate(u_hm_map(nsx, nzm), v_hm_map(nsx, nzm),  &
                t_hm_map(nsx, nzm), q_hm_map(nsx, nzm))
        allocate(w_hm_map(nsx, nz))
        

        ! master 自己这一列先放进去（约定 rank 0 → is_hm=1）
        u0_map(1,:)   = u0_local_hm(:)
        v0_map(1,:)   = v0_local_hm(:)
        t0_map(1,:)   = t0_local_hm(:)
        q0_map(1,:)   = q0_local_hm(:)
        wsub_map(1,:) = wsub_subdomain(:)

        ! 其余需要接收的消息数
        nother_hm  = max(0, nsx-1)
        nexp_zm_hm = nother_hm * 4          ! u/v/t/q, 每列 4 条，长度 nzm
        nexp_w_hm  = nother_hm * 1          ! wsub,    每列 1 条，长度 nz

        if (nexp_zm_hm > 0 .or. nexp_w_hm > 0) then

            if (nexp_zm_hm > 0) then
                allocate(rbuf_zm_hm(nzm, nexp_zm_hm), reqs_zm_hm(nexp_zm_hm), done_zm_hm(nexp_zm_hm))
                done_zm_hm = .false.
            end if
            if (nexp_w_hm > 0) then
                allocate(rbuf_w_hm(nz, nexp_w_hm), reqs_w_hm(nexp_w_hm), done_w_hm(nexp_w_hm))
                done_w_hm = .false.
            end if

            
            do is_hm = 2, nsx
                if (nexp_zm_hm > 0) then
                    jzm=jzm+1; call task_receive_float(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm))  ! 101
                    jzm=jzm+1; call task_receive_float(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm))  ! 102
                    jzm=jzm+1; call task_receive_float(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm))  ! 103
                    jzm=jzm+1; call task_receive_float(rbuf_zm_hm(1,jzm), nzm, reqs_zm_hm(jzm))  ! 104
                end if
                if (nexp_w_hm > 0) then
                    jw=jw+1;   call task_receive_float(rbuf_w_hm(1,jw), nz, reqs_w_hm(jw))       ! 105
                end if
            end do

            ! 2) 轮询直到全部完成；按 tag_hm+来源 rf_hm 分拣
            
            cnt_zm = nexp_zm_hm
            cnt_w  = nexp_w_hm

            do while (cnt_zm > 0 .or. cnt_w > 0)
                if (nexp_zm_hm > 0) then
                    do ireq_hm = 1, nexp_zm_hm
                        if (.not. done_zm_hm(ireq_hm)) then
                            call task_test(reqs_zm_hm(ireq_hm), done_zm_hm(ireq_hm), rf_hm, tag_hm)  ! 返回rf_hm(rank), tag_hm
                            if (done_zm_hm(ireq_hm)) then
                                is_hm = rf_hm + 1
                                select case (tag_hm)
                                case (101); u0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
                                case (102); v0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
                                case (103); t0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
                                case (104); q0_map(is_hm,:) = rbuf_zm_hm(:,ireq_hm)
                                case default
                                    ! 忽略与 nzm 组不匹配的 tag_hm（理论上不会来 105）
                                end select
                                cnt_zm = cnt_zm - 1
                            end if
                        end if
                    end do
                end if

                if (nexp_w_hm > 0) then
                    do ireq_hm = 1, nexp_w_hm
                        if (.not. done_w_hm(ireq_hm)) then
                            call task_test(reqs_w_hm(ireq_hm), done_w_hm(ireq_hm), rf_hm, tag_hm)
                            if (done_w_hm(ireq_hm)) then
                                is_hm = rf_hm + 1
                                if (tag_hm == 105) then
                                    wsub_map(is_hm,:) = rbuf_w_hm(:,ireq_hm)
                                end if 
                                cnt_w = cnt_w - 1
                            end if
                        end if
                    end do
                end if
            end do

            if (nexp_zm_hm > 0) deallocate(rbuf_zm_hm, reqs_zm_hm, done_zm_hm)
            if (nexp_w_hm  > 0) deallocate(rbuf_w_hm , reqs_w_hm , done_w_hm )
        end if  ! 有需要接收的

        ! ==== 至此，u0_map/v0_map/t0_map/q0_map/wsub_map 都齐了 ====



        !------------- 调用 host model -------------
        call host_model_evolve( u0_map=u0_map, v0_map=v0_map, wsub_map=wsub_map, &
                            t0_map=t0_map, q0_map=q0_map,                    &
                            rho=rho, rhow=rhow, adz=adz, adzw=adzw,          &
                            u_hm_map=u_hm_map, v_hm_map=v_hm_map,            &
                            w_hm_map=w_hm_map, t_hm_map=t_hm_map, q_hm_map=q_hm_map )

        ! ================== 发送每列 hm_map → 各子域（绝对值） ==================
        ! tag_hms: 301(u), 302(v), 303(t), 304(q), 305(w)
        

        ! master 自己的列（rank=0 → is_hm=1）直接覆盖
        if (rank .eq. 0) then
            is_hm = 1
            ug0(1:nzm)          = u_hm_map(is_hm,1:nzm)
            vg0(1:nzm)          = v_hm_map(is_hm,1:nzm)
            tg0(1:nzm)          = t_hm_map(is_hm,1:nzm)
            qg0(1:nzm)          = q_hm_map(is_hm,1:nzm)
            wsub_subdomain(1:nz)= w_hm_map(is_hm,1:nz)
        end if

        ! 其它子域：阻塞发送 5 条列向量到对应 rank
        do m = 1, nsubdomains-1  ! m=0是master本身
            call task_bsend_float(m, u_hm_map(m+1,1), nzm, 301)
            call task_bsend_float(m, v_hm_map(m+1,1), nzm, 302)
            call task_bsend_float(m, t_hm_map(m+1,1), nzm, 303)
            call task_bsend_float(m, q_hm_map(m+1,1), nzm, 304)
            call task_bsend_float(m, w_hm_map(m+1,1),  nz, 305)
        end do

        !------------- 清理分配的缓冲 -------------

        deallocate(u0_map, v0_map, t0_map, q0_map)
        deallocate(u_hm_map, v_hm_map, t_hm_map, q_hm_map)
        deallocate(wsub_map, w_hm_map)

    end if  ! masterproc

    !------------------------------------------------------------
    ! 非 master
    !------------------------------------------------------------
    if (.not. masterproc) then
        real :: ucol(nzm), vcol(nzm), tcol(nzm), qcol(nzm), wcol(nz)
        integer :: req_u, req_v, req_t, req_q, req_w
        logical :: done
        integer :: rf_hm, tag_hm
        integer :: nleft

        ! 先挂五个非阻塞接收（必须先接收再发送以避免死锁；master 端按顺序 bsend）
        call task_receive_float(ucol, nzm, req_u)  ! tag_hm=301
        call task_receive_float(vcol, nzm, req_v)  ! tag_hm=302
        call task_receive_float(tcol, nzm, req_t)  ! tag_hm=303
        call task_receive_float(qcol, nzm, req_q)  ! tag_hm=304
        call task_receive_float(wcol,  nz,  req_w)  ! tag_hm=305

        ! 轮询全部完成
        nleft = 5
        do while (nleft > 0)
            call task_test(req_u, done, rf_hm, tag_hm); if (done .and. req_u .ne. -1) then; req_u = -1; nleft = nleft - 1; end if
            call task_test(req_v, done, rf_hm, tag_hm); if (done .and. req_v .ne. -1) then; req_v = -1; nleft = nleft - 1; end if
            call task_test(req_t, done, rf_hm, tag_hm); if (done .and. req_t .ne. -1) then; req_t = -1; nleft = nleft - 1; end if
            call task_test(req_q, done, rf_hm, tag_hm); if (done .and. req_q .ne. -1) then; req_q = -1; nleft = nleft - 1; end if
            call task_test(req_w, done, rf_hm, tag_hm); if (done .and. req_w .ne. -1) then; req_w = -1; nleft = nleft - 1; end if
        end do

        
        ug0_hm(1:nzm)           = ucol(1:nzm)
        vg0_hm(1:nzm)           = vcol(1:nzm)
        tg0_hm(1:nzm)           = tcol(1:nzm)
        qg0_hm(1:nzm)           = qcol(1:nzm)
        wsub_subdomain(1:nz) = wcol(1:nz)
    end if

end subroutine hm_couple_step