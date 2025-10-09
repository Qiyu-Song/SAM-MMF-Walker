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
   
   

    real, allocatable :: u0_map(:,:),  t0_map(:,:), q0_map(:,:)
    real, allocatable :: u_out_map(:,:), v_out_map(:,:), t_out_map(:,:), q_out_map(:,:)
    real, allocatable :: w_out_map(:,:)
    real, allocatable :: dummy2d(:,:)

    real,    allocatable :: rbuf_zm_hm(:,:)      ! (nzm, nexp_zm_hm)
    integer, allocatable :: reqs_zm_hm(:)
    logical, allocatable :: done_zm_hm(:)

    real, allocatable :: tabs0_map(:,:), qv0_map(:,:), qn0_map(:,:), qp0_map(:,:)

    if (masterproc) then
        allocate(u0_map(nsx, nzm), &
                t0_map(nsx, nzm), q0_map(nsx, nzm),  &
                tabs0_map(nsx, nzm), qv0_map(nsx, nzm),  &
                qn0_map(nsx, nzm), qp0_map(nsx, nzm))
                
        allocate(u_out_map(nsx, nzm),  &
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
    
    

    if (masterproc) then

        


        !------------- 调用 host model -------------
        call host_model_evolve( u0_in=u0_map, wsub_in=wsub_map, &
                            t0_in=t0_map, q0_in=q0_map,                    &
                            tabs0_in = tabs0_map, qv0_in = qv0_map, qn0_in = qn0_map, qp0_in = qp0_map,   &
                            u_out_map=u_out_map,           &
                            w_out_map=w_out_map, t_out_map=t_out_map, q_out_map=q_out_map )
        
        wsub_map(:, :)         = w_out_map(:, :)

        

    end if  ! masterproc

    ! distribute ug0_hm
    if (masterproc) then
        call task_bscatter_float_map(0, u_out_map, nzm, nsx, ug0_hm(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, ug0_hm(1))
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

    
    !------------------------------------------------------------
    ! clean up
    !------------------------------------------------------------
    if (masterproc) then
        if (allocated(u0_map))    deallocate(u0_map)
        if (allocated(t0_map))    deallocate(t0_map)
        if (allocated(q0_map))    deallocate(q0_map)
        if (allocated(tabs0_map))    deallocate(tabs0_map)
        if (allocated(qv0_map))    deallocate(qv0_map)
        if (allocated(qn0_map))    deallocate(qn0_map)
        if (allocated(qp0_map))    deallocate(qp0_map)
        if (allocated(u_out_map))  deallocate(u_out_map)  
        if (allocated(t_out_map))  deallocate(t_out_map)  
        if (allocated(q_out_map))  deallocate(q_out_map)  
        if (allocated(w_out_map))  deallocate(w_out_map)  
    else
        if (allocated(dummy2d))  deallocate(dummy2d)
    end if

    call t_stopf('host_model')

end subroutine hm_couple_step
