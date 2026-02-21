subroutine get_domain_avg_U()
    use vars
    use params

    implicit none

    integer :: k

   
    real :: u0_local(nzm), v0_local(nzm), w0_local(nz)
    real, allocatable :: u0_map(:,:),  v0_map(:,:), w0_map(:,:)
    real, allocatable :: u0_mean_map(:,:),  v0_mean_map(:,:), w0_mean_map(:,:)
    real, allocatable :: dummy2d(:,:), dummy2dw(:,:)
   
    
    if (masterproc) then
        allocate(u0_map(nsx, nzm), v0_map(nsx, nzm), w0_map(nsx, nz),&
                   u0_mean_map(nsx, nzm), v0_mean_map(nsx, nzm), w0_mean_map(nsx, nz) )
    else
        allocate(dummy2d(nsx, nzm), dummy2dw(nsx, nz))  ! assign dummy to avoid alloc error
    end if


    do k = 1, nz
        w0_local(k) = sum( w(1:nx,1:ny,k) ) / real(nx*ny)
    end do


    u0_local = u0
    v0_local = v0
    

    
    ! gather u0
    if (masterproc) then
        call task_bgather_float_map(0, u0_local(1), nzm, nsx, u0_map)
    else
        call task_bgather_float_map(0, u0_local(1), nzm, nsx, dummy2d)
    end if

    ! gather v0
    if (masterproc) then
        call task_bgather_float_map(0, v0_local(1), nzm, nsx, v0_map)
    else
        call task_bgather_float_map(0, v0_local(1), nzm, nsx, dummy2d)
    end if
   
    ! gather w0
    if (masterproc) then
        call task_bgather_float_map(0, w0_local(1), nz, nsx, w0_map)
    else
        call task_bgather_float_map(0, w0_local(1), nz, nsx, dummy2dw)
    end if
    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    if(masterproc) then
        do k = 1, nzm
            u0_mean_map(:,k) = sum( u0_map(1:nsx,k) ) / real(nsx)
            v0_mean_map(:,k) = sum( v0_map(1:nsx,k) ) / real(nsx)
        end do

        do k = 1, nz
            w0_mean_map(:,k) = sum( w0_map(1:nsx,k) ) / real(nsx)
        end do

    end if

    !-----------------------------------------------------------------------------------------------------------------------------------------------------
    ! distribute u_domain_avg
    if (masterproc) then
        call task_bscatter_float_map(0, u0_mean_map, nzm, nsx, u_domain_avg(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, u_domain_avg(1))
    end if

    ! distribute v_domain_avg
    if (masterproc) then
        call task_bscatter_float_map(0, v0_mean_map, nzm, nsx, v_domain_avg(1))
    else
        call task_bscatter_float_map(0, dummy2d,  nzm, nsx, v_domain_avg(1))
    end if

    ! distribute w_domain_avg
    if (masterproc) then
        call task_bscatter_float_map(0, w0_mean_map, nz, nsx, w_domain_avg(1))
    else
        call task_bscatter_float_map(0, dummy2dw,  nz, nsx, w_domain_avg(1))
    end if
 

    !------------------------------------------------------------
    ! clean up   
    !------------------------------------------------------------
    if (masterproc) then
        if (allocated(u0_map))    deallocate(u0_map)
        if (allocated(v0_map))    deallocate(v0_map)
        if (allocated(w0_map))    deallocate(w0_map)
        if (allocated(u0_mean_map))    deallocate(u0_mean_map)
        if (allocated(v0_mean_map))    deallocate(v0_mean_map)
        if (allocated(w0_mean_map))    deallocate(w0_mean_map)
    else
        if (allocated(dummy2d))  deallocate(dummy2d)
        if (allocated(dummy2dw))  deallocate(dummy2dw)
    end if



end subroutine get_domain_avg_U
