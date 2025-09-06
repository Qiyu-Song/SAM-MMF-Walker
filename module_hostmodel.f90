module module_hostmodel

implicit none

private

public:: host_model_evolve

contains

subroutine host_model_evolve()
! TODO

! advect_mom_hm(u0_map, v0_map, wsub_hm_map) -> dudt_hm_map, dvdt_hm_map, dwdt_hm_map

! pressure_hm(u0_map, v0_map, w0_map, dudt_hm_map, dvdt_hm_map, dwdt_hm_map) -> dudt_hm_map, dvdt_hm_map, dwdt_hm_map

    ! update (dudt_hm_map_series, dvdt_hm_map_series, dwdt_hm_map_series) with (dudt_hm_map, dvdt_hm_map, dwdt_hm_map)

! adams_hm(dudt_hm_map_series, dvdt_hm_map_series, dwdt_hm_map_series) ->  u_hm_map, v_hm_map, wsub_hm_map     !! series存放前几个时间步的信息

! advect_scalars_hm(u_hm_map, v_hm_map, wsub_hm_map, t0, q0) -> t_hm_map, q_hm_map



end subroutine host_model_evolve

end module module_hostmodel
! external_forcing(u_hm_map, u0_map) -> ex_forcing_hm_map
! external_forcing(v_hm_map, v0_map) -> ex_forcing_hm_map
! external_forcing(t_hm_map, t0_map) -> ex_forcing_hm_map
! external_forcing(q_hm_map, q0_map) -> ex_forcing_hm_map