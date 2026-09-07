! Set the domain dimensionality, size and number of subdomains.

module domain

       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
       integer, parameter :: nx_gl = 2560 ! Number of grid points in X     # default: 2560
       integer, parameter :: ny_gl = 32 ! Number of grid points in Y
       integer, parameter :: nz_gl = 64 ! Number of pressure (scalar) levels
       integer, parameter :: nsubdomains_x  = 160 ! No of subdomains in x   # default: 80
       integer, parameter :: nsubdomains_y  = 1 ! No of subdomains in y

!--------------------------------------------------------------------
! MMF host-model grid (used only when dompimmf = .true.)
!
! One host column per subdomain (1:1). The host column need NOT be as wide
! as the CRM subdomain that samples it -- SP-CAM style, the CRM is a
! representative sample of the column, not a tile that fills it.
!
!   subdomain width = (nx_gl/nsubdomains_x) * dx      [dx set in prm]
!   host domain     = nsubdomains_x * dx_hm_km
!   CRM extent      = nx_gl * dx
!   host Nyquist    = 2 * dx_hm_km   <- the one scale the coupling cannot
!                                       transfer (null space of cell-averaging)
!
! Set dx_hm_km = (nx_gl/nsubdomains_x)*dx/1000 to make the host domain equal
! the CRM extent, which is the original 1-tile-per-column configuration.
!
!   dx_hm_km   nsubdomains_x   nx_gl    host domain   host Nyquist
!      32           320         5120      10240 km        64 km
!      64           160         2560      10240 km       128 km
!     128            80         1280      10240 km       256 km
!--------------------------------------------------------------------
       real, parameter :: dx_hm_km = 64.   ! host grid spacing [km]

       ! define # of points in x and y direction to average for 
       !   output relating to statistical moments.
       ! For example, navgmom_x = 8 means the output will be   
       !  8 times coarser grid than the original.
       ! If don't wanna such output, just set them to -1 in both directions. 
       ! See Changes_log/README.UUmods for more details.
       integer, parameter :: navgmom_x = -1 
       integer, parameter :: navgmom_y = -1 

       integer, parameter :: ntracers = 0 ! number of transported tracers (dotracers=.true.)
       
! Note:
!  * nx_gl and ny_gl should be a factor of 2,3, or 5 (see User's Guide)
!  * if 2D case, ny_gl = nsubdomains_y = 1 ;
!  * nsubdomains_x*nsubdomains_y = total number of processors
!  * if one processor is used, than  nsubdomains_x = nsubdomains_y = 1;
!  * if ntracers is > 0, don't forget to set dotracers to .true. in namelist 

end module domain
