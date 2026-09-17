	
subroutine setparm
	
!       initialize parameters:

use vars
!use micro_params
use params
use microphysics, only: micro_setparm
use sgs, only: sgs_setparm
use movies, only : irecc
use instrument_diagnostics, only: zero_instr_diag
implicit none
	
integer icondavg, ierr, ios, ios_missing_namelist, place_holder

NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
                dolongwave, doshortwave, dosgs, dz, doconstdz, &
                docoriolis, docoriolisz, dosurface, dolargescale, doradforcing, &
		            fluxt0,fluxq0,tau0,tabs_s,z0,nelapse, nelapsemin, dt, dx, dy,  &
                fcor, ug, vg, nstop, caseid, case_restart,caseid_restart, &
		            nstat, nstatfrq, nprint, nrestart, doradsimple, &
		            nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
		            donudging_uv, donudging_tq, &
                donudging_t, donudging_q, tauls,tautqls,&
                nudging_uv_z1, nudging_uv_z2, nudging_t_z1, nudging_t_z2, &
                nudging_q_z1, nudging_q_z2, dofplane, &
		            timelargescale, longitude0, latitude0, day0, nrad, &
		            OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, soil_wetness, &
                doensemble, nensemble, dowallx, dowally, &
                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, & 
                docolumn, save2Dbin, save2Davg, save3Dbin, &
                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
	              doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
                domodis, domisr, dodynamicocean, ocean_type, delta_sst, &
                depth_slab_ocean, Szero, deltaS, timesimpleocean, &
                dosolarconstant, solar_constant, zenith_angle, rundatadir, &
                dotracers, output_sep, perturb_type, &
                doSAMconditionals, dosatupdnconditionals, &
                doscamiopdata, iopfile, dozero_out_day0, &
                nstatmom, nstatmomstart, nstatmomend, savemomsep, savemombin, &
                nmovie, nmoviestart, nmovieend, nrestart_skip, &
                bubble_x0,bubble_y0,bubble_z0,bubble_radius_hor, &
                bubble_radius_ver,bubble_dtemp,bubble_dq, dosmoke, dossthomo, &
                rad3Dout, nxco2, dosimfilesout, notracegases, &
                doradlat, doradlon, ncycle_max, doseawater, SLM, LES_S

! Parameters added by Kuang Lab at Harvard
NAMELIST /KUANG_PARAMS/ dompiensemble, &
                dolayerperturb, tperturbi, qperturbi, tperturbA, qperturbA, & ! linear response perturbation: layer by layer (Song Qiyu, 2022)
                nstartperturb, nperturbstep, dorandmultisine, &
                dorandmultisineoddonly, douniformintime, dowhitenoiseforcing, &
                delt_perturbt, delt_perturbq, &
                increase_delt, nstep_increase_start, nstep_increase_end, &
                delt_perturbt_end, delt_perturbq_end, &
                iensemble, nT, nQ, maxperturbperiod, & ! perturbation signal: randmultisine (Qiyu, 2022)
                doidealizedrad, dobulksfc, icopy, &
                nstep_separate_statfile, &
                wavenumber_factor,nstartlinearwave,nsteplinearwavebg,nsteplinearwave,&
                wavedampingtime,wavetqdampingtime,&
                doadvectbg,doparameterizedwave,dointernalnoise,&
                dosavemultirestart, nrestartstart,&
                donoisywave, noiselevel, &
                dompimmf, nstephostmodel, hm_spinup_step, &
                nouvchatting, but_nudge_u, hm_only, diffuse_intensity,do_3step_adams,hm_subcycle, &
                hm_smoother, hyper_intensity, smag_cs, smag_max_diff_vel, smag_nu_max_frac, &
                CRM_damping0, CRM_dampingRM, apply_hm_u_external_nudging, large_u_profile_filename, tauls_large_scale, &
                diffuse_intensity_subdomain_large_scale, subdomain_center_at_hm_u_center, &
                suppress_k_start, &
                tau_damp_mean, do_damp_hm_mean, &
                do_remove_nyquist_u, do_remove_coupling_residual, &
                do_fix_u_halo, &
                do_hm_bubble, hm_bubble_step, hm_bubble_z_bot, hm_bubble_z_top, hm_bubble_nsubdomain_half, hm_bubble_dtemp, &
                hm_cfl_max, &
                add_initial_bubble, init_bubble_z_top, init_bubble_nsubdomain_half, init_bubble_dtemp


                

!bloss: Create dummy namelist, so that we can figure out error code
!       for a mising namelist.  This lets us differentiate between
!       missing namelists and those with an error within the namelist.
NAMELIST /BNCUIODSBJCB/ place_holder

!----------------------------------
!  Read namelist variables from the standard input:
!------------

open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
read (55,PARAMETERS,IOSTAT=ierr)
if (ierr.ne.0) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in PARAMETERS namelist'
        call task_abort()
end if
close(55)

!----------------------------------
!  Read namelist for Kuang_Lab options from same prm file:
!------------
open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

!bloss: get error code for missing namelist (by giving the name for
!       a namelist that doesn't exist in the prm file).
read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
rewind(55) !note that one must rewind before searching for new namelists

!bloss: read in UWOPTIONS namelist
read (UNIT=55,NML=KUANG_PARAMS,IOSTAT=ios)
if (ios.ne.0) then
  if(masterproc) write(*,*) 'ios_missing_namelist = ', ios_missing_namelist
  if(masterproc) write(*,*) 'ios for KUANG_PARAMS = ', ios
   !namelist error checking
   if(ios.ne.ios_missing_namelist) then
     rewind(55) !note that one must rewind before searching for new namelists
     read (UNIT=55,NML=KUANG_PARAMS)
     if(masterproc) then
       write(*,*) '****** ERROR: bad specification in KUANG_PARAMS namelist'
     end if
      call task_abort()
   elseif(masterproc) then
      write(*,*) '****************************************************'
      write(*,*) '******* No KUANG_PARAMS namelist in prm file *******'
      write(*,*) '****************************************************'
   end if
end if
close(55)

! write namelist values out to file for documentation
if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.nml',&
            form='formatted')
      write (55,nml=PARAMETERS)
      write (55,nml=KUANG_PARAMS)
      write(55,*) 
      close(55)
end if

!------------------------------------
!  Set parameters 


        ! Allow only special cases for separate output:

        output_sep = output_sep.and.RUN3D
        if(output_sep)  save2Dsep = .true.

	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif

        if(docoriolis.and..not.dofplane.or.doradlat) dowally=.true.

	if(ny.eq.1) dy=dx
	dtn = dt

	notopened2D = .true.
	notopened3D = .true.

        call zero_instr_diag() ! initialize instruments output 
        call sgs_setparm() ! read in SGS options from prm file.
        call micro_setparm() ! read in microphysical options from prm file.

        if(dosmoke) then
           epsv=0.
        else    
           epsv=0.61
        endif   

        if(navgmom_x.lt.0.or.navgmom_y.lt.0) then  
            nstatmom        = 1
            nstatmomstart    = 999999999
            nstatmomend      = 999999999
        end if

        if(doseawater) then
          salt_factor = 0.981
        else
          salt_factor = 1.
        end if

        if(tautqls.eq.99999999.) tautqls = tauls
        
        dtfactor = 1.

        !===============================================================
        ! KUANG_LAB ADDITION

        if(dompiensemble.AND.dompi) then
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using the Kuang_Lab Ensemble Run Method'
            write(*,*) '  This will turn off MPI in the model run, such that'
            write(*,*) '  each subdomain is run independently of each other.'
            write(*,*) '  However, MPI is turned on for saving of output and'
            write(*,*) '  restart files.'
            write(*,*) '*********************************************************'
          end if
        else if(dompiensemble.AND.(.NOT.dompi)) then
          dompiensemble = .false.
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Do not use the Kuang_Lab Ensemble Run Method'
            write(*,*) '  MPI is not called because number of processors = 1'
            write(*,*) '  Setting dompiensemble to FALSE'
            write(*,*) '*********************************************************'
          end if
        end if
          
        if(dolayerperturb) then
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using the Kuang_Lab layer-by-layer perturbation'
            write(*,*) '  This is meant to calculate linear response function.'
            if(tperturbi.gt.0.and.tperturbA.ne.0.) then
              write(*,*) '  Add temperature perturbation to layer ', tperturbi
            end if
            if(qperturbi.gt.0.and.qperturbA.ne.0.) then
              write(*,*) '  Add water vapor perturbation to layer ', qperturbi
            end if
            write(*,*) '*********************************************************'
          end if
        end if
        
        if(dorandmultisine) then
          if((nT.le.0).or.(nQ.le.0).or.(nT.gt.nzm).or.(nQ.gt.nzm)) then
            print*,'Check nT and nQ in randmultisine!'
            call task_abort()
          end if
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using Random Phase Multi-Sine perturbations as input'
            write(*,*) '  signal of the system. These perturbations are added to'
            write(*,*) '  field every nstat steps, while statistic field is'
            write(*,*) '  written in output at the same steps.'
            if(dorandmultisineoddonly) then
              write(*,*) '  - Only odd frequency perturbations are included.'
            else
              write(*,*) '  - All frequency perturbations are included.'
            end if
            if(increase_delt) then
              write(*,*) '  - delt_t/q are: ', delt_perturbt, delt_perturbq
              write(*,*) '  - before step number: ', nstep_increase_start
              write(*,*) '  - then linearly increase to:', delt_perturbt_end, delt_perturbq_end
              write(*,*) '  - by step number: ', nstep_increase_end
              write(*,*) '  - then stay constant.'
            else
              write(*,*) '  - delt_t/q are: ', delt_perturbt, delt_perturbq
            end if
            write(*,*) '*********************************************************'
          end if
        end if

        if(doidealizedrad) then
          if(masterproc) then
            write(*,*) '  Doing idealized radiation.'
            if(dolongwave.or.doshortwave.or.doradforcing) then
              dolongwave = .false.
              doshortwave = .false.
              doradforcing = .false.
              write(*,*) '  Setting dolongwave, doshortwave, doradforcing as false.'
            end if
          end if
        end if

        if(dompimmf) then
          dx_hm = dx_hm_km * 1000.     ! host grid spacing, set in domain.f90
          dt_hm = dt * nstephostmodel
          dt_hm_subcycle = dt_hm / hm_subcycle

          ! ---- coupling filter: validate the suppressed wavenumber range ----
          if(suppress_k_start.lt.0) suppress_k_start = nsx/2   ! default: Nyquist only
          if(suppress_k_start.lt.1 .or. suppress_k_start.gt.nsx/2) then
            if(masterproc) then
              write(*,*) '*********************************************************'
              write(*,*) '  ERROR: suppress_k_start = ', suppress_k_start
              write(*,*) '  must satisfy  1 <= suppress_k_start <= nsx/2 = ', nsx/2
              write(*,*) '  (nsx/2 is the Nyquist wavenumber of the host grid)'
              write(*,*) '*********************************************************'
            end if
            call task_abort()
          end if
          if(hm_smoother.lt.0 .or. hm_smoother.gt.3) then
            if(masterproc) then
              write(*,*) '*********************************************************'
              write(*,*) '  ERROR: hm_smoother = ', hm_smoother
              write(*,*) '  must be 0 (none), 1 (grad^2), 2 (grad^4) or 3 (Smagorinsky).'
              write(*,*) '  See the parameter block in vars.f90.'
              write(*,*) '*********************************************************'
            end if
            call task_abort()
          end if
          if(smag_nu_max_frac.le.0.) then
            if(masterproc) then
              write(*,*) '*********************************************************'
              write(*,*) '  ERROR: smag_nu_max_frac = ', smag_nu_max_frac
              write(*,*) '  must be > 0; it is the stability clamp on the Smagorinsky'
              write(*,*) '  viscosity, nu <= frac*dx_hm^2/dt_hm_subcycle (frac < 0.136).'
              write(*,*) '*********************************************************'
            end if
            call task_abort()
          end if
          if(tau_damp_mean.le.0.) then
            if(masterproc) then
              write(*,*) '*********************************************************'
              write(*,*) '  ERROR: tau_damp_mean = ', tau_damp_mean, ' s'
              write(*,*) '  must be > 0 (seconds). Use do_damp_hm_mean = .false.'
              write(*,*) '  to switch the domain-mean drag off instead.'
              write(*,*) '*********************************************************'
            end if
            call task_abort()
          end if
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using the Kuang_Lab Multi-scale Modeling Framework'
            write(*,*) '  Coupling with a host model.'
            write(*,*) '  Currently only works with 2D Walker circulation.'
            write(*,*) '  The host model is run every ', nstephostmodel, ' steps.'
            write(*,*) '  Coupling filter suppresses wavenumbers ', suppress_k_start, &
                       ' to ', nsx/2, ' (Nyquist)'
            write(*,*) '    i.e. host-scale wavelengths at or below ', &
                       nsx*dx_hm/float(suppress_k_start)/1000., ' km'
            write(*,*) '  do_fix_u_halo    = ', do_fix_u_halo, &
                       ' (refresh the subdomain u halo right after the coupling increment)'
            write(*,*) '  ----- grid geometry (from domain.f90 + prm) -----'
            write(*,*) '  dx_hm            = ', dx_hm/1000., ' km'
            write(*,*) '  nsx (host cols)  = ', nsx
            write(*,*) '  host domain      = ', nsx*dx_hm/1000., ' km'
            write(*,*) '  subdomain width  = ', nx*dx/1000., ' km'
            write(*,*) '  CRM extent       = ', nx_gl*dx/1000., ' km'
            write(*,*) '  host Nyquist     = ', 2.*dx_hm/1000., ' km'
            if(hm_cfl_max.gt.0.) then
              write(*,*) '  host CFL guard   = ', hm_cfl_max, ' (abort above; vertical term usually binds)'
            else
              write(*,*) '  host CFL guard   = OFF  <- hm_cfl_max <= 0., diagnostic only'
            end if
            ! ----- audit of every term acting on the mean wind -----------------
            ! Four separate terms can act, on two different fields, under three
            ! gates plus a spin-up branch.  Print them all so the configuration is
            ! readable from the log instead of having to be reconstructed from the
            ! prm.  This is what would have caught the fighting-targets bug that
            ! fef8700 fixed.  See CHANGES_since_tend-nudging2.md.
            write(*,*) '  ----- mean-wind forcing: every active term -----'
            if(apply_hm_u_external_nudging) then
              write(*,*) '   HOST  external-profile nudging -> u_external_profile,  tau =', &
                         tauls_large_scale, ' s   [ON]'
              write(*,*) '         profile from: ', trim(large_u_profile_filename)
            else
              write(*,*) '   HOST  external-profile nudging                          [off]'
            end if
            if(do_damp_hm_mean) then
              if(apply_hm_u_external_nudging) then
                write(*,*) '   HOST  domain-mean drag         -> u_external_profile,  tau =', &
                           tau_damp_mean, ' s   [ON]'
                write(*,*) '         *** same target as the nudging above and', &
                           tau_damp_mean/tauls_large_scale, 'x weaker:'
                write(*,*) '         *** REDUNDANT, contributes', &
                           100./(1.+tau_damp_mean/tauls_large_scale), '% of the relaxation'
              else
                write(*,*) '   HOST  domain-mean drag         -> zero,                tau =', &
                           tau_damp_mean, ' s   [ON]'
                write(*,*) '         this is the ONLY term bounding <u>; the host has no surface stress'
              end if
            else
              write(*,*) '   HOST  domain-mean drag                                  [off]'
              if(.not.apply_hm_u_external_nudging) &
                write(*,*) '         *** WARNING: nothing bounds the host domain-mean wind'
            end if
            if(dompimmf .and. .not.hm_only) then
              if(donudging_uv) then
                write(*,*) '   CRM   coupling increment ug0_hm via nudging_hm          [ON]'
              else
                write(*,*) '   CRM   coupling increment ug0_hm                         [off]'
                write(*,*) '         *** WARNING: donudging_uv = F gates the block in', &
                           ' nudging_hm, so the host increment never reaches the CRMs'
              end if
            end if
            if(hm_only) then
              write(*,*) '   CRM   nudging() -> ug0,                                tau =', &
                         tauls, ' s   [ON, hm_only takes this branch]'
              if(apply_hm_u_external_nudging) then
                write(*,*) '         ug0 is overwritten from the external profile each step'
                write(*,*) '         (set_ug0_from_external_profile, called from forcing)'
              else
                write(*,*) '         ug0 comes from snd -- whose u column is ZERO in our cases'
              end if
            end if
            write(*,*) '   CRM   nudging() -> ug0,  ALSO for the first', nstephostmodel, &
                       ' steps of every run'
            write(*,*) '         (could_hm_nudging is .false. until the first hm_couple_step)'
            write(*,*) '  ----- horizontal smoother (see vars.f90) -----'
            select case (hm_smoother)
            case (0)
              write(*,*) '  hm_smoother = 0 : none'
            case (1)
              write(*,*) '  hm_smoother = 1 : grad^2, diffuse_intensity = ', diffuse_intensity
              write(*,*) '    nu            = ', diffuse_intensity*dx_hm*dx_hm/dt_hm_subcycle, ' m2/s'
              write(*,*) '    e-fold at 2dx = ', dt_hm_subcycle/max(4.*diffuse_intensity,1.e-30)/3600., ' h'
              write(*,*) '    e-fold at 1000 km = ', &
                   1./max(diffuse_intensity*dx_hm*dx_hm/dt_hm_subcycle,1.e-30) &
                   /(2.*3.14159265/1.e6)**2/3600., ' h'
            case (2)
              write(*,*) '  hm_smoother = 2 : grad^4, hyper_intensity  = ', hyper_intensity
              write(*,*) '    nu4           = ', hyper_intensity*dx_hm**4/dt_hm_subcycle, ' m4/s'
              write(*,*) '    e-fold at 2dx = ', dt_hm_subcycle/max(16.*hyper_intensity,1.e-30)/3600., ' h'
              write(*,*) '    e-fold at 1000 km = ', &
                   1./max(hyper_intensity*dx_hm**4/dt_hm_subcycle,1.e-30) &
                   /(2.*3.14159265/1.e6)**4/3600., ' h'
            case (3)
              write(*,*) '  hm_smoother = 3 : Smagorinsky, smag_cs = ', smag_cs
              write(*,*) '    equivalent C (=pi*Cs, MITgcm/Griffies-Hallberg) = ', 3.14159265*smag_cs
              write(*,*) '    nu = (smag_cs*dx_hm)^2*|du/dx|, prefactor = ', (smag_cs*dx_hm)**2, ' m2'
              write(*,*) '    cap: min(', smag_max_diff_vel*dx_hm, ',', &
                   smag_nu_max_frac*dx_hm*dx_hm/dt_hm_subcycle, ') m2/s'
            end select
            ! AB3 real-axis stability limit is about 0.545
            if(hm_smoother.eq.1 .and. 4.*diffuse_intensity.gt.0.545) then
              write(*,*) '  *** WARNING: 4*diffuse_intensity = ', 4.*diffuse_intensity, &
                         ' exceeds the AB3 real-axis limit 0.545'
            end if
            if(hm_smoother.eq.2 .and. 16.*hyper_intensity.gt.0.545) then
              write(*,*) '  *** WARNING: 16*hyper_intensity = ', 16.*hyper_intensity, &
                         ' exceeds the AB3 real-axis limit 0.545'
            end if
            if(hm_smoother.eq.3 .and. 4.*smag_nu_max_frac.gt.0.545) then
              write(*,*) '  *** WARNING: 4*smag_nu_max_frac = ', 4.*smag_nu_max_frac, &
                         ' exceeds the AB3 real-axis limit 0.545'
            end if
            if(hm_smoother.eq.2 .and. hyper_intensity.le.0.) &
              write(*,*) '  *** NOTE: hm_smoother=2 but hyper_intensity <= 0, so no smoothing'
            if(hm_smoother.eq.3 .and. smag_cs.le.0.) &
              write(*,*) '  *** NOTE: hm_smoother=3 but smag_cs <= 0, so no smoothing'
            if(hm_smoother.ne.1 .and. diffuse_intensity.ne.0.) &
              write(*,*) '  *** NOTE: diffuse_intensity = ', diffuse_intensity, &
                         ' is IGNORED because hm_smoother /= 1'
            write(*,*) '*********************************************************'
          end if
        end if

        !===============================================================
        ! UW ADDITION

        !bloss: set up conditional averages
        ncondavg = 1 ! always output CLD conditional average
        if(doSAMconditionals) ncondavg = ncondavg + 2
        if(dosatupdnconditionals) ncondavg = ncondavg + 3
        if(allocated(condavg_factor)) then ! avoid double allocation when nrestart=2
          DEALLOCATE(condavg_factor,condavg_mask,condavgname,condavglongname)
        end if
        ALLOCATE(condavg_factor(nzm,ncondavg), & ! replaces old cloud_factor, core_factor
             condavg_mask(nx,ny,nzm,ncondavg), & ! nx x ny x nzm indicator arrays
             condavgname(ncondavg), & ! short names (e.g. CLD, COR, SATUP)
             condavglongname(ncondavg), & ! long names (e.g. cloud, core, saturated updraft)
             STAT=ierr)
        if(ierr.ne.0) then
             write(*,*) '**************************************************************************'
             write(*,*) 'ERROR: Could not allocate arrays for conditional statistics in setparm.f90'
             call task_abort()
        end if
        
        ! indicators that can be used to tell whether a particular average
        !   is present.  If >0, these give the index into the condavg arrays
        !   where this particular conditional average appears.
        icondavg_cld = -1
        icondavg_cor = -1
        icondavg_cordn = -1
        icondavg_satup = -1
        icondavg_satdn= -1
        icondavg_env = -1

        icondavg = 0
        icondavg = icondavg + 1
        condavgname(icondavg) = 'CLD'
        condavglongname(icondavg) = 'cloud'
        icondavg_cld = icondavg

        if(doSAMconditionals) then
           icondavg = icondavg + 1
           condavgname(icondavg) = 'COR'
           condavglongname(icondavg) = 'core'
           icondavg_cor = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'CDN'
           condavglongname(icondavg) = 'downdraft core'
           icondavg_cordn = icondavg
        end if
           
        if(dosatupdnconditionals) then
           icondavg = icondavg + 1
           condavgname(icondavg) = 'SUP'
           condavglongname(icondavg) = 'saturated updrafts'
           icondavg_satup = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'SDN'
           condavglongname(icondavg) = 'saturated downdrafts'
           icondavg_satdn = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'ENV'
           condavglongname(icondavg) = 'unsaturated environment'
           icondavg_env = icondavg
        end if
           
        ! END UW ADDITIONS
        !===============================================================

        irecc = 1


end
