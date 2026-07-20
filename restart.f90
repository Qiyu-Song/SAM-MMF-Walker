	subroutine write_all()
	
  use vars
  use params, only: dompimmf
	implicit none
	character *4 rankchar
        character *8 nrestartfile
	character *256 filename
	integer irank
        integer, external :: lenstr

        call t_startf ('restart_out')
        
        ! number of restart file (Qiyu, 2024)
        if(dosavemultirestart) then
          write(nrestartfile,'(i8)') nstep/(nstat*(1+nrestart_skip))
        end if

        if(masterproc) then
         print*,'Writing restart file ...'
         if(dosavemultirestart) then
           filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//trim(adjustl(nrestartfile))//'_misc_restart.bin'
         else
           filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
         end if
         open(66,file=trim(filename), status='unknown',form='unformatted')
        end if


	if(restart_sep) then

          write(rankchar,'(i4)') rank
          
          if(dosavemultirestart) then
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_'//trim(adjustl(nrestartfile))//'_restart.bin'
          else
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          end if

          open(65,file=trim(filename), status='unknown',form='unformatted')
          write(65) nsubdomains, nsubdomains_x, nsubdomains_y

	  call write_statement


	else
	  write(rankchar,'(i4)') nsubdomains
	  !filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
          !      rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          if(dosavemultirestart) then
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_'//trim(adjustl(nrestartfile))//'_restart.bin'
          else
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          end if

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       if(masterproc) then
	      
	        open(65,file=trim(filename), status='unknown',form='unformatted')
	        write(65) nsubdomains, nsubdomains_x, nsubdomains_y

	       else

                open(65,file=trim(filename), status='unknown',form='unformatted',&
                   position='append')

	       end if

               call write_statement

             end if
	  end do

	end if ! restart_sep

        call task_barrier()

        if(dompimmf) call write_hm_restart()

        call t_stopf ('restart_out')

        return
        end
 
 
 
 
     
	subroutine read_all()
	
	use vars
        use params, only: dompiensemble, dompimmf
  implicit none
	character *4 rankchar
	character *256 filename
	integer irank, ii
        integer, external :: lenstr

        if(masterproc) print*,'Reading restart file ...'

        if(nrestart.ne.2) then
          filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
        else
          filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_misc_restart.bin'
        end if
        open(66,file=trim(filename), status='unknown',form='unformatted')


	if(restart_sep) then

           write(rankchar,'(i4)') rank

           if(nrestart.ne.2) then
             filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
           else
             filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
           end if


           open(65,file=trim(filename), status='unknown',form='unformatted')
           read(65)

	   call read_statement


	else

	  write(rankchar,'(i4)') nsubdomains

          if(nrestart.ne.2) then
	    filename='./RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          else
	    filename='./RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          end if
          open(65,file=trim(filename), status='unknown',form='unformatted')

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       read (65)
 
               do ii=0,irank-1 ! skip records
                 read(65)
	       end do

	       call read_statement

             end if

	  end do

	end if ! restart_sep

	call task_barrier()

        if(dompimmf) call read_hm_restart()

        dtfactor = -1.

! update the boundaries 
! (just in case when some parameterization initializes and needs boundary points)
        
        ! Kuang Ensemble run: turn off mpi for boundaries (Song Qiyu, 2022)
        if(dompiensemble.or.dompimmf) dompi = .false.

        call boundaries(1)
        call boundaries(4)

        ! Kuang Ensemble run: turn on mpi after boundaries (Song Qiyu, 2022)
        if(dompiensemble.or.dompimmf) dompi = .true.

        return
        end



        subroutine write_hm_restart()

        use vars
        use params, only: dompimmf
        implicit none
        character *4 rankchar
        character *8 nrestartfile
        character *256 filename
        integer irank
        integer, external :: lenstr

        if(.not.dompimmf) return

        if(dosavemultirestart) then
          write(nrestartfile,'(i8)') nstep/(nstat*(1+nrestart_skip))
        end if

        if(masterproc) print*,'Writing host-model restart file ...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank
          if(dosavemultirestart) then
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_'//trim(adjustl(nrestartfile))//'_restart_hm.bin'
          else
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_hm.bin'
          end if

          open(67,file=trim(filename), status='unknown',form='unformatted')
          write(67) nsubdomains, nsubdomains_x, nsubdomains_y
          call write_hm_statement

        else

          write(rankchar,'(i4)') nsubdomains
          if(dosavemultirestart) then
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_'//trim(adjustl(nrestartfile))//'_restart_hm.bin'
          else
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_hm.bin'
          end if

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then

                open(67,file=trim(filename), status='unknown',form='unformatted')
                write(67) nsubdomains, nsubdomains_x, nsubdomains_y

               else

                open(67,file=trim(filename), status='unknown',form='unformatted',&
                   position='append')

               end if

               call write_hm_statement

             end if
          end do

        end if

        call task_barrier()

        return
        end



        subroutine read_hm_restart()

        use vars
        use params, only: dompimmf
        implicit none
        character *4 rankchar
        character *256 filename
        integer irank, ii, ios
        integer, external :: lenstr

        if(.not.dompimmf) return

        if(masterproc) print*,'Reading host-model restart file ...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          if(nrestart.ne.2) then
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_hm.bin'
          else
            filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_hm.bin'
          end if

          open(67,file=trim(filename), status='old',form='unformatted',iostat=ios)
          if(ios.ne.0) call missing_hm_restart(filename)
          read(67)

          call read_hm_statement

        else

          write(rankchar,'(i4)') nsubdomains

          if(nrestart.ne.2) then
            filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_hm.bin'
          else
            filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart_hm.bin'
          end if

          open(67,file=trim(filename), status='old',form='unformatted',iostat=ios)
          if(ios.ne.0) call missing_hm_restart(filename)

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               read (67)

               do ii=0,irank-1
                 call skip_hm_statement
               end do

               call read_hm_statement

             end if
          end do

        end if

        call task_barrier()

        return
        end



        subroutine skip_hm_statement()

        implicit none
        integer rank1, nx1, ny1, nz1, nsx1, hm_step1
        logical wsub_inited1, wsub_alloc1, could_hm_nudging1

        read(67) rank1, nx1, ny1, nz1, nsx1, wsub_inited1, wsub_alloc1, &
             could_hm_nudging1, hm_step1
        if(wsub_alloc1) read(67)
        read(67)

        return
        end



        subroutine write_hm_statement()

        use vars
        implicit none
        logical wsub_alloc

        wsub_alloc = allocated(wsub_map)

        write(67) rank, nx, ny, nz, nsx, wsub_inited, wsub_alloc, &
             could_hm_nudging, hm_step
        if(wsub_alloc) write(67) wsub_map
        write(67) &
             u0_local_hm, t0_local_hm, q0_local_hm, tabs0_local_hm, &
             qn0_local_hm, qp0_local_hm, qni0_local_hm, qnl0_local_hm, &
             qpi0_local_hm, qpl0_local_hm, prec_flx_local_hm, &
             ug0_hm, tg0_hm, qg0_hm, ug0_press_modify, &
             u_hm_map_save, u_hm_updated_map_save, u_sub_map_save, &
             t_hm_map_save, t_hm_updated_map_save, t_sub_map_save, &
             q_hm_map_save, q_hm_updated_map_save, q_sub_map_save, &
             dudt_hm_hist, dwdt_hm_hist, dudt_subdomain_diffuse, &
             u_external_profile, u_external_profile_loaded, &
             prec_xy_crm, prec_xy_save, shf_xy_crm, shf_xy_save, &
             lhf_xy_crm, lhf_xy_save
        close(67)

        if(rank.eq.nsubdomains-1) then
          print *,'Host-model restart file was saved. nstep=',nstep
        endif

        return
        end



        subroutine read_hm_statement()

        use vars
        implicit none
        integer rank1, nx1, ny1, nz1, nsx1
        logical wsub_alloc

        read(67) rank1, nx1, ny1, nz1, nsx1, wsub_inited, wsub_alloc, &
             could_hm_nudging, hm_step

        if(rank.ne.rank1) then
           print *,'Error: rank of host-model restart data is not the same as rank of the process'
           print *,'rank1=',rank1,'   rank=',rank
           call task_abort()
        endif
        if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1.or.nsx.ne.nsx1) then
           print *,'Error: host-model restart grid does not match current executable.'
           print *,'in executable:   nx, ny, nz, nsx:',nx,ny,nz,nsx
           print *,'in restart file: nx, ny, nz, nsx:',nx1,ny1,nz1,nsx1
           call task_abort()
        endif

        if(wsub_alloc) then
          if(.not.allocated(wsub_map)) allocate(wsub_map(nsx,nz))
          read(67) wsub_map
        else
          if(allocated(wsub_map)) deallocate(wsub_map)
        end if

        if(wsub_inited.and..not.allocated(wsub_map)) then
           print *,'Error: host-model restart says MMF is initialized, but wsub_map was not saved.'
           call task_abort()
        endif

        read(67) &
             u0_local_hm, t0_local_hm, q0_local_hm, tabs0_local_hm, &
             qn0_local_hm, qp0_local_hm, qni0_local_hm, qnl0_local_hm, &
             qpi0_local_hm, qpl0_local_hm, prec_flx_local_hm, &
             ug0_hm, tg0_hm, qg0_hm, ug0_press_modify, &
             u_hm_map_save, u_hm_updated_map_save, u_sub_map_save, &
             t_hm_map_save, t_hm_updated_map_save, t_sub_map_save, &
             q_hm_map_save, q_hm_updated_map_save, q_sub_map_save, &
             dudt_hm_hist, dwdt_hm_hist, dudt_subdomain_diffuse, &
             u_external_profile, u_external_profile_loaded, &
             prec_xy_crm, prec_xy_save, shf_xy_crm, shf_xy_save, &
             lhf_xy_crm, lhf_xy_save
        close(67)

        if(rank.eq.nsubdomains-1) then
           print *,'Restarting host model at SAM step:',nstep
           print *,'Host-model step:',hm_step
        endif

        return
        end



        subroutine missing_hm_restart(filename)

        use vars
        implicit none
        character(len=*), intent(in) :: filename

        if(masterproc) then
          print *,'Error: MMF restart requested, but host-model restart file was not found.'
          print *,'Missing file:',trim(filename)
          print *,'This restart cannot preserve MMF state without the _restart_hm.bin file.'
        end if
        call task_abort()

        return
        end
 
 

        subroutine write_statement()

        use vars
        use microphysics, only: micro_field, nmicro_fields
        use sgs, only: sgs_field, nsgs_fields, sgs_field_diag, nsgs_fields_diag
        use tracers
        use params
        use movies, only: irecc
        implicit none

        write(65)  &
         u, v, w, t, p, qv, qcl, qci, qpl, qpi, dudt, dvdt, dwdt, &
         tracer, micro_field, sgs_field, sgs_field_diag, z, pres, prespot, presi, prespoti, &
         rho, rhow, bet, sstxy, precinst, rank, nx, ny, nz, irecc
        close(65)
        if(masterproc) then
           write(66) version, &
            at, bt, ct, dt, dtn, dt3, time, dx, dy, dz, doconstdz,&
            day, day0, nstep, na, nb, nc, caseid, case, &
            dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
            dolongwave, doshortwave, dosgs, dosubsidence, dotracers,  dosmoke, &
            docoriolis, dosurface, dolargescale,doradforcing, dossthomo, &
            dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
            dowallx, dowally, doperpetual, doseasons, &
            docup, docolumn, soil_wetness, dodynamicocean, ocean_type, &
            delta_sst, depth_slab_ocean, Szero, deltaS, timesimpleocean, &
            pres0, ug, vg, fcor, fcorz, tabs_s, z0, fluxt0, fluxq0, tau0, &
            tauls, tautqls, timelargescale, epsv, nudging_uv_z1, nudging_uv_z2, &
            donudging_t, donudging_q, doisccp, domodis, domisr, dosimfilesout, & 
            dosolarconstant, solar_constant, zenith_angle, notracegases, &
            doSAMconditionals, dosatupdnconditionals, LES_S, &
            nudging_t_z1, nudging_t_z2, nudging_q_z1, nudging_q_z2, &
            ocean, land, sfc_flx_fxd, sfc_tau_fxd, &
            nrad, nxco2, latitude0, longitude0, dofplane, SLM, &
            docoriolisz, doradlon, doradlat, doseawater, salt_factor, &
            ntracers, nmicro_fields, nsgs_fields, nsgs_fields_diag
            close(66)
        end if
        if(rank.eq.nsubdomains-1) then
            print *,'Restart file was saved. nstep=',nstep
        endif

        return
        end




        subroutine read_statement()

        use vars
        use microphysics, only: micro_field, nmicro_fields
        use sgs, only: sgs_field, nsgs_fields, sgs_field_diag, nsgs_fields_diag
        use tracers
        use params
        use movies, only: irecc
        implicit none
        integer  nx1, ny1, nz1, rank1, ntr, nmic, nsgs, nsgsd
        character(100) case1,caseid1
        character(7) version1

        read(65)  &
         u, v, w, t, p, qv, qcl, qci, qpl, qpi, dudt, dvdt, dwdt, &
         tracer, micro_field, sgs_field, sgs_field_diag, z, pres, prespot, presi, prespoti, &
         rho, rhow, bet, sstxy, precinst, rank1, nx1, ny1, nz1, irecc
        close(65)
        read(66) version1, &
            at, bt, ct, dt, dtn, dt3, time, dx, dy, dz, doconstdz, &
            day, day0, nstep, na, nb, nc, caseid1(1:sizeof(caseid)), case1(1:sizeof(case)), &
            dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
            dolongwave, doshortwave, dosgs, dosubsidence, dotracers,  dosmoke, &
            docoriolis, dosurface, dolargescale,doradforcing, dossthomo, &
            dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
            dowallx, dowally, doperpetual, doseasons, &
            docup, docolumn, soil_wetness, dodynamicocean, ocean_type,&
            delta_sst, depth_slab_ocean, Szero, deltaS, timesimpleocean, &
            pres0, ug, vg, fcor, fcorz, tabs_s, z0, fluxt0, fluxq0, tau0, &
            tauls, tautqls, timelargescale, epsv, nudging_uv_z1, nudging_uv_z2, &
            donudging_t, donudging_q, doisccp, domodis, domisr, dosimfilesout,  &
            dosolarconstant, solar_constant, zenith_angle, notracegases, &
            doSAMconditionals, dosatupdnconditionals, LES_S, &
            nudging_t_z1, nudging_t_z2, nudging_q_z1, nudging_q_z2, &
            ocean, land, sfc_flx_fxd, sfc_tau_fxd, &
            nrad, nxco2, latitude0, longitude0, dofplane, SLM, &
            docoriolisz, doradlon, doradlat, doseawater, salt_factor, &
            ntr, nmic, nsgs, nsgsd
        close(66)

        if(version1.ne.version) then
          print *,'Wrong restart file!'
          print *,'Version of SAM that wrote the restart files:',version1
          print *,'Current version of SAM',version
          call task_abort()
        end if
        if(nrestart.ne.3) then
          if(rank.ne.rank1) then
             print *,'Error: rank of restart data is not the same as rank of the process'
             print *,'rank1=',rank1,'   rank=',rank
             call task_abort()
          endif
          if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then
             print *,'Error: domain dims (nx,ny,nz) set by grid.f'
             print *,' not correspond to ones in the restart file.'
             print *,'in executable:   nx, ny, nz:',nx,ny,nz
             print *,'in restart file: nx, ny, nz:',nx1,ny1,nz1
             print *,'Exiting...'
             call task_abort()
          endif
        end if
        if(nmic.ne.nmicro_fields) then
           print*,'Error: number of micro_field in restart file is not the same as nmicro_fields'
           print*,'nmicro_fields=',nmicro_fields,'   in file=',nmic
           print*,'Exiting...'
           call task_abort()
        end if
        if(nsgs.ne.nsgs_fields.or.nsgsd.ne.nsgs_fields_diag) then
           print*,'Error: number of sgs_field in restart file is not the same as nsgs_fields'
           print*,'nsgs_fields=',nsgs_fields,'   in file=',nsgs
           print*,'nsgs_fields_diag=',nsgs_fields_diag,'   in file=',nsgsd
           print*,'Exiting...'
           call task_abort()
        end if
        if(ntr.ne.ntracers) then
           print*,'Error: number of tracers in restart file is not the same as ntracers.'
           print*,'ntracers=',ntracers,'   ntracers(in file)=',ntr
           print*,'Exiting...'
           call task_abort()
        end if
        close(65)
        if(rank.eq.nsubdomains-1) then
           print *,'Case:',caseid
           print *,'Restarting at step:',nstep
           print *,'Time(s):',nstep*dt
        endif

        return
        end
 
