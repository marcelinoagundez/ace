
       subroutine compute_chemical_equilibrium
       
       include 'ace.common'       
       integer i,j,k,l,js,iwk,nspec_save,ncondensed_save,
     # idcondensed(nmaxcondensed),ibad
       real*8 rwk,dg,dgmin,abuntot_save
       logical usd(nmaxcondensed),bad(nmaxcondensed)

c      ...initialize output files
       iwk=0                           ! to avoid warning
       init=.true.
       call write_results
       init=.false.

c      ...compute terms b0(i) of initial elemental abundances (+charge)
       rwk=0.0d0
       do i=1,nelem
          rwk=rwk+elfab(i)*mat(i)*amu
       enddo
       abuntot=0.0d0
       b0max=0.0d0
       b0min=huge
       do i=1,nelem
          b0(i)=elfab(i)/navg/rwk      ! [=] mole [g of mixture]-1
          abuntot=abuntot+b0(i)        ! n_tot(gas) for a fully atomic gas
          if(b0(i).gt.b0max)b0max=b0(i)
          if(b0(i).lt.b0min)b0min=b0(i)
       enddo
       abuntot_save=abuntot
       if(ion)then                     ! charge
          b0(nelem+1)=0.0d0
          do j=1,nspec                 ! add charge to nat(i,j) variable
             nat(j,nelem+1)=charge(j)
          enddo
       endif

c      ...compute the maximum abundance reachable by each species
       do j=1,nspec
          abunmax(j)=huge
          do i=1,nelem
             iwk=nat(j,i)
             if(iwk.ne.0.and.abunmax(j).gt.(b0(i)/dble(iwk)))
     #       abunmax(j)=b0(i)/dble(iwk)
          enddo
       enddo
       if(ion.and.id_electr.ne.0)abunmax(id_electr)=abuntot

c      ...assign initial guess for species abundances
       do j=1,nspec
          abun(j)=abuntot/(dble(ngas))
          if(condensed(j))abun(j)=0.0d0
          if(abun(j).gt.abunmax(j))abun(j)=abunmax(j)
       enddo

c Loop over (p,T) grid
       do iapt=1,napt

          altitude=a_apt(iapt)
          pressure=p_apt(iapt)
          tk=t_apt(iapt)

c ...if only gas phase species
          if(ncondensed.eq.0)then
             call minimize_gibbs_energy
             if(.not.converge)then
                if(iapt.eq.1)then
                   abuntot=abuntot_save
                   call initialize_abun
                else
                   stop' E - Convergence not reached'
                endif
             endif
             call write_results
             cycle
          endif

c ...if condensed species are included
          nspec_save=nspec                       ! save condensed species data on js=nspec+1,...
          ncondensed_save=ncondensed
          call compute_chemical_potential
          do i=1,ncondensed_save
             j=ngas+i
             js=ngas+ncondensed_save+i
             call assign_species(j,js)
             abun(js)=0.0d0
          enddo
          nspec=ngas                             ! solve gas phase chemical equilibrium
          ncondensed=0
          call minimize_gibbs_energy

c      ...loop: one calculation per new condensed species
          usd(1:nmaxcondensed)=.false.           ! add (one by one) condensed species and
          bad(1:nmaxcondensed)=.false.           ! solve gas+condensed chemical equilibrium
          do l=1,ncondensed_save

             dgmin=huge                          ! look for new condensed species that reduces Gibbs energy
             do i=1,ncondensed_save
                js=ngas+ncondensed_save+i
                if(usd(i).or.bad(i))cycle
                if(spec_tkout(js))cycle
                rwk=0.0d0
                do k=1,npilag
                   rwk=rwk+nat(js,k)*pilag(k)
                enddo
                dg=mu(js)-rwk
                if(dg.lt.dgmin)then
                   dgmin=dg
                   iwk=i
                endif
             enddo
             if(dgmin.ge.0.0d0)exit

             usd(iwk)=.true.                      ! if condensed species reduces Gibbs energy
10           ncondensed=0                         ! include it and solve gas+condensed chemical equilibrium
             do i=1,ncondensed_save
                if(.not.usd(i).or.bad(i))cycle
                ncondensed=ncondensed+1
                idcondensed(ncondensed)=ngas+ncondensed_save+i
                j=ngas+ncondensed
                js=idcondensed(ncondensed)
                call assign_species(js,j)
                abun(j)=abun(js)
             enddo
             nspec=ngas+ncondensed
             call minimize_gibbs_energy
             do i=1,ncondensed                   ! assign computed abundances of condensed species
                abun(idcondensed(i))=abun(ngas+i)! to js=nspec+1,...
             enddo
             ibad=0
             do i=1,ncondensed
                if(abun(idcondensed(i)).ge.0.0d0)cycle
                ibad=ibad+1
                bad(idcondensed(i)-ngas-ncondensed_save)=.true.
                abun(idcondensed(i))=0.0d0
             enddo
             if(ibad.gt.0)goto 10

c      ...end of loop: one calculation per new condensed species
          enddo

          do i=1,ncondensed_save                 ! assign back condensed species data and abundances
             j=ngas+i
             js=ngas+ncondensed_save+i
             call assign_species(js,j)
             abun(j)=abun(js)
          enddo
          nspec=nspec_save
          ncondensed=ncondensed_save

          call write_results                     ! write out computed abundances

c End of loop over (p,T) grid
       enddo

       return
       end




c_______________________________________________________________________

       subroutine assign_species(ji,jo)

       include 'ace.common'       
       integer ji,jo

c      ...Intrinsic data
       spec(jo)(1:nmaxcharspec)=spec(ji)(1:nmaxcharspec)
       ilenspec(jo)=ilenspec(ji)
       condensed(jo)=condensed(ji)
       nattot(jo)=nattot(ji)
       nat(jo,:)=nat(ji,:)
       charge(jo)=charge(ji)
       type_therm(jo)=type_therm(ji)
       natherm(jo)=natherm(ji)
       ntemp_therm(jo)=ntemp_therm(ji)
       temp_therm(jo,:)=temp_therm(ji,:)
       atherm(jo,:,:)=atherm(ji,:,:)

c      ...T-dependent data (also p-dependent for gas species)
       spec_tkout(jo)=spec_tkout(ji)
       mu(jo)=mu(ji)

       return
       end




c_______________________________________________________________________

       subroutine initialize_abun

       include 'ace.common'       
       integer j

c      ...assign initial guess for species abundances
       do j=1,nspec
          abun(j)=abuntot/(dble(ngas))
          if(abun(j).gt.abunmax(j))abun(j)=abunmax(j)
       enddo

c      ...initialize abun values at 0.1 mbar and 2000 K
       pressure=1.0d-4                 ! [bar]
       tk=2000.0d0                     ! [K]
       call minimize_gibbs_energy

c      ...minimize Gibbs energy at first (p,T) grid point
       pressure=p_apt(1)
       tk=t_apt(1)
       call minimize_gibbs_energy
       if(.not.converge)stop' E- Gibbs minimization not converged'

       return
       end
