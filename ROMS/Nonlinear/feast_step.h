      !==================================================================
      !  FEAST
      !==================================================================
      
      !---------------------------------
      ! Setup
      !---------------------------------
      
      ! NOTE the following nfeast signifies do the biology every nfeast 
      ! time steps, i.e. every day.

      ! number of feast calls per day
      feast_calls = 1

      ! number of dt-sized timesteps per feast call
      ! dt(ng) is 288 (or whatever) steps per day
      nfeast=int((3600*24/feast_calls)/dt(ng))  

      ! Do this the first time through only, syncs modulus of feast_time 
      ! with whatever the starting time-moment of the model was.
      if (feast_time==0) then
        feast_time = MOD(iic(ng),nfeast) + 1
      end if
      
      !---------------------------------
      ! Daily FEAST calculations
      !---------------------------------
      
      ! Note: the calculations are mostly confined to the feast_biol 
      ! subroutine (in mod_feast).  This bit of code just prepares all 
      ! the necessary input and output arguments.
       
      IF ((MOD(iic(ng),nfeast)+1).eq.feast_time) THEN

        ! Print to log at the beginning of each daily calculation

        if (tile==0) then
          print *,"Feast Step ",yday,feast_time,iic(ng),nfeast
        end if
          
        !KYA fish-indices moved to feast_initialize in mod_biology
        !call fish_indices(tspal,tspl,tsp,kspal,kspl,ksp,UBk)

        !KYA no longer a subroutine      
        !subroutine tarray_to_fish(fal,fl,fsp,&
        !          &t,nnew,nstp,LBi,UBi,LBj,UBj,Istr,Iend,Jstr,Jend)

        ! Initialize values
        
        GF % outcatch = 0.0
        GF % ozm      = 0.0
        GF % zoop_death = 0.0   
        ROMS_zoop     = 0.0
        ROMS_temp     = 0.0
        ROMS_depth    = 0.0

        ictr=LBi+(UBi-LBi)/2
        jctr=LBj+(UBj-LBj)/2

        do i=LBi,UBi
          do j=LBj,UBj

            !do fvaral=1,nfvaral
            do sp=1,NUM_AGED_SPECIES
              do lc=1,NUM_AGED_LENGTHS
                do ac=1,NUM_AGES
                  NNstart = t(i,j,kspal(1,sp,lc,ac),nstp,tspal(1,sp,lc,ac))
                  CFstart = t(i,j,kspal(2,sp,lc,ac),nstp,tspal(2,sp,lc,ac))
                  CAstart = t(i,j,kspal(3,sp,lc,ac),nstp,tspal(3,sp,lc,ac))
                  if (isnan(NNstart).or.isnan(CFstart).or.isnan(CAstart)) then
                    NNstart = 0.0
                    CFstart = 1.0
                    CAstart = CALbase_all(rlook_fal(sp,ac,lc))
                  end if      
                  fal(i,j,1,sp,lc,ac)=NNstart
                  fal(i,j,2,sp,lc,ac)=CFstart
                  fal(i,j,3,sp,lc,ac)=CAstart
                end do
              end do
            end do
            !end do

            do fvarl=1,nfvarl
              do spl=1,NUM_LENGTHED_SPECIES
                do lc=1,NUM_NOAGE_LENGTHS 
                  klev=kspl(fvarl,spl,lc)
                  itrc=tspl(fvarl,spl,lc)
                  fl(i,j,fvarl,spl,lc)=t(i,j,klev,nstp,itrc)
                  if (fl(i,j,fvarl,spl,lc)<F_EPSILON) then
                    if (fvarl==1) fl(i,j,fvarl,spl,lc) = F_EPSILON
                    if (fvarl==2) fl(i,j,fvarl,spl,lc) = 1.
                  end if
                end do
              end do
            end do

            do fvar=1,nfvar
              do sp=1,NUM_SIMPLE_SPECIES
                klev=ksp(fvar,sp)
                itrc=tsp(fvar,sp)
                fsp(i,j,fvar,sp)=t(i,j,klev,nstp,itrc)
                if (fsp(i,j,fvar,sp)<F_EPSILON) then
                  if (fvar==1) fsp(i,j,fvar,sp) = F_EPSILON
                end if
              end do
            end do
          end do  ! j loop
        end do   ! i loop

        !END SUBROUTINE tarray_to_fish

        ! Promotion at beginning of year
!! NOPROMOTE if (feast_promotion>0) then
!! NOPROMOTE 
!! NOPROMOTE   if (yday < 1.01) then
!! NOPROMOTE      if (tile==0) then
!! NOPROMOTE       print *,"Beginning year, promoting for day ",yday
!! NOPROMOTE      end if    
!! NOPROMOTE 
!! NOPROMOTE   do i=LBi,UBi
!! NOPROMOTE    do j=LBj,UBj
!! NOPROMOTE 
!! NOPROMOTE      do spal=1,NUM_AGED_SPECIES
!! NOPROMOTE       ac = NUM_AGES ! Plus Group
!! NOPROMOTE        do lc=1,NUM_AGED_LENGTHS 
!! NOPROMOTE         elder = labs_all(rlook_fal(spal,ac,lc) )
!! NOPROMOTE         do lcp = 1,NUM_AGED_LENGTHS
!! NOPROMOTE          younger = labs_all(rlook_fal(spal,ac-1,lcp))
!! NOPROMOTE          if (elder == younger) then
!! NOPROMOTE           ! average CF
!! NOPROMOTE           fal(i,j,2,spal,lc,ac)=&
!! NOPROMOTE            &(fal(i,j,2,spal,lc,ac) * fal(i,j,1,spal,lc,ac) +&
!! NOPROMOTE            &fal(i,j,2,spal,lcp,ac-1) * fal(i,j,1,spal,lcp,ac-1))/&
!! NOPROMOTE            &(fal(i,j,1,spal,lc,ac) + fal(i,j,1,spal,lcp,ac-1)+1.0e-30)
!! NOPROMOTE           ! average CAL
!! NOPROMOTE           fal(i,j,3,spal,lc,ac) =&
!! NOPROMOTE            &(fal(i,j,3,spal,lc,ac) * fal(i,j,1,spal,lc,ac) +&
!! NOPROMOTE            &fal(i,j,3,spal,lcp,ac-1) * fal(i,j,1,spal,lcp,ac-1))/&
!! NOPROMOTE            &(fal(i,j,1,spal,lc,ac) + fal(i,j,1,spal,lcp,ac-1)+1.0e-30)
!! NOPROMOTE                     
!! NOPROMOTE           ! change N                    
!! NOPROMOTE           fal(i,j,1,spal,lc,ac) = fal(i,j,1,spal,lc,ac) +&
!! NOPROMOTE                                  &fal(i,j,1,spal,lcp,ac-1)
!! NOPROMOTE          end if
!! NOPROMOTE         end do
!! NOPROMOTE        end do
!! NOPROMOTE        do lc=1,NUM_AGED_LENGTHS
!! NOPROMOTE         fal(i,j,1,spal,lc,ac-1) = 0.0 
!! NOPROMOTE         fal(i,j,2,spal,lc,ac-1) = 1.0                    
!! NOPROMOTE         fal(i,j,3,spal,lc,ac-1) = 1.0
!! NOPROMOTE        end do
!! NOPROMOTE             ! end plus group                 
!! NOPROMOTE             do ac = NUM_AGES-1,2,-1
!! NOPROMOTE               do lc=1,NUM_AGED_LENGTHS 
!! NOPROMOTE                 elder = labs_all(rlook_fal(spal,ac,lc) )
!! NOPROMOTE                 do lcp = 1,NUM_AGED_LENGTHS
!! NOPROMOTE                    younger = labs_all(rlook_fal(spal,ac-1,lcp))          
!! NOPROMOTE                    if (elder == younger) then
!! NOPROMOTE                       !if (tile==0) then
!! NOPROMOTE                       !   print *,"promoting",spal,ac-1,lcp," to ",ac,lc
!! NOPROMOTE                       !end if
!! NOPROMOTE                       fal(i,j,1,spal,lc,ac) = fal(i,j,1,spal,lcp,ac-1)
!! NOPROMOTE                       fal(i,j,2,spal,lc,ac) = fal(i,j,2,spal,lcp,ac-1)
!! NOPROMOTE                       fal(i,j,3,spal,lc,ac) = fal(i,j,3,spal,lcp,ac-1)
!! NOPROMOTE                    end if
!! NOPROMOTE                 end do          
!! NOPROMOTE               end do
!! NOPROMOTE              do lc=1,NUM_AGED_LENGTHS
!! NOPROMOTE               fal(i,j,1,spal,lc,ac-1) = 0.0 
!! NOPROMOTE               fal(i,j,2,spal,lc,ac-1) = 1.0                    
!! NOPROMOTE               fal(i,j,3,spal,lc,ac-1) = 1.0
!! NOPROMOTE              end do  
!! NOPROMOTE             end do ! end age loop
!! NOPROMOTE           end do !end spal loop         
!! NOPROMOTE 
!! NOPROMOTE   ! Now copy back (CHANGES DON'T HAPPEN HERE!)
!! NOPROMOTE     do fvaral=1,nfvaral
!! NOPROMOTE      do spal=1,NUM_AGED_SPECIES
!! NOPROMOTE       do lc=1,NUM_AGED_LENGTHS
!! NOPROMOTE        do ac=1,NUM_AGES
!! NOPROMOTE         klev=kspal(fvaral,spal,lc,ac)
!! NOPROMOTE         itrc=tspal(fvaral,spal,lc,ac)
!! NOPROMOTE         !if ((tile==0).and.(yday < 1.1)) then
!! NOPROMOTE           !if (fvaral==1) then
!! NOPROMOTE           ! print *,"setting",t(i,j,klev,nnew,itrc),fal(i,j,fvaral,spal,lc,ac)
!! NOPROMOTE           !end if
!! NOPROMOTE         !end if   
!! NOPROMOTE         t(i,j,klev,nstp,itrc) = fal(i,j,fvaral,spal,lc,ac)
!! NOPROMOTE         t(i,j,klev,nnew,itrc) = fal(i,j,fvaral,spal,lc,ac)
!! NOPROMOTE        end do
!! NOPROMOTE       end do
!! NOPROMOTE      end do
!! NOPROMOTE     end do
!! NOPROMOTE    end do  ! j loop
!! NOPROMOTE   end do   ! i loop
!! NOPROMOTE   end if   ! end promotion
!! NOPROMOTE end if     ! end promotion flag

        do i=LBi,UBi
          do j=LBj,UBj
        
            ! copy layer thickness over
            ! ROMS_edges stores depths as negative (down from zero)
            !CurD = 0.0
            !DO k=N(ng),1,-1
            !  ROMS_depth(i,j,k)= 2 * ((-z_r(i,j,k) - CurD))
            !  ROMS_edges(i,j,k+1) = -z_w(i,j,k)
            !  ROMS_edges(i,j,k+1) = -CurD 
            !  CurD = CurD + ROMS_depth(i,j,k)
            !end DO
            !ROMS_edges(i,j,1) = -z_w(i,j,0)
            !ROMS_edges(i,j,1) = -CurD
       
            !if ((i.eq.ictr).and.(j.eq.jctr)) then
            ROMS_edges(i,j,N(ng)+1) = 0.0
            do k=1,N(ng)
              ROMS_edges(i,j,k) = z_w(i,j,k-1) - z_w(i,j,N(ng))
              ROMS_depth(i,j,k) = z_w(i,j,k)   - z_w(i,j,k-1) 
              !write(*,'(a6,3I6,2ES14.4)')"EDGE",i,j,k,ROMS_edges(i,j,k),z_w(i,j,k)
            end do
        
            !if ((i.eq.100).and.(j.eq.100)) then
            !  do k=1,N(ng)
            !    write(*,'(a6,3I6,2ES14.4)')"EDGE",i,j,k,ROMS_edges(i,j,k),ROMS_depth(i,j,k)
            !  end do
            !end if
        
            ! copy other needed values
            DO k=1,N(ng)
              ! Temperature
              ROMS_temp(i,j,k)   = MAX(0.0,t(i,j,k,nstp,itemp))
              ! copepods  
              ROMS_zoop(i,j,k,1) = MAX(0.0,t(i,j,k,nstp,iCop))
              ! neocalanus
              ROMS_zoop(i,j,k,2) = MAX(0.0,t(i,j,k,nstp,iNCaS))  
              ROMS_zoop(i,j,k,3) = MAX(0.0,t(i,j,k,nstp,iNCaO))
              ! euphausiids 
              ROMS_zoop(i,j,k,4) = MAX(0.0,t(i,j,k,nstp,iEupS))  
              ROMS_zoop(i,j,k,5) = MAX(0.0,t(i,j,k,nstp,iEupO))
              ! zero out benthos
              ROMS_zoop(i,j,k,6) = 0.0
              do sp=1,6
                if (isnan(ROMS_zoop(i,j,k,sp))) then
                  !write(*,'(A6,4i4)'),"ZNAN",i,j,k,sp
                  ROMS_zoop(i,j,k,sp) = 0.0
                end if
              end do
            END DO ! k loop
            ! benthos in bottom layer, converted to m^3
            ROMS_zoop(i,j,1,6) = MAX(0.0,                                 &
     &                           bt(i,j,1,nstp,iBen))/ROMS_depth(i,j,1)  
        
          end do ! j loop
        end do ! i loop
        
        !Now recruitment copy tracer first if, then copy tracer back second if
        !if (feast_recruitment>0) then    
        !do sp=1,TOT_LENGTHED
        !  GF % eggs(sp,:,:) = t(:,:,keggs(sp),nstp,teggs(sp))
        !end do
        !end if            
        
        do ac = 1,11
          GF%sp_ibm(ac) = 1
          GF%lc_ibm(ac) = 7
          GF%ac_ibm(ac) = ac
        end do
          
        !call the main feast procedure
        
        call feast_biol(N(ng),rmask,ROMS_depth,ROMS_edges,ROMS_temp,      &
     &                  ROMS_zoop,fal,fl,fsp,GF,yday,LBi,UBi,LBj,UBj)
             
        if (tile==0) then
          print *,"Feast Biol called ",iic(ng)
        end if
        
        !if (feast_recruitment>0) then
        !  do sp=1,TOT_LENGTHED
        !         t(:,:,keggs(sp),nstp,teggs(sp)) = GF % eggs(sp,:,:)  
        !         t(:,:,keggs(sp),nnew,teggs(sp)) = GF % eggs(sp,:,:)
        !  end do
        !end if   !end recruitment flag
        
        ! 11/6/11 The happiness calculations that used to be right here
        ! were moved to the end of mod_feast 

      END IF ! conditional execution of daily feast stuff 

      !---------------------------------
      ! Place all newly-calculated 
      ! fish variables back in tracer 
      ! arrays
      !---------------------------------
      ! Note: Is HZ needed?

      !ftstp = dtdays*feast_calls*2.
      ftstp = dtdays ! *feast_calls*1.

      ! first replace the zooplankton
      ictr=LBi+(UBi-LBi)/2
      jctr=LBj+(UBj-LBj)/2
      DO k=1,N(ng)
        do i=Istr,Iend
          do j=Jstr,Jend
            !DO i=Istr,Iend
            ! DO j=Jstr,Jend
            ! March 21, 2012 changed ozm to a survival rate scaled to fstep in mod_feast
            ! ozm is guaranteed in mod feast to be between 1e-6 and 1.0 
            GF%zoop_death(i,j,k,1) = GF%zoop_death(i,j,k,1) + t(i,j,k,nnew,iCop  ) * (1.0-GF%ozm(i,j,k,1))
            GF%zoop_death(i,j,k,2) = GF%zoop_death(i,j,k,2) + t(i,j,k,nnew,iNCaS ) * (1.0-GF%ozm(i,j,k,2))    
            GF%zoop_death(i,j,k,3) = GF%zoop_death(i,j,k,3) + t(i,j,k,nnew,iNCaO ) * (1.0-GF%ozm(i,j,k,3))
            GF%zoop_death(i,j,k,4) = GF%zoop_death(i,j,k,4) + t(i,j,k,nnew,iEupS ) * (1.0-GF%ozm(i,j,k,4))
            GF%zoop_death(i,j,k,5) = GF%zoop_death(i,j,k,5) + t(i,j,k,nnew,iEupO ) * (1.0-GF%ozm(i,j,k,5))
            if (k==1) GF%zoop_death(i,j,k,6) = GF%zoop_death(i,j,k,6) + t(i,j,k,nnew,iBen ) * (1.0-GF%ozm(i,j,k,6))
#ifdef STATIONARY
            st(i,j,k,nnew,114) = (1.0-GF%ozm(i,j,k,1))!GF%zoop_death(i,j,k,1)
            st(i,j,k,nnew,115) = (1.0-GF%ozm(i,j,k,2))!GF%zoop_death(i,j,k,2) 
            st(i,j,k,nnew,116) = (1.0-GF%ozm(i,j,k,3))!GF%zoop_death(i,j,k,3)
            st(i,j,k,nnew,117) = (1.0-GF%ozm(i,j,k,4))!GF%zoop_death(i,j,k,4)
            st(i,j,k,nnew,118) = (1.0-GF%ozm(i,j,k,5))!GF%zoop_death(i,j,k,5)

!             pt3(i,j,k,nnew,iFCopMort)  = (1.0-GF%ozm(i,j,k,1))!GF%zoop_death(i,j,k,1)
!             pt3(i,j,k,nnew,iFNCaSMort) = (1.0-GF%ozm(i,j,k,2))!GF%zoop_death(i,j,k,2)
!             pt3(i,j,k,nnew,iFNCaOMort) = (1.0-GF%ozm(i,j,k,3))!GF%zoop_death(i,j,k,3)
!             pt3(i,j,k,nnew,iFEupSMort) = (1.0-GF%ozm(i,j,k,4))!GF%zoop_death(i,j,k,4)
!             pt3(i,j,k,nnew,iFEupOMort) = (1.0-GF%ozm(i,j,k,5))!GF%zoop_death(i,j,k,5)
            !ictr=129-1
            !jctr=67-1
            !if ((i.eq.ictr).and.(j.eq.jctr)) then
            !   write (*,'(A8,i6,4i4,F8.2,6ES18.8)'),"DEATH",iic(ng),tile,k,i,j, &
            !  & rmask(i,j), pt3(i,j,k,nnew,iFCopMort), &
            !  & pt3(i,j,k,nnew,iFNCaSMort),pt3(i,j,k,nnew,iFNCaOMort),  &
            !  & pt3(i,j,k,nnew,iFEupSMort),pt3(i,j,k,nnew,iFEupOMort)
            !end if
#endif

            if (feast_coupled>0) then     
              t(i,j,k,nnew,iCop )=max(0.0,t(i,j,k,nnew,iCop ) * GF%ozm(i,j,k,1))
              t(i,j,k,nnew,iNCaS)=max(0.0,t(i,j,k,nnew,iNCaS) * GF%ozm(i,j,k,2))   
              t(i,j,k,nnew,iNCaO)=max(0.0,t(i,j,k,nnew,iNCaO) * GF%ozm(i,j,k,3))           
              t(i,j,k,nnew,iEupS)=max(0.0,t(i,j,k,nnew,iEupS) * GF%ozm(i,j,k,4))  
              t(i,j,k,nnew,iEupO)=max(0.0,t(i,j,k,nnew,iEupO) * GF%ozm(i,j,k,5))
              if (k==1) bt(i,j,1,nnew,iBen)=max(0.0,bt(i,j,1,nnew,iBen) * GF%ozm(i,j,1,6))         
              if (isnan(t(i,j,k,nnew, iCop ))) write(*,'(A6,4i4,3ES14.4)'),"ONAN",1,i,j,k,t(i,j,k,nstp, iCop),t(i,j,k,nnew, iCop),GF%ozm(i,j,k,1)
              if (isnan(t(i,j,k,nnew,iNCaS ))) write(*,'(A6,4i4,3ES14.4)'),"ONAN",2,i,j,k,t(i,j,k,nstp,iNCaS),t(i,j,k,nnew,iNCaS),GF%ozm(i,j,k,2)
              if (isnan(t(i,j,k,nnew,iNCaO ))) write(*,'(A6,4i4,3ES14.4)'),"ONAN",3,i,j,k,t(i,j,k,nstp,iNCaO),t(i,j,k,nnew,iNCaO),GF%ozm(i,j,k,3)
              if (isnan(t(i,j,k,nnew,iEupS ))) write(*,'(A6,4i4,3ES14.4)'),"ONAN",4,i,j,k,t(i,j,k,nstp,iEupS),t(i,j,k,nnew,iEupS),GF%ozm(i,j,k,4)
              if (isnan(t(i,j,k,nnew,iEupO ))) write(*,'(A6,4i4,3ES14.4)'),"ONAN",5,i,j,k,t(i,j,k,nstp,iEupO),t(i,j,k,nnew,iEupO),GF%ozm(i,j,k,5)
            end if
     
        
          END DO
        END DO
      END DO

      ! FISH AND FISHERIES, then GROWTH    
      do i=LBi,UBi
        do j=LBj,UBj

          ! FISHING AND FISHERIES
          if (feast_fishing == 1) then
            !do sp=1,NUM_AGED_SPECIES
            !  CC = 0.0
            !  do lc=1,NUM_AGED_LENGTHS
            !    do ac=1,NUM_AGES

            do spn=1,ALL_LINKS
              sp = sp_all(spn)
              lc = lc_all(spn)
              ac = age_all(spn)
      
              ! First set starting values
              SELECT CASE (type_all(spn))
                CASE (1)   ! aged
                  NNstart = t(i,j,kspal(1,sp,lc,ac),nstp,tspal(1,sp,lc,ac))
                  CFstart = t(i,j,kspal(2,sp,lc,ac),nstp,tspal(2,sp,lc,ac))           
                  CAstart = t(i,j,kspal(3,sp,lc,ac),nstp,tspal(3,sp,lc,ac))
                  ratepar = GF%ghfal(i,j,:,sp,lc,ac) 
                CASE (2)   ! lengthed
                  NNstart = t(i,j,kspl(1,sp-NUM_AGED_SPECIES,lc),nstp,tspl(1,sp-NUM_AGED_SPECIES,lc))
                  CFstart = t(i,j,kspl(2,sp-NUM_AGED_SPECIES,lc),nstp,tspl(2,sp-NUM_AGED_SPECIES,lc))           
                  CAstart = t(i,j,kspl(3,sp-NUM_AGED_SPECIES,lc),nstp,tspl(3,sp-NUM_AGED_SPECIES,lc))           
                  ratepar = GF%ghfl(i,j,:,sp-NUM_AGED_SPECIES,lc)
                CASE (3)   ! simple       
              END SELECT
            
              if (isnan(NNstart).or.isnan(CFstart).or.isnan(CAstart)) then
                NNstart = 0.0
                CFstart = 1.0
                CAstart = CALbase_all(spn)
              end if  

              spy = spn  !spy     = rlook_fal(sp,ac,lc)
              WW      = all_AL_LL_BL(spy) * CFstart  
              !BB      = BB + NNstart * WW 
              !delN   = 0.0
              Nloss  = 0.0
              Ftot   = 0.0
              do gr=1,NUM_GEARS
                !delN = delN + GF%catch_all(gr,spy,i,j) !* ftstp
                Ftot = Ftot + GF%frate_all(gr,spy,i,j)
              end do                 
              if (Ftot > 0.0) then
                !Navail = min(delN,NNstart*0.99)
                Fmax = min(0.99,Ftot) * ftstp
                do gr=1,NUM_GEARS
                  !Ngloss  =  Navail * GF%catch_all(gr,spy,i,j)/delN
                  Ngloss = NNstart * Fmax *  GF%frate_all(gr,spy,i,j) / Ftot
                  Nloss  = Nloss + Ngloss
                  !NNstart = NNstart - Nloss
                  GF%outcatch(i,j,gr,sp) = GF%outcatch(i,j,gr,sp) + Ngloss * WW
                  !CC = CC + Ngloss * WW
                end do
              end if  

              !if ((i .eq. (88-1)).and.(j .eq. (60-1))) then
              !   if (sp .eq. 1) then
              !       write(*,'(a10,2I4,3E14.4)'), "Nstart",&
              !      & lc,ac,real(NNstart),real(Nloss),real(Fmax)
              !   end if
              !end if
      
              ! REMOVE FISHING FROM NNStart
              !NNstart = NNstart - Nloss    
              SELECT CASE (type_all(spn))
                CASE (1)   ! aged
                  t(i,j,kspal(1,sp,lc,ac),nnew,tspal(1,sp,lc,ac)) = max(NNstart-Nloss ,   0.0)
                CASE (2)
                  ! Turned off fishing update 10/26/12
                  ! t(i,j,kspl(1,sp-NUM_AGED_SPECIES,lc),nnew,tspl(1,sp-NUM_AGED_SPECIES,lc)) = max(NNstart-Nloss ,   0.0)
              END SELECT
              !    end do ! end ac loop
              !  end do  ! end lc loop        
              !end do ! END OF AGED SPECIES LOOP
            end do
          end if ! END OF FISHING AGED SPECIES

          ! GROWTH, NATURAL MORTALITY, and RECRUITMENT
          Npromote = 0.0

          do sp=1,TOT_LENGTHED
            if (yday.lt.fsh_sp_sday(sp) .or. yday.gt.fsh_z_eday(sp)) then
              eggs(sp) = 0.0
            else
              eggs(sp) = t(i,j,keggs(sp),nstp,teggs(sp))
            end if
            !if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
            ! if (sp.eq.1) then
            !  write (*,'(A14,2I6,F14.2,4I4,ES14.4)'),    & 
            !  & "FTILE EGGBEGIN",tile,iic(ng),yday,i,j,ac,lc, &
            !  & eggs(sp)     
            ! end if
            !end if    
          end do

          !do sp=1,NUM_AGED_SPECIES
          !    do lc=1,NUM_AGED_LENGTHS
          !     do ac=1,NUM_AGES
          do spn=1,ALL_LINKS
            sp = sp_all(spn)
            lc = lc_all(spn)
            ac = age_all(spn)
            SELECT CASE (type_all(spn))
              CASE (1)   ! aged
                NNstart = t(i,j,kspal(1,sp,lc,ac),nnew,tspal(1,sp,lc,ac))
                CFstart = t(i,j,kspal(2,sp,lc,ac),nstp,tspal(2,sp,lc,ac))           
                CAstart = t(i,j,kspal(3,sp,lc,ac),nstp,tspal(3,sp,lc,ac))
                ratepar = GF%ghfal(i,j,:,sp,lc,ac) 
              CASE (2)   ! lengthed
                NNstart = t(i,j,kspl(1,sp-NUM_AGED_SPECIES,lc),nnew,tspl(1,sp-NUM_AGED_SPECIES,lc))
                CFstart = t(i,j,kspl(2,sp-NUM_AGED_SPECIES,lc),nstp,tspl(2,sp-NUM_AGED_SPECIES,lc))           
                CAstart = t(i,j,kspl(3,sp-NUM_AGED_SPECIES,lc),nstp,tspl(3,sp-NUM_AGED_SPECIES,lc))           
                ratepar = GF%ghfl(i,j,:,sp-NUM_AGED_SPECIES,lc)
              CASE (3)   ! simple       
            END SELECT

            !spn = rlook_fal(sp,ac,lc)      
            !NNstart  = t(i,j,kspal(1,sp,lc,ac),nnew,tspal(1,sp,lc,ac))  
            !CFstart  = t(i,j,kspal(2,sp,lc,ac),nstp,tspal(2,sp,lc,ac))
            !CAstart  = t(i,j,kspal(3,sp,lc,ac),nstp,tspal(3,sp,lc,ac))
            if (isnan(NNstart).or.isnan(CFstart).or.isnan(CAstart)) then
              NNstart = 0.0
              CFstart = 1.0
              CAstart = CALbase_all(spn)
            end if  
            ! ftstep replaced with exprate in mod_feast
            ! promote in delG is stored as proportion remaining, e.g. delG=0.95 means 
            ! that 5% promote.  To get change, multiply NNstart by (1-delG)
            if (next_len(spn).gt.0) then
              Npromote(spn)  = NNstart * (1.0 - ratepar(id_delG) ) * feast_growth
            else
              Npromote(spn) = 0.0
            end if
            ! not a rate so use ftstep
            !rec = GF%ghfal(i,j,id_delRec,sp,lc,ac) * feast_recruitment * ftstp
      
            if (yday.ge.fsh_z_sday(sp) .and. yday.le.fsh_z_eday(sp) ) then
              rec      = eggs(sp) * 1.0/(fsh_z_eday(sp)-yday+1) *       &
     &                   all_zeros(spn) * feast_recruitment * ftstp
              eggs(sp) = MAX(eggs(sp)-rec,0.0)
            end if   
        
            NNbase(spn) = max(0.0,NNstart + rec)  
            if (NNbase(spn).gt.0) then
              CFbase(spn) = (1.0             *rec + CFstart*NNstart) / NNbase(spn)
              CAbase(spn) = (CALbase_all(spn)*rec + CAstart*NNstart) / NNbase(spn)      
            else
              CFbase(spn) = 1.0
              CAbase(spn) = CALbase_all(spn)    !# CALbase_all(rlook_fal(sp,ac,lc)) 
            end if
            !IF ((MOD(iic(ng),nfeast)+1).eq.feast_time) THEN      
            !  if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
            !    if (sp.eq.1) then
            !   write (*,'(A14,2I6,F14.2,4I4,8ES14.4)'),    & 
            !    & "FTILE PROM",tile,iic(ng),yday,i,j,ac,lc, &
            !   & NNstart,rec,NNbase(spn),CFstart,CFbase(spn),CAstart,CAbase(spn),Npromote(spn)       
            !    end if
            !    end if
            !END IF
            !    end do
            !   end do
            !end do    
          end do

          do spn=1,ALL_LINKS
            sp = sp_all(spn)
            lc = lc_all(spn)
            ac = age_all(spn)
            SELECT CASE (type_all(spn))
              CASE (1)   ! aged
                ratepar = GF%ghfal(i,j,:,sp,lc,ac) 
              CASE (2)   ! lengthed      
                ratepar = GF%ghfl(i,j,:,sp-NUM_AGED_SPECIES,lc)
              CASE (3)   ! simple       
            END SELECT
          
            if (prev_len(spn).gt.0) then
              Nloss   = Npromote(prev_len(spn))
            else
              Nloss   = 0.0
            end if
     
            NNstart = MAX(0.0, NNbase(spn) - Npromote(spn) + Nloss)                            
            if (NNstart.gt.0) then
              if (prev_len(spn).gt.0) then
                CFstart  = (CFbase(spn) * (NNbase(spn) -              &
     &                       Npromote(spn)) + CFbase(prev_len(spn)) *  &
     &                       Nloss) / NNstart 
                CAstart  = (CAbase(spn) * (NNbase(spn) -              &
     &                       Npromote(spn)) + CAbase(prev_len(spn)) *  &
     &                       Nloss) / NNstart 
              else
                CFstart  = CFbase(spn)
                CAstart  = CAbase(spn)
              end if   
            else
              CFstart = 1.0
              CAstart = CALbase_all(spn)
            end if

            !IF ((MOD(iic(ng),nfeast)+1).eq.feast_time) THEN 
            !     if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
            !       if (sp.eq.1) then
            !       write (*,'(A14,2I6,F14.2,4I4,9ES14.4)'),    & 
            !       & "FTILE GROW",tile,iic(ng),yday,i,j,ac,lc, &
            !       & NNstart,NNbase(spn),CFstart,CFbase(spn),CAstart,CAbase(spn),Npromote(spn),Npromote(prev_len(spn))       
            !      end if
            !    end if
            !END IF
  
            ! ftstp replaced with exprate in mod feast    
            !if (sp .eq. 1) then
            ! These are all stored as proportion change; converted here for addition
            delN   = NNstart * (ratepar(id_delTZ) - 1.0 ) * feast_mort               
            delCF  = CFstart * (ratepar(id_delCF) - 1.0 ) * feast_growth ! rate
            delCAL = CAstart * (ratepar(id_delCAL) -1.0 ) * feast_growth ! rate    
            !end if
                
            !IF ((MOD(iic(ng),nfeast)+1).eq.feast_time) THEN          
            !     if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
            !      if (sp.eq.1) then
            !       write (*,'(A14,2I6,F14.2,4I4,6ES14.4)'),    & 
            !       & "FTILE BAL",tile,iic(ng),yday,i,j,ac,lc, &
            !       & NNstart,delN,CFstart,delCF,CAstart,delCAL
            !        end if
            !      end if
            !END IF
            eggs(sp) = eggs(sp) + NNstart * ratepar(id_delSpawn) * feast_recruitment
      
            ! UPDATE ALL STATE VARIABLES 
            SELECT CASE (type_all(spn))
              CASE (1)   ! aged
                t(i,j,kspal(1,sp,lc,ac),nnew,tspal(1,sp,lc,ac)) = min(max(NNstart+delN,   0.0),1e35)
                t(i,j,kspal(2,sp,lc,ac),nnew,tspal(2,sp,lc,ac)) = min(max(CFstart+delCF, 0.01),1e35)
                t(i,j,kspal(3,sp,lc,ac),nnew,tspal(3,sp,lc,ac)) = min(max(CAstart+delCAL, 1.0),1e35)
              CASE (2)  
                ! 10/26/12 NO MORE GROWTH OR DEATH UPDATES
                !   t(i,j,kspl(1,sp-NUM_AGED_SPECIES,lc),nnew,tspl(1,sp-NUM_AGED_SPECIES,lc)) = min(max(NNstart+delN,   0.0),1e35)
                !   t(i,j,kspl(2,sp-NUM_AGED_SPECIES,lc),nnew,tspl(2,sp-NUM_AGED_SPECIES,lc)) = min(max(CFstart+delCF, 0.01),1e35)
                !   t(i,j,kspl(3,sp-NUM_AGED_SPECIES,lc),nnew,tspl(3,sp-NUM_AGED_SPECIES,lc)) = min(max(CAstart+delCAL, 1.0),1e35)  
            END SELECT        
            !end do ! end ac loop
            !end do  ! end lc loop    

            !if (((i .eq. (97-1)).and.(j .eq. (101-1)))) then
            ! if (sp.eq.1) then
            !  write (*,'(A14,2I6,F14.2,4I4,2ES14.4)'),    & 
            !  & "FTILE EGGEND",tile,iic(ng),yday,i,j,ac,lc, &
            !  & eggs(sp),t(i,j,keggs(sp),nnew,teggs(sp)) 
            !  end if
            !end if
    
          end do ! END OF AGED OR LENGTHED SPECIES LOOP

          do sp=1,TOT_LENGTHED
            t(i,j,keggs(sp),nnew,teggs(sp)) = min(max(eggs(sp), 0.0),1e35)
          end do
!------------------------------------------------------------------------

          ! SIMPLE UPDATE LOOP
          do fvar=1,nfvar
            do sp=1,NUM_SIMPLE_SPECIES
              klev=ksp(fvar,sp)
              itrc=tsp(fvar,sp)
              SELECT CASE (fvar)                    
                CASE (1)   ! Biomass
                  !t(i,j,klev,nnew,itrc)= t(i,j,klev,nnew,itrc)
                CASE (2)   ! Calories

                !t(i,j,klev,nnew,itrc)= t(i,j,klev,nnew,itrc)
                !t(i,j,klev,nnew,itrc)= t(i,j,klev,nnew,itrc)
              END SELECT
    
              t(i,j,klev,nnew,itrc) = min(max(t(i,j,klev,nnew,itrc),0.0),1e35)
            end do
          end do

        end do ! j loop
      end do ! i loop

      !end if
      
      !---------------------------------
      ! Fish movement
      !---------------------------------

      !test a bit of advection..... 
      if (feast_movement>0) then
        !write(*,'(A6,5i6)'),"FHEY",tile,LBi,LBj,UBi,UBj         
        do spn=1,ALL_LINKS
          ! sp is the species numbers of aged fish followed by lengthed fish
          ! e.g. 1=POL, 2=COD, 3=ATF, 4=HER, 5=CAP etc.
          sp = sp_all(spn)
          ! Length class (1 through 14 for aged, 1 through 20 for lengthed)
          lc = lc_all(spn)
          ! Age class 1=age 0, 2=age 1, 0 if not an aged fish
          ac = age_all(spn)
          !write(4,"(I5) (I5) (I5) (I5)") spn, sp, lc, ac
          SELECT CASE (type_all(spn))
            CASE (1)   ! aged
              !klev=kspal(1,sp,lc,ac)
              !itrc=tspal(1,sp,lc,ac)
              t_new = t(:,:,kspal(1,sp,lc,ac),nnew,tspal(1,sp,lc,ac))
              CFmat = t(:,:,kspal(2,sp,lc,ac),nnew,tspal(2,sp,lc,ac))           
              CAmat = t(:,:,kspal(3,sp,lc,ac),nnew,tspal(3,sp,lc,ac)) 
            CASE (2)   ! lengthed
              !klev=kspl(1,sp-NUM_AGED_SPECIES,lc)
              !itrc=tspl(1,sp-NUM_AGED_SPECIES,lc)
              t_new = t(:,:,kspl(1,sp-NUM_AGED_SPECIES,lc),nnew,tspl(1,sp-NUM_AGED_SPECIES,lc))
              CFmat = t(:,:,kspl(2,sp-NUM_AGED_SPECIES,lc),nnew,tspl(2,sp-NUM_AGED_SPECIES,lc))           
              CAmat = t(:,:,kspl(3,sp-NUM_AGED_SPECIES,lc),nnew,tspl(3,sp-NUM_AGED_SPECIES,lc))           
            CASE (3)   ! simple       
          END SELECT

          do i=LBi,UBi
            do j=LBj,UBj 
              if ((i.le.1).or.(i.ge.182).or.(j.le.1).or.(j.ge.258)) then
                t_new(i,j)=0.0
                CFmat(i,j)=1.0
                CAmat(i,j)=CALbase_all(spn)
              end if
            end do
          end do
                
          !if ((sp .eq. 1).and.(((ac.eq.2).and.(lc.eq.7)).or.((ac.eq.11).and.(lc.eq.14)))) then
          !   do i=LBi,UBi
          !    do j=LBj,UBj  
          !       if ((i.eq.1).or.(i.eq.182).or.(j.eq.1).or.(j.eq.258)) then
          !          write(*,'(A6,2i6,3ES18.8)'),"FSID",i,j,t_new(i,j),CFmat(i,j),CAmat(i,j)
          !       end if
          !    end do
          !   end do
          !end if
          ! swim_all in m/sec; divide by approx 10km cell edge and mult by
          ! number of seconds
          base_speed = swim_all(spn) * (1./10000.) * dt(ng)     
          diff_coeff = 0.4
       
          ! diffusion      
          diff_xiup = 0.0 
          diff_xidn = 0.0
          diff_etup = 0.0 
          diff_etdn = 0.0    
          
          diff_xiup(LBi:(UBi-1),:) = (t_new(LBi:(UBi-1),:) - t_new((LBi+1):UBi,:))
          diff_xidn((LBi+1):UBi,:) = (t_new((LBi+1):UBi,:) - t_new(LBi:(UBi-1),:))
          diff_etup(:,LBj:(UBj-1)) = (t_new(:,LBj:(UBj-1)) - t_new(:,(LBj+1):UBj))
          diff_etdn(:,(LBj+1):UBj) = (t_new(:,(LBj+1):UBj) - t_new(:,LBj:(UBj-1)))     
          
          diff_xiup(LBi:(UBi-1),:) = diff_xiup(LBi:(UBi-1),:) * rmask((LBi+1):UBi,:)
          diff_xidn((LBi+1):UBi,:) = diff_xidn((LBi+1):UBi,:) * rmask(LBi:(UBi-1),:)
          diff_etup(:,LBj:(UBj-1)) = diff_etup(:,LBj:(UBj-1)) * rmask(:,(LBj+1):UBj)
          diff_etdn(:,(LBj+1):UBj) = diff_etdn(:,(LBj+1):UBj) * rmask(:,LBj:(UBj-1))
          
          diff_xiup = diff_xiup * diff_coeff * base_speed * rmask
          diff_xidn = diff_xidn * diff_coeff * base_speed * rmask
          diff_etup = diff_etup * diff_coeff * base_speed * rmask
          diff_etdn = diff_etdn * diff_coeff * base_speed * rmask
      
          ! directed
          u_move = 0.0
          v_move = 0.0
          u_up   = 0.0
          u_dn   = 0.0
          v_up   = 0.0
          v_dn   = 0.0

          u_move    = GF%u_swim(spn,:,:) * base_speed * t_new * rmask
          v_move    = GF%v_swim(spn,:,:) * base_speed * t_new * rmask       
       
          u_up((LBi+1):UBi,:) = max(0.0, -u_move((LBi+1):UBi,:)) * rmask(LBi:(UBi-1),:)  
          u_dn(LBi:(UBi-1),:) = max(0.0, +u_move(LBi:(UBi-1),:)) * rmask((LBi+1):UBi,:)
          v_up(:,(LBj+1):UBj) = max(0.0, -v_move(:,(LBj+1):UBj)) * rmask(:,LBj:(UBj-1))
          v_dn(:,LBj:(UBj-1)) = max(0.0, +v_move(:,LBj:(UBj-1))) * rmask(:,(LBj+1):UBj)   
                
          ! bringing together
       
          t_new  = t_new - u_up - u_dn - v_up - v_dn       

          t_left  = t_new
          CF_left = CFmat * t_left 
          CA_left = CAmat * t_left
       
          t_new(LBi:(UBi-1),:) = t_new(LBi:(UBi-1),:) + u_up((LBi+1):UBi,:)
          t_new((LBi+1):UBi,:) = t_new((LBi+1):UBi,:) + u_dn(LBi:(UBi-1),:)
          t_new(:,LBj:(UBj-1)) = t_new(:,LBj:(UBj-1)) + v_up(:,(LBj+1):UBj)
          t_new(:,(LBj+1):UBj) = t_new(:,(LBj+1):UBj) + v_dn(:,LBj:(UBj-1))
           
          gT_u_up(LBi:(UBi-1),:)   = u_up((LBi+1):UBi,:)
          gT_u_dn((LBi+1):UBi,:)   = u_dn(LBi:(UBi-1),:)
          gT_v_up(:,LBj:(UBj-1))   = v_up(:,(LBj+1):UBj)
          gT_v_dn(:,(LBj+1):UBj)   = v_dn(:,LBj:(UBj-1))
                                         
          gCF_u_up(LBi:(UBi-1),:)  = u_up((LBi+1):UBi,:) * CFmat((LBi+1):UBi,:)
          gCF_u_dn((LBi+1):UBi,:)  = u_dn(LBi:(UBi-1),:) * CFmat(LBi:(UBi-1),:)
          gCF_v_up(:,LBj:(UBj-1))  = v_up(:,(LBj+1):UBj) * CFmat(:,(LBj+1):UBj)
          gCF_v_dn(:,(LBj+1):UBj)  = v_dn(:,LBj:(UBj-1)) * CFmat(:,LBj:(UBj-1))
          
          gCA_u_up(LBi:(UBi-1),:)  = u_up((LBi+1):UBi,:) * CAmat((LBi+1):UBi,:)
          gCA_u_dn((LBi+1):UBi,:)  = u_dn(LBi:(UBi-1),:) * CAmat(LBi:(UBi-1),:)
          gCA_v_up(:,LBj:(UBj-1))  = v_up(:,(LBj+1):UBj) * CAmat(:,(LBj+1):UBj)
          gCA_v_dn(:,(LBj+1):UBj)  = v_dn(:,LBj:(UBj-1)) * CAmat(:,LBj:(UBj-1))       
                 
          t_left = t_left + gT_u_up +  gT_u_dn +  gT_v_up +  gT_v_dn
          where (t_left .gt. 0)
            CFnew = (CF_left + gCF_u_up + gCF_u_dn + gCF_v_up + gCF_v_dn) / t_left
            CAnew = (CA_left + gCA_u_up + gCA_u_dn + gCA_v_up + gCA_v_dn) / t_left
          elsewhere
            CFnew = 1.0
            CAnew = CALbase_all(spn)
          end where
       
          !if ((sp .eq. 1).and.(((ac.eq.2).and.(lc.eq.7)).or.((ac.eq.11).and.(lc.eq.14)))) then
          ! do i=LBi,UBi
          !  do j=LBj,UBj   
          !   if (((i .eq. (182-1)).and.(j .eq. (165-1)))) then
          !     write(*,'(A5,6ES18.8)'),"CM", CFmat(i,j),CFmat(i+1,j),CFmat(i-1,j),CFmat(i,j+1),CFmat(i,j-1)
          !     write(*,'(A5,6ES18.8)'),"TT", t_left(i,j),gT_u_up(i,j),gT_u_dn(i,j),gT_v_up(i,j),gT_v_dn(i,j),t_new(i,j)
          !     write(*,'(A5,6ES18.8)'),"CF", CF_left(i,j),gCF_u_up(i,j),gCF_u_dn(i,j),gCF_v_up(i,j),gCF_v_dn(i,j),CFnew(i,j)  
          !   end if 
          !  end do
          ! end do         
          ! end if       
         
          CFmat = CFnew
          CAmat = CAnew

          ! diffusion part
          t_left = t_new              
          t_new = t_new - diff_xiup - diff_xidn - diff_etup - diff_etdn
      
          ! take out the positive losses only
          t_left = t_left - max(0.0,diff_xiup) - max(0.0,diff_xidn) - max(0.0,diff_etup) - max(0.0,diff_etdn)
          CF_left = CFmat * t_left 
          CA_left = CAmat * t_left
            
          ! now change these to represent the negative (gains)
          diff_xiup = abs(min(0.0,diff_xiup))
          diff_xidn = abs(min(0.0,diff_xidn))
          diff_etup = abs(min(0.0,diff_etup))
          diff_etdn = abs(min(0.0,diff_etdn))
       
          gCF_u_up(LBi:(UBi-1),:)  = diff_xiup(LBi:(UBi-1),:) * CFmat((LBi+1):UBi,:)
          gCF_u_dn((LBi+1):UBi,:)  = diff_xidn((LBi+1):UBi,:) * CFmat(LBi:(UBi-1),:)
          gCF_v_up(:,LBj:(UBj-1))  = diff_etup(:,LBj:(UBj-1)) * CFmat(:,(LBj+1):UBj)
          gCF_v_dn(:,(LBj+1):UBj)  = diff_etdn(:,(LBj+1):UBj) * CFmat(:,LBj:(UBj-1))
          
          gCA_u_up(LBi:(UBi-1),:)  = diff_xiup(LBi:(UBi-1),:) * CAmat((LBi+1):UBi,:)
          gCA_u_dn((LBi+1):UBi,:)  = diff_xidn((LBi+1):UBi,:) * CAmat(LBi:(UBi-1),:)
          gCA_v_up(:,LBj:(UBj-1))  = diff_etup(:,LBj:(UBj-1)) * CAmat(:,(LBj+1):UBj)
          gCA_v_dn(:,(LBj+1):UBj)  = diff_etdn(:,(LBj+1):UBj) * CAmat(:,LBj:(UBj-1))
       
          t_left = t_left + diff_xiup + diff_xidn + diff_etup +  diff_etdn              
          where (t_left .gt. 0)
            CFnew = (CF_left + gCF_u_up + gCF_u_dn + gCF_v_up + gCF_v_dn) / t_left
            CAnew = (CA_left + gCA_u_up + gCA_u_dn + gCA_v_up + gCA_v_dn) / t_left
          elsewhere
            CFnew = 1.0
            CAnew = CALbase_all(spn)
          end where

          !if ((sp .eq. 1).and.(((ac.eq.2).and.(lc.eq.7)).or.((ac.eq.11).and.(lc.eq.14)))) then
          ! do i=LBi,UBi
          !  do j=LBj,UBj   
          !   if (((i .eq. (182-1)).and.(j .eq. (165-1)))) then
          !     write(*,'(A5,6ES18.8)'),"CMD", CFmat(i,j),CFmat(i+1,j),CFmat(i-1,j),CFmat(i,j+1),CFmat(i,j-1)
          !     write(*,'(A5,6ES18.8)'),"TTD", t_left(i,j),diff_xiup(i,j),diff_xidn(i,j),diff_etup(i,j),diff_etdn(i,j),t_new(i,j)
          !     write(*,'(A5,6ES18.8)'),"CFD", CF_left(i,j),gCF_u_up(i,j),gCF_u_dn(i,j),gCF_v_up(i,j),gCF_v_dn(i,j),CFnew(i,j)  
          !   end if 
          !  end do
          ! end do         
          ! end if   

          do i=LBi,UBi
            do j=LBj,UBj 
              if ((i.le.1).or.(i.ge.181).or.(j.le.1).or.(j.ge.257)) then
                t_new(i,j)=0.0
                CFnew(i,j)=1.0
                CAnew(i,j)=CALbase_all(spn)
              end if
              if (isnan(t_new(i,j)).or.isnan(CFnew(i,j)).or.isnan(CAnew(i,j))) then
                !write (*,"(A6,5i4,3ES14.4)"),"FNAN",i,j, &
                !& sp_all(spn),age_all(spn),lc_all(spn),  &
                !& t_new(i,j),CFnew(i,j),CAnew(i,j)
                t_new(i,j)=0.0
                CFnew(i,j)=1.0
                CAnew(i,j)=CALbase_all(spn)
              end if  
            end do
          end do
     
     
                 
          SELECT CASE (type_all(spn))
            CASE (1)   ! aged
              t(:,:,kspal(1,sp,lc,ac),nnew,tspal(1,sp,lc,ac)) = min(max(t_new,0.00),1e35) 
              t(:,:,kspal(2,sp,lc,ac),nnew,tspal(2,sp,lc,ac)) = min(max(CFnew,0.01),1e35)          
              t(:,:,kspal(3,sp,lc,ac),nnew,tspal(3,sp,lc,ac)) = min(max(CAnew,1.00),1e35)
            CASE (2)   ! lengthed
              t(:,:,kspl(1,sp-NUM_AGED_SPECIES,lc),nnew,tspl(1,sp-NUM_AGED_SPECIES,lc)) = min(max(t_new,0.00),1e35)
              t(:,:,kspl(2,sp-NUM_AGED_SPECIES,lc),nnew,tspl(2,sp-NUM_AGED_SPECIES,lc)) = min(max(CFnew,0.01),1e35)          
              t(:,:,kspl(3,sp-NUM_AGED_SPECIES,lc),nnew,tspl(3,sp-NUM_AGED_SPECIES,lc)) = min(max(CAnew,1.00),1e35)          
            CASE (3)   ! simple       
          END SELECT
           
        end do ! END OF MOVING ALL_LINKS N, CF, and CAL ------------------------------    

        ! Move some eggs
        do sp=1,TOT_LENGTHED
          !if (sp == 1) then
          !   print *,sp,sum(t_new)
          !end if

          !t_new = GF % eggs(sp,:,:)
          t_new = t(:,:,keggs(sp),nnew,teggs(sp))
          base_speed = 1.0 * (1./10000.) * dt(ng)     
          diff_coeff = 0.4    ! started at 0.1 KYA 2/18/2012 
          ! diffusion     
          diff_xiup = 0.0 
          diff_xidn = 0.0
          diff_etup = 0.0 
          diff_etdn = 0.0    
       
          diff_xiup(LBi:(UBi-1),:) = (t_new(LBi:(UBi-1),:) - t_new((LBi+1):UBi,:))
          diff_xidn((LBi+1):UBi,:) = (t_new((LBi+1):UBi,:) - t_new(LBi:(UBi-1),:))
          diff_etup(:,LBj:(UBj-1)) = (t_new(:,LBj:(UBj-1)) - t_new(:,(LBj+1):UBj))
          diff_etdn(:,(LBj+1):UBj) = (t_new(:,(LBj+1):UBj) - t_new(:,LBj:(UBj-1)))     
          
          diff_xiup(LBi:(UBi-1),:) = diff_xiup(LBi:(UBi-1),:) * rmask((LBi+1):UBi,:)
          diff_xidn((LBi+1):UBi,:) = diff_xidn((LBi+1):UBi,:) * rmask(LBi:(UBi-1),:)
          diff_etup(:,LBj:(UBj-1)) = diff_etup(:,LBj:(UBj-1)) * rmask(:,(LBj+1):UBj)
          diff_etdn(:,(LBj+1):UBj) = diff_etdn(:,(LBj+1):UBj) * rmask(:,LBj:(UBj-1))
          
          diff_xiup = diff_xiup * diff_coeff * base_speed * rmask
          diff_xidn = diff_xidn * diff_coeff * base_speed * rmask
          diff_etup = diff_etup * diff_coeff * base_speed * rmask
          diff_etdn = diff_etdn * diff_coeff * base_speed * rmask        
          
          ! directed
          u_move = 0.0
          v_move = 0.0
          u_up   = 0.0
          u_dn   = 0.0              
          v_up   = 0.0
          v_dn   = 0.0
          
          u_move    = u(:,:,UBk,nstp) * base_speed * t_new * rmask
          v_move    = v(:,:,UBk,nstp) * base_speed * t_new * rmask
          
          u_up((LBi+1):UBi,:) = max(0.0, -u_move((LBi+1):UBi,:)) * rmask(LBi:(UBi-1),:)  
          u_dn(LBi:(UBi-1),:) = max(0.0, +u_move(LBi:(UBi-1),:)) * rmask((LBi+1):UBi,:)
          v_up(:,(LBj+1):UBj) = max(0.0, -v_move(:,(LBj+1):UBj)) * rmask(:,LBj:(UBj-1))
          v_dn(:,LBj:(UBj-1)) = max(0.0, +v_move(:,LBj:(UBj-1))) * rmask(:,(LBj+1):UBj)   

          !if (sp == 1) then
          !if ((LBi<100).and.(UBi>100)) then
          !  if ((LBj<100).and.(UBj>100)) then
          !    do i = LBi,UBi
          !       do j = LBj,UBj
          !         ! print something                     
          !       end do
          !    end do
          !  end if
          !end if                         
          !end if
      
          ! bringing together EGGS
          t_new = t_new - u_up - u_dn - v_up - v_dn 
          t_new(LBi:(UBi-1),:) = t_new(LBi:(UBi-1),:) + u_up((LBi+1):UBi,:)
          t_new((LBi+1):UBi,:) = t_new((LBi+1):UBi,:) + u_dn(LBi:(UBi-1),:)
          t_new(:,LBj:(UBj-1)) = t_new(:,LBj:(UBj-1)) + v_up(:,(LBj+1):UBj)
          t_new(:,(LBj+1):UBj) = t_new(:,(LBj+1):UBj) + v_dn(:,LBj:(UBj-1))
                   
          t_new = t_new - diff_xiup - diff_xidn - diff_etup - diff_etdn
                                                                                                                                         
          ! FINAL EGGS
          t(:,:,keggs(sp),nnew,teggs(sp)) = min(max(t_new,0.0),1e35)
     
        end do

      end if  ! END FEAST MOVEMENT CONDITIONAL
 

              