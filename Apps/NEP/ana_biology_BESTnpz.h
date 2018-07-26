!
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions for the BEST NPZ model .                                       !
!                                                                      !
!=======================================================================
!



#ifdef BARRENSTART

!
!-----------------------------------------------------------------------
! This set of initial conditions starts with an empty Bering Sea shelf; 
! NO3 is only present at depth, with no NH4 or detritus anywhere.  
! Living critters are seeded with only a tiny amount to allow for future 
! growth.  This start condition is intended to allow the model domain to 
! move towards its own internally-regulated nutrient/benthos steady-
! state
!-----------------------------------------------------------------------
!
      

#else
!
!-----------------------------------------------------------------------
! This set of initial conditions reflects those used in the original 
! Bering10K runs.  
!-----------------------------------------------------------------------
!
      
      DO i=IstrR,IendR
        DO j=JstrR,JendR
          
# ifdef IRON_LIMIT
          ! Iron top/bottom values set based on bottom (matches nudging)
          FeSurf = CalcLinearCapped(Feinh, Feinlo, Feoffh, Feofflo, GRID(ng)%h(i,j))
          FeDeep = CalcLinearCapped(Feinh, Feinhi, Feoffh, Feoffhi, GRID(ng)%h(i,j))
# endif
          
          DO k=1,N(ng)
# ifdef IRON_LIMIT

            ! Iron: depth-dependant linear-capped profile (same as nudging)
            
            t(i,j,k,1,iFe) = CalcLinearCapped(-50.0_r8, FeSurf, -300.0_r8, FeDeep, GRID(ng)%z_r(i,j,k))
# endif
            ! Nitrate: depth-dependant linear-capped profile
            
            t(i,j,k,1,iNO3) = CalcLinearCapped(-100.0_r8, 18.0_r8, -300.0_r8, 30.0_r8, GRID(ng)%z_r(i,j,k))
            
            ! NH4: 0
            
            t(i,j,k,1,iNH4) = 0.0_r8
            
            ! Phytoplankton and small zooplankton: constant profiles
            
            t(i,j,k,1,iPhS) = 5.0_r8
            t(i,j,k,1,iPhL) = 1.0_r8
            t(i,j,k,1,iMZS) = 0.0_r8
            t(i,j,k,1,iMZL) = 1.0_r8
            t(i,j,k,1,iCop) = 0.1_r8
            
            ! Large copepods: Offshore Neocalanus are only added to grid cells deeper than 
            ! 400m, with a concentration at mid-depth under the assumption that we're 
            ! starting in January and they would be in diapause at that time.  Onshore 
            ! Calanus marshallae are added only to grid cells shallower than 200m, and 
            ! placed in the top 5 layers (regardless of model resolution).

	     		  if (z_r(i,j,1) .le. -400.0_r8) then 
	            if (z_r(i,j,k).le.-400.and.z_r(i,j,k).ge.-800) THEN
                t(i,j,k,1,iNCaO) = 0.1_r8
	            endif
	     		  elseif (z_r(i,j,1) .ge. -200.0_r8) then 
  	          if(k.le.5)THEN
                t(i,j,k,1,iNCaS) =  8.0_r8
	            endif
	          endif

            ! Euphausiids: Constant concentration, split into onshore (<100m) and offshore 
            ! (<100m) groups.
            
	          if (-z_r(i,j,1) .gt. 100.0_r8 ) then 
              t(i,j,k,1,iEupO) =  0.1_r8
	            t(i,j,k,1,iEupS) =  0_r8
	          else
              t(i,j,k,1,iEupS) =  0.1_r8
	            t(i,j,k,1,iEupO) =  0_r8
	          endif
            
# ifdef JELLY
            ! Jellyfish: depth-dependant linear-capped profile
            
            t(i,j,k,1,iJel) = CalcLinearCapped(-100.0_r8, 0.1_r8, -300.0_r8, eps, GRID(ng)%z_r(i,j,k))
# endif

            ! Detritus: nearly nothing
            
            t(i,j,k,1,iDet) =  eps
	          t(i,j,k,1,iDetF) =  eps
            

          END DO
# ifdef BENTHIC
          DO k = 1,NBL(ng)
            
            ! Benthos: start high (note: due to a missing check for analytical conditions, 
            ! these were overrided by input file initial conditions in most earlier 
            ! Bering10K runs)
            
            bt(i,j,k,1,iBen) = 8000.0_r8 
            
            ! Benthic detritus
            
		        bt(i,j,k,1,iDetBen) = 500.0_r8
            
          END DO
# endif        
# ifdef ICE_BIO

          ! Ice: start with no-ice conditions
          
#  ifdef CLIM_ICE_1D
          it(i,j,1,iIcePhL)  =  0.0_r8      
          it(i,j,1,iIceNO3)  =  0.0_r8      
          it(i,j,1,iIceNH4)  =  0.0_r8      
		      itL(i,j,1,iIceLog) = -1.0_r8    
#  elif defined BERING_10K
          IcePhL(i,j,1) =  0.0_r8
          IceNO3(i,j,1) =  0.0_r8
          IceNH4(i,j,1) =  0.0_r8
          IceLog(i,j,1) = -1.0_r8
#  endif
# endif

        END DO
      END DO

#endif

!
!-----------------------------------------------------------------------
! Additional setup, applicable to all compilation options above
!-----------------------------------------------------------------------
!

      ! Set other time index, and add a minimum eps value for all 3D variables

      DO i=IstrR,IendR
        DO j=JstrR,JendR
          DO k=1,N(ng)
            DO is=1,NBT
              itrc=idbio(is)
              t(i,j,k,1,itrc) = MAX(t(i,j,k,1,itrc),eps)
              t(i,j,k,2,itrc) = t(i,j,k,1,itrc)
            END DO
          END DO
#ifdef BENTHIC
          DO k=1,NBL(ng)
		        bt(i,j,k,2,iBen)    = bt(i,j,k,1,iBen)
		        bt(i,j,k,2,iDetBen) = bt(i,j,k,1,iDetBen)
          END DO
#endif
#ifdef ICE_BIO
# ifdef CLIM_ICE_1D
          it(i,j,2,iIcePhL)  = it(i,j,1,iIcePhL)    
          it(i,j,2,iIceNO3)  = it(i,j,1,iIceNO3)    
          it(i,j,2,iIceNH4)  = it(i,j,1,iIceNH4)    
		      itL(i,j,2,iIceLog) = itL(i,j,1,iIceLog) 
# elif defined BERING_10K
          IcePhL(i,j,2) = IcePhL(i,j,1)
          IceNO3(i,j,2) = IceNO3(i,j,1)
          IceNH4(i,j,2) = IceNH4(i,j,1)
          IceLog(i,j,2) = IceLog(i,j,1)
# endif
#endif
        END DO
      END DO
      
#ifdef STATIONARY
      ! Initialize all diagnostics to 0
      
      st = 0.0_r8
#endif
            

!
!
!
!
!
!
!
!
!
!
!
! #undef DEPAVG      /* TEST case: set all equal to depth averaged values */
!
! #ifdef DEPAVG
! !
! ! Set all points equal to the same number
! !
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             t(i,j,k,1,iNO3) = 18.0_r8
!             t(i,j,k,1,iNH4) = 0.20_r8
!             t(i,j,k,1,iPhS) = 0.10_r8
!             t(i,j,k,1,iPhL) = 0.10_r8
!             t(i,j,k,1,iMZS) = 0.00_r8
!             t(i,j,k,1,iMZL) = 0.10_r8
!             t(i,j,k,1,iCop) = 2.041_r8
!             t(i,j,k,1,iNCa) = 2.417_r8
!             t(i,j,k,1,iEup) = 1.337_r8
!             t(i,j,k,1,iDet) = eps
!             t(i,j,k,1,iDetF) = eps
! # ifdef IRON_LIMIT
!             t(i,j,k,1,iFe) = 0.6_r8
! # endif
!           enddo
!         enddo
!       enddo
! #else
! !
! ! Make the curves.  First come parameters that are constant through
! ! entire water column.
! !
!       deepval(iNO3)  = 30.0_r8
!       deepval(iNH4)  = eps
!       deepval(iPhS)  = eps
!       deepval(iPhL)  = eps
!       deepval(iMZS)  = 0.0_r8
!       deepval(iMZL)  = 0.0_r8
!       deepval(iCop)  = 0.0_r8
!       deepval(iNCaS) = 0.0_r8
!       deepval(iEupS) = 0.0_r8
!       deepval(iNCaO) = 0.0_r8
!       deepval(iEupO) = 0.0_r8
!       deepval(iDet)  = eps
!       deepval(iDetF) = eps
!       deepval(iJel)  = eps
!
!       do k=1,N(ng)
!         do j=JstrR,JendR
!           do i=IstrR,IendR
!             t(i,j,k,1,iNO3) = 18.0_r8 !  18.0_r8
!             t(i,j,k,1,iNH4) =  0.0_r8 !2.0_r8
!             t(i,j,k,1,iDet) =  eps
!             t(i,j,k,1,iDetF) =  eps
!             biod(i,j,k) = -1 * ( z_r(i,j,k) + 2.5_r8 )
! # ifdef JELLY
!             t(i,j,k,1,iJel) =  0.1_r8
! # endif
!             enddo
!          enddo
!       enddo
! !
! ! PS - a combination of 2 curves: 2nd order polynomial curve to create
! ! a subsurface maximum at ~20 m and an exponentially decreasing curve to
! ! assure values aren't negative at depth.
! !
!       var1 = 25.58_r8
!       var2 = -0.250_r8 / (5.0_r8**2)
!       var3 =  0.008_r8 / (5.0_r8**3)
!       var4 = 3.82_r8 * 5.0_r8
!       var5 = 75.0_r8
!       var6 = var5 - var4
!       var7 = var1 + var2*(var6**2) + var3*(var6**3)
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             if ( biod(i,j,k) .le. var5 ) then
!               var6 = biod(i,j,k) - var4
!               t(i,j,k,1,iPhS) = var1 + var2*(var6**2) + var3*(var6**3)
!             else
!               t(i,j,k,1,iPhS) =                                         &
!      &          var7 * exp( ( -1.0_r8*biod(i,j,k) + var5 ) / 5.0_r8 )
!             endif
!           enddo
!         enddo
!       enddo
! !
! ! PL - exponentially decreasing with depth.
! !
!       var1 = 8.25_r8
!       var2 = 0.322_r8 / 5.0_r8
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             t(i,j,k,1,iPhL) = var1 * exp( -1.0_r8 * var2 * biod(i,j,k) )
!           enddo
!         enddo
!       enddo
! !
! ! Microzooplankton, from Howell-Kubler, 1996
! ! approximated with a straight line and an exponentially decreasing
! ! curve to assure values aren't negative at depth.  Curves meet ~60m
! !
!       var1 = 3.1714_r8
!       var2 = -0.1865_r8 / 5.0_r8
!       var3 = 60.0_r8
!       var4 = 0.5_r8
!       var5 = var1 + var2 * var3
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             if ( biod(i,j,k) .le. var3 ) then
!               t(i,j,k,1,iMZL) = var1 + var2 * biod(i,j,k)
!             else
!               t(i,j,k,1,iMZL) = var4 +                                  &
!      &          ( var5 - var4 ) * exp( ( var3 - biod(i,j,k) ) / 5.0_r8 )
!             endif
! !           t(i,j,k,1,iMZS) = t(i,j,k,1,iMZL)
!             t(i,j,k,1,iMZS) = 0.0_r8
!           enddo
!         enddo
!       enddo
! !
! ! Pseudocalanus Copepods - from Shelikof data via Shelikof NPZ
! ! Step at ~ 50 m
! !
!       var1 = pi / 2.0_r8
!       var2 = 2.876_r8 / pi
!       var3 = 5.0_r8 / 5.0_r8
!       var4 = 9.151_r8 * 5.0_r8
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             var5 = var3 * ( biod(i,j,k) - var4 )
!             t(i,j,k,1,iCop) = var1 - var2 * atan( var5 )
!           enddo
!         enddo
!       enddo
!
!
!
! !
! ! Neocalanus, from Shelikof NPZ. Step at ~ 30m
! !
!
!       var1 = 4.0_r8
!       var2 = 1.3_r8
!       var3 = 5.0_r8 / 5.0_r8
!       var4 = 5.2_r8 * 5.0_r8
!
! !On and off shelf neocalanus - at depth for january
!
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             var5 = var3 * ( biod(i,j,k) - var4 )
!              if (z_r(i,j,1) .le. -400.0_r8) then
!               if(z_r(i,j,k).le.-400.and.z_r(i,j,k).ge.-800)THEN
!                 t(i,j,k,1,iNCaO) = 0.1_r8
!               endif
!              elseif (z_r(i,j,1) .ge. -200.0_r8) then
!               if(k.le.5)THEN
!                 t(i,j,k,1,iNCaS) =  8.0_r8
!               endif
!             endif
!           enddo
!         enddo
!       enddo
! !
! !
! ! Euphausiids, Wild guesses. Step at ~ 30m
! !
!       var1 = 1.78_r8
!       var2 = 0.8_r8
!       var3 = 5.0_r8 / 5.0_r8
!       var4 = 5.2_r8 * 5.0_r8
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             var5 = var3 * ( biod(i,j,k) - var4)
!             if (-z_r(i,j,1) .gt. 100.0_r8 ) then
! !             t(i,j,k,1,iEupO) =  var1 - var2 * atan( var5 )
!               t(i,j,k,1,iEupO) =  0.1_r8
!               t(i,j,k,1,iEupS) =  0_r8
!             else
! !              t(i,j,k,1,iEupS) =  var1 - var2 * atan( var5 )
!               t(i,j,k,1,iEupS) =  0.1_r8
!               t(i,j,k,1,iEupO) =  0_r8
!             endif
!           enddo
!         enddo
!       enddo
!
!
! # ifdef IRON_LIMIT
! ! Iron - linear from surface value to value at 100m and increase onshore
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           var1 = MAX(0._r8,MIN(1._r8,                                   &
!      &             (GRID(ng)%h(i,j)-Feinh)/(Feoffh-Feinh)))
!           FeSurf = Feinlo + var1*(Feofflo-Feinlo)
!           FeDeep = Feinhi + var1*(Feoffhi-Feinhi)
!           var1 = (FeDeep-FeSurf) / 200._r8
!           do k=1,N(ng)
!             t(i,j,k,1,iFe) = MIN(FeDeep, FeSurf - z_r(i,j,k)*var1)
!           enddo
!         enddo
!       enddo
! # endif
! !
! ! Concentrations of everything below 100m - i.e. below
! ! depths where calculations are performed.  Have linear slope
! ! between values above and below.
! ! Iron deep values have already been determined.
! !
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=N(ng),1,-1
!             if ( biod(i,j,k) .gt. 300.0_r8 ) then       !was 120 for GOA
!               t(i,j,k,1,iNO3) = deepval(iNO3)
!               t(i,j,k,1,iNH4) = deepval(iNH4)
!               t(i,j,k,1,iPhS) = deepval(iPhS)
!               t(i,j,k,1,iPhL) = deepval(iPhL)
!               t(i,j,k,1,iMZS) = deepval(iMZS)
!               t(i,j,k,1,iMZL) = deepval(iMZL)
!               t(i,j,k,1,iCop) = deepval(iCop)
! !             t(i,j,k,1,iNCaS) = deepval(iNCaS)
!               t(i,j,k,1,iEupS) = deepval(iEupS)
! !              t(i,j,k,1,iNCaO) = deepval(iNCaO)
!               t(i,j,k,1,iEupO) = deepval(iEupO)
!               t(i,j,k,1,iDet) = deepval(iDet)
!               t(i,j,k,1,iDetF) = deepval(iDetF)
! # ifdef JELLY
!               t(i,j,k,1,iJel) = deepval(iJel)
! # endif
!             else if ( biod(i,j,k) .gt. 100.0_r8 .and.                   &
!      &                biod(i,j,k) .le. 300.0_r8) then
!               var1 = ( 100.0_r8 - biod(i,j,k) ) / ( 100.0_r8-300.0_r8 )
!               t(i,j,k,1,iNO3) = loval(iNO3) +                           &
!      &                          ( deepval(iNO3) - loval(iNO3) ) * var1
!               t(i,j,k,1,iNH4) = loval(iNH4) +                           &
!      &                          ( deepval(iNH4) - loval(iNH4) ) * var1
!               t(i,j,k,1,iPhS) = loval(iPhS) +                           &
!      &                          ( deepval(iPhS) - loval(iPhS) ) * var1
!               t(i,j,k,1,iPhL) = loval(iPhL) +                           &
!      &                          ( deepval(iPhL) - loval(iPhL) ) * var1
!               t(i,j,k,1,iMZS) = loval(iMZS) +                           &
!      &                          ( deepval(iMZS) - loval(iMZS) ) * var1
!               t(i,j,k,1,iMZL) = loval(iMZL) +                           &
!      &                          ( deepval(iMZL) - loval(iMZL) ) * var1
!               t(i,j,k,1,iCop) = loval(iCop) +                           &
!      &                          ( deepval(iCop) - loval(iCop) ) * var1
! !             t(i,j,k,1,iNCaS) = loval(iNCaS) +                         &
! !     &                          ( deepval(iNCaS) - loval(iNCaS) ) * var1
!               t(i,j,k,1,iEupS) = loval(iEupS) +                         &
!      &                          ( deepval(iEupS) - loval(iEupS) ) * var1
! !             t(i,j,k,1,iNCaO) = loval(iNCaO) +                         &
! !     &                          ( deepval(iNCaO) - loval(iNCaO) ) * var1
!               t(i,j,k,1,iEupO) = loval(iEupO) +                         &
!      &                          ( deepval(iEupO) - loval(iEupO) ) * var1
!               t(i,j,k,1,iDet) = loval(iDet) +                           &
!      &                          ( deepval(iDet) - loval(iDet) ) * var1
!               t(i,j,k,1,iDetF) = loval(iDetF) +                         &
!      &                          ( deepval(iDetF) - loval(iDetF) ) * var1
! # ifdef JELLY
!               t(i,j,k,1,iJel) = loval(iJel) +                           &
!                                 ( deepval(iJel) - loval(iJel) ) * var1
! # endif
!             else
!               loval(iNO3) = t(i,j,k,1,iNO3)
!               loval(iNH4) = t(i,j,k,1,iNH4)
!               loval(iPhS) = t(i,j,k,1,iPhS)
!               loval(iPhL) = t(i,j,k,1,iPhL)
!               loval(iMZS) = t(i,j,k,1,iMZS)
!               loval(iMZL) = t(i,j,k,1,iMZL)
!               loval(iCop) = t(i,j,k,1,iCop)
! !             loval(iNCaS) = t(i,j,k,1,iNCaS)
!               loval(iEupS) = t(i,j,k,1,iEupS)
! !              loval(iNCaO) = t(i,j,k,1,iNCaO)
!               loval(iEupO) = t(i,j,k,1,iEupO)
!               loval(iDet) = t(i,j,k,1,iDet)
!               loval(iDetF) = t(i,j,k,1,iDetF)
!               loval(iJel) = t(i,j,k,1,iJel)
!             endif
!           enddo
!         enddo
!       enddo
! #endif /* DEPAVG */
! #ifdef GAK1D
! !
! !  This is a hack for sensitivity studies - to test warmer or cooler water
! !
! !      do i=IstrR,IendR
! !        do j=JstrR,JendR
! !          do k=1,N(ng)
! !            t(i,j,k,1,itemp) = t(i,j,k,1,itemp) + 2.0_r8
! !          enddo
! !        enddo
! !      enddo
! #endif
!
!
!
!
! #ifdef ICE_BIO
! # ifdef CLIM_ICE_1D
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           it(i,j,1,iIcePhL) =0.0_r8       !was eps
!           it(i,j,1,iIceNO3) =0.0_r8        !was eps
!           it(i,j,1,iIceNH4) =0.0_r8       !was eps
!           itL(i,j,1,iIceLog) =-1.0_r8       !was eps
!         enddo
!       enddo
! # elif defined BERING_10K
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           IcePhL(i,j,1) =0.0_r8       !was eps
!           IceNO3(i,j,1) =0.0_r8        !was eps
!           IceNH4(i,j,1) =0.0_r8       !was eps
!           IceLog(i,j,1) =-1.0_r8       !was eps
!           IcePhL(i,j,2) =0.0_r8       !was eps
!           IceNO3(i,j,2) =0.0_r8        !was eps
!           IceNH4(i,j,2) =0.0_r8       !was eps
!           IceLog(i,j,2) =-1.0_r8       !was eps
!         enddo
!       enddo
! # endif
! #endif
!
!
! !reduce the initial concs by 2 orders of magnitude
! !test for Bering Sea - initialize in Jan
!
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             t(i,j,k,1,iPhS) = 5_r8 !t(i,j,k,1,iPhS)*.01_r8
!             t(i,j,k,1,iPhL) = 1_r8 !t(i,j,k,1,iPhL)*.01_r8
!             t(i,j,k,1,iMZS) = 0.0_r8
!             t(i,j,k,1,iMZL) = 1_r8
!             t(i,j,k,1,iCop) = 0.1_r8
! !              t(i,j,k,1,iNCaS) = t(i,j,k,1,iNCaS)*.0001_r8
! !              t(i,j,k,1,iEupS) = t(i,j,k,1,iEupS)*.0001_r8
! !        t(i,j,k,1,iNCaO) = t(i,j,k,1,iNCaO)*.0001_r8
! !              t(i,j,k,1,iEupO) = t(i,j,k,1,iEupO)*.0001_r8
!           enddo
!         enddo
!       enddo
!
! ! Check for size, set other time index, and periodic BC's
! !
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             DO is=1,NBT
!               itrc=idbio(is)
!               t(i,j,k,1,itrc) = MAX(t(i,j,k,1,itrc),eps)
!               t(i,j,k,2,itrc) = t(i,j,k,1,itrc)
!             enddo
!           enddo
!         enddo
!       enddo
!
! #ifdef BENTHIC
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,NBL(ng)
!             bt(i,j,k,1,iBen) = 8000.0_r8    !was eps
!             bt(i,j,k,1,iDetBen) = 500.0_r8! 00.0_r8    !was eps  500
!             bt(i,j,k,2,iBen) = bt(i,j,k,1,iBen)
!             bt(i,j,k,2,iDetBen) = bt(i,j,k,1,iDetBen)
!           enddo
!         enddo
!       enddo
! #endif
! !g#ifdef BENTHIC
! !g         do i=IstrR,IendR
! !g           do j=JstrR,JendR
! !g             do k=1,NBL(ng)
! !g                 bt(i,j,k,1,iBen) = MAX(bt(i,j,k,1,iBen),eps)
! !g     bt(i,j,k,1,iDetBen) = MAX(bt(i,j,k,1,iDetBen),eps)
! !g     bt(i,j,k,2,iBen) = bt(i,j,k,1,iBen)
! !g     bt(i,j,k,2,iDetBen) = bt(i,j,k,1,iDetBen)
! !g      enddo
! !g     enddo
! !g            enddo
! !g#endif
!
! #ifdef STATIONARY
!       do i=IstrR,IendR
!         do j=JstrR,JendR
!           do k=1,N(ng)
!             st(i,j,k,1,1) = 0.0_r8
!             st(i,j,k,1,2) = 0.0_r8
!             st(i,j,k,1,3) = 0.0_r8
!             st(i,j,k,1,4) = 0.0_r8
!             st(i,j,k,1,5) = 0.0_r8
!             st(i,j,k,1,6) = 0.0_r8
!             st(i,j,k,1,7) = 0.0_r8
!             st(i,j,k,1,8) = 0.0_r8
!             st(i,j,k,2,1) = 0.0_r8
!             st(i,j,k,2,2) = 0.0_r8
!             st(i,j,k,2,3) = 0.0_r8
!             st(i,j,k,2,4) = 0.0_r8
!             st(i,j,k,2,5) = 0.0_r8
!             st(i,j,k,2,6) = 0.0_r8
!             st(i,j,k,2,7) = 0.0_r8
!             st(i,j,k,2,8) = 0.0_r8
!
!
!             st(i,j,k,3,1) = 0.0_r8
!             st(i,j,k,3,2) = 0.0_r8
!             st(i,j,k,3,3) = 0.0_r8
!             st(i,j,k,3,4) = 0.0_r8
!             st(i,j,k,3,5) = 0.0_r8
!             st(i,j,k,3,6) = 0.0_r8
!             st(i,j,k,3,7) = 0.0_r8
!             st(i,j,k,3,8) = 0.0_r8
!           enddo
!         enddo
!       enddo
! #endif
!
