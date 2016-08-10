      MODULE biology_mod
!
!svn $Id: biology.F 1076 2009-09-25 23:18:43Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the source and sink terms for selected        !
!   biology model.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: biology
      CONTAINS
      SUBROUTINE biology (ng,tile)
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Georgina Gibsons BESTNPZ Code August 2016  
!
!  Originally modified from Sarah Hinckleys GOANPZ code which had been
!   implemented by Craig Lewis (CVL) and modified by Liz Dobbins and Sarah Hinckley   
!                       !
!                                                                        !
!=======================================================================
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_mixing
      integer, intent(in) :: ng, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private storage
!  arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J)- and MAX(I,J)-directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
!  Set header file name.
! 
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)="ROMS/Nonlinear/bestnpz.h"
      END IF
      CALL wclock_on (ng, iNLM, 15)
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   N(ng), NT(ng),                                 &
     &                   nnew(ng), nstp(ng),                            &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % t                                  &
     &                   ,OCEAN(ng) % bt                                &
     &                  ,MIXING(ng) % Akt                               &
     &                              )     
      CALL wclock_off (ng, iNLM, 15)
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         UBk, UBt,                                &
     &                         nnew, nstp,                              &
     &                         Hz, z_r, z_w, srflx, t                   & 
     &                            ,bt                                   &
     &                          ,Akt                                    & 
     &                                  )    
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_scalars
      USE mod_ocean
      USE mod_grid
      USE exchange_3d_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: UBk, UBt
      integer, intent(in) :: nnew, nstp
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bt(LBi:,LBj:,:,:,:)
      real(r8) :: predSumCop, predSumNCaS, predSumEupS, predSumNCaO, predSumEupO
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!   
      integer :: i, j, k, ibio, ibio2,itr, itrmx, itrc, itrc2
      real(r8) :: cff5,cff6,cff6b,cff7,cff8,cff9,cff10,cff11
      integer :: ibioB     
      real(r8) :: bf,fbase,TSS,Ifs,atss,btss,SF
      real(r8) ::avgD,avgDF,avgPS,avgPL,dw,wcPS,wcPL,wcD,wcDF,PSsum  
      real(r8) ::sumD,sumDF,sumPL
      integer :: Iter,is
      integer :: iday, month, year
      real(r8) :: cff0,cff1, cff1b,cff2,cff3,cff4,pmaxs,dz
      real(r8) :: TFMZS,TFMZL,TFCop,TFNCa,TFEup,TFJel
      real(r8) :: Drate, Pmax,NOup, NHup,offset
      real(r8) :: dtdays,Ra,Rf
      real(r8) :: LightLim,NOLim,NHLim,IronLim
      real(r8) :: hour,yday,lat,k_phy,Dl,Par1,k_extV,k_chlV
      real(r8) :: Sal1,Temp1
      real(r8) :: ParMax,BasalMetMZL, BasalMetCop
      real(r8) :: BasalMetCM,BasalMetNC,BasalMetEup     
      real(r8) :: Iron1,kfePh,respPh,BasalMet,BasalMetJel
      real(r8) :: PON,Dep1,Nitrif,NH4R
      real(r8) :: NitrifMax,DLNitrif
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: DBio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: BioB
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: DBioB
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: Bio_bakB
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS) :: PARs
      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dens
      real(r8), dimension(N(ng)) :: DensV
      real(r8), dimension(N(ng)) :: ZW_V
      real(r8), dimension(IminS:ImaxS) :: StabParam
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TestVal
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncCop
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncNeo
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncEup
      real(r8),dimension(IminS:ImaxS,N(ng))::TempFuncJel
      real(r8), dimension(IminS:ImaxS,N(ng)) :: HzL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: z_wL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: sinkIN,sinkOUT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: riseIN,riseOUT      
!
      logical :: downward = .false., upward = .false.
      logical :: downwardCM = .false., upwardCM = .false.
      real(r8), parameter :: eps  = 1.0E-20_r8
      real(r8), parameter :: minv = 0.0E-20_r8
      real(r8) :: Alpha
      real(r8) :: ALPHA_N,ALPHA_P, kN, kP, alphaPhSv, alphaPhLv
      real(r8) ::respNC, respCM 
			! Vertical movement
			real(r8), dimension(N(ng)) :: Btmp, Hztmp
			real(r8), dimension(0:N(ng)) :: zwtmp
			real(r8) :: sinkout2   
			real(r8) :: RSNC, RENC, SSNC, SENC, RSCM, RECM, SSCM, SECM
!----------------------
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
      dtdays = dt(ng)*sec2day/REAL(BioIter(ng),r8)
      k_phy = k_chl / ccr
		 ! For diapause, set movement direction flags for large copepods, 
	 	 ! and lower respiration rates if they're in the diapause phase.
	 	 ! NCaS = CM = mostly C. melanaster, on-shelf
	 	 ! NCaO = NC = mostly Neocalanus, off-shelf
	 	 RSNC = MOD(RiseStart, 366.0_r8)
	 	 RENC = MOD(RiseEnd,   366.0_r8)
	 	 SSNC = MOD(SinkStart, 366.0_r8)
	 	 SENC = MOD(SinkEnd,   366.0_r8)
	 	 RSCM = MOD(RiseStart + 30, 366.0_r8)
	 	 RECM = MOD(RiseEnd   + 30, 366.0_r8)
	 	 SSCM = MOD(SinkStart + 30, 366.0_r8)
	 	 SECM = MOD(SinkEnd   + 30, 366.0_r8)
	 	 upward =     ((RSNC.lt.RENC) .and.                                  &
	 	&              (yday.ge.RSNC .and. yday.le.RENC))                    &
	 	&             .or.                                                   &
	 	&             ((RSNC.gt.RENC) .and.                                  &
	 	&              (yday.ge.RSNC .or.  yday.le.RENC))
	 	 upwardCM =   ((RSCM.lt.RECM) .and.                                  &
	 	&              (yday.ge.RSCM .and. yday.le.RECM))                    &
	 	&             .or.                                                   &
	 	&             ((RSCM.gt.RECM) .and.                                  &
	 	&              (yday.ge.RSCM .or.  yday.le.RECM))
	 	 downward   = ((SSNC.lt.SENC) .and.                                  &
	 	&              (yday.ge.SSNC .and. yday.le.SENC))                    &
	 	&             .or.                                                   &
	 	&             ((SSNC.gt.SENC) .and.                                  &
	 	&              (yday.ge.SSNC .or.  yday.le.SENC))
	 	 downwardCM = ((SSCM.lt.SECM) .and.                                  &
	 	&              (yday.ge.SSCM .and. yday.le.SECM))                    &
	 	&             .or.                                                   &
	 	&             ((SSCM.gt.SECM) .and.                                  &
	 	&              (yday.ge.SSCM .or.  yday.le.SECM))
	 	 if (upward .or. downward) then
	 	   respNC = respNCa * 0.3_r8
	 	 else
	 	   respNC = respNCa
	 	 end if
	 	 if (upwardCM .or. downwardCM) then
	 	   respCM = respNCa * 0.3_r8
	 	 else
	 	   respCM = respNCa
	 	 end if
!
!-----------------------------------------------------------------------
! Begin HORIZONTAL INDEX LOOPING
!-----------------------------------------------------------------------
!
      J_LOOP : DO j=Jstr,Jend
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &        (MOD(iic(ng)-1,nHIS(ng)).eq.0)) THEN
          IF (nrrec(ng).eq.0.or.iic(ng).ne.ntstart(ng)) THEN
          END IF
        END IF
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, and
!  restrict their values to be positive definite.  Removed CVL''s
!  conservation of mass correction because conflicted with SPLINES.
!  For ROMS 2.2+, convert from "flux form" to concentrations by
!  dividing by grid cell thickness.
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_bak(i,k,ibio)=max(t(i,j,k,nstp,ibio),0.0_r8)
              Bio(i,k,ibio)=Bio_bak(i,k,ibio)
              DBio(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
        DO itrc=1,NBEN
           ibioB=idben(itrc)
          DO k=1,NBL(ng)
            DO i=Istr,Iend
              Bio_bakB(i,k,ibioB)=max(bt(i,j,k,nstp,ibioB),0.0_r8)
              BioB(i,k,ibioB)=bt(i,j,k,nstp,ibioB)   
              DBioB(i,k,ibioB)=0.0_r8
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
            Bio(i,k,isalt)=t(i,j,k,nstp,isalt)
            IF (Bio(i,k,itemp) .gt. 35._r8) THEN
              print *, 'Temperature: ',                                 &
     &             Bio(i,k,itemp),i, j, k,ng, yday
              print *,'Tracer: ',t(i,j,k,1,itemp),t(i,j,k,2,itemp),     &
     &              t(i,j,k,3,itemp),Hz(i,j,k),nnew
              print *,'Others: ', z_w(i,j,N(ng)),                       &
     &                  GRID(ng) % h(i,j)
            END IF
            IF ((grid(ng) % h(i,j) + z_w(i,j,N(ng))) .lt.0.0_r8) THEN
              print *, 'zeta & h: ',z_w(i,j,N(ng)),'  ',grid(ng) % h(i,j)
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Calculate Day Length and Surface PAR
!-----------------------------------------------------------------------
!
!  Calculate Day Length
        DO i=Istr,Iend
!  Day Length is already accounted for in ANA_SWRAD so disable correction
          Dl = 24.0_r8
!
!  Calculate PAR at the surface
!
!
!  For PAR, use Shortwave radiation ( = surface solar irradiance)
!  converted from deg C m/s to E/m2/day   rho0=1025 Cp=3985
!             
          PARs(i) =PARfrac(ng) * srflx(i,j) * rho0 * Cp * 0.394848_r8
        END DO
! 
!  Calculate light decay in the water column    
!   
!----------------------------------------------------------------------- 
!  George Gibsons version after Morel 1988 (in Loukos 1997) 
!-----------------------------------------------------------------------
!  Version from Sarah Hinckley old C code 
        DO k=N(ng),1,-1
          DO i=Istr,Iend
            cff3 = z_r(i,j,k)+2.5_r8
            IF ( cff3 .gt. -71.0_r8 ) THEN
              cff1 = k_ext + k_chl *                                    &
     &                  ( Bio(i,k,iPhS) + Bio(i,k,iPhL) ) / ccr
              ELSE
                cff1 = 0.077_r8
              END IF
              PAR(i,k) = PARfrac(ng) * cff2 * exp( cff1 * cff3 )
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              HzL(i,k) = Hz(i,j,k)
              Sal1 = Bio(i,k,isalt)
              Temp1 = Bio(i,k,itemp)
!         
!-----------------------------------------------------------------------
!Compute sigma-t for each depth
!-----------------------------------------------------------------------
! 
              Dens(i,k) = ComputeDensity(Temp1,Sal1)
            END DO
          END DO
        DO i=Istr,Iend
!           
!-----------------------------------------------------------------------
!  Compute dsigmat/dz for each depth
!  Return the maximum value
!-----------------------------------------------------------------------
! 
          Do k=1,N(ng)
            DensV(k)=Dens(i,k)
            ZW_V(k)=z_wL(i,k)
          end do
          StabParam(i) = ComputeStability(ng,ZW_V,DensV)
        END DO
        DO k=0,N(ng)
          DO i=Istr,Iend
            z_wL(i,k) = z_w(i,j,k)
          END DO
        END DO
!************************************************************************
!************************************************************************
! Begin BIOITER LOOP
!************************************************************************
!************************************************************************
        ITER_LOOP: DO Iter=1,BioIter(ng)
!         if ( .not. downward) then
! 
!-----------------------------------------------------------------------
!  Make Neocalanus go down if temp > 12
!-----------------------------------------------------------------------
!                   
!           if (Bio(i,k,itemp) .gt. 12._r8 .and. NCa(k) .gt. 0.2_r8) then
!             downward = .true.
!             goto 111
!           end if
!         end if
! 111     Continue
          LightLim = 1.0_r8
          NOLim = 1.0_r8
          NHLim = 1.0_r8
          IronLim = 1.0_r8
!=======================================================================
!  Nutrient uptake by Small Phytoplankton
!=======================================================================
! Options for alpha
          DO k=1,N(ng)
            DO i=Istr,Iend
! Variable -depends on light/nutrient - gradual not step       
!             ALPHA_N = (Bio(i,N(ng),iNO3) +Bio(i,N(ng),iNH4))/ ( kN +Bio(i,N(ng),iNO3) +Bio(i,N(ng),iNH4));
!             ALPHA_P =  PAR(i,N(ng))/ ( kP +  PAR(i,N(ng)));
!             cff1= (ALPHA_N* ALPHA_P);
!             alphaPhSv=max(1.0_r8,cff1*alphaPhS);
! Constant
              alphaPhSv=alphaPhS
! 
!-----------------------------------------------------------------------
!  Growth rate computations
!-----------------------------------------------------------------------
! 
              if(DiS.gt.0_r8) THEN
                Drate = DiS * 10.0_r8 ** (DpS * Bio(i,k,itemp) )
                Pmax = (2.0_r8 ** Drate - 1.0_r8 )   !maximum daily mass specific growth rate FROST (1987)
                Pmaxs=Pmax*ccr !max chla specific growth rate from FROST (1987)
! day length fraction scalar
                Pmax = Pmax * Dl / 24.0_r8      
!-----------------------------------------------------------------------
!  Nitrate limitation
!-----------------------------------------------------------------------
!
!      
! limitation folowing Wroblewski (JMR 1977)    
! 
!               NOLim = Bio(i,k,iNO3) * EXP( -psiPhS * Bio(i,k,iNH4) )  &
!    &               / ( k1PhS + Bio(i,k,iNO3) )
! 
! limitation following Lomas (marine Biology 1999)    
! 
                NOLim = (Bio(i,k,iNO3)/ ( k1PhS +Bio(i,k,iNO3)))        &
     &                *(1-(0.8_r8*Bio(i,k,iNH4)/(k2PhS + Bio(i,k,iNH4))));  
! 
! Limitation following Vallina and Quere (Ecological Modelling 2008)
!               NOLim = Bio(i,k,iNO3) *                                 &
!     &                (1-((Bio(i,k,iNH4)/(k2PhS + Bio(i,k,iNH4))))) /  &
!     &                ( k1PhS + Bio(i,k,iNO3) 
!-----------------------------------------------------------------------
!  Light limitation function
!-----------------------------------------------------------------------
                Par1 = PAR(i,k)
                LightLim = TANH( alphaPhSv * MAX((Par1 - OffSet),0.0_r8)&
     &                         / Pmaxs)
!-----------------------------------------------------------------------
!  Nitrate uptake
!-----------------------------------------------------------------------
!  Multiplicative limitation
!               NOup = max(0.0_r8,                                      &
!    &               (Bio(i,k,iPhS)/ccr) * Pmaxs * LightLim * NOLim * IronLim)
!  Denman limitation 
                NOup = MAX(0.0_r8,(Bio(i,k,iPhS)/ccr) *                 &
     &               Pmaxs*MIN(LightLim, NOLim, IronLim))
!-----------------------------------------------------------------------
! Ammonium limitation
!-----------------------------------------------------------------------
                NHLim = Bio(i,k,iNH4) / ( k2PhS + Bio(i,k,iNH4) )
                if((NOLim+NHLim).gt.1.0_r8) then
                  NHLim= 1.0_r8-NOLim
                endif
!-----------------------------------------------------------------------
!  Light limitation for ammonium
!-----------------------------------------------------------------------
                LightLim = TANH( alphaPhSv * MAX((Par1 - OffSet),0.0_r8)&
     &                        / Pmaxs) 
!-----------------------------------------------------------------------
!  Ammonium uptake
!-----------------------------------------------------------------------
!               NHup = max(0.0_r8,                                      &
!    &                 (Bio(i,k,iPhS)/ccr) * Pmaxs * LightLim * NHLim)
                NHup = MAX(0.0_r8,(Bio(i,k,iPhS)/ccr) * Pmaxs * MIN(LightLim, NHLim))
              else
                NOup=0.0_r8
                NHup=0.0_r8
                cff3=0.0_r8
              endif
!  ajh limit NOup and NHup to amount of NO and NH present
              NOup=MIN(NOup,Bio(i,k,iNO3)/(xi*dtdays))
              NHup=MIN(NHup,Bio(i,k,iNH4)/(xi*dtdays))
              NOup=MAX(0.0_r8,NOup)
              NHup=MAX(0.0_r8,NHup)
!-----------------------------------------------------------------------
!  Change in nitrate concentration
!-----------------------------------------------------------------------
              DBio(i,k,iNO3) = DBio(i,k,iNO3) - xi * NOup * dtdays
!-----------------------------------------------------------------------
!  Change in ammonium concentration
!-----------------------------------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) - xi * NHup * dtdays
!-----------------------------------------------------------------------
!  Change in concentration of small phytoplankton
!-----------------------------------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS)+(NOup+NHup)*dtdays
!-----------------------------------------------------------------------
!  Primary production of small phytoplankton
!-----------------------------------------------------------------------
            END DO
          END DO
!=========================================================================
!  Nutrient uptake by Large Phytoplankton
!=========================================================================
          kN=0;
          kP=48;
!         
! Coyle step alpha function   
! 
!         if (ParMax.lt.48) then
!           Alpha = 1.0_r8
!         else
!           Alpha = 4.0_r8
!         endif
! 
! Depends on light/nutrient - gradual not step  
! 
!         ALPHA_N = (Bio(i,N(ng),iNO3) +Bio(i,N(ng),iNH4))/ ( kN +Bio(i,N(ng),iNO3) +Bio(i,N(ng),iNH4));
!         ALPHA_P =  PAR(i,k)/ ( kP +  PAR(i,k));
!         ALPHA_P =  PAR(i,N(ng))/ ( kP +  PAR(i,N(ng)));
!         cff1= (ALPHA_N* ALPHA_P);
!         alphaPhLv=max(1.0_r8,cff1*alphaPhL);
          alphaPhLv=alphaPhL          
          DO k=1,N(ng)
            DO i=Istr,Iend
              LightLim=1.0_r8 
              IronLim=1.0_r8    
              NOLim=1.0_r8
!  
!-----------------------------------------------------------------------
! Growth rate computations
!-----------------------------------------------------------------------
! 
              if(DiL.gt.0_r8) THEN
                Drate = DiL * 10.0_r8 ** (DpL * Bio(i,k,itemp) )
                Pmax = (2.0_r8 ** Drate - 1.0_r8 )
                Pmaxs=Pmax*ccrPhL
! 
!  day length fraction scalar
!             
!               Pmax = (2.0_r8 ** Drate - 1.0_r8 ) * Dl / 24.0_r8
!-----------------------------------------------------------------------
!  Nitrate limitation
!-----------------------------------------------------------------------
! 
! limitation following Wroblewski (JMR 1977)    
!               NOLim = Bio(i,k,iNO3) * EXP( -psiPhL * Bio(i,k,iNH4) )  &
!     &                  / ( k1PhL + Bio(i,k,iNO3) )
! limitation following Lomas (marine Biology 1999) 
                NOLim = (Bio(i,k,iNO3)/ ( k1PhL +Bio(i,k,iNO3)))        &
     &                *(1-(0.8_r8*Bio(i,k,iNH4)/(k2PhL + Bio(i,k,iNH4))));
! Limitation following Vallina and Quere (Ecological Modelling 2008)
!               NOLim = Bio(i,k,iNO3) *                                 &
!    &                 (1-(Bio(i,k,iNH4)/(k2PhL + Bio(i,k,iNH4)))))     &
!    &                 / ( k1PhL + Bio(i,k,iNO3)       
! 
!-----------------------------------------------------------------------
!  Light limitation for nitrate uptake function
!-----------------------------------------------------------------------
! 
                Par1 = Par(i,k)
                ParMax = Par(i,N(ng))
                OffSet=0.0_r8
                LightLim = TANH( alphaPhLv * MAX((PAR1 - OffSet),0.0_r8)&
     &                         / Pmaxs)
! 
!-----------------------------------------------------------------------
!  Nitrate uptake
!-----------------------------------------------------------------------
! 
!               NOup = MAX(0.0_r8,(Bio(i,k,iPhL)/ccrPhL)                &
!     &              * Pmaxs * LightLim * NOLim * IronLim)
                NOup=MAX(0.0_r8,(Bio(i,k,iPhL)/ccrPhL)*Pmaxs*MIN(LightLim ,NOLim,IronLim))
!               
!-----------------------------------------------------------------------
!  Ammonium limitation
!-----------------------------------------------------------------------
!  
                NHLim = Bio(i,k,iNH4) / ( k2PhL + Bio(i,k,iNH4) )
                if((NOLim+NHLim).gt.1.0_r8) then
                  NHLim= 1.0_r8-NOLim
                endif
!-----------------------------------------------------------------------
!  light limitation for ammonium uptake 
!-----------------------------------------------------------------------
                OffSet = 0.0_r8
                LightLim = TANH( alphaPhLv * MAX((PAR1 - OffSet),0.0_r8)&
     &                     / Pmaxs )
!-----------------------------------------------------------------------
!  Ammonium uptake
!-----------------------------------------------------------------------
!               NHup = max(0.0_r8,                                      &
!     &                  (Bio(i,k,iPhL)/ccrPhL) * Pmaxs * LightLim * NHLim)
                NHup = MAX(0.0_r8,(Bio(i,k,iPhL)/ccrPhL) * Pmaxs * MIN(LightLim ,NHLim))
              else
                NOup=0.0_r8
                NHup=0.0_r8
                cff3=0.0_r8
              endif
! 
!  ajh limit NOup and NHup to amount of NO and NH present
! 
              NOup=MIN(NOup,Bio(i,k,iNO3)/(xi*dtdays))
              NHup=MIN(NHup,Bio(i,k,iNH4)/(xi*dtdays))
              NOup=MAX(0.0_r8,NOup)
              NHup=MAX(0.0_r8,NHup)
!-----------------------------------------------------------------------
!  Change in nitrate concentration
!-----------------------------------------------------------------------
              DBio(i,k,iNO3) = DBio(i,k,iNO3) - xi * NOup * dtdays
!-----------------------------------------------------------------------
!  Change in ammonium concentration
!-----------------------------------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) - xi * NHup * dtdays
!-----------------------------------------------------------------------
!  Change in concentration of large phytoplankton
!-----------------------------------------------------------------------
! 
              DBio(i,k,iPhL) = DBio(i,k,iPhL)                           &
     &                       + ( NOup + NHup ) * dtdays
! 
!-----------------------------------------------------------------------
!  Primary production of large phytoplankton
!-----------------------------------------------------------------------
! 
            END DO
          END DO
!=======================================================================
! Grazing by MZS
!=======================================================================
!         DO k=1,N(ng)
!           DO i=Istr,Iend
!       
!-----------------------------------------------------------------------
!  Food preferences
!-----------------------------------------------------------------------
!             cff1 = fpPhSMZS * Bio(i,k,iPhS)**2                        &
!    &             + fpPhLMZS * Bio(i,k,iPhL)**2
!     
!-----------------------------------------------------------------------
!  Food consumption
!-----------------------------------------------------------------------
!             cff2 = eMZS * Bio(i,k,iMZS) / (fMZS**2 + cff1)
!             cff2 = eMZS * Bio(i,k,iMZS) / (fMZS + cff1)
!
!-----------------------------------------------------------------------
!  Temperature correction
!-----------------------------------------------------------------------
!             cff3 = Q10MZS ** ( (Bio(i,k,itemp)-Q10MZST)/ 10.0_r8)    
!             cff3 =1.0_r8
!
!                 
!-----------------------------------------------------------------------
!  Change in small and large phytoplankton due to predation
!-----------------------------------------------------------------------
!             DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSMZS *              &
!    &                        (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
!             DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLMZS *              &
!    &                        (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
!                     
!-----------------------------------------------------------------------  
!  Growth of small microzooplankton due to consumption
!-----------------------------------------------------------------------
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) +                         &
!    &                         gammaMZS * cff1 * cff2 * cff3 * dtdays
!     
!      
!-----------------------------------------------------------------------
!  Production for small microzooplankton
!-----------------------------------------------------------------------
!             Prod(i,k,iMZSprd) = Prod(i,k,iMZSprd) + DBio(i,k,iMZS)
!-----------------------------------------------------------------------
!  Additions to detritus pool - unassimilated food
!-----------------------------------------------------------------------
!             DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
!    &              (1.0_r8 - gammaMZS) * cff1 * cff2 * cff3 * dtdays
!     
!             IF (i.eq.3.and.j.eq.3) THEN
!
!               bflx(iPhS,iMZS) = bflx(iPhS,iMZS) +                     &            
!    &          fpPhSMZS * (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays*xi
!     
!               bflx(iPhL,iMZS) = bflx(iPhL,iMZS) +                     &
!    &          fpPhLMZS * (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays*xi
!       
!     
!               bflx(iMZS,iDet) = bflx(iMZS,iDet)  +                    &
!    &          ( 1.0_r8-gammaMZS )*cff1*cff2* cff3 * dtdays*xi
!             END IF
!           END DO
!         END DO
!========================================================================
! Grazing by MZL
!========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!  Food preferences
!-----------------------------------------------------------------------
              cff1 = fpPhSMZL * Bio(i,k,iPhS)**2                        &
     &             + fpPhLMZL * Bio(i,k,iPhL)**2                   
!    &             + fpMZSMZL * Bio(i,k,iMZS)**2
!-----------------------------------------------------------------------
! Food consumption
!-----------------------------------------------------------------------
!             cff2 = eMZL * Bio(i,k,iMZL) / (fMZL**2 + cff1)
              cff2 = Bio(i,k,iMZL) / (fMZL + cff1)
!               
!-----------------------------------------------------------------------
!  Temperature correction
!-----------------------------------------------------------------------
              cff3= eMZL *Q10MZL**((Bio(i,k,itemp)-Q10MZLT)/10.0_r8)
!             cff3= max(0.6_r8,cff3)
!             cff3= 1.0_r8
!-----------------------------------------------------------------------
!  Change in small and large phytoplankton due to predation
!-----------------------------------------------------------------------
! 
              DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSMZL *              &
     &                        (Bio(i,k,iPhS)**2) * cff2 * cff3* dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLMZL *              &
     &                        (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) - fpMZSMZL *              &
!    &                        (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Growth of large microzooplankton
!-----------------------------------------------------------------------
              DBio(i,k,iMZL) = DBio(i,k,iMZL) +                         &
     &                         gammaMZL * cff1 * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Additions to detritus pool - unassimilated food
!-----------------------------------------------------------------------
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &          (1.0_r8 - gammaMZL) * cff1 * cff2 * cff3 * dtdays
            END DO
          END DO
!===========================================================================
! Grazing and Predation by Copepods
!===========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!Food preferences
!-----------------------------------------------------------------------
! 
              cff1 = fpPhSCop * Bio(i,k,iPhS)**2                        &
     &               + fpPhLCop * Bio(i,k,iPhL)**2                      &
!    &               + fpMZSCop * Bio(i,k,iMZS)**2                      &
     &               + fpMZLCop * Bio(i,k,iMZL)**2     
!-----------------------------------------------------------------------
!  Food consumption
!-----------------------------------------------------------------------
!             cff2 = eCop * Bio(i,k,iCop) / (fCop**2 + cff1)
              cff2 = eCop * Bio(i,k,iCop) / (fCop + cff1)
!         
!-----------------------------------------------------------------------
!  Temperature correction
!-----------------------------------------------------------------------
! 
              cff3 = Q10Cop ** ( (Bio(i,k,itemp)-Q10CopT)/10.0_r8)
!             cff3=1.0_r8
! 
!-----------------------------------------------------------------------
!  Growth of small copepods
!-----------------------------------------------------------------------
              DBio(i,k,iCop) = DBio(i,k,iCop) +                         &
     &                        gammaCop * cff1 * cff2 * cff3 * dtdays
! 
!-----------------------------------------------------------------------
!  Copepod production  iCop=9
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!Changes in prey concentration due to predation
!-----------------------------------------------------------------------
! 
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -  fpPhSCop               &
     &                * (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -  fpPhLCop               &
     &                * (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -  fpMZSCop               &
!    &                * (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -  fpMZLCop               &
     &                * (Bio(i,k,iMZL)**2) * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Additions to detritus pool - unassimilated food
!-----------------------------------------------------------------------
! 
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &            (1.0_r8 - gammaCop) * cff1 * cff2 * cff3 * dtdays
            END DO
          END DO
!g#ifdef ICE_BIO
!g        DO i=Istr,Iend
!g          DBioBI(i,iIcePhL) = DBioBI(i,iIcePhL)                       &
!g   &         - Hz(i,j,N(ng))*(fpPhLCop*                               &
!g   &         (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays)
!g
!g#if defined BIOFLUX && defined 
!g         bflx(NT(ng)+3,iCop)= bflx(NT(ng)+3,iCop)                     &
!g   &         + Hz(i,j,N(ng))*(fpPhLCop*                               &
!g   &       (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays)
!g#endif
!g      
!g        END DO
!g
!g#endif    
!========================================================================
! Grazing and Predation by NCa initiated ON the shelf
!========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1 = fpPhSNCa * Bio(i,k,iPhS)**2                        &
     &                 + fpPhLNCa * Bio(i,k,iPhL)**2                    &
!    &                 + fpMZSNCa * Bio(i,k,iMZS)**2                    &
     &                 + fpMZLNCa * Bio(i,k,iMZL)**2
!-----------------------------------------------------------------------
!Food consumption
!-----------------------------------------------------------------------
!      
!             cff2 = eNCa * Bio(i,k,iNCaS) / (fNCa**2 + cff1)
              cff2 = eNCa * Bio(i,k,iNCaS) / (fNCa + cff1)
!               
!-----------------------------------------------------------------------
!Temperature correction 
!-----------------------------------------------------------------------
! 
              cff3 = Q10NCa ** ( (Bio(i,k,itemp)-Q10NCaT)/10.0_r8) 
!-----------------------------------------------------------------------
!  Growth of Neocalanus
!-----------------------------------------------------------------------
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) +                       &
     &                     gammaNCa * cff1 * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!Neocalanus production  iNCaS=10
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!  Changes in prey concentration due to predation
!-----------------------------------------------------------------------
! 
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &           fpPhSNCa * (Bio(i,k,iPhS)**2) * cff2 * cff3 *dtdays
                 DBio(i,k,iPhL) = DBio(i,k,iPhL) -                      &
     &           fpPhLNCa * (Bio(i,k,iPhL)**2) * cff2 * cff3 *dtdays     
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!    &           fpMZSNCa * (Bio(i,k,iMZS)**2) * cff2 * cff3 *dtdays
                 DBio(i,k,iMZL) = DBio(i,k,iMZL) -                      &
     &           fpMZLNCa * (Bio(i,k,iMZL)**2) * cff2 * cff3 *dtdays
!-----------------------------------------------------------------------
!  Additions to Fast Sinking detritus pool - unassimilated food
!-----------------------------------------------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &          (1.0_r8 - gammaNCa) * cff1 * cff2 * cff3 * dtdays
            END DO
          END DO
!g#ifdef ICE_BIO
!g        DO i=Istr,Iend
!g       
!g          DBioBI(i,iIcePhL) = DBioBI(i,iIcePhL) - Hz(i,j,N(ng)) *     &
!g            (fpPhLNCa*(BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 *      &
!g            cff3 * dtdays)
!g
!g#if defined BIOFLUX && defined            
!g          bflx(NT(ng)+3,iNCaS)= bflx(NT(ng)+3,iNCaS)                  &
!g   &        + Hz(i,j,N(ng))*(fpPhLNCa*                                &
!g   &        (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays*xi)
!g#endif
!g      
!g
!g
!g        END DO
!g#endif
!=========================================================================
! Grazing and Predation by Euphuasiids initiated ON the shelf
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!Food preferences
!-----------------------------------------------------------------------
              cff1 = fpPhSEup * Bio(i,k,iPhS)**2                        &

     &               + fpPhLEup * Bio(i,k,iPhL)**2                      &

!    &               + fpMZSEup * Bio(i,k,iMZS)**2                      &

     &               + fpMZLEup * Bio(i,k,iMZL)**2                      &

     &               + fpCopEup * Bio(i,k,iCop)**2  
!-----------------------------------------------------------------------
!  Food consumption
!-----------------------------------------------------------------------
! 
!             cff2 = eEup * Bio(i,k,iEupS) / (fEup**2 + cff1)
              cff2 = eEup * Bio(i,k,iEupS) / (fEup + cff1)
!-----------------------------------------------------------------------
!  Temperature correction 
!-----------------------------------------------------------------------
! 
              cff3 = Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8 )
!-----------------------------------------------------------------------
!  Growth of Euphausiids
!-----------------------------------------------------------------------
              DBio(i,k,iEupS) = DBio(i,k,iEupS) +                       &
     &          gammaEup * cff1 * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Euphausiid production   iEupS=11
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Changes in prey concentration due to predation
!-----------------------------------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSEup *              &
     &          (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLEup *              &
     &          (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) - fpMZSEup *              &
!    &          (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) - fpMZLEup *              &
     &          (Bio(i,k,iMZL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iCop) = DBio(i,k,iCop) - fpCopEup *              &
     &          (Bio(i,k,iCop)**2) * cff2 * cff3 * dtdays
!                     
!-----------------------------------------------------------------------
! Additions to Fast Sinking detritus pool- unassimilated food
!-----------------------------------------------------------------------
! 
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &          (1.0_r8 - gammaEup) * cff1 * cff2 * cff3 * dtdays
            END DO
          END DO
!g#ifdef ICE_BIO
!g        DO i=Istr,Iend
!g          DBioBI(i,iIcePhL) = DBioBI(i,iIcePhL) -                     &
!g   &        Hz(i,j,N(ng))*(fpPhLEup*                                  &
!g   &        (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays)
!g
!g#if defined BIOFLUX && defined            
!g          bflx(NT(ng)+3,iEupO)= bflx(NT(ng)+3,iEupS)                  &
!g   &        + Hz(i,j,N(ng))*(fpPhLEup*                                &
!g   &        (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays*xi)
!g#endif
!g      
!g        END DO
!g#endif     
!========================================================================
! Grazing and Predation by NCa initiated OFF the shelf
!========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!  Food preferences
!-----------------------------------------------------------------------
              cff1 = fpPhSNCa * Bio(i,k,iPhS)**2                        &
     &                 + fpPhLNCa * Bio(i,k,iPhL)**2                    &
!    &                 + fpMZSNCa * Bio(i,k,iMZS)**2                    &
     &                 + fpMZLNCa * Bio(i,k,iMZL)**2
!-----------------------------------------------------------------------
!  Food consumption
!-----------------------------------------------------------------------
!             cff2 = eNCa * Bio(i,k,iNCaO) / (fNCa**2 + cff1)
              cff2 = eNCa * Bio(i,k,iNCaO) / (fNCa + cff1)
!-----------------------------------------------------------------------
!  Temperature correction 
!-----------------------------------------------------------------------
              cff3 = Q10NCa ** ( (Bio(i,k,itemp)-Q10NCaT) / 10.0_r8 ) 
!-----------------------------------------------------------------------
!  Growth of Neocalanus
!-----------------------------------------------------------------------
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) +                       &
     &                        gammaNCa * cff1 * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Changes in prey concentration due to predation
!-----------------------------------------------------------------------
!   
              DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSNCa *              &
     &          (Bio(i,k,iPhS)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLNCa *              &
     &          (Bio(i,k,iPhL)**2) * cff2 * cff3 *dtdays
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) - fpMZSNCa *              &
!    &          (Bio(i,k,iMZS)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) - fpMZLNCa *              &
     &          (Bio(i,k,iMZL)**2) * cff2 * cff3 *dtdays
! 
!-----------------------------------------------------------------------
!  Additions to detritus pool - unassimilated food
!-----------------------------------------------------------------------
! 
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &          (1.0_r8 - gammaNCa) * cff1 * cff2 * cff3 * dtdays
            END DO
          END DO
!g#ifdef ICE_BIO         
!g        DO i=Istr,Iend
!g          DBioBI(i,iIcePhL) = DBioBI(i,iIcePhL) -                     &
!g   &        Hz(i,j,N(ng))*(fpPhLNCa*                                  &
!g   &        (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays)
!g
!g
!g#if defined BIOFLUX && defined            
!g          bflx(NT(ng)+3,iNCaO)= bflx(NT(ng)+3,iNCaO)                  &
!g   &        + Hz(i,j,N(ng))*(fpPhLNCa*                                &
!g   &        (BioBI(i,iIcePhL)/Hz(i,j,N(ng)))**2*cff2 * cff3 * dtdays*xi)
!g#endif
!g 
!g        END DO
!g#endif
!=========================================================================
! Grazing and Predation by Euphuasiids initiated OFF the shelf
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!  Food preferences
!-----------------------------------------------------------------------
! 
              cff1 = fpPhSEup * Bio(i,k,iPhS)**2                        &
     &               + fpPhLEup * Bio(i,k,iPhL)**2                      &
!    &               + fpMZSEup * Bio(i,k,iMZS)**2                      &
     &               + fpMZLEup * Bio(i,k,iMZL)**2                      &
     &               + fpCopEup * Bio(i,k,iCop)**2 
!-----------------------------------------------------------------------
!Food consumption
!-----------------------------------------------------------------------
! 
!             cff2 = eEup * Bio(i,k,iEupO) / (fEup**2 + cff1)
              cff2 = eEup * Bio(i,k,iEupO) / (fEup + cff1)
!-----------------------------------------------------------------------
!  Temperature correction 
!-----------------------------------------------------------------------
! 
              cff3 = Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8 )
!-----------------------------------------------------------------------
!Growth of Euphausiids
!-----------------------------------------------------------------------
              DBio(i,k,iEupO) = DBio(i,k,iEupO) +                       &
     &          gammaEup * cff1 * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Euphausiid production
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Changes in prey concentration due to predation
!-----------------------------------------------------------------------
! 
              DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSEup *              &
     &          (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLEup *              &
     &          (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) - fpMZSEup *              &
!    &          (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) - fpMZLEup *              &
     &          (Bio(i,k,iMZL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iCop) = DBio(i,k,iCop) - fpCopEup *              &
     &          (Bio(i,k,iCop)**2) * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Additions to detritus pool- unassimilated food
!-----------------------------------------------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &          (1.0_r8 - gammaEup) * cff1 * cff2 * cff3 * dtdays
            END DO
          END DO    
!=========================================================================
! Grazing and Predation by Jellyfish
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!  Food preferences
!-----------------------------------------------------------------------
              cff1 = fpCopJel * Bio(i,k,iCop)**2 +                      &
     &               fpNCaJel * Bio(i,k,iNCaS)**2 +                     &
     &               fpNCaJel * Bio(i,k,iNCaO)**2 +                     &
     &               fpEupJel * Bio(i,k,iEupS)**2 +                     &
     &               fpEupJel * Bio(i,k,iEupO)**2
!-----------------------------------------------------------------------
!  Food consumption 
!-----------------------------------------------------------------------
! 
!  Linear
! 
!             cff2 = eJel 
! 
!  MM
! 
              cff2 = eJel * Bio(i,k,iJel) / (fJel + cff1)
!-----------------------------------------------------------------------
!  Temperature correction 
!-----------------------------------------------------------------------
! 
              cff3= Q10Jele ** ((Bio(i,k,itemp)-Q10JelTe) / 10.0_r8)  
!             cff3=1.0_r8
!-----------------------------------------------------------------------
!  Growth of Jellies
!-----------------------------------------------------------------------  
! 
              DBio(i,k,iJel) =  DBio(i,k,iJel) +                        &
     &          gammaJel * cff1 * cff2 * cff3 * dtdays
!-----------------------------------------------------------------------
!  Jellyfish production
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Changes in prey concentration due to predation
!-----------------------------------------------------------------------
              DBio(i,k,iCop) = DBio(i,k,iCop) - fpCopJel *              &
     &          Bio(i,k,iCop)**2 * cff2 * cff3 * dtdays
              DBio(i,k,iEupS) = DBio(i,k,iEupS) - fpEupJel *            &
     &          Bio(i,k,iEupS)**2* cff2 * cff3 * dtdays
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) - fpNCaJel *            &
     &          Bio(i,k,iNCaS)**2* cff2 * cff3 * dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupO) - fpEupJel *            &
     &          Bio(i,k,iEupO)**2* cff2 * cff3 * dtdays
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) - fpNCaJel *            &
     &          Bio(i,k,iNCaO)**2* cff2 * cff3 * dtdays
              DBio(i,k,iDetF)=DBio(i,k,iDetF) +(1-gammaJel)             &
     &             * cff1 * cff2 * cff3 * Bio(i,k,iJel)*dtdays
            END DO
          END DO
!=======================================================================
! Phytoplankton Linear Mortality and Senescence Terms
!=======================================================================
! 
          DO k=1,N(ng)
            DO i=Istr,Iend
!             cff1 = MAX( minmPhS , maxmPhS -                           &
!    &              ( maxmPhS - minmPhS ) * Bio(i,k,iNO3) / NcritPhS)
!             cff2 = MAX( minmPhL , maxmPhL -                           &
!    &              ( maxmPhL - minmPhL ) * Bio(i,k,iNO3) / NcritPhL)
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &                         mPhS* Bio(i,k,iPhS) * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &                         mPhL* Bio(i,k,iPhL) * dtdays
!-----------------------------------------------------------------------
!  Additions to detritus pool - phytoplankton mort
!-----------------------------------------------------------------------
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &                         ( mPhS * Bio(i,k,iPhS)) * dtdays
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &                         ( mPhL * Bio(i,k,iPhL)) * dtdays
            END DO
          END DO
!           
!=======================================================================
!  Microzooplankton Mortality - use only linear OR QUADRATIC
!=======================================================================
! 
          DO k=1,N(ng)
            DO i=Istr,Iend
! 
!  Linear   
! 
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!     &                        mMZS * Bio(i,k,iMZS) * dtdays
!             DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
!     &                         mMZL * Bio(i,k,iMZL) * dtdays
!
!  Quadratic 
!  
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!    &                         mpredMZS*dtdays*Bio(i,k,iMZS)**2
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &                         mpredMZL*dtdays*Bio(i,k,iMZL)**2
!-----------------------------------------------------------------------
!  Additions to detritus pool - natural microzoo mortality
!-----------------------------------------------------------------------
! 
!  if linear (George)
! 
!             Bio(i,k,iDet) = DBio(i,k,iDet) +                          &
!    &                        mMZL * Bio(i,k,iMZL) * dtdays   
! 
!  if quadratic (Ken)
! 
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
      &                        (mpredMZL * Bio(i,k,iMZL)**2 ) * dtdays 
            END DO
          END DO
!=======================================================================
!  Mesozooplankton Mortality (Closure terms)
!=======================================================================
! 
          DO k=1,N(ng)
            DO i=Istr,Iend
              TFEup = Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8)
              predSumCop  = TFEup*(mpredCop)*Bio(i,k,iCop)**2
              predSumNCaS = TFEup*(mpredNca)*Bio(i,k,iNCaS)**2     
              predSumEupS = TFEup*(mpredEup)*Bio(i,k,iEupS)**2
              predSumNCaO = TFEup*(mpredNca)*Bio(i,k,iNCaO)**2
              predSumEupO = TFEup*(mpredEup)*Bio(i,k,iEupO)**2  
              DBio(i,k,iCop)  = DBio(i,k,iCop)  - predSumCop  * dtdays
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) - predSumNCaS * dtdays               
              DBio(i,k,iEupS) = DBio(i,k,iEupS) - predSumEupS * dtdays                   
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) - predSumNCaO * dtdays                   
              DBio(i,k,iEupO) = DBio(i,k,iEupO) - predSumEupO * dtdays                    
!-----------------------------------------------------------------------
!  Detritus from nonlinear mortality
!-----------------------------------------------------------------------  
! 
              DBio(i,k,iDetF) = DBio(i,k,iDetF) + dtdays *              &
     &                          (predSumCop + predSumNCaS + predSumEupS &
     &                          + predSumNCaO + predSumEupO)
              DBio(i,k,iJel) = DBio(i,k,iJel) - mpredJel *              &
     &                         Bio(i,k,iJel)**2  * dtdays
              DBio(i,k,iDetF) = DBio(i,k,iDetF)                         &
     &                   + mpredJel * Bio(i,k,iJel)**2  * dtdays
            END DO
          END DO
!=======================================================================
!Phytoplankton respiration losses
!=======================================================================
! 
          DO k=1,N(ng)
            DO i=Istr,Iend
              BasalMet = respPhS
!-----------------------------------------------------------------------
!  Arhonditsis temperature functions
!-----------------------------------------------------------------------
! 
              TempFuncPhS(i,k) = GetPhytoResp2(Temp1,TmaxPhS,KtBm_PhS)
!-----------------------------------------------------------------------
!  Change in concentration of Small Phytoplankton
!----------------------------------------------------------------------- 
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &          TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &          xi * TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
!-----------------------------------------------------------------------
!  Primary production of Small phytoplankton
!-----------------------------------------------------------------------
              BasalMet = respPhL
 !-----------------------------------------------------------------------
 !  Arhonditsis temperature functions
 !-----------------------------------------------------------------------
              TempFuncPhL(i,k) = GetPhytoResp2(Temp1,TmaxPhL,KtBm_PhL)
!-----------------------------------------------------------------------
!  Change in concentration of Large Phytoplankton
!-----------------------------------------------------------------------
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &          TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
               DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
     &           xi * TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
!-----------------------------------------------------------------------
!  Primary production of Large phytoplankton
!-----------------------------------------------------------------------
! 
            END DO
          END DO
!=======================================================================
!  Microzooplankton respiration losses
!=======================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!             cff1 = fpPhSMZL * Bio(i,k,iPhS)+ fpPhLMZL * Bio(i,k,iPhL)
!             if(cff1.lt.1.0_r8)THEN
!               BasalMetMZL=0.0_r8
!             else
!               BasalMetMZL=respMZL 
!             end if
!-----------------------------------------------------------------------      
!  Small Microzooplankton 
!-----------------------------------------------------------------------
! 
!  Arhonditsis temperature functions
! 
!             TempFuncMZS(i,k) = GetPhytoResp2(Temp1,TmaxMZS,           &
!    &                           KtBm_MZS)
!             BasalMet = respMZS
!             TFMZL = Q10MZS ** ( (Bio(i,k,itemp)-Q10MZST) / 10.0_r8 )
! 
!-----------------------------------------------------------------------
!  Change in concentration of small microzooplankton
!-----------------------------------------------------------------------
!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!     &         TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)          
!     &         TFMZS*respMZS*dtdays*Bio(i,k,iMZS)
! 
!-----------------------------------------------------------------------
!  Small Microzooplankton production
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!  Add ammonium to correct for excretion related to metabolism
!-----------------------------------------------------------------------
!             DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
!     &         xi*(TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS))
!     &         xi*(TFMZS*respMZS*dtdays*Bio(i,k,iMZS))
!-----------------------------------------------------------------------      
!  Large Microzooplankton 
!-----------------------------------------------------------------------
! 
!  Arhonditsis temperature functions
! 
              TFMZL = exp(KtBm_MZL * (Bio(i,k,itemp) - TmaxMZL))
              BasalMetMZL = respMZL
!-----------------------------------------------------------------------
!  Change in concentration of large microzooplankton
!-----------------------------------------------------------------------
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
!    &          TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)
     &          TFMZL* BasalMetMZL*dtdays*Bio(i,k,iMZL)
! 
!-----------------------------------------------------------------------
!  Large Microzooplankton net production
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!  Add ammonium to correct for excretion related to metabolism
!-----------------------------------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &

!    &          xi*(TempFuncMZL(i,k)*BasalMetMZL*dtdays*Bio(i,k,iMZL))
     &          xi*(TFMZL* BasalMetMZL*dtdays*Bio(i,k,iMZL))
            END DO
          END DO
!=======================================================================
!  Mesozooplankton respiration losses
!=======================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1 = fpPhSCop * Bio(i,k,iPhS)+ fpPhLCop * Bio(i,k,iPhL) &
     &             + fpMZLCop * Bio(i,k,iMZL)
              if(cff1.lt.0.01_r8)THEN  !0.05
                BasalMetCop= respCop*cff1/0.01_r8 
              else
                BasalMetCop= respCop
              endif 
!-----------------------------------------------------------------------
!  Copepod respiration correction
!-----------------------------------------------------------------------
              TFCop = exp(ktbmC * (Bio(i,k,itemp) - TrefC))
!-----------------------------------------------------------------------
!  Neocalanus respiration correction
!-----------------------------------------------------------------------
              TFNCa = exp(ktbmN * (Bio(i,k,itemp) - TrefN))
!-----------------------------------------------------------------------
!  Euphausiid respiration correction
!-----------------------------------------------------------------------
              TFEup = exp(ktbmE * (Bio(i,k,itemp) - TrefE))
! 
!  starvation response  
! 
              cff1 = fpPhSEup * Bio(i,k,iPhS)                           &
     &               + fpPhLEup * Bio(i,k,iPhL)                         &
!    &               + fpMZSEup * Bio(i,k,iMZS)                         &
     &               + fpMZLEup * Bio(i,k,iMZL)                         &
     &               + fpCopEup * Bio(i,k,iCop)                    
              cff2 = fpPhSNCa * Bio(i,k,iPhS)+ fpPhLNCa * Bio(i,k,iPhL) &
!    &               + fpMZSNCa * Bio(i,k,iMZS)                         &
     &               + fpMZLNCa * Bio(i,k,iMZL)
              if(cff1.lt.0.01_r8)THEN
                BasalMetEup=respEup*cff1/0.01_r8
              else
                BasalMetEup= respEup 
              endif 
              if(cff2.lt.0.01_r8)THEN
                BasalMetNC= respNC*cff2/0.01_r8 
                BasalMetCM= respCM*cff2/0.01_r8  
              else
                BasalMetNC=respNC
                BasalMetCM=respCM
              end if 
!-----------------------------------------------------------------------
!  Change in concentration from small copepod respiration
!-----------------------------------------------------------------------
! 
              DBio(i,k,iCop) = DBio(i,k,iCop) -                         &
     &          BasalMetCop*TFCop*Bio(i,k,iCop)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &          xi*(TFCop*BasalMetCop*dtdays*Bio(i,k,iCop))
!-----------------------------------------------------------------------
!  Change in concentration from large copepods respiration 
!-----------------------------------------------------------------------
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) -                       &
     &            TFNCa*BasalMetCM*Bio(i,k,iNCaS)*dtdays
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) -                       &
     &            TFNCa*BasalMetNC*Bio(i,k,iNCaO)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFNCa*BasalMetCM*dtdays*Bio(i,k,iNCaS))
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFNCa*BasalMetNC*dtdays*Bio(i,k,iNCaO))
!-----------------------------------------------------------------------
!  Change in concentration from euphausiid respiration
!-----------------------------------------------------------------------
! 
              DBio(i,k,iEupS) = DBio(i,k,iEupS) -                       &
     &            TFEup*BasalMetEup*Bio(i,k,iEupS)*dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupO) -                       &
     &            TFEup*BasalMetEup*Bio(i,k,iEupO)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFEup*BasalMetEup* dtdays*Bio(i,k,iEupS))
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFEup*BasalMetEup*dtdays*Bio(i,k,iEupO))
              BasalMetJel = respJel
              TFJel = Q10Jelr ** ( (Bio(i,k,itemp)-Q10JelTr)/10.0_r8)
!             TFJel =1.0_r8
              DBio(i,k,iJel) = DBio(i,k,iJel) -                         &
     &            TFJel *BasalMetJel* Bio(i,k,iJel)*dtdays   
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFJel*BasalMetJel*dtdays*Bio(i,k,iJel))
            END DO
          END DO
!=========================================================================
! Molting: 
!=========================================================================
!-----------------------------------------------------------------------
!  NOTE: It is unclear where molting equation came from.
!  This is present only for euphausiids, not copepods
!-----------------------------------------------------------------------
!g        DO k=1,N(ng)
!g          DO i=Istr,Iend
!g            cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*     &
!g   &               Bio(i,k,iEupS)             
!g            DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
!g            DBio(i,k,iEupS) = DBio(i,k,iEupS) - cff1 * dtdays
!g
!g            cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*     &
!g   &               Bio(i,k,iEupO)                 
!g            DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
!g            DBio(i,k,iEupO) = DBio(i,k,iEupO) - cff1 * dtdays
!g     
!g#if defined BIOFLUX && defined 
!g            IF (i.eq.3.and.j.eq.3) THEN
!g              bflx(iEupS,iDet)= bflx(iEupS,iDet) + cff1* dtdays*xi
!g              bflx(iEupO,iDet)= bflx(iEupO,iDet) + cff1* dtdays*xi
!g            END IF
!g#endif
!g          END DO
!g        END DO
!=========================================================================
! Detrital Remineralization   (Det -> NH4)
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
              Temp1 = Bio(i,k,itemp)
!-----------------------------------------------------------------------
!  From Frost (1993).
!-----------------------------------------------------------------------
!             cff1 = regen * dgrad * Bio(i,k,iDet)
!             DBio(i,k,iNH4) = DBio(i,k,iNH4) + xi * cff1 * dtdays
!             DBio(i,k,iDet) = DBio(i,k,iDet) - cff1 * dtdays
!-----------------------------------------------------------------------
!  From Kawamiya(2000)
!-----------------------------------------------------------------------
              PON = Bio(i,k,iDet)*xi  !Particulate organic nitrogen
              DBio(i,k,iDet) = DBio(i,k,iDet) -                         &
     &            ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            ((Pv0*exp(PvT*Temp1)*PON))*dtdays
              PON = Bio(i,k,iDetF)*xi  !Particulate organic nitrogen
              DBio(i,k,iDetF) = DBio(i,k,iDetF) -                       &
     &            ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            ((Pv0*exp(PvT*Temp1)*PON))*dtdays
            END DO
          END DO
!=========================================================================
!Nitrification  (NH4 -> NO3)
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
!-----------------------------------------------------------------------
!Temperature dependance
!-----------------------------------------------------------------------
! 
!  Kawamiya 2000  NitrMax=Nitr0*exp(KnT*Temp1) - Ken
!   
!             NitrifMax=GetNitrifMaxK(Temp1,KnT,Nitr0)
! 
!  Arhonditsis NitrMax=Nitr0*exp(-ktntr*(Temp1 - ToptNtr)**2)
! 
              NitrifMax=GetNitrifMaxA(Nitr0, ktntr,Temp1,ToptNtr)
! 
!  No temperaure effects - NitrMax is constant
! 
!             NitrifMax=Nitr0
!-----------------------------------------------------------------------
!  Light/Depth dependance
!-----------------------------------------------------------------------
! 
!  Fennel  
!             DLNitrif=GetNitrifLight(Par1,tI0,KI)
! 
! Denman
!             DLNitrif = (z_wL(k)**10_r8)/( (20_r8**10_r8) +            &
!    &                    z_wL(k)**10_r8)
! 
!  No Depth/Ligt effects
! 
              DLNitrif= 1.0_r8
!-----------------------------------------------------------------------
!  Saturation
!-----------------------------------------------------------------------
! 
!  Arhonditsis
!         
!             cff1 = Bio(i,k,iNH4)/(KNH4Nit +Bio(i,k,iNH4))
! 
!  No saturation -ken
! 
              cff1 =1.0_r8
              DBio(i,k,iNH4) = DBio(i,k,iNH4)  - NitrifMax *            &
     &            Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays 
              DBio(i,k,iNO3) = DBio(i,k,iNO3) + NitrifMax *             &
     &            Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays    
            END DO
          END DO
!======================================================================
!Benthic Sub Model
!======================================================================
!         DO k=1,NBL(ng)
          DO i=Istr,Iend
!-----------------------------------------------------------------------
!  Growth of Benthic Infauna
!-----------------------------------------------------------------------
! 
!  calculate average water column concentrations
!   
            cff1=0.0_r8
            cff2=0.0_r8
            cff3=0.0_r8
            cff4=0.0_r8
            do k = 1,N(ng)
              cff1 = cff1 + Bio(i,k,iDet)  * Hz(i,j,k)
              cff2 = cff2 + Bio(i,k,iDetF) * Hz(i,j,k)
              cff3 = cff3 + Bio(i,k,iPhS)  * Hz(i,j,k)
              cff4 = cff4 + Bio(i,k,iPhL)  * Hz(i,j,k)
            end do
            avgD  = cff1/grid(ng) % h(i,j)
            avgDF = cff2/grid(ng) % h(i,j)
            avgPS = cff3/grid(ng) % h(i,j)
            avgPL = cff4/grid(ng) % h(i,j)
! use  food concentration in bottom layer of WC only
!           avgD  =  Bio(i,k,iDet)
!           avgDF =  Bio(i,k,iDetF)
!           avgPS =  Bio(i,k,iPhS)
!           avgPL =  Bio(i,k,iPhL)
            k = 1
            Temp1 = Bio(i,k,itemp)
            cff0 = q10r**((Temp1-T0benr)/10.0_r8) 
! 
! Potential food available from water column
! 
!           dw=3.0_r8 !assume bottom 3ms available
            dw=1.0_r8 !assume only bottom meter/layer available
            cff1=(prefD *avgD *dw/((prefD *avgD *dw)+LupP))*prefD *avgD *dw
            cff2=(prefD *avgDF*dw/((prefD *avgDF*dw)+LupP))*prefD *avgDF*dw
            cff3=(prefPS*avgPS*dw/((prefPS*avgPS*dw)+LupP))*prefPS*avgPS*dw
            cff4=(prefPL*avgPL*dw/((prefPL*avgPL*dw)+LupP))*prefPL*avgPL*dw
!           print*, 'Food=',cff1, cff2, cff3, cff4
! 
!Potential food available from  sea floor
! 
            cff5 = (prefD * BioB(i,k,iBenDet) / (prefD *                &
     &             BioB(i,k,iBenDet) + LupD)) * prefD * BioB(i,k,iBenDet)
! 
!  Total pelagic food available       
! 
            cff6 = cff1+cff2+cff3+cff4
!  uptake of each food category
            cff7  = min(cff1,(cff0*cff1*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff8  = min(cff2,(cff0*cff2*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff9  = min(cff3,(cff0*cff3*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff10 = min(cff4,(cff0*cff4*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff11 = min(cff5,(cff0*cff5*BioB(i,k,iBen)*Rup/(cff5+KupD)))
            DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)-dtdays*cff11
! 
!  if just removing from base layer 
! 
            k=1     
            DBio(i,k,iDet) = DBio(i,k,iDet)  - dtdays*cff7/ Hz(i,j,k)
            DBio(i,k,iDetF)= DBio(i,k,iDetF) - dtdays*cff8/ Hz(i,j,k)
            DBio(i,k,iPhS) = DBio(i,k,iPhS)  - dtdays*cff9/ Hz(i,j,k)
            DBio(i,k,iPhL) = DBio(i,k,iPhL)  - dtdays*cff10/Hz(i,j,k)
!  Spread removal evenly over water column 
!  ajh note I have changed the following section of benthic removal    
!           DO k=1,N(ng)
!             DBio(i,k,iDet) =DBio(i,k,iDet)-dtdays*cff7/grid(ng) % h(i,j)
!             DBio(i,k,iDetF)=DBio(i,k,iDetF)-dtdays*cff8/grid(ng) % h(i,j)
!             DBio(i,k,iPhS) =DBio(i,k,iPhS)-dtdays*cff9/grid(ng) % h(i,j)
!             DBio(i,k,iPhL) =DBio(i,k,iPhL)-dtdays*cff10/grid(ng)% h(i,j)
!           END DO
            k=1
            DBioB(i,k,iBen)=DBioB(i,k,iBen) + cff11*dtdays
            DBioB(i,k,iBen)=DBioB(i,k,iBen) + cff7 *dtdays
            DBioB(i,k,iBen)=DBioB(i,k,iBen) + cff8 *dtdays
            DBioB(i,k,iBen)=DBioB(i,k,iBen) + cff9 *dtdays
            DBioB(i,k,iBen)=DBioB(i,k,iBen) + cff10*dtdays
!-----------------------------------------------------------------------
! Benthic Production
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!  Excretion
!-----------------------------------------------------------------------
            cff1=cff7 *eexD 
            cff2=cff8 *eexD
            cff3=cff9 *eex
            cff4=cff10*eex
            cff5=cff11*eexD
            DBioB(i,k,iBen)=DBioB(i,k,iBen)                             &
     &           -(cff1+cff2+cff3+cff4+cff5)*dtdays   
!  Material excreted is considered inorganic and not available
!  for further secondary production
            DBio(i,k,iNH4)= DBio(i,k,iNH4)+xi*dtdays*0.5_r8             &
     &         *(cff1+cff2+cff3+cff4+cff5)/Hz(i,j,k)
            DBioB(i,k,iBenDet)= DBioB(i,k,iBenDet)+dtdays               &
     &         *(cff1+cff2+cff3+cff4+cff5)*0.5_r8
!----------------------------------------------------------------------- 
!  Respiration
!-----------------------------------------------------------------------
            cff3=cff0*BioB(i,k,iBen)*Rres 
            cff4=Qres*(((1_r8-eexD)*cff1) + ((1_r8-eexD)*cff2) +        &
    &            ((1_r8-eex)*cff3) + ((1_r8-eex)*cff4) +                &
    &            ((1_r8-eexD)*cff5))
            cff6=cff3+cff4
            DBioB(i,k,iBen)=DBioB(i,k,iBen) -cff6*dtdays
            DBio(i,k,iNH4)= DBio(i,k,iNH4)+ xi*dtdays*cff6/Hz(i,j,k) 
!-----------------------------------------------------------------------
!  Benthic Production
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Mortality
!-----------------------------------------------------------------------
            cff1= rmort*BioB(i,k,iBen)*cff0
            DBioB(i,k,iBen) = DBioB(i,k,iBen)- cff1 * dtdays
            DBioB(i,k,iBenDet)= DBioB(i,k,iBenDet)+cff1*dtdays
!-----------------------------------------------------------------------
!  Predation
!-----------------------------------------------------------------------
            DBioB(i,k,iBen) =DBioB(i,k,iBen)                            &
     &           -cff0*BenPred*dtdays*BioB(i,k,iBen)**2
            DBioB(i,k,iBenDet) =DBioB(i,k,iBenDet)                      &
     &           +cff0*BenPred*dtdays*BioB(i,k,iBen)**2
!-----------------------------------------------------------------------
!  (Det -> NH4) temperature dependent  
!-----------------------------------------------------------------------
            PON = BioB(i,k,iBenDet)*xi/Hz(i,j,k)  !Benthic Particulate organic nitrogen
!           cff5= q10**((Temp1-T0ben)/10.0_r8)
            cff1=Pv0*exp(PvT*Temp1)*PON       !-Kawamiya 2000 
!           cff1=Pv0*cff5*PON
! 
!  Arhonditsis
! 
            cff2 = Bio(i,k,iNH4)/(KNH4Nit +Bio(i,k,iNH4))
            DBioB(i,k,iBenDet) =DBioB(i,k,iBenDet)                      &
     &           - (cff1/xi)*dtdays*Hz(i,j,k)
            DBio(i,k,iNH4)= DBio(i,k,iNH4)                              &
     &           + (cff1)*dtdays
          END DO
!
!==================================================================
! Vertical Sinking of Particles(Phyto and Det) 
!==================================================================
!
! Use Sinking Code adapted from J. Warner ROMS sediment code 
!   - this is similar to the approach taken in Fasham 
!     biology code but adapted to be used in subroutine
!   - Also addapted to return particle flux in and out of each
!     level - used in  sub model
!   - Incorporated Liz Dobbins zlimit to determine sinking rate
!     zlimit =1 for constant sinking - use for Phyto and Det
!     zlimit =-1*NcmaxZ for Neocalanus- sink rate is then attenuated
!     as max sinking depth is approached
! 
! G.Gibson  July 2008        
! K.Kearney March 2016: updated to correct mass-accumulation bug
					! Note: all sinking and rising uses the subroutine BioSink, 
					! which modifies Btmp and sinkout
					! Initialize temporary arrays to 0
					Btmp = 0    ! tracer profile, (1:n)
					Hztmp = 0   ! layer thickness, (1:n)
					zwtmp = 0   ! layer edges, (0:n)
					sinkout2 = 0 ! loss out of bottom cell (1)
					! Small phytoplankton: sinks, and 79% of what sinks out of the 
					! bottom goes to benthic detritus 
					DO i=Istr,Iend
					  DO k = 1,N(ng)
							Btmp(k) = Bio(i,k,iPhS)
					    Hztmp(k) = Hz(i,j,k)
					  END DO
					  DO k = 0,N(ng)
					    zwtmp(k) = z_w(i,j,k)
					  END DO
					  call BioSink(N(ng), Btmp, wPhS, Hztmp, dtdays, zwtmp, 1.0_r8, sinkout2)
					  DO k = 1,N(ng)
							DBio(i,k,iPhS) = DBio(i,k,iPhS) + (Btmp(k) - Bio(i,k,iPhS))
					  END DO
						DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8
					END DO
          ! Large phytoplankton: sinks, and 79% of what sinks out of the 
          ! bottom goes to benthic detritus 
					DO i=Istr,Iend
					  DO k = 1,N(ng)
							Btmp(k) = Bio(i,k,iPhL)
					    Hztmp(k) = Hz(i,j,k)
					  END DO
					  DO k = 0,N(ng)
					    zwtmp(k) = z_w(i,j,k)
					  END DO
					  call BioSink(N(ng), Btmp, wPhL, Hztmp, dtdays, zwtmp, 1.0_r8, sinkout2)
					  DO k = 1,N(ng)
							DBio(i,k,iPhL) = DBio(i,k,iPhL) + (Btmp(k) - Bio(i,k,iPhL))
					  END DO
						DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8
					END DO
          ! Slow-sinking detritus: sinks, and 79% of what sinks out of the 
          ! bottom goes to benthic detritus 
					DO i=Istr,Iend
					  DO k = 1,N(ng)
							Btmp(k) = Bio(i,k,iDet)
					    Hztmp(k) = Hz(i,j,k)
					  END DO
					  DO k = 0,N(ng)
					    zwtmp(k) = z_w(i,j,k)
					  END DO
					  call BioSink(N(ng), Btmp, wDet, Hztmp, dtdays, zwtmp, 1.0_r8, sinkout2)
					  DO k = 1,N(ng)
							DBio(i,k,iDet) = DBio(i,k,iDet) + (Btmp(k) - Bio(i,k,iDet))
					  END DO
						DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8
					END DO
          ! Fast-sinking detritus: sinks, and 79% of what sinks out of the 
          ! bottom goes to benthic detritus 
					DO i=Istr,Iend
					  DO k = 1,N(ng)
							Btmp(k) = Bio(i,k,iDetF)
					    Hztmp(k) = Hz(i,j,k)
					  END DO
					  DO k = 0,N(ng)
					    zwtmp(k) = z_w(i,j,k)
					  END DO
					  call BioSink(N(ng), Btmp, wDetF, Hztmp, dtdays, zwtmp, 1.0_r8, sinkout2)
					  DO k = 1,N(ng)
							DBio(i,k,iDetF) = DBio(i,k,iDetF) + (Btmp(k) - Bio(i,k,iDetF))
					  END DO
						DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8
					END DO
          ! On-shelf large copepods (NCaS i.e. CM): Move up and down 
          ! based on dates set in input file.  If water is deeper than 
          ! 400 m, stop at 400m.  If shallower, when they hit the bottom, 
          ! the biomass is transferred to benthic detritus.  Uses zlimit 
          ! to enforce 400-m limit when going down, and to prevent rising 
          ! out of the surface layer on the way up
          DO i=Istr,Iend
            if (downwardCM) then
              DO k = 1,N(ng)
								Btmp(k) = Bio(i,k,iNCaS)
                Hztmp(k) = Hz(i,j,k)
              END DO
              DO k = 0,N(ng)
                zwtmp(k) = z_w(i,j,k)
              END DO
              call BioSink(N(ng), Btmp, wNCsink, Hztmp, dtdays, zwtmp, -400.0_r8, sinkout2)
              DO k = 1,N(ng)
								DBio(i,k,iNCaS) = DBio(i,k,iNCaS) + (Btmp(k) - Bio(i,k,iNCaS))
              END DO
							DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8
            else if (upwardCM) then
              DO k = 1,N(ng)
                Btmp(k) = Bio(i,N(ng)+1-k,iNCaS) ! flip
                Hztmp(k) = Hz(i,j,N(ng)+1-k)     ! flip
              END DO
              DO k = 0,N(ng)
                zwtmp(k) = z_w(i,j,0) - z_w(i,j,N(ng)-k) ! make surface the bottom
              END DO
              call BioSink(N(ng), Btmp, wNCrise, Hztmp, dtdays, zwtmp, z_w(i,j,0)-z_w(i,j,N(ng))+eps, sinkout2)
              DO k = 1,N(ng)
								DBio(i,k,iNCaS) = DBio(i,k,iNCaS) + (Btmp(N(ng)+1-k) - Bio(i,k,iNCaS)) ! flip back
              END DO
            end if
          END DO
          ! Off-shelf large copepods (NCaO i.e. NC): Same as above, but
          ! with different dates
          DO i=Istr,Iend
            if (downward) then
              DO k = 1,N(ng)
                Btmp(k) = Bio(i,k,iNCaO)
                Hztmp(k) = Hz(i,j,k)
              END DO
              DO k = 0,N(ng)
                zwtmp(k) = z_w(i,j,k)
              END DO
              call BioSink(N(ng), Btmp, wNCsink, Hztmp, dtdays, zwtmp, -400.0_r8, sinkout2)
              DO k = 1,N(ng)
                DBio(i,k,iNCaO) = DBio(i,k,iNCaO) + (Btmp(k) - Bio(i,k,iNCaO))
              END DO
              DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8
            else if (upward) then
              DO k = 1,N(ng)
                Btmp(k) = Bio(i,N(ng)+1-k,iNCaO) ! flip
                Hztmp(k) = Hz(i,j,N(ng)+1-k)     ! flip
              END DO
              DO k = 0,N(ng)
                zwtmp(k) = z_w(i,j,0) - z_w(i,j,N(ng)-k) ! make surface the bottom
              END DO
              call BioSink(N(ng), Btmp, wNCrise, Hztmp, dtdays, zwtmp, z_w(i,j,0)-z_w(i,j,N(ng))+eps, sinkout2)
              DO k = 1,N(ng)
                DBio(i,k,iNCaO) = DBio(i,k,iNCaS) + (Btmp(N(ng)+1-k) - Bio(i,k,iNCaO)) ! flip back
              END DO
            end if
          END DO
!=======================================================================
!Ice Sub Model
!=======================================================================
! 
!=======================================================================
! Update Bio array
!=======================================================================
          DO i=Istr,Iend
            DO itrc=1,NBT
              ibio=idbio(itrc)
              DO k=1,N(ng)
                Bio(i,k,ibio)=Bio(i,k,ibio)+DBio(i,k,ibio)
              END DO
            END DO
          END DO
          DO itrc=1,NBEN
            ibioB=idben(itrc)
            DO k=1,NBL(ng)
              DO i=Istr,Iend
                BioB(i,k,ibioB)=BioB(i,k,ibioB)+DBioB(i,k,ibioB)
              END DO
            END DO
          END DO
        END DO ITER_LOOP
!=======================================================================
!  Update global tracer variables (m Tunits).
!=======================================================================
        DO i=Istr,Iend
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              t(i,j,k,nnew,ibio)=MAX(t(i,j,k,nnew,ibio)+                &
      &                               (Bio(i,k,ibio)-Bio_bak(i,k,ibio)) &
      &                               *Hz(i,j,k)                        &
      &                               ,0.0_r8)
              t(i,j,k,nnew,iMZS)=0.0_r8
            END DO
          END DO
        END DO
        DO itrc=1,NBEN
          ibioB=idben(itrc)
          DO k=1,NBL(ng)
            DO i=Istr,Iend
!  check indexing here. think it ok.  
              bt(i,j,k,nnew,ibioB)=MAX(bt(i,j,k,nstp,ibioB)+            &
     &                           (BioB(i,k,ibioB)-Bio_bakB(i,k,ibioB))  &
     &                           ,0.0_r8) 
            END DO
          END DO
        END DO
      END DO J_LOOP
!
!  Apply periodic boundary conditions.
!
      DO itrc=1,NBT
        ibio=idbio(itrc)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,ibio))
      END DO
      RETURN
      END SUBROUTINE biology_tile
!=====================================================================
! BIOSINK_1  particle sinking subroutine After J. Warner sed sink code
! G. Gibson July 2008
!=====================================================================
      subroutine BIOSINK_1(ng,wBio,Bio,sinkIN,sinkOUT,HzL,dtdays,z_wL,  &
     &                       zlimit,LBi,UBi,IminS, ImaxS)
      USE mod_param
!
      implicit none
!
      integer, intent(in)     :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in)    :: wBio
      real(r8), intent(in)    :: zlimit
      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(inout) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
      real(r8), intent(in)    :: dtdays
      real(r8), intent(out) :: sinkIN(IminS:ImaxS,N(ng))
      real(r8), intent(out) :: sinkOUT(IminS:ImaxS,N(ng))
      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff
      real(r8):: dBio(0:N(ng)), wBiod(LBi:UBi,0:N(ng))
      real(r8) :: FC(IminS:ImaxS,0:N(ng))
      real(r8) :: Hz_inv(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv2(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv3(IminS:ImaxS,N(ng))
      integer :: ksource(IminS:ImaxS,N(ng))
      real(r8) :: qR(IminS:ImaxS,N(ng))
      real(r8) :: qL(IminS:ImaxS,N(ng))
      real(r8) :: WL(IminS:ImaxS,N(ng))
      real(r8) :: WR(IminS:ImaxS,N(ng))
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      IF ( zlimit .lt. 0 ) THEN
       DO k=0,N(ng)
         DO i=LBi, UBi
          IF ( z_wL(i,k) .ge. zlimit ) THEN
            wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2/        &
     &        (zlimit/2)**2 )
          ELSE
            wBiod(i,k) = 0.0_r8
          END IF
        END DO
       END DO
      ELSE
        DO k=0,N(ng)
          DO i=LBi, UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
! 
!
!  Compute inverse thickness to avoid repeated divisions.
!
      DO k=1,N(ng)
       DO i=LBi,UBi
         !if (HzL(i,k)==0) then
         !   print *,i,k, HzL(i,k)
         !end if
         Hz_inv(i,k)=1.0_r8/HzL(i,k)
       END DO
      END DO
      DO k=1,N(ng)-1
       DO i=LBi,UBi
         Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=2,N(ng)-1
       DO i=LBi,UBi
         Hz_inv3(i,k)=1.0_r8/(HzL(i,k-1)+HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=1,N(ng)
         DO i=LBi,UBi
           qc(i,k)=Bio(i,k)
         END DO
      END DO
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
          FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
        END DO
      END DO
!     DO k=N(ng)-1,1,-1
!       DO i=LBi,UBi
!         print*,'LBi=',LBi,'UBi=',UBi
!         print*,'i=',i,'k=',k
!         if(i.le.UBi)THEN
!           FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
!         endif
!       END DO
!     END DO
      DO k=2,N(ng)-1
        DO i=LBi,UBi
          dltR=HzL(i,k)*FC(i,k)
          dltL=HzL(i,k)*FC(i,k-1)
          cff=HzL(i,k-1)+2.0_r8*HzL(i,k)+HzL(i,k+1)
          cffR=cff*FC(i,k)
          cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
          cff=(dltR-dltL)*Hz_inv3(i,k)
          dltR=dltR-cff*HzL(i,k+1)
          dltL=dltL+cff*HzL(i,k-1)
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
          WR(i,k)=(2.0_r8*dltR-dltL)**2
          WL(i,k)=(dltR-2.0_r8*dltL)**2
        END DO
      END DO
      cff=1.0E-14_r8
      DO k=2,N(ng)-2
         DO i=LBi,UBi
          dltL=MAX(cff,WL(i,k  ))
          dltR=MAX(cff,WR(i,k+1))
          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
          qL(i,k+1)=qR(i,k)
        END DO
      END DO
      DO i=LBi,UBi
        FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
        qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
        qL(i,N(ng))=qc(i,N(ng))         ! conditions
        qR(i,N(ng)-1)=qc(i,N(ng))
        qL(i,2)=qc(i,1)                 ! bottom grid boxes are
        qR(i,1)=qc(i,1)                 ! re-assumed to be
        qL(i,1)=qc(i,1)                 ! piecewise constant.
      END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
          dltR=qR(i,k)-qc(i,k)
          dltL=qc(i,k)-qL(i,k)
          cffR=2.0_r8*dltR
          cffL=2.0_r8*dltL
          IF ((dltR*dltL).lt.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
        END DO
      END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)). 
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
          cff=dtdays*ABS(wBiod(i,k))
          FC(i,k-1)=0.0_r8
          WL(i,k)=z_wL(i,k-1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO
      DO k=1,N(ng)
        DO ks=k,N(ng)-1
          DO i=LBi,UBi
            IF (WL(i,k).gt.z_wL(i,ks)) THEN
              ksource(i,k)=ks+1
              FC(i,k-1)=FC(i,k-1)+WR(i,ks)
            END IF
          END DO
        END DO
      END DO
!
!  Finalize computation of flux: add fractional part.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
          ks=ksource(i,k)
          cu=MIN(1.0_r8,(WL(i,k)-z_wL(i,ks-1))*Hz_inv(i,ks))
          FC(i,k-1)=FC(i,k-1)+                                          &
     &                  HzL(i,ks)*cu*                                   &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
! 
!  G.Gibson  - FC is the flux into the level 
!            - should be 0 at the surface
          if (k.eq.N(ng)) then
            FC(i,k)=0.0_r8
          endif
        END DO
      END DO
      DO k=1,N(ng)
        DO i=LBi,UBi
! 					
!  The Bio variables are now updated in the main subroutine
! 
!         Bio(i,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
          sinkIN(i,k)=FC(i,k) 
          sinkOUT(i,k)=FC(i,k-1) 
          if(k.eq.N(ng))THEN
            sinkIN(i,k)=0.0_r8
          endif
        END DO
      END DO
      RETURN
      END SUBROUTINE BIOSINK_1
!=====================================================================
! BIORISE_1  rising particle subroutine a reversal of J. Warner sed sink code
!      plus and attenuation of rise rate based on closeness to max sink depth
! G.Gibson July 2008
!=====================================================================
     subroutine BIORISE_1(ng,wBio,Bio,rOUT,rIN,HzL,dtdays,z_wL,         &
     &                       LBi,UBi,zlimit,IminS, ImaxS)
      USE mod_param
!
      implicit none
!
!     real(r8), intent(in)    :: z_w(0:N(ng))
      integer, intent(in)     :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in)    :: wBio
      real(r8), intent(in)    :: zlimit
      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(in) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
!     real(r8), intent(in)    :: HzL(N(ng))
      real(r8), intent(in)    :: dtdays
      real(r8), intent(out) :: rIN(IminS:ImaxS,N(ng))
      real(r8), intent(out) :: rOUT(IminS:ImaxS,N(ng))
      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff
      real(r8):: dBio(0:N(ng)), wBiod(LBi:UBi,0:N(ng))
      real(r8) :: FC(IminS:ImaxS,0:N(ng))
      real(r8) :: Hz_inv(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv2(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv3(IminS:ImaxS,N(ng))
      integer :: ksource(IminS:ImaxS,N(ng))
      real(r8) :: qR(IminS:ImaxS,N(ng))
      real(r8) :: qL(IminS:ImaxS,N(ng))
      real(r8) :: WL(IminS:ImaxS,N(ng))
      real(r8) :: WR(IminS:ImaxS,N(ng))
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      IF ( zlimit .lt. 0 ) THEN
       DO k=0,N(ng)
         DO i=LBi, UBi
          IF ( z_wL(i,k) .ge. zlimit ) THEN
            wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2/        &
     &                   (zlimit/2)**2 )
          ELSE
            wBiod(i,k) = 0.0_r8
          END IF
        END DO
       END DO
      ELSE
        DO k=0,N(ng)
          DO i=LBi, UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
! 
!
!  Compute inverse thickness to avoid repeated divisions.
!
      DO k=1,N(ng)
       DO i=LBi,UBi
         Hz_inv(i,k)=1.0_r8/HzL(i,k)
       END DO
      END DO
      DO k=1,N(ng)-1
       DO i=LBi,UBi
         Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=2,N(ng)-1
       DO i=LBi,UBi
         Hz_inv3(i,k)=1.0_r8/(HzL(i,k-1)+HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=1,N(ng)
         DO i=LBi,UBi
           qc(i,k)=Bio(i,k)
         END DO
      END DO
! !
! !  Reconstruct vertical profile of suspended sediment "qc" in terms
! !  of a set of parabolic segments within each grid box. Then, compute
! !  semi-Lagrangian flux due to sinking.
! !
      DO k=2,N(ng),1
        DO i=LBi,UBi
          FC(i,k)=(qc(i,k-1)-qc(i,k))*Hz_inv2(i,k)
        END DO
      END DO
      DO k=N(ng)-1,2,-1
         DO i=LBi,UBi
          dltR=HzL(i,k)*FC(i,k)
          dltL=HzL(i,k)*FC(i,k+1)
          cff=HzL(i,k+1)+2.0_r8*HzL(i,k)+HzL(i,k-1)
          cffR=cff*FC(i,k)
          cffL=cff*FC(i,k+1)
! !
! !  Apply PPM monotonicity constraint to prevent oscillations within the
! !  grid box.
! !
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
! !
! !  Compute right and left side values (qR,qL) of parabolic segments
! !  within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
! !
! !  NOTE: Although each parabolic segment is monotonic within its grid
! !        box, monotonicity of the whole profile is not guaranteed,
! !        because qL(k+1)-qR(k) may still have different sign than
! !        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
! !        are reconciled using WENO procedure.
! !
          cff=(dltR-dltL)*Hz_inv3(i,k)
          dltR=dltR-cff*HzL(i,k-1)
          dltL=dltL+cff*HzL(i,k+1)
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
          WR(i,k)=(2.0_r8*dltR-dltL)**2
          WL(i,k)=(dltR-2.0_r8*dltL)**2
        END DO
      END DO
      cff=1.0E-14_r8
      DO k=N(ng),2,-1
        DO i=LBi,UBi
          dltL=MAX(cff,WL(i,k  ))
          dltR=MAX(cff,WR(i,k-1))
          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k-1))/(dltR+dltL)
          qL(i,k-1)=qR(i,k)
        END DO
      END DO
      DO i=LBi,UBi
             FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
             qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
             qL(i,N(ng))=qc(i,N(ng))         ! conditions
             qR(i,N(ng)-1)=qc(i,N(ng))
             qL(i,2)=qc(i,1)                 ! bottom grid boxes are
             qR(i,1)=qc(i,1)                 ! re-assumed to be
             qL(i,1)=qc(i,1)                 ! piecewise constant.
      END DO
! !
! !  Apply monotonicity constraint again, since the reconciled interfacial
! !  values may cause a non-monotonic behavior of the parabolic segments
! !  inside the grid box.
! !
! 
      DO k=N(ng),1,-1
        DO i=LBi,UBi
          dltR=qR(i,k)-qc(i,k)
          dltL=qc(i,k)-qL(i,k)
          cffR=2.0_r8*dltR
          cffL=2.0_r8*dltL
          IF ((dltR*dltL).lt.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
         END IF
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
        END DO
      END DO
!       
! !   
! !  After this moment reconstruction is considered complete. The next
! !  stage is to compute vertical advective fluxes, FC. It is expected
! !  that sinking may occurs relatively fast, the algorithm is designed
! !  to be free of CFL criterion, which is achieved by allowing
! !  integration bounds for semi-Lagrangian advective flux to use as
! !  many grid boxes in upstream direction as necessary.
! !
! !  In the two code segments below, WL is the z-coordinate of the
! !  departure point for grid box interface z_w with the same indices;
! !  FC is the finite volume flux; ksource(:,k) is index of vertical
! !  grid box which contains the departure point (restricted by N(ng)). 
! !  During the search: also add in content of whole grid boxes
! !  participating in FC.
! !
!            
!       
      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
          cff=dtdays*ABS(wBiod(i,k))
          FC(i,k)=0.0_r8
          WL(i,k)=z_wL(i,k+1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO
      DO k=N(ng)-1,1,-1
        DO ks=N(ng),k,1
          DO i=LBi,UBi
            IF (WL(i,k).gt.z_wL(i,ks)) THEN
              ksource(i,k)=ks-1
              FC(i,k+1)=FC(i,k+1)+WR(i,ks)
            END IF
          END DO
        END DO
      END DO
! 
!        
! !
! !  Finalize computation of flux: add fractional part.
! !
      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
          ks=ksource(i,k)
          cu=MIN(1.0_r8,(WL(i,k)-z_wL(i,ks+1))*Hz_inv(i,ks))
          FC(i,k+1)=FC(i,k+1)+                                          &
      &                  HzL(i,ks)*cu*                                  &
      &                  (qL(i,ks)+                                     &
      &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-               &
      &                       (1.5_r8-cu)*                              &
      &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
        END DO
     END DO
     DO k=1,N(ng),1
       DO i=LBi,UBi
         rIN(i,k)=MAX(0.0_r8,FC(i,k)) 
         if (k.eq.1) THEN
           rIN(i,1)=0.0_r8
         endif
       END DO
     END DO
     DO k=1,N(ng)-1,1
       DO i=LBi,UBi
         rOUT(i,k)=MAX(0.0_r8,FC(i,k+1))
         rOUT(i,N(ng))=0.0_r8
       END DO
     END DO
      RETURN
      END SUBROUTINE BIORISE_1
!=====================================================================
! BIOSINK_2      -original goanpz sink code (C.Lewis+ E. Dobbins)
!=====================================================================
!=====================================================================
!     BioSink: New subroutine to replace BIOSINK_1, BIORISE_1
!=====================================================================
      Subroutine BioSink(n, Btmp, wBio, HzL, dtdays, z_wL, zlimit, sinkout)
      !------------------------------------------------------------------
      ! Computes redistribution of a tracer in the water column due to 
      ! sinking.  After J. Warner sed sink code.
      !
      !   n       = # layers
      !   Btmp    = n x 1 array, tracer concentration (bottom to top)
      !   wBio    = sinking rate (m d^-1)
      !   HzL     = n x 1 array, thickness of layers (m)
      !   dtdays  = time step (d)
      !   z_wL    = n+1 x 1 array, depth of layer edges (m, negative)
      !   zlimit  = maximum depth for sinking (m, negative)
      !   sinkout = amount lost out of bottom cell (concentration)
      !
      ! Modifies Btmp and sinkout in the calling program
      !------------------------------------------------------------------ 
      USE mod_param
      implicit none
      integer,  intent(in) :: n
      real(r8), intent(in) :: wBio
      real(r8), intent(in) :: zlimit
      real(r8), intent(in) :: dtdays
      real(r8), intent(in) :: z_wL(0:n)
      real(r8), intent(inout) :: Btmp(n)
      real(r8), intent(in) :: HzL(n)
      real(r8), intent(inout) :: sinkout
      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff
      real(r8):: dBio(0:n), wBiod(0:n)
      real(r8) :: FC(0:n)
      real(r8) :: Hz_inv(n)
      real(r8) :: Hz_inv2(n)
      real(r8) :: Hz_inv3(n)
      integer :: ksource(n)
      real(r8) :: qR(n)
      real(r8) :: qL(n)
      real(r8) :: WL(n)
      real(r8) :: WR(n)
      real(r8), dimension(n) :: qc
      ! Modify sinking rates as necessary 
      IF ( zlimit .lt. 0 ) THEN
        DO k=0,n
          IF ( z_wL(k) .ge. zlimit ) THEN
            wBiod(k) = wBio*exp(-1*(z_wL(k)-(zlimit/2))**2/(zlimit/2)**2)
          ELSE
            wBiod(k) = 0.0_r8
          END IF
        END DO
      ELSE
        DO k=0,n
          wBiod(k) = wBio
        END DO
      END IF
      ! Compute inverse thickness to avoid repeated divisions.
      DO k=1,n
        Hz_inv(k)=1.0_r8/HzL(k)
      END DO
      DO k=1,n-1
        Hz_inv2(k)=1.0_r8/(HzL(k)+HzL(k+1))
      END DO
      DO k=2,n-1
        Hz_inv3(k)=1.0_r8/(HzL(k-1)+HzL(k)+HzL(k+1))
      END DO
      DO k=1,n
        qc(k)=Btmp(k)
      END DO
      ! Reconstruct vertical profile of suspended sediment "qc" in terms
      ! of a set of parabolic segments within each grid box. Then, 
      ! compute semi-Lagrangian flux due to sinking.
      DO k=n-1,1,-1
          FC(k)=(qc(k+1)-qc(k))*Hz_inv2(k)
      END DO
      DO k=2,n-1
          dltR=HzL(k)*FC(k)
          dltL=HzL(k)*FC(k-1)
          cff=HzL(k-1)+2.0_r8*HzL(k)+HzL(k+1)
          cffR=cff*FC(k)
          cffL=cff*FC(k-1)
          ! Apply PPM monotonicity constraint to prevent oscillations 
          ! within the grid box.
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
          ! Compute right and left side values (qR,qL) of parabolic 
          ! segments within grid box Hz(k); (WR,WL) are measures of 
          ! quadratic variations.
          !
          !  NOTE: Although each parabolic segment is monotonic within 
          !        its grid box, monotonicity of the whole profile is not 
          !        guaranteed, because qL(k+1)-qR(k) may still have 
          !        different sign than qc(k+1)-qc(k).  This possibility 
          !        is excluded, after qL and qR are reconciled using WENO 
          !        procedure.
          cff=(dltR-dltL)*Hz_inv3(k)
          dltR=dltR-cff*HzL(k+1)
          dltL=dltL+cff*HzL(k-1)
          qR(k)=qc(k)+dltR
          qL(k)=qc(k)-dltL
          WR(k)=(2.0_r8*dltR-dltL)**2
          WL(k)=(dltR-2.0_r8*dltL)**2
      END DO
      cff=1.0E-14_r8
      DO k=2,n-2
        dltL=MAX(cff,WL(k  ))
        dltR=MAX(cff,WR(k+1))
        qR(k)=(dltR*qR(k)+dltL*qL(k+1))/(dltR+dltL)
        qL(k+1)=qR(k)
      END DO
      FC(n)=0.0_r8                  ! no-flux boundary condition
      qR(n)=qc(n)                   ! default strictly monotonic
      qL(n)=qc(n)                   ! conditions
      qR(n-1)=qc(n)
      qL(2)=qc(1)                   ! bottom grid boxes are
      qR(1)=qc(1)                   ! re-assumed to be
      qL(1)=qc(1)                   ! piecewise constant.
      ! Apply monotonicity constraint again, since the reconciled 
      ! interfacial values may cause a non-monotonic behavior of the 
      ! parabolic segments inside the grid box.
      DO k=1,n
          dltR=qR(k)-qc(k)
          dltL=qc(k)-qL(k)
          cffR=2.0_r8*dltR
          cffL=2.0_r8*dltL
          IF ((dltR*dltL).lt.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
          qR(k)=qc(k)+dltR
          qL(k)=qc(k)-dltL
      END DO
      ! After this moment reconstruction is considered complete. The next
      ! stage is to compute vertical advective fluxes, FC. It is expected
      ! that sinking may occurs relatively fast, the algorithm is designed
      ! to be free of CFL criterion, which is achieved by allowing
      ! integration bounds for semi-Lagrangian advective flux to use as
      ! many grid boxes in upstream direction as necessary.
      !
      ! In the two code segments below, WL is the z-coordinate of the
      ! departure point for grid box interface z_w with the same indices;
      ! FC is the finite volume flux; ksource(:,k) is index of vertical
      ! grid box which contains the departure point (restricted by N(ng)). 
      ! During the search: also add in content of whole grid boxes
      ! participating in FC.
      DO k=1,n
        cff=dtdays*ABS(wBiod(k))
        FC(k-1)=0.0_r8
        WL(k)=z_wL(k-1)+cff
        WR(k)=HzL(k)*qc(k)
        ksource(k)=k
      END DO
      DO k=1,n
        DO ks=k,n-1
          IF (WL(k).gt.z_wL(ks)) THEN
            ksource(k)=ks+1
            FC(k-1)=FC(k-1)+WR(ks)
          END IF
        END DO
      END DO
      ! Finalize computation of flux: add fractional part.
      DO k=1,n
        ks=ksource(k)
        cu=MIN(1.0_r8,(WL(k)-z_wL(ks-1))*Hz_inv(ks))
        FC(k-1)=FC(k-1)+                                                &
     &                HzL(ks)*cu*                                       &
     &                (qL(ks)+                                          &
     &                 cu*(0.5_r8*(qR(ks)-qL(ks))-                      &
     &                     (1.5_r8-cu)*                                 &
     &                     (qR(ks)+qL(ks)-2.0_r8*qc(ks))))
      END DO
      ! New profile of tracer Btmp (mass per volume)
      DO k=1,n
        Btmp(k)=qc(k)+(FC(k)-FC(k-1))*Hz_inv(k)
      END DO
      ! Amount lost out of bottom cell (mass per area)
      sinkout = FC(0)
      RETURN
      END SUBROUTINE BioSink
!=======================================================================
!    DetSINK 
!----------------
!====================================================================
        Function ComputeDensity(Temp1,Sal1)
!----------------------------------------------------------------
! Computes the water column density from salinity and temperature
! Returns sigma-t
!----------------------------------------------------------------
        USE mod_kinds
        Real(r8) ComputeDensity
        Real(r8) Temp1, Sal1
        Real(r8) Sig
        Sig = 999.842594 + 0.06793952 * Temp1
        Sig = Sig - 0.00909529 * Temp1 ** 2 +                           &
     &          0.0001001685 * Temp1 ** 3
        Sig = Sig - 0.000001120083 * Temp1 ** 4 +                       &
     &          0.000000006536332 * Temp1 ** 5
        Sig = Sig + 0.824493 * Sal1 - 0.0040899 * Temp1 * Sal1
        Sig = Sig + 0.000076438 * Temp1 ** 2 * Sal1 -                   &
     &          0.00000082467 * Temp1 ** 3 * Sal1
        Sig = Sig + 0.0000000053875 * Temp1 ** 4 * Sal1 -               &
     &          0.00572466 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00010227 * Temp1 * Sal1 ** (3 / 2) -              &
     &          0.0000016546 * Temp1 ** 2 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00048314 * Sal1 ** 2
         ComputeDensity = Sig - 1000
        End Function ComputeDensity
!===============================================================
        Function GetPhytoResp2(Temp1, Tref, KbmPh)
!------------------------------------------------------
! Computes the temperature correction for phytoplankton
! respiration according to Arhonditsis 2005.
!------------------------------------------------------
        USE mod_kinds
        Real(r8) GetPhytoResp2
        Real(r8) Temp1      !Temperature, passed
        Real(r8) Tref       !Reference temperature
        Real(r8) KbmPh      !Half saturation, temperature
        Real(r8) Resp       !Returned variable
        Resp = exp(KbmPh * (Temp1 - Tref))
        GetPhytoResp2 = Resp
        Return
        End Function GetPhytoResp2
!=====================================================================
      FUNCTION GetLightLimIronSml(alphaPh, PAR1, Pmax1,                 &
     &     CrChlRatio1,IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none
        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIronSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)        &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIronSml = LightLim
      END FUNCTION GetLightLimIronSml
!=====================================================================
      FUNCTION GetLightLimSml(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none
        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)        &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLimSml = LightLim
      END FUNCTION GetLightLimSml
!=====================================================================
      FUNCTION GetLightLimIron(alphaPh, PAR1, Pmax1, CrChlRatio1,       &
     &     IronLim1, ParMax)
!------------------------------------------------------------------
! Light lim with varying alpha. Works with iron limitation. Alph is
! a function of the surface light intensity.
!------------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none
       Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
       Real(r8) GetLightLimIron
       Real(r8) Alpha,LightLim,OffSet,ParMax
       !Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
       !if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = alphaPh
       end if
       LightLim = TANH( Alpha * Par1/Pmax1/CrChlRatio1/IronLim1)
       GetLightLimIron = LightLim
      END FUNCTION GetLightLimIron
!=======================================================================
      FUNCTION GetLightLim(alphaPh, PAR1, Pmax1, CrChlRatio1, ParMax)
!-----------------------------------------------------------------
! Generates a light lim with varying alphaPh without iron
!-----------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none
       Real(r8) :: alphaPh,Pmax1,PAR1,CrChlRatio1
       Real(r8) GetLightLim
       Real(r8) Alpha,LightLim,OffSet,ParMax
!       Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
!       if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = alphaPh
       end if
       LightLim = TANH(alphaPh * PAR1/Pmax1/CrChlRatio1)
       GetLightLim = LightLim
      END FUNCTION GetLightLim
!==============================================================
        Function GetCopepodResp(Temp1,respVal,ktbm,Tref)
!--------------------------------------------------------------
! Computes copepod respiration according to Arhonditsis (2005).
!--------------------------------------------------------------
        USE mod_kinds
        USE mod_param
        Real(r8) GetCopepodResp
        real(r8) :: respVal
        Real(r8) Temp1       !Passed variable
!        Real(r8) :: bm  = 0.04   !Basal metabolic rate day**-1
        Real(r8) :: ktbm  != 0.05 !Temperature response degrees C**-1
        Real(r8) :: Tref  != 20   !Reference temperature degrees C
        Real(r8) Resp        !Returned variable
        Resp = respVal * exp(ktbm * (Temp1 - Tref))
        GetCopepodResp = Resp
        Return
        End Function GetCopepodResp
!==============================================================
        Function GetJelResp(Temp1,bmJ,ktbmJ,TrefJ)
        USE mod_param
        implicit none
!
        real(r8) :: GetJelResp
        real(r8) :: Temp1          !Passed variable
        real(r8) :: bmJ   != 0.04   !Basal metabolic rate day**-1
        real(r8) :: ktbmJ != 0.05 !Temperature response degrees C**-1
        real(r8) :: TrefJ != 20   !Reference temperature degrees C
        real(r8) :: Resp           !Returned variable
        Resp = bmJ * exp(ktbmJ * (Temp1 - TrefJ))
        GetJelResp = Resp
        Return
        End Function GetJelResp 
!=================================================================
        Function GetBasalMetabolism(respPh,kfePh,Iron1)
!---------------------------------------------------------
! Computes an iron correction for the basal metabolism for
! the phytoplankton respiration calculation
!---------------------------------------------------------
        USE mod_kinds
        Real(r8) GetBasalMetabolism
        Real(r8) Iron1     !Iron concentration
        Real(r8) kfePh     !Half saturation for iron
        Real(r8) respPh    !Phytoplankton uncorrected basal metabolism
        Real(r8) BaseMet   !Phytoplankton basal metabolism
        BaseMet = Iron1/(kfePh + Iron1)
        BaseMet = BaseMet * ((kfePh +2)/2)
        GetBasalMetabolism = BaseMet * respPh
        Return
        End Function GetBasalMetabolism
!=========================================================
        Function GetNitrif(Temp1,Dep1,NH4R)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrif
                               !--------------------------------
        Real(r8) Temp1, Dep1, NH4R !Passed variables
        Real(r8) NH4conv  /14/     !mg N/mol N
        Real(r8) KNH4Nit  /0.08/   !Half Sat Con mg N/m3/day
        Real(r8) KTNitr   /0.002/  !Temperature responce dec C^2
        Real(r8) ToptNtr  /28/     !Optimum nitrification temp
        Real(r8) Zox      /20/     !50% nitrification depth
        Real(r8) Nexp     /6/      !Exponent to adjust profile shape
        Real(r8) NitrMax  /0.011/  !Maximum nitrification (mM/m3/d
                               !--------------------------------
        Real(r8) Nitr, DepCor
        NH4R = NH4R * NH4conv
        Nitr = NH4R/(KNH4Nit + NH4R)
        Nitr = Nitr * exp(-KTNitr*(Temp1 - ToptNtr)**2)
        DepCor = (Dep1**Nexp)/( (Zox**Nexp) + Dep1**Nexp)
        Nitr = (Nitr * DepCor) * NitrMax
        GetNitrif = Nitr
        Return
        End Function GetNitrif
!========================================================================
        Function GetNitrif2(Temp1,PAR1,NH4R)
        !---------------------------------------------------------
        !Computes nitrificaton from Kawamiya with light correction
        !from Fennel; Kawamiya (2000), Fennel (2006)
        !---------------------------------------------------------
        USE mod_kinds
        USE mod_param
         Real(r8) GetNitrif2
                               !--------------------------------
        Real(r8) :: Temp1, NH4R       !Passed variables
        Real(r8) :: I0 = 0.0095     !Threshold,light inhibition, W m-2
        Real(r8) :: KI = 4.0        !Half Saturation light intensity, W m-2
        Real(r8) :: KN0 = 0.03       !Nitrification at 0 deg C, day-1
        Real(r8) :: KNT = 0.0693     !Temperature coefficient
        Real(r8) :: ParW              !Par in watts
        Real(r8) :: NitrMax           !Maximum nitrification
        Real(r8) :: Nitr              !Nitrification
                               !---------------------------------
        Real(r8) :: cff1, PAR1
        !-----------------------------------
        !Temperature dependent nitrification
        !-----------------------------------
        KI = 1.5
        KN0 = 0.15
        !KNT = 0.07
        NitrMax = (KN0*Exp(KNT*Temp1))*NH4R
        !-----------------------------------
        !Convert PAR in E m-2 d-1 to W day-1
        !-----------------------------------
        ParW = PAR1/0.394848_r8
        !---------------------------------
        !Light correction of nitrification
        !---------------------------------
        cff1 = (ParW-I0)/(KI+ParW-I0)
        Nitr = NitrMax*(1-MAX(0.0_r8,cff1))
        GetNitrif2 = Nitr
        Return
        End Function GetNitrif2
!-----------------------------------------------------------  
  Function GetNitrifLight(Par1,tI0,KI)
   USE mod_kinds
        Real(r8) GetNitrifLight
        Real(r8) :: Par1               !--------------------------------
        Real(r8) :: tI0                 !Threshold,light inhibition, W m-2
        Real(r8) :: KI                 !Half Saturation light intensity, W m-2
        Real(r8) :: ParW               !Par in watts
        Real(r8) :: cff1
      ParW = Par1/0.394848_r8 !convert PAR back to watts
            cff1 = (ParW-tI0)/(KI+ParW-tI0)
      GetNitrifLight=(1-MAX(0.0_r8,cff1))
        Return
        End Function GetNitrifLight
!-----------------------------------------------
   Function GetNitrifMaxK(Temp1,KnT,Nitr0)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrifMaxK
        Real(r8) :: Temp1         !--------------------------------
        Real(r8) :: KnT           !Passed variables
        Real(r8) :: Nitr0         !--------------------------------
  GetNitrifMaxK=Nitr0*exp(KnT*Temp1)
        Return
        End Function GetNitrifMaxK
!-----------------------------------------------------------  
  Function GetNitrifMaxA(Nitr0, ktntr,Temp1,ToptNtr)
   USE mod_kinds
        Real(r8) GetNitrifMaxA
        Real(r8) :: Temp1         
        Real(r8) :: Ktntr           !--------------------------------
        Real(r8) :: Nitr0           !Passed variables
        Real(r8) :: ToptNtr         !--------------------------------
  GetNitrifMaxA=Nitr0*exp(-ktntr*(Temp1 - ToptNtr)**2)
        Return
        End Function GetNitrifMaxA
!-----------------------------------------------------------    
!=====================================================================
      FUNCTION GetLightLimIron2(alphaPh, PAR1, Pmax1, CrChlRatio1,      &
     &     IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none
        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIron2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)        &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIron2 = LightLim
      END FUNCTION GetLightLimIron2
!=====================================================================
      FUNCTION GetLightLim2(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none
        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLim2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)        &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLim2 = LightLim
      END FUNCTION GetLightLim2
!===================================================================
     FUNCTION ComputeStability(ng,z_wL,Dens)
!--------------------------------------------------------
! Computes the stability parameter from density and depth
! Simpson stability parameter.
!--------------------------------------------------------
     USE mod_param
!
     implicit none
          integer, intent(in) ::  ng
     real(r8), intent(in) :: z_wL(0:N(ng))
     real(r8), intent(in) :: Dens(1:N(ng))
     Integer k,indx
     real ComputeStability
     real dZ,SumSigma,SumDep,RMean,Z1,Dep1,Dep2
     !SumSigma = Dens(30)
     !SumSigma = 0
     SumDep = 0
     indx = N(ng)     !30
     SumDep = z_wL(indx) * (-1.0_r8)
     SumSigma = Dens(indx) * SumDep
     indx = indx - 1
    DO k = 1, N(ng) - 1
       IF (z_wL(indx) .ge. -140) THEN
          dZ = z_wL(indx + 1) - z_wL(indx)
          !Print *, '**',Dens(k),z_wL(k),SumSigma,SumDep,dZ,k
          SumSigma = SumSigma + 0.5 *                                   &
    &                (Dens(indx) + Dens(indx+1)) * dZ
          SumDep = SumDep + dZ
       END IF
       indx = indx - 1
     END DO
     !if (SumDep == 0) then
     !    print *,ng,z_wL(indx),Dens(indx)
     !end if
     RMean = SumSigma/SumDep
     !print *, '&&&&%%%%',RMean, SumDep, SumSigma,z_wL(1)
     !Print *, z_wL(30), Dens(1),Dens(30)
     indx = N(ng)-1          !29
     Z1 = 0.5 * (z_wL(indx))
     dZ = z_wL(indx)
     SumSigma = Z1 * (RMean - Dens(indx - 1)) * dZ
     !print *, 'xxx',SumSigma
     DO k = 1, N(ng) - 1
       IF (z_wL(indx) .ge. -140) THEN
          Dep1 = z_wL(indx) * (-1.0_r8)
          Dep2 = z_wL(indx -1) * (-1.0_r8)
          Z1 = 0.5 * (Dep1 + Dep2)
          dZ = (Dep2 - Dep1)
!         print *, '&&**',z_wL(indx),z_wL(indx + 1),dZ,Z1
!        print *, '&&**',SumSigma,z_wL(k),dens(k),Z1*(Rmean -           &
!     &   0.5*(Dens(k) + Dens(k+1))) * dZ, dZ,rmean
          SumSigma = SumSigma + Z1*(Rmean -                             &
    &         0.5*(Dens(indx) + Dens(indx + 1))) * dZ
       END IF
       indx = indx - 1
     END DO
     !Z1 = 0.5 * z_wL(N(ng))
     !dZ = z_wL(N(Ng))
     !SumSigma = SumSigma + Z1*(RMean - Dens(N(ng))) * dZ
     ComputeStability = -9.8 * SumSigma / SumDep
     END FUNCTION ComputeStability
      END MODULE biology_mod
