      MODULE set_avg2_mod
!
!svn $Id: set_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC :: set_avg2
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg2 (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_average2
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
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
      CALL wclock_on (ng, iNLM, 5)
      CALL set_avg2_tile (ng, tile,                                     &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp(ng),                                      &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % rho,                               &
     &                   MIXING(ng) % hsbl,                             &
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % lhflx,                            &
     &                   FORCES(ng) % shflx,                            &
     &                   FORCES(ng) % lrflx,                            &
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
     &                   FORCES(ng) % evap,                             &
     &                   FORCES(ng) % rain,                             &
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   OCEAN(ng) % zeta,                              &
     &                   AVERAGE2(ng) % avgu3d,                         &
     &                   AVERAGE2(ng) % avgv3d,                         &
     &                   AVERAGE2(ng) % avgt,                           &
     &                   AVERAGE2(ng) % avgrho,                         &
     &                   AVERAGE2(ng) % avghsbl,                        &
     &                   AVERAGE2(ng) % avgstf,                         &
     &                   AVERAGE2(ng) % avgswf,                         &
     &                   AVERAGE2(ng) % avglhf,                         &
     &                   AVERAGE2(ng) % avgshf,                         &
     &                   AVERAGE2(ng) % avglrf,                         &
     &                   AVERAGE2(ng) % avguwind,                       &
     &                   AVERAGE2(ng) % avgvwind,                       &
     &                   AVERAGE2(ng) % avgevap,                        &
     &                   AVERAGE2(ng) % avgrain,                        &
     &                   AVERAGE2(ng) % avgsrf,                         &
     &                   AVERAGE2(ng) % avgsus,                         &
     &                   AVERAGE2(ng) % avgsvs,                         &
     &                   AVERAGE2(ng) % avgbus,                         &
     &                   AVERAGE2(ng) % avgbvs,                         &
     &                   AVERAGE2(ng) % avgu2d,                         &
     &                   AVERAGE2(ng) % avgv2d,                         &
     &                   AVERAGE2(ng) % avgzeta)
      CALL wclock_off (ng, iNLM, 5)
      RETURN
      END SUBROUTINE set_avg2
!
!***********************************************************************
      SUBROUTINE set_avg2_tile (ng, tile,                               &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Kout,                                    &
     &                         Nout,                                    &
     &                         pm, pn,                                  &
     &                         u, v,                                    &
     &                         t,                                       &
     &                         rho,                                     &
     &                         hsbl,                                    &
     &                         stflx,                                   &
     &                         lhflx, shflx, lrflx, Uwind, Vwind,       &
     &                         evap, rain,                              &
     &                         srflx,                                   &
     &                         sustr, svstr, bustr, bvstr,              &
     &                         ubar, vbar,                              &
     &                         zeta,                                    &
     &                         avgu3d, avgv3d,                          &
     &                         avgt,                                    &
     &                         avgrho,                                  &
     &                         avghsbl,                                 &
     &                         avgstf, avgswf,                          &
     &                         avglhf, avgshf, avglrf,                  &
     &                         avguwind, avgvwind,                      &
     &                         avgevap, avgrain,                        &
     &                         avgsrf,                                  &
     &                         avgsus, avgsvs, avgbus, avgbvs,          &
     &                         avgu2d, avgv2d,                          &
     &                         avgzeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
      integer, intent(in) :: Nout
!
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: lhflx(LBi:,LBj:)
      real(r8), intent(in) :: shflx(LBi:,LBj:)
      real(r8), intent(in) :: lrflx(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(in) :: evap(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgu3d(LBi:,LBj:)
      real(r8), intent(inout) :: avgv3d(LBi:,LBj:)
      real(r8), intent(inout) :: avgt(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgrho(LBi:,LBj:)
      real(r8), intent(inout) :: avghsbl(LBi:,LBj:)
      real(r8), intent(inout) :: avgstf(LBi:,LBj:)
      real(r8), intent(inout) :: avgswf(LBi:,LBj:)
      real(r8), intent(inout) :: avglhf(LBi:,LBj:)
      real(r8), intent(inout) :: avgshf(LBi:,LBj:)
      real(r8), intent(inout) :: avglrf(LBi:,LBj:)
      real(r8), intent(inout) :: avguwind(LBi:,LBj:)
      real(r8), intent(inout) :: avgvwind(LBi:,LBj:)
      real(r8), intent(inout) :: avgevap(LBi:,LBj:)
      real(r8), intent(inout) :: avgrain(LBi:,LBj:)
      real(r8), intent(inout) :: avgsrf(LBi:,LBj:)
      real(r8), intent(inout) :: avgsus(LBi:,LBj:)
      real(r8), intent(inout) :: avgsvs(LBi:,LBj:)
      real(r8), intent(inout) :: avgbus(LBi:,LBj:)
      real(r8), intent(inout) :: avgbvs(LBi:,LBj:)
      real(r8), intent(inout) :: avgu2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgv2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgzeta(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, itrc, j
      real(r8) :: fac, fac1
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
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG2(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG2(ng)).and.                                &
     &     (MOD(iic(ng)-1,nAVG2(ng)).eq.1)).or.                         &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=zeta(i,j,Kout)
            avgu2d (i,j)=ubar(i,j,Kout)
            avgv2d (i,j)=vbar(i,j,Kout)
            avghsbl(i,j)=hsbl(i,j)
            avgstf(i,j)=stflx(i,j,itemp)
            avgswf(i,j)=stflx(i,j,isalt)
            avglhf(i,j)=lhflx(i,j)
            avgshf(i,j)=shflx(i,j)
            avglrf(i,j)=lrflx(i,j)
            avguwind(i,j)=Uwind(i,j)
            avgvwind(i,j)=Vwind(i,j)
            avgevap(i,j)=evap(i,j)
            avgrain(i,j)=rain(i,j)
            avgsrf(i,j)=srflx(i,j)
            avgsus(i,j)=sustr(i,j)
            avgsvs(i,j)=svstr(i,j)
            avgbus(i,j)=bustr(i,j)
            avgbvs(i,j)=bvstr(i,j)
          END DO
        END DO
!
!  Initialize fields associated with 3D horizontal momentum.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgu3d(i,j)=u(i,j,n(ng),Nout)
            avgv3d(i,j)=v(i,j,n(ng),Nout)
            avgrho(i,j)=rho(i,j,n(ng))
          END DO
        END DO
!
!  Initialized fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgt(i,j,itrc)=t(i,j,N(ng),Nout,itrc)
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG2(ng)) THEN
!
!  Accumulate 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=avgzeta(i,j)+zeta(i,j,Kout)
            avgu2d (i,j)=avgu2d (i,j)+ubar(i,j,Kout)
            avgv2d (i,j)=avgv2d (i,j)+vbar(i,j,Kout)
            avghsbl(i,j)=avghsbl(i,j)+hsbl(i,j)
            avgstf(i,j)=avgstf(i,j)+stflx(i,j,itemp)
            avgswf(i,j)=avgswf(i,j)+stflx(i,j,isalt)
            avglhf(i,j)=avglhf(i,j)+lhflx(i,j)
            avgshf(i,j)=avgshf(i,j)+shflx(i,j)
            avglrf(i,j)=avglrf(i,j)+lrflx(i,j)
            avguwind(i,j)=avguwind(i,j)+Uwind(i,j)
            avgvwind(i,j)=avgvwind(i,j)+Vwind(i,j)
            avgevap(i,j)=avgevap(i,j)+evap(i,j)
            avgrain(i,j)=avgrain(i,j)+rain(i,j)
            avgsrf(i,j)=avgsrf(i,j)+srflx(i,j)
            avgsus(i,j)=avgsus(i,j)+sustr(i,j)
            avgsvs(i,j)=avgsvs(i,j)+svstr(i,j)
            avgbus(i,j)=avgbus(i,j)+bustr(i,j)
            avgbvs(i,j)=avgbvs(i,j)+bvstr(i,j)
          END DO
        END DO
!
!  Accumulate fields associated with 3D horizontal momentum.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgu3d(i,j)=avgu3d(i,j)+u(i,j,N(ng),Nout)
            avgv3d(i,j)=avgv3d(i,j)+v(i,j,N(ng),Nout)
            avgrho(i,j)=avgrho(i,j)+rho(i,j,N(ng))
          END DO
        END DO
!
!  Accumulate fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgt(i,j,itrc)=avgt(i,j,itrc)+t(i,j,N(ng),Nout,itrc)
            END DO
          END DO
        END DO
!
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG2(ng)).and.                                 &
     &    (MOD(iic(ng)-1,nAVG2(ng)).eq.0).and.                          &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG2(ng),r8)
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          AVG2time(ng)=AVG2time(ng)+REAL(nAVG2(ng),r8)*dt(ng)
        END IF
!
!  Process 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=fac*avgzeta(i,j)
            avgu2d (i,j)=fac*avgu2d (i,j)
            avgv2d (i,j)=fac*avgv2d (i,j)
            avghsbl(i,j)=fac*avghsbl(i,j)
            avgstf(i,j)=fac*avgstf(i,j)
            avgswf(i,j)=fac*avgswf(i,j)
            avglhf(i,j)=fac*avglhf(i,j)
            avgshf(i,j)=fac*avgshf(i,j)
            avglrf(i,j)=fac*avglrf(i,j)
            avguwind(i,j)=fac*avguwind(i,j)
            avgvwind(i,j)=fac*avgvwind(i,j)
            avgevap(i,j)=fac*avgevap(i,j)
            avgrain(i,j)=fac*avgrain(i,j)
            avgsrf(i,j)=fac*avgsrf(i,j)
            avgsus(i,j)=fac*avgsus(i,j)
            avgsvs(i,j)=fac*avgsvs(i,j)
            avgbus(i,j)=fac*avgbus(i,j)
            avgbvs(i,j)=fac*avgbvs(i,j)
          END DO
        END DO
!
!  Process fields associated with 3D horizontal momentum.
!
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgu3d(i,j)=fac*avgu3d(i,j)
              avgv3d(i,j)=fac*avgv3d(i,j)
              avgrho(i,j)=fac*avgrho(i,j)
          END DO
        END DO
!
!  Process fields associated with tracers.
!
        DO itrc=1,NT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,itrc)=fac*avgt(i,j,itrc)
            END DO
          END DO
        END DO
!-----------------------------------------------------------
! Computes average of time series
! For production, totals are needed so this is commented out
!-----------------------------------------------------------
!        DO itrc=1,NTS(ng)
!	    DO j-JstrR,JendR
!	      DO i-IstrR,IendR
!	        anvst(i,j,itrc)=fac*avgt(i,j,itrc)
!	    END DO
!	  END DO
!	END DO
!
      END IF
      RETURN
      END SUBROUTINE set_avg2_tile
      END MODULE set_avg2_mod
