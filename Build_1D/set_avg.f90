      MODULE set_avg_mod
!
!svn $Id: set_avg.F 956 2009-03-19 22:44:49Z kate $
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
      PUBLIC :: set_avg
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_average
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
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp(ng),                                      &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % pm, GRID(ng) % pn,                  &
     &                   OCEAN(ng) % u, OCEAN(ng) % v,                  &
     &                   OCEAN(ng) % W, OCEAN(ng) % wvel,               &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % bt,                                &
     &                   OCEAN(ng) % rho,                               &
     &                   MIXING(ng) % hsbl,                             &
     &                   MIXING(ng) % Akv,                              &
     &                   MIXING(ng) % Akt,                              &
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % lhflx,                            &
     &                   FORCES(ng) % shflx,                            &
     &                   FORCES(ng) % lrflx,                            &
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
     &                   FORCES(ng) % evap, FORCES(ng) % rain,          &
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % sustr, FORCES(ng) % svstr,        &
     &                   FORCES(ng) % bustr, FORCES(ng) % bvstr,        &
     &                   OCEAN(ng) % ubar, OCEAN(ng) % vbar,            &
     &                   OCEAN(ng) % zeta,                              &
     &                   AVERAGE(ng) % avgu3d,                          &
     &                   AVERAGE(ng) % avgv3d,                          &
     &                   AVERAGE(ng) % avgw3d,                          &
     &                   AVERAGE(ng) % avgwvel,                         &
     &                   AVERAGE(ng) % avgt,                            &
     &                   AVERAGE(ng) % avgbt,                           &
     &                   AVERAGE(ng) % avgrho,                          &
     &                   AVERAGE(ng) % avghsbl,                         &
     &                   AVERAGE(ng) % avgAKv,                          &
     &                   AVERAGE(ng) % avgAKt,                          &
     &                   AVERAGE(ng) % avgAKs,                          &
     &                   AVERAGE(ng) % avgstf,                          &
     &                   AVERAGE(ng) % avgswf,                          &
     &                   AVERAGE(ng) % avglhf,                          &
     &                   AVERAGE(ng) % avgshf,                          &
     &                   AVERAGE(ng) % avglrf,                          &
     &                   AVERAGE(ng) % avguwind,                        &
     &                   AVERAGE(ng) % avgvwind,                        &
     &                   AVERAGE(ng) % avgevap,                         &
     &                   AVERAGE(ng) % avgrain,                         &
     &                   AVERAGE(ng) % avgsrf,                          &
     &                   AVERAGE(ng) % avgsus,                          &
     &                   AVERAGE(ng) % avgsvs,                          &
     &                   AVERAGE(ng) % avgbus,                          &
     &                   AVERAGE(ng) % avgbvs,                          &
     &                   AVERAGE(ng) % avgu2d,                          &
     &                   AVERAGE(ng) % avgv2d,                          &
     &                   AVERAGE(ng) % avgzeta)
      CALL wclock_off (ng, iNLM, 5)
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Kout,                                    &
     &                         Nout,                                    &
     &                         pm, pn,                                  &
     &                         u, v,                                    &
     &                         W, wvel, t,                              &
     &                         bt,                                      &
     &                         rho,                                     &
     &                         hsbl,                                    &
     &                         Akv,                                     &
     &                         Akt,                                     &
     &                         stflx,                                   &
     &                         lhflx, shflx, lrflx, Uwind, Vwind,       &
     &                         evap, rain,                              &
     &                         srflx,                                   &
     &                         sustr, svstr, bustr, bvstr,              &
     &                         ubar, vbar,                              &
     &                         zeta,                                    &
     &                         avgu3d, avgv3d,                          &
     &                         avgw3d, avgwvel,                         &
     &                         avgt,                                    &
     &                         avgbt,                                   &
     &                         avgrho,                                  &
     &                         avghsbl,                                 &
     &                         avgAKv,                                  &
     &                         avgAKt,                                  &
     &                         avgAKs,                                  &
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
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(in) :: wvel(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: bt(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
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
      real(r8), intent(inout) :: avgu3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgv3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgw3d(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgwvel(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgt(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgbt(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgrho(LBi:,LBj:,:)
      real(r8), intent(inout) :: avghsbl(LBi:,LBj:)
      real(r8), intent(inout) :: avgAKv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgAKt(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgAKs(LBi:,LBj:,0:)
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
      integer :: i, itrc, j, k
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
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
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
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgu3d(i,j,k)=u(i,j,k,Nout)
              avgv3d(i,j,k)=v(i,j,k,Nout)
              avgrho(i,j,k)=rho(i,j,k)
            END DO
          END DO
        END DO
!
!  Initialized fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,k,itrc)=t(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
        END DO
           DO itrc=1,NBeT(ng)
             DO k=1,NBL(ng)
               DO j=JstrR,JendR
                 DO i=IstrR,IendR
                   avgbt(i,j,k,itrc)=bt(i,j,k,Nout,itrc)
                 END DO
               END DO
             END DO
          END DO
!
!  Initialize fields at W-points.
!
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgw3d(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
              avgwvel(i,j,k)=wvel(i,j,k)
              avgAKv(i,j,k)=Akv(i,j,k)
              avgAKt(i,j,k)=Akt(i,j,k,itemp)
              avgAKs(i,j,k)=Akt(i,j,k,isalt)
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
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
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgu3d(i,j,k)=avgu3d(i,j,k)+u(i,j,k,Nout)
              avgv3d(i,j,k)=avgv3d(i,j,k)+v(i,j,k,Nout)
              avgrho(i,j,k)=avgrho(i,j,k)+rho(i,j,k)
            END DO
          END DO
        END DO
!
!  Accumulate fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,k,itrc)=avgt(i,j,k,itrc)+t(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
        END DO
!--------------------------------------------
! values are accumulated in the bestnpz.h file
! so do not sum them here
!--------------------------------------------
          DO itrc=1,NBeT(ng)
            DO k=1,NBL(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgbt(i,j,k,itrc)=avgbt(i,j,k,itrc)+bt(i,j,k,Nout,itrc)
!         if(i.eq.144.and.j.eq.144) THEN
!          print*,'setAvg'
!          print*,'avgbt(i,j,k,itrc)=', avgbt(i,j,k,itrc)
!          
!          print*,'bt(i,j,k,Nout,1)=', bt(i,j,k,Nout,1)
!          print*,''
!          end if 
                END DO
              END DO
            END DO
          END DO
!
!  Accumulate fields at W-points.
!
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgw3d(i,j,k)=avgw3d(i,j,k)+W(i,j,k)*pm(i,j)*pn(i,j)
              avgwvel(i,j,k)=avgwvel(i,j,k)+wvel(i,j,k)
              avgAKv(i,j,k)=avgAKv(i,j,k)+Akv(i,j,k)
              avgAKt(i,j,k)=avgAKt(i,j,k)+Akt(i,j,k,itemp)
              avgAKs(i,j,k)=avgAKs(i,j,k)+Akt(i,j,k,isalt)
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG(ng),r8)
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
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
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgu3d(i,j,k)=fac*avgu3d(i,j,k)
              avgv3d(i,j,k)=fac*avgv3d(i,j,k)
              avgrho(i,j,k)=fac*avgrho(i,j,k)
            END DO
          END DO
        END DO
!
!  Process fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,k,itrc)=fac*avgt(i,j,k,itrc)
              END DO
            END DO
          END DO
        END DO
!-----------------------------------------------------------
! Computes average of time series
! For production, totals are needed so this is commented out
!-----------------------------------------------------------
!        DO itrc=1,NTS(ng)
!          DO k=1,N(ng)
!            DO j-JstrR,JendR
!              DO i-IstrR,IendR
!                anvst(i,j,k,itrc)=fac*avgt(i,j,k,itrc)
!              END DO
!            END DO
!          END DO
!        END DO
         DO itrc=1,NBeT(ng)
            DO k=1,NBL(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgbt(i,j,k,itrc)=fac*avgbt(i,j,k,itrc)
                END DO
              END DO
            END DO
         END DO
!        DO itrc=1,NTS(ng)
!          DO k=1,N(ng)
!            DO j-JstrR,JendR
!              DO i-IstrR,IendR
!                anvst(i,j,k,itrc)=fac*avgt(i,j,k,itrc)
!              END DO
!            END DO
!          END DO
!        END DO
!
!  Process fields at W-points.
!
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgw3d(i,j,k)=fac*avgw3d(i,j,k)
              avgwvel(i,j,k)=fac*avgwvel(i,j,k)
              avgAKv(i,j,k)=fac*avgAKv(i,j,k)
              avgAKt(i,j,k)=fac*avgAKt(i,j,k)
              avgAKs(i,j,k)=fac*avgAKs(i,j,k)
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE set_avg_tile
      END MODULE set_avg_mod
