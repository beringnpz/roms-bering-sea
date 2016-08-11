      MODULE metrics_mod
!
!svn $Id: metrics.F 999 2009-06-09 23:48:31Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes various horizontal metric terms.              !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: metrics
      CONTAINS
!
!***********************************************************************
      SUBROUTINE metrics (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL metrics_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
     &                   GRID(ng) % f,                                  &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % grdscl,                             &
     &                   MIXING(ng) % Hviscosity,                       &
     &                   GRID(ng) % om_p,                               &
     &                   GRID(ng) % om_r,                               &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % om_v,                               &
     &                   GRID(ng) % on_p,                               &
     &                   GRID(ng) % on_r,                               &
     &                   GRID(ng) % on_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   GRID(ng) % fomn,                               &
     &                   GRID(ng) % omn,                                &
     &                   GRID(ng) % pnom_p,                             &
     &                   GRID(ng) % pnom_r,                             &
     &                   GRID(ng) % pnom_u,                             &
     &                   GRID(ng) % pnom_v,                             &
     &                   GRID(ng) % pmon_p,                             &
     &                   GRID(ng) % pmon_r,                             &
     &                   GRID(ng) % pmon_u,                             &
     &                   GRID(ng) % pmon_v)
      RETURN
      END SUBROUTINE metrics
!
!***********************************************************************
      SUBROUTINE metrics_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
     &                         f, h, pm, pn,                            &
     &                         Hz, z_r, z_w,                            &
     &                         grdscl,                                  &
     &                         Hviscosity,                              &
     &                         om_p, om_r, om_u, om_v,                  &
     &                         on_p, on_r, on_u, on_v,                  &
     &                         fomn, omn,                               &
     &                         pnom_p, pnom_r, pnom_u, pnom_v,          &
     &                         pmon_p, pmon_r, pmon_u, pmon_v)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_iounits
!
      USE exchange_2d_mod
      USE set_depth_mod, ONLY : set_depth_tile
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(out) :: grdscl(LBi:,LBj:)
      real(r8), intent(out) :: Hviscosity(LBi:,LBj:)
      real(r8), intent(out) :: om_p(LBi:,LBj:)
      real(r8), intent(out) :: om_r(LBi:,LBj:)
      real(r8), intent(out) :: om_u(LBi:,LBj:)
      real(r8), intent(out) :: om_v(LBi:,LBj:)
      real(r8), intent(out) :: on_p(LBi:,LBj:)
      real(r8), intent(out) :: on_r(LBi:,LBj:)
      real(r8), intent(out) :: on_u(LBi:,LBj:)
      real(r8), intent(out) :: on_v(LBi:,LBj:)
      real(r8), intent(out) :: fomn(LBi:,LBj:)
      real(r8), intent(out) :: omn(LBi:,LBj:)
      real(r8), intent(out) :: pnom_p(LBi:,LBj:)
      real(r8), intent(out) :: pnom_r(LBi:,LBj:)
      real(r8), intent(out) :: pnom_u(LBi:,LBj:)
      real(r8), intent(out) :: pnom_v(LBi:,LBj:)
      real(r8), intent(out) :: pmon_p(LBi:,LBj:)
      real(r8), intent(out) :: pmon_r(LBi:,LBj:)
      real(r8), intent(out) :: pmon_u(LBi:,LBj:)
      real(r8), intent(out) :: pmon_v(LBi:,LBj:)
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: NSUB, bry, i, is, j, k, rec
      real(r8), parameter :: Large = 1.0E+20_r8
      real(r8), parameter :: PecletCoef = 1.0_r8 / 12.0_r8
      real(r8), parameter :: Uscale = 0.1_r8
      real(r8) :: cff, cff1, cff2
      real(r8) :: my_DXmax, my_DXmin, my_DYmax, my_DYmin
      real(r8) :: my_DZmax, my_DZmin
      real(r8) :: my_Cu_Cor, my_Cu_max, my_Cu_min, my_grdmax
      real(r8) :: my_ViscMax, my_ViscMin
      character (len=4) :: units
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
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
!  Compute 1/m, 1/n, 1/mn, and f/mn at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          om_r(i,j)=1.0_r8/pm(i,j)
          on_r(i,j)=1.0_r8/pn(i,j)
          omn(i,j)=1.0_r8/(pm(i,j)*pn(i,j))
          fomn(i,j)=f(i,j)*omn(i,j)
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_r)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_r)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        omn)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        fomn)
!
!-----------------------------------------------------------------------
!  Compute n/m, and m/n at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          pnom_r(i,j)=pn(i,j)/pm(i,j)
          pmon_r(i,j)=pm(i,j)/pn(i,j)
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_r)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_r)
!
!-----------------------------------------------------------------------
!  Compute m/n, 1/m, and 1/n at horizontal U-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          pmon_u(i,j)=(pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
          pnom_u(i,j)=(pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
          om_u(i,j)=2.0_r8/(pm(i-1,j)+pm(i,j))
          on_u(i,j)=2.0_r8/(pn(i-1,j)+pn(i,j))
        END DO
      END DO
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_u)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_u)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_u)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_u)
!
!-----------------------------------------------------------------------
!  Compute n/m, 1/m, and 1/m at horizontal V-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          pmon_v(i,j)=(pm(i,j-1)+pm(i,j))/(pn(i,j-1)+pn(i,j))
          pnom_v(i,j)=(pn(i,j-1)+pn(i,j))/(pm(i,j-1)+pm(i,j))
          om_v(i,j)=2.0_r8/(pm(i,j-1)+pm(i,j))
          on_v(i,j)=2.0_r8/(pn(i,j-1)+pn(i,j))
        END DO
      END DO
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_v)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_v)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_v)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_v)
!
!-----------------------------------------------------------------------
!  Compute n/m and m/n at horizontal PSI-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          pnom_p(i,j)=(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))/        &
     &                (pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          pmon_p(i,j)=(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))/        &
     &                (pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          om_p(i,j)=4.0_r8/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j)=4.0_r8/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
        END DO
      END DO
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pnom_p)
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        pmon_p)
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        om_p)
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        on_p)
!
!-----------------------------------------------------------------------
! Determine the squared root of the area of each grid cell used to
! rescale horizontal mixing by the grid size.
!-----------------------------------------------------------------------
!
      cff=0.0_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          grdscl(i,j)=SQRT(om_r(i,j)*on_r(i,j))
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        grdscl)
!
!-----------------------------------------------------------------------
! Compute time invariant, horizontal mixing coefficient using grid
! scale. Following Holland et (1998), Webb et al. (1998), Griffies
! and Hallberg (2000), and Lee et al. (2002), the horizontal mixing
! coefficient can be estimated as:
!
!   Hmixing = 1/12 * Uscale * grdscl           (Harmonic)
!   Bmixing = 1/12 * Uscale * grdscl**3        (Biharmonic)
!-----------------------------------------------------------------------
!
      my_ViscMin= Large
      my_ViscMax=-Large
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Hviscosity(i,j)=PecletCoef*Uscale*grdscl(i,j)
          my_ViscMin=MIN(my_ViscMin, Hviscosity(i,j))
          my_ViscMax=MAX(my_ViscMax, Hviscosity(i,j))
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Hviscosity)
!
!-----------------------------------------------------------------------
!  Compute minimum and maximum grid spacing.
!-----------------------------------------------------------------------
!
!  Compute time invariant depths (use zero free-surface).
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8
        END DO
      END DO
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
     &                     A2d,                                         &
     &                     Hz, z_r, z_w)
!
!  Compute grid spacing range.
!
      my_DXmin= Large
      my_DXmax=-Large
      my_DYmin= Large
      my_DYmax=-Large
      my_DZmin= Large
      my_DZmax=-Large
      my_grdmax=-Large
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=grdscl(i,j)
          my_DXmin=MIN(my_DXmin,om_r(i,j))
          my_DXmax=MAX(my_DXmax,om_r(i,j))
          my_DYmin=MIN(my_DYmin,on_r(i,j))
          my_DYmax=MAX(my_DYmax,on_r(i,j))
          my_grdmax=MAX(my_grdmax,cff)
          DO k=1,N(ng)
            my_DZmin=MIN(my_DZmin,Hz(i,j,k))
            my_DZmax=MAX(my_DZmax,Hz(i,j,k))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute Courant number.
!-----------------------------------------------------------------------
!
!  The Courant number is defined as:
!
!     Cu = c * dt * SQRT (1/dx^2 + 1/dy^2)
!
!  where c=SQRT(g*h) is phase speed for barotropic mode, and dx, dy
!  are grid spacing in each direction.
!
      my_Cu_min= Large
      my_Cu_max=-Large
      my_Cu_Cor=-Large
      DO j=JstrR,JendR
        DO i=IstrR,IendR
            cff=dtfast(ng)*                                             &
     &          SQRT(g*ABS(h(i,j))*(pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
            my_Cu_min=MIN(my_Cu_min,cff)
            my_Cu_max=MAX(my_Cu_max,cff)
            cff=dt(ng)*ABS(f(i,j))
            my_Cu_Cor=MAX(my_Cu_Cor,cff)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform global reductions.
!-----------------------------------------------------------------------
!
      IF ((Istr.eq.1).and.(Jstr.eq.1).and.                              &
     &    (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
      IF (tile_count.eq.0) THEN
        Cu_min=my_Cu_min
        Cu_max=my_Cu_max
        Cu_Cor=my_Cu_Cor
        grdmax(ng)=my_grdmax
        DXmin(ng)=my_DXmin
        DXmax(ng)=my_DXmax
        DYmin(ng)=my_DYmin
        DYmax(ng)=my_DYmax
        DZmin(ng)=my_DZmin
        DZmax(ng)=my_DZmax
        ViscMin(ng)=my_ViscMin
        ViscMax(ng)=my_ViscMax
      ELSE
        Cu_min=MIN(Cu_min,my_Cu_min)
        Cu_max=MAX(Cu_max,my_Cu_max)
        Cu_Cor=MAX(Cu_Cor,my_Cu_Cor)
        grdmax(ng)=MAX(grdmax(ng),my_grdmax)
        DXmin(ng)=MIN(DXmin(ng),my_DXmin)
        DXmax(ng)=MAX(DXmax(ng),my_DXmax)
        DYmin(ng)=MIN(DYmin(ng),my_DYmin)
        DYmax(ng)=MAX(DYmax(ng),my_DYmax)
        DZmin(ng)=MIN(DZmin(ng),my_DZmin)
        DZmax(ng)=MAX(DZmax(ng),my_DZmax)
        ViscMin(ng)=MIN(ViscMin(ng),my_ViscMin)
        ViscMax(ng)=MAX(ViscMax(ng),my_ViscMax)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
        IF (Master.and.LwrtInfo(ng)) THEN
          WRITE(stdout,10) DXmin(ng)/1000.0_r8, DXmax(ng)/1000.0_r8,    &
     &                     DYmin(ng)/1000.0_r8, DYmax(ng)/1000.0_r8
  10      FORMAT (/,' Minimum X-grid spacing, DXmin = ',1pe15.8,' km',  &
     &            /,' Maximum X-grid spacing, DXmax = ',1pe15.8,' km',  &
     &            /,' Minimum Y-grid spacing, DYmin = ',1pe15.8,' km',  &
     &            /,' Maximum Y-grid spacing, DYmax = ',1pe15.8,' km')
          WRITE(stdout,20) DZmin(ng), DZmax(ng)
  20      FORMAT (' Minimum Z-grid spacing, DZmin = ',1pe15.8,' m',/,   &
     &            ' Maximum Z-grid spacing, DZmax = ',1pe15.8,' m')
          WRITE (stdout,30) Cu_min, Cu_max, Cu_Cor
  30      FORMAT (/,' Minimum barotropic Courant Number = ', 1pe15.8,/, &
     &              ' Maximum barotropic Courant Number = ', 1pe15.8,/, &
     &              ' Maximum Coriolis   Courant Number = ', 1pe15.8,/)
          WRITE (stdout,40) grdmax(ng)/1000.0_r8
  40      FORMAT (' Horizontal mixing scaled by grid size,',            &
     &            ' GRDMAX = ',1pe15.8,' km')
          units='m2/s'
          WRITE (stdout,60) ViscMin(ng), TRIM(units),                   &
     &                      ViscMax(ng), TRIM(units)
  60      FORMAT (/,' Minimum horizontal viscosity coefficient = ',     &
     &             1pe15.8,1x,a,                                        &
     &            /,' Maximum horizontal viscosity coefficient = ',     &
     &             1pe15.8,1x,a)
        END IF
      END IF
      RETURN
      END SUBROUTINE metrics_tile
      END MODULE metrics_mod
