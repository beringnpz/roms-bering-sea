      MODULE mod_average
!
!svn $Id: mod_average.F 990 2009-05-28 00:54:01Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  2D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgu2d     2D velocity component (m/s) in the XI-direction.         !
!  avgv2d     2D velocity component (m/s) in the ETA-direction.        !
!  avgzeta    Free surface (m).                                        !
!                                                                      !
!  3D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgAKs     Vertical diffusion of Salinity (m2/s).                   !
!  avgAKt     Vertical diffusion of temperature (m2/s).                !
!  avgAKv     Vertical viscosity (m2/s).                               !
!  avgbedldu  Bed load flux u-direction (kg/m2/s).                     !
!  avgbedldv  Bed load flux v-direction (kg/m2/s).                     !
!  avglhf     Latent heat flux (W/m2).                                 !
!  avglrf     Longwave radiation flux (W/m2).                          !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                         !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).             !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                         !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).             !
!  avgbus     Bottom u-momentum stress (N/m2).                         !
!  avgbvs     Bottom v-momentum stress (N/m2).                         !
!  avghbbl    Depth of oceanic bottom boundary layer (m).              !
!  avghsbl    Depth of oceanic surface boundary layer (m).             !
!  avgrho     Density anomaly (kg/m3).                                 !
!  avgsssflx  Sea surface salinity flux correction.                    !
!  avgshf     Sensible heat flux (W/m2).                               !
!  avgsrf     Shortwave radiation flux (W/m2).                         !
!  avgstf     Surface net heat flux (W/m2).                            !
!  avgswf     Surface net salt flux (kg/m2/s).                         !
!  avgevap    Surface net evaporation (kg/m2/s).                       !
!  avgrain    Surface net rain fall (kg/m2/s).                         !
!  avgsus     Surface u-momentum stress (N/m2).                        !
!  avgsvs     Surface v-momentum stress (N/m2).                        !
!  avgt       Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  avgUT      Quadratic term <u*t> for potential temperature and       !
!               salinity at U-points.                                  !
!  avgVT      Quadratic term <v*t> for potential temperature and       !
!               salinity at V-points.                                  !
!  avgTT      Quadratic term <t*t> for tracers.                        !
!  avgUU      Quadratic term <u*u> for 3D momentum at U-points.        !
!  avgUV      Quadratic term <u*v> for 3D momentum at RHO-points.      !
!  avgVV      Quadratic term <v*v> for 3D momentum at V-points.        !
!  avgU2      Quadratic term <ubar*ubar> for 2D momentum at U-points.  !
!  avgV2      Quadratic term <vbar*vbar> for 2D momentum at V-points.  !
!  avgZZ      Quadratic term <zeta*zeta> for free-surface.             !
!  avguwind   2D wind velocity component (m/s) in the XI-direction.    !
!  avgvwind   2D wind velocity component (m/s) in the ETA-direction.   !
!  avgu3d     3D velocity component (m/s) in the XI-direction.         !
!  avgv3d     3D velocity component (m/s) in the ETA-direction.        !
!  avgw3d     S-coordinate [omega*Hz/mn] vertical velocity (m3/s).     !
!  avgwvel    3D "true" vertical velocity (m/s).                       !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_AVERAGE
          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)
          real(r8), pointer :: avgzeta(:,:)
          real(r8), pointer :: avgrho(:,:,:)
          real(r8), pointer :: avgt(:,:,:,:)
          real(r8), pointer :: avgbt(:,:,:,:)
          real(r8), pointer :: avgu3d(:,:,:)
          real(r8), pointer :: avgv3d(:,:,:)
          real(r8), pointer :: avgw3d(:,:,:)
          real(r8), pointer :: avgwvel(:,:,:)
          real(r8), pointer :: avgAKs(:,:,:)
          real(r8), pointer :: avgAKt(:,:,:)
          real(r8), pointer :: avgAKv(:,:,:)
          real(r8), pointer :: avgstf(:,:)
          real(r8), pointer :: avgswf(:,:)
          real(r8), pointer :: avglhf(:,:)
          real(r8), pointer :: avglrf(:,:)
          real(r8), pointer :: avgshf(:,:)
          real(r8), pointer :: avguwind(:,:)
          real(r8), pointer :: avgvwind(:,:)
          real(r8), pointer :: avgevap(:,:)
          real(r8), pointer :: avgrain(:,:)
          real(r8), pointer :: avgsrf(:,:)
          real(r8), pointer :: avghsbl(:,:)
          real(r8), pointer :: avgsus(:,:)
          real(r8), pointer :: avgsvs(:,:)
          real(r8), pointer :: avgbus(:,:)
          real(r8), pointer :: avgbvs(:,:)
        END TYPE T_AVERAGE
        TYPE (T_AVERAGE), allocatable :: AVERAGE(:)
      CONTAINS
      SUBROUTINE allocate_average (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1 ) allocate ( AVERAGE(Ngrids) )
!
      allocate ( AVERAGE(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgzeta(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgrho(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
       allocate ( AVERAGE(ng) % avgbt(LBi:UBi,LBj:UBj,NBL(ng),NBeT(ng)) )
      allocate ( AVERAGE(ng) % avgu3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgv3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( AVERAGE(ng) % avgw3d(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( AVERAGE(ng) % avgwvel(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( AVERAGE(ng) % avgAKs(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( AVERAGE(ng) % avgAKt(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( AVERAGE(ng) % avgAKv(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( AVERAGE(ng) % avgstf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgswf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avglhf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avglrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgshf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avguwind(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgvwind(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgevap(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgrain(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgsrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avghsbl(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgsus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgbus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE(ng) % avgbvs(LBi:UBi,LBj:UBj) )
      RETURN
      END SUBROUTINE allocate_average
      SUBROUTINE initialize_average (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      integer :: itrc, itrc2, k
      real(r8), parameter :: IniVal = 0.0_r8
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
!  Set array initialization range.
!
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          AVERAGE(ng) % avgu2d(i,j) = IniVal
          AVERAGE(ng) % avgv2d(i,j) = IniVal
          AVERAGE(ng) % avgzeta(i,j) = IniVal
          AVERAGE(ng) % avgsus(i,j) = IniVal
          AVERAGE(ng) % avgsvs(i,j) = IniVal
          AVERAGE(ng) % avgbus(i,j) = IniVal
          AVERAGE(ng) % avgbvs(i,j) = IniVal
          AVERAGE(ng) % avgstf(i,j) = IniVal
          AVERAGE(ng) % avgswf(i,j) = IniVal
          AVERAGE(ng) % avglhf(i,j) = IniVal
          AVERAGE(ng) % avglrf(i,j) = IniVal
          AVERAGE(ng) % avgshf(i,j) = IniVal
          AVERAGE(ng) % avguwind(i,j) = IniVal
          AVERAGE(ng) % avgvwind(i,j) = IniVal
          AVERAGE(ng) % avgevap(i,j) = IniVal
          AVERAGE(ng) % avgrain(i,j) = IniVal
          AVERAGE(ng) % avgsrf(i,j) = IniVal
          AVERAGE(ng) % avghsbl(i,j) = IniVal
        END DO
        DO k=1,N(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgrho(i,j,k) = IniVal
            AVERAGE(ng) % avgu3d(i,j,k) = IniVal
            AVERAGE(ng) % avgv3d(i,j,k) = IniVal
          END DO
        END DO
        DO k=0,N(ng)
          DO i=Imin,Imax
            AVERAGE(ng) % avgw3d(i,j,k) = IniVal
            AVERAGE(ng) % avgwvel(i,j,k) = IniVal
            AVERAGE(ng) % avgAKs(i,j,k) = IniVal
            AVERAGE(ng) % avgAKt(i,j,k) = IniVal
            AVERAGE(ng) % avgAKv(i,j,k) = IniVal
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgt(i,j,k,itrc) = IniVal
            END DO
          END DO
        END DO
       DO itrc=1,NBeT(ng)
          DO k=1,NBL(ng)
            DO i=Imin,Imax
              AVERAGE(ng) % avgbt(i,j,k,itrc) = IniVal
              END DO
          END DO
        END DO
        DO i=Imin,Imax
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_average
      END MODULE mod_average
