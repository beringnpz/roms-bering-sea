      MODULE mod_average2
!
!svn $Id: mod_average.F 707 2008-08-19 17:58:55Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! FOR the AVERAGE2 option, these become surface only fields.           !
!                                                                      !
!  2D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avgu2d     2D velocity component (m/s) in the XI-direction.        !
!  avgv2d     2D velocity component (m/s) in the ETA-direction.       !
!  avgzeta    Free surface (m).                                       !
!                                                                      !
!  3D Time-averaged fields for output purposes.                        !
!                                                                      !
!  avglhf     Latent heat flux (W/m2).                                !
!  avglrf     Longwave radiation flux (W/m2).                         !
!  avgHuon    U-momentum flux, Hz*u/pn (m3/s).                        !
!  avgHuonT   Tracer u-transport, Hz*u*t/pn (Tunits m3/s).            !
!  avgHvom    V-momentum flux, Hz*v/pm (m3/s).                        !
!  avgHvomT   Tracer v-transport, Hz*v*t/pn (Tunits m3/s).            !
!  avgbus     Bottom u-momentum stress (N/m2).                        !
!  avgbvs     Bottom v-momentum stress (N/m2).                        !
!  avghbbl    Depth of oceanic bottom boundary layer (m).             !
!  avghsbl    Depth of oceanic surface boundary layer (m).            !
!  avgrho     Density anomaly (kg/m3).                                !
!  avgsssflx  Sea surface salinity flux correction.                   !
!  avgshf     Sensible heat flux (W/m2).                              !
!  avgsrf     Shortwave radiation flux (W/m2).                        !
!  avgstf     Surface net heat flux (W/m2).                           !
!  avgswf     Surface net salt flux (kg/m2/s).                        !
!  avgevap    Surface net evaporation (kg/m2/s).                      !
!  avgrain    Surface net rain fall (kg/m2/s).                        !
!  avgsus     Surface u-momentum stress (N/m2).                       !
!  avgsvs     Surface v-momentum stress (N/m2).                       !
!  avgt       Tracer type variables (usually, potential temperature   !
!               and salinity).                                         !
!  avguwind   2D wind velocity component (m/s) in the XI-direction.   !
!  avgvwind   2D wind velocity component (m/s) in the ETA-direction.  !
!  avgu3d     3D velocity component (m/s) in the XI-direction.        !
!  avgv3d     3D velocity component (m/s) in the ETA-direction.       !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_AVERAGE2
          real(r8), pointer :: avgu2d(:,:)
          real(r8), pointer :: avgv2d(:,:)
          real(r8), pointer :: avgzeta(:,:)
          real(r8), pointer :: avgrho(:,:)
          real(r8), pointer :: avgt(:,:,:)
          real(r8), pointer :: avgu3d(:,:)
          real(r8), pointer :: avgv3d(:,:)
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
        END TYPE T_AVERAGE2
        TYPE (T_AVERAGE2), allocatable :: AVERAGE2(:)
      CONTAINS
      SUBROUTINE allocate_average2 (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1 ) allocate ( AVERAGE2(Ngrids) )
!
      allocate ( AVERAGE2(ng) % avgu2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgv2d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgzeta(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsvs(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgbus(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgbvs(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgrho(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgt(LBi:UBi,LBj:UBj,NT(ng)) )
      allocate ( AVERAGE2(ng) % avgu3d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgv3d(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgstf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgswf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avglhf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avglrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgshf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avguwind(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgvwind(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgevap(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgrain(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avgsrf(LBi:UBi,LBj:UBj) )
      allocate ( AVERAGE2(ng) % avghsbl(LBi:UBi,LBj:UBj) )
      RETURN
      END SUBROUTINE allocate_average2
      SUBROUTINE initialize_average2 (ng, tile)
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
      integer :: itrc, k
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
          AVERAGE2(ng) % avgu2d(i,j) = IniVal
          AVERAGE2(ng) % avgv2d(i,j) = IniVal
          AVERAGE2(ng) % avgzeta(i,j) = IniVal
          AVERAGE2(ng) % avgsus(i,j) = IniVal
          AVERAGE2(ng) % avgsvs(i,j) = IniVal
          AVERAGE2(ng) % avgbus(i,j) = IniVal
          AVERAGE2(ng) % avgbvs(i,j) = IniVal
          AVERAGE2(ng) % avgstf(i,j) = IniVal
          AVERAGE2(ng) % avgswf(i,j) = IniVal
          AVERAGE2(ng) % avglhf(i,j) = IniVal
          AVERAGE2(ng) % avglrf(i,j) = IniVal
          AVERAGE2(ng) % avgshf(i,j) = IniVal
          AVERAGE2(ng) % avguwind(i,j) = IniVal
          AVERAGE2(ng) % avgvwind(i,j) = IniVal
          AVERAGE2(ng) % avgevap(i,j) = IniVal
          AVERAGE2(ng) % avgrain(i,j) = IniVal
          AVERAGE2(ng) % avgsrf(i,j) = IniVal
          AVERAGE2(ng) % avghsbl(i,j) = IniVal
        END DO
        DO i=Imin,Imax
            AVERAGE2(ng) % avgrho(i,j) = IniVal
            AVERAGE2(ng) % avgu3d(i,j) = IniVal
            AVERAGE2(ng) % avgv3d(i,j) = IniVal
        END DO
        DO itrc=1,NT(ng)
          DO i=Imin,Imax
            AVERAGE2(ng) % avgt(i,j,itrc) = IniVal
          END DO
        END DO
        DO i=Imin,Imax
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_average2
      END MODULE mod_average2
