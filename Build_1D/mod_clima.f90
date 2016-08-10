      MODULE mod_clima
!
!svn $Id: mod_clima.F 988 2009-05-26 19:40:39Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sea surface height fields.                                          !
!                                                                      !
!   ssh         Climatology for sea surface height (m).                !
!   sshG        Latest two-time snapshots of input "ssh" grided        !
!                 data used for interpolation.                         !
!   zeta_ads    Sensitivity functional for sea surface height.         !
!   zeta_adsF   Latest two-time snapshots of input "zeta_ads" grided   !
!                 data used fot interpolation.                         !
!                                                                      !
!  2D momentum fields.                                                 !
!                                                                      !
!   ubarclm     Vertically integrated U-momentum climatology (m/s).    !
!   ubarclmG    Latest two-time snapshots of input "ubarclm" grided    !
!                 data used for interpolation.                         !
!   ubar_ads    Sensitivity functional for vertically integrated       !
!                 U-momentum.                                          !
!   ubar_adsG   Latest two-time snapshots of input "ubar_ads" grided   !
!                 data used for interpolation.                         !
!   vbarclm     Vertically integrated V-momentum climatology (m/s).    !
!   vbarclmG    Latest two-time snapshots of input "vbarclm" grided    !
!                 data used for interpolation.                         !
!   vbar_ads    Sensitivity functional for vertically integrated       !
!                 V-momentum.                                          !
!   vbar_adsG   Latest two-time snapshots of input "vbar_ads" grided   !
!                 data used for interpolation.                         !
!                                                                      !
!  Tracer fields.                                                      !
!                                                                      !
!   tclm        Climatology for tracer type variables (usually,        !
!                 temperature: degC; salinity: PSU).                   !
!   tclmG       Latest two-time snapshots of input "tclm" grided       !
!                 data used for interpolation.                         !
!   t_ads       Sensitivity functional for tracer type variables.      !
!   t_adsG      Latest two-time snapshots of input "t_ads" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  3D momentum climatology.                                            !
!                                                                      !
!   uclm        3D U-momentum climatology (m/s).                       !
!   uclmG       Latest two-time snapshots of input "uclm" grided       !
!                 data used for interpolation.                         !
!   u_ads       Sensitivity functional for 3D U-momentum.              !
!   u_adsG      Latest two-time snapshots of input "u_ads" grided      !
!                 data used for interpolation.                         !
!   vclm        3D V-momentum climatology (m/s).                       !
!   vclmG       Latest two-time snapshots of input "vclm" grided       !
!                 data used for interpolation.                         !
!   v_ads       Sensitivity functional for 3D V-momentum.              !
!   v_adsG      Latest two-time snapshots of input "v_ads" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Nudging variables.                                                  !
!                                                                      !
!   M2nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 2D momentum data.                                    !
!   M3nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 3D momentum data.                                    !
!   Tnudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 tracer data.                                         !
!   Znudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 sea surface height data.                             !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        USE mod_scalars
        implicit none
        TYPE T_CLIMA
          real(r8), pointer :: ubarclm(:,:)
          real(r8), pointer :: vbarclm(:,:)
          real(r8), pointer :: M2nudgcof(:,:)
          real(r8), pointer :: tclm(:,:,:,:)
          real(r8), pointer :: tclmG(:,:,:,:,:)
          real(r8), pointer :: Tnudgcof(:,:,:)
        END TYPE T_CLIMA
        TYPE (T_CLIMA), allocatable :: CLIMA(:)
      CONTAINS
      SUBROUTINE allocate_clima (ng, LBi, UBi, LBj, UBj,ic)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj,ic
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( CLIMA(Ngrids) )
!
      allocate ( CLIMA(ng) % ubarclm(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % vbarclm(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % M2nudgcof(LBi:UBi,LBj:UBj) )
!      allocate ( CLIMA(ng) % tclm(LBi:UBi,LBj:UBj,N(ng),ic))
      allocate ( CLIMA(ng) % tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng)+2))
!      allocate ( CLIMA(ng) % tclmG(LBi:UBi,LBj:UBj,N(ng),2,ic))
      allocate ( CLIMA(ng) % tclmG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)+2))
      allocate ( CLIMA(ng) % Tnudgcof(LBi:UBi,LBj:UBj,NT(ng)) )
      RETURN
      END SUBROUTINE allocate_clima
      SUBROUTINE initialize_clima (ng, tile)
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
          CLIMA(ng) % ubarclm(i,j) = IniVal
          CLIMA(ng) % vbarclm(i,j) = IniVal
          CLIMA(ng) % M2nudgcof(i,j) = IniVal
        END DO
        DO k=1,N(ng)
          DO i=Imin,Imax
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=Imin,Imax
              CLIMA(ng) % tclm(i,j,k,itrc) = IniVal
              CLIMA(ng) % tclmG(i,j,k,1,itrc) = IniVal
              CLIMA(ng) % tclmG(i,j,k,2,itrc) = IniVal
            END DO
          END DO
          DO i=Imin,Imax
            CLIMA(ng) % Tnudgcof(i,j,itrc) = IniVal
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_clima
      END MODULE mod_clima
