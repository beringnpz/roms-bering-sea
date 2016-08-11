      MODULE mod_mixing
!
!svn $Id: mod_mixing.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Horizontal and vertical mixing coefficients:                        !
!                                                                      !
!  Akt          Vertical mixing coefficient (m2/s) for tracers.        !
!  Akv          Vertical mixing coefficient (m2/s) for momentum.       !
!  Hviscosity   Time invariant, horizontal viscosity (harmonic or      !
!                 (biharmonic) term at RHO-points.                     !
!  dAktdz       Vertical gradient in mixing coefficient (m/s) for      !
!                 tracer 1, used in float random walk calculations     !
!  diff2        Horizontal, time invariant harmonic coefficient        !
!                 (m2/s) for tracers.                                  !
!  diff4        Horizontal, time invariant biharmonic coefficient      !
!                 SQRT(m4/s) for tracers.                              !
!  visc2_r      Horizontal, time invariant harmonic viscosity          !
!                 coefficient (m2/s) at RHO-points.                    !
!  visc2_p      Horizontal, time invariant harmonic viscosity          !
!                 coefficient (m2/s) at PSI-points.                    !
!  visc4_r      Horizontal, time invariant harmonic viscosity          !
!                 coefficient SQRT(m4/s) at RHO-points.                !
!  visc4_p      Horizontal, time invariant harmonic viscosity          !
!                 coefficient SQRT(m4/s) at RHO-points.                !
!                                                                      !
!  Variables associated with the equation of state:                    !
!                                                                      !
!  alpha        Surface thermal expansion coefficient (1/Celsius).     !
!  beta         Surface saline contraction coefficient (1/PSU).        !
!  bvf          Brunt-Vaisala frequency squared (1/s2).                !
!  neutral      Coefficient to convert "in situ" density to neutral    !
!                 surface.                                             !
!                                                                      !
!  tke          Turbulent energy squared (m2/s2) at horizontal         !
!                 at W-points.                                         !
!  gls          Turbulent energy squared times turbulent length        !
!                 scale (m3/s2) at W-points.                           !
!                                                                      !
!  Large/McWilliams/Doney interior vertical mixing variables:          !
!                                                                      !
!  alfaobeta    Ratio of thermal expansion and saline contraction      !
!                 coefficients (Celsius/PSU) used in double            !
!                 diffusion.                                           !
!                                                                      !
!  Water clarity parameters:                                           !
!                                                                      !
!  Jwtype       Water clarity (Jerlov water type classification).      !
!                                                                      !
!  Large/McWilliams/Doney oceanic boundary layer variables:            !
!                                                                      !
!  ghats        Boundary layer nonlocal transport (T units/m).         !
!  hbbl         Depth of bottom oceanic boundary layer (m).            !
!  hsbl         Depth of surface oceanic boundary layer (m).           !
!  kbbl         Index of grid level above bottom  boundary layer.      !
!  ksbl         Index of grid level below surface boundary layer.      !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_MIXING
!
!  Nonlinear model state.
!
          integer,  pointer :: Jwtype(:,:)
          integer,  pointer :: ksbl(:,:)
          real(r8), pointer :: Hviscosity(:,:)
          real(r8), pointer :: visc2_p(:,:)
          real(r8), pointer :: visc2_r(:,:)
          real(r8), pointer :: visc3d_r(:,:,:)
          real(r8), pointer :: diff2(:,:,:)
          real(r8), pointer :: Akv(:,:,:)
          real(r8), pointer :: Akt(:,:,:,:)
          real(r8), pointer :: alpha(:,:)
          real(r8), pointer :: beta(:,:)
          real(r8), pointer :: bvf(:,:,:)
          real(r8), pointer :: neutral(:,:,:)
          real(r8), pointer :: hsbl(:,:)
          real(r8), pointer :: ghats(:,:,:,:)
        END TYPE T_MIXING
        TYPE (T_MIXING), allocatable :: MIXING(:)
      CONTAINS
      SUBROUTINE allocate_mixing (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1) allocate ( MIXING(Ngrids) )
! 
!  Nonlinear model state.
!
      allocate ( MIXING(ng) % visc2_p(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc2_r(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % Hviscosity(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc3d_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % diff2(LBi:UBi,LBj:UBj,NT(ng)) )
      allocate ( MIXING(ng) % Akv(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT) )
      allocate ( MIXING(ng) % alpha(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % beta(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % bvf(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % neutral(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % Jwtype(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % ksbl(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % hsbl(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % ghats(LBi:UBi,LBj:UBj,0:N(ng),NAT) )
      RETURN
      END SUBROUTINE allocate_mixing
      SUBROUTINE initialize_mixing (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_mixing" for all nested grids.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      integer :: itrc, k
      real(r8), parameter :: IniVal = 0.0_r8
      real(r8) :: cff1, cff2, cff3, cff4
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
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            MIXING(ng) % visc2_p(i,j) = visc2(ng)
            MIXING(ng) % visc2_r(i,j) = visc2(ng)
          END DO
          DO k=1,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % Hviscosity(i,j) = IniVal
              MIXING(ng) % visc3d_r(i,j,k) = IniVal
            END DO
          END DO
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
              MIXING(ng) % diff2(i,j,itrc) = tnu2(itrc,ng)
            END DO
          END DO
          DO i=Imin,Imax
            MIXING(ng) % Akv(i,j,0) = IniVal
            MIXING(ng) % Akv(i,j,N(ng)) = IniVal
          END DO
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              MIXING(ng) % Akv(i,j,k) = Akv_bak(ng)
            END DO
          END DO
          DO itrc=1,NAT
            DO i=Imin,Imax
              MIXING(ng) % Akt(i,j,0,itrc) = IniVal
              MIXING(ng) % Akt(i,j,N(ng),itrc) = IniVal
            END DO
            DO k=1,N(ng)-1
              DO i=Imin,Imax
                MIXING(ng) % Akt(i,j,k,itrc) = Akt_bak(itrc,ng)
              END DO
            END DO
          END DO
          DO i=Imin,Imax
            MIXING(ng) % alpha(i,j) = IniVal
            MIXING(ng) % beta(i,j) = IniVal
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % bvf(i,j,k) = IniVal
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % neutral(i,j,k) = IniVal
            END DO
          END DO
          DO i=Imin,Imax
            MIXING(ng) % Jwtype(i,j) = lmd_Jwt(ng)
          END DO
          DO i=Imin,Imax
            MIXING(ng) % ksbl(i,j) = 0
            MIXING(ng) % hsbl(i,j) = IniVal
          END DO
          DO itrc=1,NAT
            DO k=0,N(ng)
              DO i=Imin,Imax
                MIXING(ng) % ghats(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_mixing
      END MODULE mod_mixing
