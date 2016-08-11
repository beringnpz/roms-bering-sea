      MODULE mod_coupling
!
!svn $Id: mod_coupling.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  DU_avg1   Time averaged U-flux for 2D equations (m3/s).             !
!  DU_avg2   Time averaged U-flux for 3D equations coupling (m3/s).    !
!  DV_avg1   Time averaged V-flux for 2D equations (m3/s).             !
!  DV_avg2   Time averaged V-flux for 3D equations coupling (m3/s).    !
!  Zt_avg1   Free-surface averaged over all short time-steps (m).      !
!  rhoA      Normalized vertical averaged density.                     !
!  rhoS      Normalized vertical averaged density perturbation.        !
!  rufrc     Right-hand-side forcing term for 2D U-momentum (m4/s2)    !
!  rvfrc     Right-hand-side forcing term for 2D V-momentum (m4/s2)    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_COUPLING
!
!  Nonlinear model state.
!
          real(r8), pointer :: DU_avg1(:,:)
          real(r8), pointer :: DU_avg2(:,:)
          real(r8), pointer :: DV_avg1(:,:)
          real(r8), pointer :: DV_avg2(:,:)
          real(r8), pointer :: Zt_avg1(:,:)
          real(r8), pointer :: rufrc(:,:)
          real(r8), pointer :: rvfrc(:,:)
          real(r8), pointer :: rhoA(:,:)
          real(r8), pointer :: rhoS(:,:)
        END TYPE T_COUPLING
        TYPE (T_COUPLING), allocatable :: COUPLING(:)
        CONTAINS
      SUBROUTINE allocate_coupling (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!   grids.                                                             !
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
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( COUPLING(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( COUPLING(ng) % DU_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DV_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % DV_avg2(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % Zt_avg1(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rufrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rvfrc(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rhoA(LBi:UBi,LBj:UBj) )
      allocate ( COUPLING(ng) % rhoS(LBi:UBi,LBj:UBj) )
      RETURN
      END SUBROUTINE allocate_coupling
      SUBROUTINE initialize_coupling (ng, tile, model)
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
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
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
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            COUPLING(ng) % DU_avg1(i,j) = IniVal
            COUPLING(ng) % DU_avg2(i,j) = IniVal
            COUPLING(ng) % DV_avg1(i,j) = IniVal
            COUPLING(ng) % DV_avg2(i,j) = IniVal
            COUPLING(ng) % Zt_avg1(i,j) = IniVal
            COUPLING(ng) % rufrc(i,j) = IniVal
            COUPLING(ng) % rvfrc(i,j) = IniVal
            COUPLING(ng) % rhoA(i,j) = IniVal
            COUPLING(ng) % rhoS(i,j) = IniVal
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_coupling
      END MODULE mod_coupling
