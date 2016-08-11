      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 1076 2009-09-25 23:18:43Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and assimilation data    !
!  from input NetCDF files.  If there is more than one time-record,    !
!  data  is loaded  into global two-time record arrays.  The actual    !
!  interpolation is carried elsewhere.                                 !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
      integer :: LBi, UBi, LBj, UBj
      integer :: i, sp, gr, ii, jj
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
!integer  :: idUnms         ! surface U-momentum stress (NCEP)
!integer  :: ncFRCid(NV,Ngrids)  ! input forcing
!ncFRCid(i,ng)=-1
!
!
!
!-----------------------------------------------------------------------
!  Surface wind components.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idUair, ncFRCid(idUair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % UwindG(LBi,LBj,1))
      CALL get_2dfld (ng , iNLM, idVair, ncFRCid(idVair,ng),            &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % VwindG(LBi,LBj,1))
!
!-----------------------------------------------------------------------
!  Surface air pressure.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idPair, ncFRCid(idPair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % PairG(LBi,LBj,1))
!
!-----------------------------------------------------------------------
!  Surface solar shortwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idSrad, ncFRCid(idSrad,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % srflxG(LBi,LBj,1))
!
!-----------------------------------------------------------------------
!  Surface downwelling longwave radiation.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idLdwn, ncFRCid(idLdwn,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % lrflxG(LBi,LBj,1))
!
!-----------------------------------------------------------------------
!  Surface air temperature.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idTair, ncFRCid(idTair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % TairG(LBi,LBj,1))
!
!-----------------------------------------------------------------------
!  Surface air humidity.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idQair, ncFRCid(idQair,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % HairG(LBi,LBj,1))
!
!-----------------------------------------------------------------------
!  Rain fall rate.
!-----------------------------------------------------------------------
!
      CALL get_2dfld (ng, iNLM, idrain, ncFRCid(idrain,ng),             &
     &                nFfiles(ng), FRCname(1,ng), update(1),            &
     &                LBi, UBi, LBj, UBj, 2, 1,                         &
     &                FORCES(ng) % rainG(LBi,LBj,1))
!
!=======================================================================
!  Read in climatology data from  NetCDF file.
!=======================================================================
!
! Only nudge the active tracers if doing a 1D simulation
      DO i=1,NAT
        CALL get_3dfld (ng, iNLM, idTclm(i), ncCLMid(ng), 1,            &
     &                  CLMname(ng), update(1),                         &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  CLIMA(ng) % tclmG(LBi,LBj,1,1,i))
      END DO
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
      RETURN
      END SUBROUTINE get_data
