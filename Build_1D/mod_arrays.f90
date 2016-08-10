      SUBROUTINE mod_arrays (allocate_vars)
!
!svn $Id: mod_arrays.F 999 2009-06-09 23:48:31Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine routine allocates and initializa model state arrays    !
!  for each nested and/or multiple connected grids.                    !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
      USE mod_param
      USE mod_parallel
!
      USE mod_average, ONLY : allocate_average, initialize_average
      USE mod_average2, ONLY : allocate_average2, initialize_average2
      USE mod_clima, ONLY : allocate_clima, initialize_clima
      USE mod_coupling, ONLY : allocate_coupling, initialize_coupling
      USE mod_forces, ONLY : allocate_forces, initialize_forces
      USE mod_grid, ONLY : allocate_grid, initialize_grid
      USE mod_mixing, ONLY : allocate_mixing, initialize_mixing
      USE mod_ocean, ONLY : allocate_ocean, initialize_ocean
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: allocate_vars
!
!  Local variable declarations.
!
      integer :: ng,ic
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
      integer :: tile, subs, thread
      integer, parameter :: model = 0
!
!-----------------------------------------------------------------------
!  Allocate model structures.
!-----------------------------------------------------------------------
!
      IF (allocate_vars) then
        tile=0
        DO ng=1,Ngrids
          LBi=BOUNDS(ng)%LBi(tile)
          UBi=BOUNDS(ng)%UBi(tile)
          LBj=BOUNDS(ng)%LBj(tile)
          UBj=BOUNDS(ng)%UBj(tile)
          LBij=BOUNDS(ng)%LBij
          UBij=BOUNDS(ng)%UBij
          CALL allocate_average (ng, LBi, UBi, LBj, UBj)
          CALL allocate_average2 (ng, LBi, UBi, LBj, UBj)
          CALL allocate_clima (ng, LBi, UBi, LBj, UBj,ic)
          CALL allocate_coupling (ng, LBi, UBi, LBj, UBj)
          CALL allocate_forces (ng, LBi, UBi, LBj, UBj)
          CALL allocate_grid (ng, LBi, UBi, LBj, UBj, LBij, UBij)
          CALL allocate_mixing (ng, LBi, UBi, LBj, UBj)
          CALL allocate_ocean  (ng, LBi, UBi, LBj, UBj)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Allocate and intialize variables within structures for each grid.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        DO ng=1,Ngrids
          subs=1
          DO tile=subs*thread,subs*(thread+1)-1
            CALL initialize_average (ng, tile)
            CALL initialize_average2 (ng, tile)
            CALL initialize_clima (ng, tile)
            CALL initialize_coupling (ng, tile, model)
            CALL initialize_forces (ng, tile, model)
            CALL initialize_grid (ng, tile, model)
            CALL initialize_mixing (ng, tile, model)
            CALL initialize_ocean (ng, tile, model)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE mod_arrays
