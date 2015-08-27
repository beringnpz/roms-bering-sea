#include "cppdefs.h"
      MODULE IceNH4bc_mod
#ifdef ICE_MODEL
!
!***********************************************************************
!  Compute lateral boundary conditions for ice ammonium
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC IceNH4bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE IceNH4bc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL IceNH4bc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 liold(ng), linew(ng),                            &
     &                 ICE(ng) % ui,                                    &
     &                 ICE(ng) % vi,                                    &
     &                 ICE(ng) % IceNH4)
      RETURN
      END SUBROUTINE IceNH4bc

!
!***********************************************************************
      SUBROUTINE IceNH4bc_tile (ng, tile,                                  &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           ui, vi, IceNH4)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: liold, linew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNH4(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNH4(LBi:UBi,LBj:UBj,2)
# endif

!
!  Local variable declarations.
!
      integer :: i, j, know

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liold

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# ifdef WEST_IceNH4CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             IceNH4(0,j,linew)=BOUNDARY(ng)%IceNH4_west(j)
#  ifdef MASKING
             IceNH4(0,j,linew)=IceNH4(0,j,linew)*                             &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             IceNH4(0,j,linew)=IceNH4(0,j,linew)*                             &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             IceNH4(0,j,linew)=IceNH4(1,j,liold)
#  ifdef MASKING
             IceNH4(0,j,linew)=IceNH4(0,j,linew)*                             &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             IceNH4(0,j,linew)=IceNH4(0,j,linew)*                             &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ENDIF
        END DO
# elif defined WEST_IceNH4GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          IceNH4(0,j,linew)=IceNH4(1,j,linew)
#  ifdef MASKING
          IceNH4(0,j,linew)=IceNH4(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          IceNH4(0,j,linew)=IceNH4(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          IceNH4(0,j,linew)=IceNH4(1,j,linew)
#  ifdef MASKING
          IceNH4(0,j,linew)=IceNH4(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          IceNH4(0,j,linew)=IceNH4(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# ifdef EAST_IceNH4CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             IceNH4(Lm(ng)+1,j,linew)=BOUNDARY(ng)%IceNH4_east(j)
#  ifdef MASKING
             IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng),j,liold)
#  ifdef MASKING
             IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ENDIF
        END DO
# elif defined EAST_IceNH4GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng),j,linew)
#  ifdef MASKING
          IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng),j,linew)
#  ifdef MASKING
          IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          IceNH4(Lm(ng)+1,j,linew)=IceNH4(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# endif
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

# ifdef SOUTH_IceNH4CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,1,linew).ge.0._r8) THEN
             IceNH4(i,0,linew)=BOUNDARY(ng)%IceNH4_south(i)
#  ifdef MASKING
             IceNH4(i,0,linew)=IceNH4(i,0,linew)*                             &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             IceNH4(i,0,linew)=IceNH4(i,0,linew)*                             &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             IceNH4(i,0,linew)=IceNH4(i,1,liold)
#  ifdef MASKING
             IceNH4(i,0,linew)=IceNH4(i,0,linew)*                             &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             IceNH4(i,0,linew)=IceNH4(i,0,linew)*                             &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ENDIF
        END DO
# elif defined SOUTH_IceNH4GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          IceNH4(i,0,linew)=IceNH4(i,1,linew)
#  ifdef MASKING
          IceNH4(i,0,linew)=IceNH4(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          IceNH4(i,0,linew)=IceNH4(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          IceNH4(i,0,linew)=IceNH4(i,1,linew)
#  ifdef MASKING
          IceNH4(i,0,linew)=IceNH4(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          IceNH4(i,0,linew)=IceNH4(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# ifdef NORTH_IceNH4CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             IceNH4(i,Mm(ng)+1,linew)=BOUNDARY(ng)%IceNH4_north(i)
#  ifdef MASKING
             IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
             IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng),liold)
#  ifdef MASKING
             IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ENDIF
        END DO
# elif defined NORTH_IceNH4GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng),linew)
#  ifdef MASKING
          IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng),linew)
#  ifdef MASKING
          IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          IceNH4(i,Mm(ng)+1,linew)=IceNH4(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# endif
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        IceNH4(0,0,linew)=0.5_r8*(IceNH4(1,0,linew)+                          &
     &                         IceNH4(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        IceNH4(Lm(ng)+1,0,linew)=0.5_r8*(IceNH4(Lm(ng)+1,1,linew)+            &
     &                                IceNH4(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        IceNH4(0,Mm(ng)+1,linew)=0.5_r8*(IceNH4(0,Mm(ng)  ,linew)+            &
     &                                IceNH4(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        IceNH4(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                            &
     &             (IceNH4(Lm(ng)+1,Mm(ng)  ,linew)+                       &
     &              IceNH4(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE IceNH4bc_tile
#endif

      END MODULE IceNH4bc_mod
