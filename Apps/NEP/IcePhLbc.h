#include "cppdefs.h"
      MODULE IcePhLbc_mod
#ifdef ICE_MODEL
!
!***********************************************************************
!  Compute lateral boundary conditions for ice algae
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC IcePhLbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE IcePhLbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL IcePhLbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 liold(ng), linew(ng),                            &
     &                 ICE(ng) % ui,                                    &
     &                 ICE(ng) % vi,                                    &
     &                 ICE(ng) % IcePhL)
      RETURN
      END SUBROUTINE IcePhLbc

!
!***********************************************************************
      SUBROUTINE IcePhLbc_tile (ng, tile,                                  &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           ui, vi, IcePhL)
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
      real(r8), intent(inout) :: IcePhL(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IcePhL(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_IcePhLCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             IcePhL(0,j,linew)=BOUNDARY(ng)%IcePhL_west(j)
#  ifdef MASKING
             IcePhL(0,j,linew)=IcePhL(0,j,linew)*                             &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             IcePhL(0,j,linew)=IcePhL(0,j,linew)*                             &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             IcePhL(0,j,linew)=IcePhL(1,j,liold)
#  ifdef MASKING
             IcePhL(0,j,linew)=IcePhL(0,j,linew)*                             &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             IcePhL(0,j,linew)=IcePhL(0,j,linew)*                             &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ENDIF
        END DO
# elif defined WEST_IcePhLGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          IcePhL(0,j,linew)=IcePhL(1,j,linew)
#  ifdef MASKING
          IcePhL(0,j,linew)=IcePhL(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          IcePhL(0,j,linew)=IcePhL(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          IcePhL(0,j,linew)=IcePhL(1,j,linew)
#  ifdef MASKING
          IcePhL(0,j,linew)=IcePhL(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          IcePhL(0,j,linew)=IcePhL(0,j,linew)*                                &
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

# ifdef EAST_IcePhLCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             IcePhL(Lm(ng)+1,j,linew)=BOUNDARY(ng)%IcePhL_east(j)
#  ifdef MASKING
             IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng),j,liold)
#  ifdef MASKING
             IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ENDIF
        END DO
# elif defined EAST_IcePhLGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng),j,linew)
#  ifdef MASKING
          IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng),j,linew)
#  ifdef MASKING
          IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          IcePhL(Lm(ng)+1,j,linew)=IcePhL(Lm(ng)+1,j,linew)*                  &
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

# ifdef SOUTH_IcePhLCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,1,linew).ge.0._r8) THEN
             IcePhL(i,0,linew)=BOUNDARY(ng)%IcePhL_south(i)
#  ifdef MASKING
             IcePhL(i,0,linew)=IcePhL(i,0,linew)*                             &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             IcePhL(i,0,linew)=IcePhL(i,0,linew)*                             &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             IcePhL(i,0,linew)=IcePhL(i,1,liold)
#  ifdef MASKING
             IcePhL(i,0,linew)=IcePhL(i,0,linew)*                             &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             IcePhL(i,0,linew)=IcePhL(i,0,linew)*                             &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ENDIF
        END DO
# elif defined SOUTH_IcePhLGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          IcePhL(i,0,linew)=IcePhL(i,1,linew)
#  ifdef MASKING
          IcePhL(i,0,linew)=IcePhL(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          IcePhL(i,0,linew)=IcePhL(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          IcePhL(i,0,linew)=IcePhL(i,1,linew)
#  ifdef MASKING
          IcePhL(i,0,linew)=IcePhL(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          IcePhL(i,0,linew)=IcePhL(i,0,linew)*                                &
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

# ifdef NORTH_IcePhLCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             IcePhL(i,Mm(ng)+1,linew)=BOUNDARY(ng)%IcePhL_north(i)
#  ifdef MASKING
             IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
             IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng),liold)
#  ifdef MASKING
             IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*               &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ENDIF
        END DO
# elif defined NORTH_IcePhLGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng),linew)
#  ifdef MASKING
          IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng),linew)
#  ifdef MASKING
          IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          IcePhL(i,Mm(ng)+1,linew)=IcePhL(i,Mm(ng)+1,linew)*                  &
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
        IcePhL(0,0,linew)=0.5_r8*(IcePhL(1,0,linew)+                          &
     &                         IcePhL(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        IcePhL(Lm(ng)+1,0,linew)=0.5_r8*(IcePhL(Lm(ng)+1,1,linew)+            &
     &                                IcePhL(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        IcePhL(0,Mm(ng)+1,linew)=0.5_r8*(IcePhL(0,Mm(ng)  ,linew)+            &
     &                                IcePhL(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        IcePhL(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                            &
     &             (IcePhL(Lm(ng)+1,Mm(ng)  ,linew)+                       &
     &              IcePhL(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE IcePhLbc_tile
#endif

      END MODULE IcePhLbc_mod
