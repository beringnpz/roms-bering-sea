#include "cppdefs.h"
      MODULE v3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: v3dbc_ex.F 895 2009-01-12 21:06:20Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  V-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: v3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE v3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng)                        &
     &                 nstp(ng), nnew(ng),                              &
     &                 OCEAN(ng) % v)
      RETURN
      END SUBROUTINE v3dbc
!
!***********************************************************************
      SUBROUTINE v3dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       nstp, nnew,                                &
     &                       v)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: nstp, nnew

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dVde, dVdt, dVdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_M3RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jstr  )=v(i  ,Jstr  ,k,nstp)-                        &
     &                     v(i-1,Jstr  ,k,nstp)
            grad(i,Jstr+1)=v(i  ,Jstr+1,k,nstp)-                        &
     &                     v(i-1,Jstr+1,k,nstp))
          END DO
          DO i=Istr,Iend
            dVdt=v(i,Jstr+1,k,nstp)-v(i,Jstr+1,k,nnew)
            dVde=v(i,Jstr+1,k,nstp)-v(i,Jstr+1,k,nstp)
#   ifdef SOUTH_M3NUDGING
            tau=M3obc_out(ng,isouth)
            IF ((dVdt*dVde).lt.0.0_r8) tau=M3obc_in(ng,isouth)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.0.0_r8) THEN
              dVdx=grad(i  ,Jstr+1)
            ELSE
              dVdx=grad(i+1,Jstr+1)
            END IF
            cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dVdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dVde)
            v(i,Jstr,k,nnew)=(1.0_r8-Ce)*v(i,Jstr,k,nstp)+              &
     &                       Ce*v(i,Jstr+1,k,nstp)-                     &
     &                       MAX(Cx,0.0_r8)*grad(i  ,Jstr)-             &
     &                       MIN(Cx,0.0_r8)*grad(i+1,Jstr)
#   ifdef SOUTH_M3NUDGING
            v(i,Jstr,k,nnew)=v(i,Jstr,k,nnew)+                          &
     &                       tau*(BOUNDARY(ng)%v_south(i,k)-            &
     &                            v(i,Jstr,k,nstp))
#   endif
#   ifdef MASKING
            v(i,Jstr,k,nnew)=v(i,Jstr,k,nnew)*                          &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jstr,k,nnew)=BOUNDARY(ng)%v_south(i,k)
#   ifdef MASKING
            v(i,Jstr,k,nnew)=v(i,Jstr,k,nnew)*                          &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jstr,k,nnew)=v(i,Jstr+1,k,nnew)
#   ifdef MASKING
            v(i,Jstr,k,nnew)=v(i,Jstr,k,nnew)*                          &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jstr,k,nnew)=0.0_r8
          END DO
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M3RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jend  )=v(i  ,Jend  ,k,nstp)-                        &
     &                     v(i-1,Jend  ,k,nstp)
            grad(i,Jend+1)=v(i  ,Jend+1,k,nstp)-                        &
     &                     v(i-1,Jend+1,k,nstp)
          END DO
          DO i=Istr,Iend
            dVdt=v(i,Jend,k,nstp)-v(i,Jend  ,k,nnew)
            dVde=v(i,Jend,k,nstp)-v(i,Jend-1,k,nstp)
#   ifdef NORTH_M3NUDGING
            tau=M3obc_out(ng,inorth)
            IF ((dVdt*dVde).lt.0.0_r8) tau=M3obc_in(ng,inorth)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
              dVdx=grad(i  ,Jend)
            ELSE
              dVdx=grad(i+1,Jend)
            END IF
            cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dVdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dVde)
            v(i,Jend+1,k,nnew)=(1.0_r8-Ce)*v(i,Jend+1,k,nstp)+          &
     &                         Ce*v(i,Jend,k,nstp)-                     &
     &                         MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i+1,Jend+1)
#   ifdef NORTH_M3NUDGING
            v(i,Jend+1,k,nnew)=v(i,Jend+1,k,nnew)+                      &
     &                          tau*(BOUNDARY(ng)%v_north(i,k)-         &
     &                               v(i,Jend+1,k,nstp))
#   endif
#   ifdef MASKING
            v(i,Jend+1,k,nnew)=v(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jend+1,k,nnew)=BOUNDARY(ng)%v_north(i,k)
#   ifdef MASKING
            v(i,Jend+1,k,nnew)=v(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jend+1,k,nnew)=v(i,Jend,k,nnew)
#   ifdef MASKING
            v(i,Jend+1,k,nnew)=v(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jend+1,k,nnew)=0.0_r8
          END DO
        END DO
#  endif
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_M3RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend
            grad(Istr-1,j)=v(Istr-1,j+1,k,nstp)-                         &
     &                     v(Istr-1,j  ,k,nstp)
            grad(Istr  ,j)=v(Istr  ,j+1,k,nstp)-                         &
     &                     v(Istr  ,j  ,k,nstp)
          END DO
          DO j=JstrV,Jend
            dVdt=v(Istr,j,k,nstp)-v(Istr  ,j,k,nnew)
            dVdx=v(Istr,j,k,nstp)-v(Istr+1,j,k,nstp)
#   ifdef WEST_M3NUDGING
            tau=M3obc_out(ng,iwest)
            IF ((dVdt*dVdx).lt.0.0_r8) tau=M3obc_in(ng,iwest)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(Istr,j-1)+grad(Istr,j))).gt.0.0_r8) THEN
              dVde=grad(Istr,j-1)
            ELSE
              dVde=grad(Istr,j  )
            END IF
            cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
            Cx=MIN(1.0_r8,cff*dVdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dVde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            v(Istr-1,j,k,nnew)=(1.0_r8-Cx)*v(Istr-1,j,k,nstp)+          &
     &                         Cx*v(Istr,j,k,nstp)-                     &
     &                         MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-         &
     &                         MIN(Ce,0.0_r8)*grad(Istr-1,j  )
#   ifdef WEST_M3NUDGING
            v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)+                      &
     &                         tau*(BOUNDARY(ng)%v_west(j,k)-           &
     &                              v(Istr-1,j,k,nstp))
#    endif
#   ifdef MASKING
            v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Istr-1,j,k,nnew)=BOUNDARY(ng)%v_west(j,k)
#   ifdef MASKING
            v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Istr-1,j,k,nnew)=v(Istr,j,k,nnew)
#   ifdef MASKING
            v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
            v(Istr-1,j,k,nnew)=gamma2(ng)*v(Istr,j,k,nnew)
#   ifdef MASKING
            v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_M3RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend
            grad(Iend  ,j)=v(Iend  ,j+1,k,nstp)-                        &
     &                     v(Iend  ,j  ,k,nstp)
            grad(Iend+1,j)=v(Iend+1,j+1,k,nstp)-                        &
     &                     v(Iend+1,j  ,k,nstp)
          END DO
          DO j=JstrV,Jend
            dVdt=v(Iend,j,k,nstp)-v(Iend  ,j,k,nnew)
            dVdx=v(Iend,j,k,nstp)-v(Iend-1,j,k,nstp)
#   ifdef EAST_M3NUDGING
            tau=M3obc_out(ng,ieast)
            IF ((dVdt*dVdx).lt.0.0_r8) tau=M3obc_in(ng,ieast)
            tau=tau*dt(ng)
#   endif
            IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
            IF ((dVdt*(grad(Iend,j-1)+grad(Iend,j))).gt.0.0_r8) THEN
              dVde=grad(Iend,j-1)
            ELSE
              dVde=grad(Iend,j  )
            END IF
            cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
            Cx=MIN(1.0_r8,cff*dVdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dVde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            v(Iend+1,j,k,nnew)=(1.0_r8-Cx)*v(Iend+1,j,k,nstp)+          &
     &                         Cx*v(Iend,j,k,nstp)-                     &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j  )
#   ifdef EAST_M3NUDGING
            v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)+                      &
     &                         tau*(BOUNDARY(ng)%v_east(j,k)-           &
     &                              v(Iend+1,j,k,nstp))
#   endif
#   ifdef MASKING
            v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Iend+1,j,k,nnew)=BOUNDARY(ng)%v_east(j,k)
#   ifdef MASKING
            v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Iend+1,j,k,nnew)=v(Iend,j,k,nnew)
#   ifdef MASKING
            v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
            v(Iend+1,j,k,nnew)=gamma2(ng)*v(Iend,j,k,nnew)
#   ifdef MASKING
            v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END DO
#   undef J_RANGE
#  endif
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Istr-1,Jstr,k,nnew)=0.5_r8*(v(Istr-1,Jstr+1,k,nnew)+        &
     &                                  v(Istr  ,Jstr  ,k,nnew))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Iend+1,Jstr,k,nnew)=0.5_r8*(v(Iend  ,Jstr  ,k,nnew)+        &
     &                                  v(Iend+1,Jstr+1,k,nnew))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Istr-1,Jend+1,k,nnew)=0.5_r8*(v(Istr-1,Jend  ,k,nnew)+      &
     &                                    v(Istr  ,Jend+1,k,nnew))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Iend+1,Jend+1,k,nnew)=0.5_r8*(v(Iend+1,Jend  ,k,nnew)+      &
     &                                    v(Iend  ,Jend+1,k,nnew))
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE v3dbc_tile
#endif
      END MODULE v3dbc_mod
