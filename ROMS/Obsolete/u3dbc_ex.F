#include "cppdefs.h"
      MODULE u3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: u3dbc_ex.F 895 2009-01-12 21:06:20Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  U-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: u3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE u3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                      &
     &                 LBi, UBi, LBj, UBj, N(ng)                        &
     &                 nstp(ng), nnew(ng),                              &
     &                 OCEAN(ng) % u)
      RETURN
      END SUBROUTINE u3dbc
!
!***********************************************************************
      SUBROUTINE u3dbc_tile (ng, Istr, Iend, Jstr, Jend,                &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       nstp, nnew,                                &
     &                       u)
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
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dUde, dUdt, dUdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"

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
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=u(Istr  ,j  ,k,nstp)-                        &
     &                     u(Istr  ,j-1,k,nstp)
            grad(Istr+1,j)=u(Istr+1,j  ,k,nstp)-                        &
     &                     u(Istr+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Istr+1,j,k,nstp)-u(Istr+1,j,k,nnew)
            dUdx=u(Istr+1,j,k,nstp)-u(Istr+2,j,k,nstp)
#   ifdef WEST_M3NUDGING
            tau=M3obc_out(ng,iwest)
            IF ((dUdt*dUdx).lt.0.0_r8) tau=M3obc_in(ng,iwest)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
              dUde=grad(Istr+1,j  )
            ELSE
              dUde=grad(Istr+1,j+1)
            END IF
            cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=MIN(1.0_r8,cff*dUdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dUde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            u(Istr,j,k,nnew)=(1.0_r8-Cx)*u(Istr,j,k,nstp)+              &
     &                       Cx*u(Istr+1,j,k,nstp)-                     &
     &                       MAX(Ce,0.0_r8)*grad(Istr,j  )-             &
     &                       MIN(Ce,0.0_r8)*grad(Istr,j+1)
#   ifdef WEST_M3NUDGING
            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)+                          &
     &                       tau*(BOUNDARY(ng)%u_west(j,k)-             &
     &                            u(Istr,j,k,nstp))
#   endif
#   ifdef MASKING
            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                          &
     &                         GRID(ng)%umask(Istr,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Istr,j,k,nnew)=BOUNDARY(ng)%u_west(j,k)
#   ifdef MASKING
            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                          &
     &                       GRID(ng)%umask(Istr,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Istr,j,k,nnew)=u(Istr+1,j,k,nnew)
#   ifdef MASKING
            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                          &
     &                       GRID(ng)%umask(Istr,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Istr,j,k,nnew)=0.0_r8
          END DO
        END DO
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
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=u(Iend  ,j  ,k,nstp)-                        &
     &                     u(Iend  ,j-1,k,nstp)
            grad(Iend+1,j)=u(Iend+1,j  ,k,nstp)-                        &
     &                     u(Iend+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Iend,j,k,nstp)-u(Iend  ,j,k,nnew)
            dUdx=u(Iend,j,k,nstp)-u(Iend-1,j,k,nstp)
#   ifdef EAST_M3NUDGING
            tau=M3obc_out(ng,ieast)
            IF ((dUdt*dUdx).lt.0.0_r8) tau=M3obc_in(ng,ieast)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dUde=grad(Iend,j  )
            ELSE
              dUde=grad(Iend,j+1)
            END IF
            cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=MIN(1.0_r8,cff*dUdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dUde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            u(Iend+1,j,k,nnew)=(1.0_r8-Cx)*u(Iend+1,j,k,nstp)+          &
     &                         Cx*u(Iend,j,k,nstp)-                     &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j  )-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j+1)
#   ifdef EAST_M3NUDGING
            u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)+                      &
     &                         tau*(BOUNDARY(ng)%u_east(j,k)-           &
     &                              u(Iend+1,j,k,nstp))
#   endif
#   ifdef MASKING
            u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          u(Iend+1,j,k,nnew)=BOUNDARY(ng)%u_east(j,k)
#   ifdef MASKING
          u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                        &
     &                       GRID(ng)%umask(Iend+1,j)
#   endif
        END DO
      END DO

#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          u(Iend+1,j,k,nnew)=u(Iend,j,k,nnew)
#   ifdef MASKING
          u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                        &
     &                       GRID(ng)%umask(Iend+1,j)
#   endif
        END DO
      END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Iend+1,j,k,nnew)=0.0_r8
          END DO
        END DO
#  endif
      END IF
# endif

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
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-                        &
     &                     u(i  ,Jstr-1,k,nstp)
            grad(i,Jstr  )=u(i+1,Jstr  ,k,nstp)-                        &
     &                     u(i  ,Jstr  ,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Jstr,k,nstp)-u(i,Jstr  ,k,nnew)
            dUde=u(i,Jstr,k,nstp)-u(i,Jstr+1,k,nstp)
#   ifdef SOUTH_M3NUDGING
            tau=M3obc_out(ng,isouth)
            IF ((dUdt*dUde).lt.0.0_r8) tau=M3obc_in(ng,isouth)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Jstr)
            ELSE
              dUdx=grad(i  ,Jstr)
            END IF
            cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dUdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dUde)
            u(i,Jstr-1,k,nnew)=(1.0_r8-Ce)*u(i,Jstr-1,k,nstp)+          &
     &                         Ce*u(i,Jstr,k,nstp)-                     &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jstr-1)
#   ifdef SOUTH_M3NUDGING
            u(i,Jstr-1,k,nnew)=u(i,Jstr-1,k,nnew)+                      &
     &                         tau*(BOUNDARY(ng)%u_south(i,k)-          &
     &                              u(i,Jstr-1,k,nstp))
#   endif
#   ifdef MASKING
            u(i,Jstr-1,k,nnew)=u(i,Jstr-1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jstr-1,k,nnew)=BOUNDARY(ng)%u_south(i,k)
#   ifdef MASKING
            u(i,Jstr-1,k,nnew)=u(i,Jstr-1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jstr-1,k,nnew)=u(i,Jstr,k,nnew)
#   ifdef MASKING
            u(i,Jstr-1,k,nnew)=u(i,Jstr-1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
            u(i,Jstr-1,k,nnew)=gamma2(ng)*u(i,Jstr,k,nnew)
#   ifdef MASKING
            u(i,Jstr-1,k,nnew)=u(i,Jstr-1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
        END DO
#   undef I_RANGE
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
          DO i=IstrU-1,Iend
            grad(i,Jend  )=u(i+1,Jend  ,k,nstp)-                        &
     &                     u(i  ,Jend  ,k,nstp)
            grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-                        &
     &                     u(i  ,Jend+1,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Jend,k,nstp)-u(i,Jend  ,k,nnew)
            dUde=u(i,Jend,k,nstp)-u(i,Jend-1,k,nstp)
#   ifdef NORTH_M3NUDGING
            tau=M3obc_out(ng,inorth)
            IF ((dUdt*dUde).lt.0.0_r8) tau=M3obc_in(ng,inorth)
            tau=tau*dt(ng)
#   endif
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Jend)
            ELSE
              dUdx=grad(i  ,Jend)
            END IF
            cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dUdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dUde)
            u(i,Jend+1,k,nnew)=(1.0_r8-Ce)*u(i,Jend+1,k,nstp)+          &
     &                         Ce*u(i,Jend,k,nstp)-                     &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jend+1)
#   ifdef NORTH_M3NUDGING
            u(i,Jend+1,k,nnew)=u(i,Jend+1,k,nnew)+                      &
     &                         tau*(BOUNDARY(ng)%u_north(i,k)-          &
     &                              u(i,Jend+1,k,nstp))
#    endif
#   ifdef MASKING
            u(i,Jend+1,k,nnew)=u(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jend+1,k,nnew)=BOUNDARY(ng)%u_north(i,k)
#   ifdef MASKING
            u(i,Jend+1,k,nnew)=u(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jend+1,k,nnew)=u(i,Jend,k,nnew)
#   ifdef MASKING
            u(i,Jend+1,k,nnew)=u(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
            u(i,Jend+1,k,nnew)=gamma2(ng)*u(i,Jend,k,nnew)
#   ifdef MASKING
            u(i,Jend+1,k,nnew)=u(i,Jend+1,k,nnew)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
        END DO
#   undef I_RANGE
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
          u(Istr,Jstr-1,k,nnew)=0.5_r8*(u(Istr+1,Jstr-1,k,nnew)+        &
     &                                  u(Istr  ,Jstr  ,k,nnew))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          u(Iend+1,Jstr-1,k,nnew)=0.5_r8*(u(Iend  ,Jstr-1,k,nnew)+      &
     &                                    u(Iend+1,Jstr  ,k,nnew))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          u(Istr,Jend+1,k,nnew)=0.5_r8*(u(Istr+1,Jend+1,k,nnew)+        &
     &                                  u(Istr  ,Jend  ,k,nnew))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          u(Iend+1,Jend+1,k,nnew)=0.5_r8*(u(Iend  ,Jend+1,k,nnew)+      &
     &                                    u(Iend+1,Jend  ,k,nnew))
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE u3dbc_tile
#endif
      END MODULE u3dbc_mod
