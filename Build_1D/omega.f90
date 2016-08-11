      MODULE omega_mod
!
!svn $Id: omega.F 1039 2009-08-11 22:52:28Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: omega, scale_omega
      CONTAINS
!
!***********************************************************************
      SUBROUTINE omega (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private storage
!  arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J)- and MAX(I,J)-directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 13)
      CALL omega_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
     &                 GRID(ng) % z_w,                                  &
     &                 OCEAN(ng) % W)
      CALL wclock_off (ng, iNLM, 13)
      RETURN
      END SUBROUTINE omega
!
!***********************************************************************
      SUBROUTINE omega_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Huon, Hvom,                                &
     &                       z_w,                                       &
     &                       W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(out) :: W(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), dimension(IminS:ImaxS) :: wrk
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
!------------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence.
!------------------------------------------------------------------------
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          W(i,j,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)-                                        &
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))
          END DO
        END DO
!
        DO i=Istr,Iend
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionally coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-                                          &
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
      RETURN
      END SUBROUTINE omega_tile
!
!***********************************************************************
      SUBROUTINE scale_omega (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &                        pm, pn, W, Wscl)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: W(LBi:,LBj:,LBk:)
      real(r8), intent(out) :: Wscl(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
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
!-----------------------------------------------------------------------
!  Scale omega vertical velocity to m/s.
!-----------------------------------------------------------------------
!
      DO k=LBk,UBk
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Wscl(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE scale_omega
      END MODULE omega_mod
