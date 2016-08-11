      MODULE hmixing_mod
!
!svn $Id: hmixing.F 1023 2009-07-20 21:45:24Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes time-dependent 3D horizontal mixing           !
!  coefficients.                                                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Smagorinsky, J, 1963: General circulation experiments with        !
!      the primitive equations: I. The basic experiment, Mon.          !
!      Wea. Rev., 91, 99-164.                                          !
!                                                                      !
!    Holland, W.R., J.C. Chow, and F.O. Bryan, 1998: Application       !
!      of a Third-Order Upwind Scheme in the NCAR Ocean Model, J.      !
!      Climate, 11, 1487-1493.                                         !
!                                                                      !
!    Webb, D.J., B.A. De Cuevas, and C.S. Richmond, 1998: Improved     !
!      Advection Schemes for Ocean Models, J. Atmos. Oceanic           !
!      Technol., 15, 1171-1187.                                        !
!                                                                      !
!    Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic Friction       !
!      with a Smagorinsky-like Viscosity for Use in Large-Scale        !
!      Eddy-Permitting Ocean Models, Monthly Weather Review, 128,      !
!      2935-2946.                                                      !
!                                                                      !
!    Marchesiello, P., L. Debreu, and Xavien Couvelard, 2008:          !
!      Spurious diapycnal mixing in terrain-following coordinate       !
!      models" advection problem and solutions, DRAFT.                 !
!                                                                      !
!  This routine was adapted from a routine provided by Patrick         !
!  Marchiesello (April 2008).                                          !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: hmixing
      CONTAINS
!
!***********************************************************************
      SUBROUTINE hmixing (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
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
      CALL wclock_on (ng, iNLM, 28)
      CALL hmixing_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   GRID(ng) % omn,                                &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   MIXING(ng) % Hviscosity,                       &
     &                   MIXING(ng) % visc3d_r,                         &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v)
      CALL wclock_off (ng, iNLM, 28)
      RETURN
      END SUBROUTINE hmixing
!
!***********************************************************************
      SUBROUTINE hmixing_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         pm, pn, omn, om_u, on_v,                 &
     &                         Hz, z_r,                                 &
     &                         Hviscosity,                              &
     &                         visc3d_r,                                &
     &                         u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: Hviscosity(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: visc3d_r(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: SmagorCoef = 0.1_r8
      real(r8), parameter :: PecletCoef = 1.0_r8 / 12.0_r8
      real(r8) :: DefRate, cff, clip_diff, clip_scale
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
!  Compute velocity dependent, Smagorinsky (1963) horizontal
!  viscosity coefficients. This are based on the local deformation
!  rate, DefRate, which has a horizontal tension and horizontal
!  shearing strain terms:
!
!     DefRate = SQRT [ (du/dx)^2 + (dvdy)^2 + 0.5 * (dvdx + dudy)^2 ]
!                            tension                shearing strain
!
!  The harmonic viscosity coefficient is computed as:
!
!     Asmag = SmagorCoef * dx * dy * DefRate
!  
!  The biharmonic viscosity coefficient follows Griffies and Hallberg
!  (2000) formulation:
!
!     Bsmag = PecletCoef * (dx * dy)^2 * DefRate
!
!-----------------------------------------------------------------------
!
!  Compute viscosity and clipping scale.
!
      clip_scale=0.01_r8*grdmax(ng)**3
      DO k=1,N(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!
!  Compute local deformation rate at RHO-points.
!
            DefRate=SQRT(((u(i+1,j,k,nrhs)-                             &
                           u(i  ,j,k,nrhs))*pm(i,j))**2+                &
     &                   ((v(i,j+1,k,nrhs)-                             &
     &                     v(i,j  ,k,nrhs))*pn(i,j))**2+                &
     &                   0.5_r8*(0.25_r8*pn(i,j)*                       &
     &                           (u(i  ,j+1,k,nrhs)+                    &
     &                            u(i+1,j+1,k,nrhs)-                    &
     &                            u(i  ,j-1,k,nrhs)-                    &
     &                            u(i+1,j-1,k,nrhs))+                   &
     &                           0.25_r8*pm(i,j)*                       &
     &                           (v(i+1,j  ,k,nrhs)+                    &
     &                            v(i+1,j+1,k,nrhs)-                    &
     &                            v(i-1,j  ,k,nrhs)-                    &
     &                            v(i-1,j+1,k,nrhs)))**2)
!
!  Smagorinsky viscosity.
!
            visc3d_r(i,j,k)=Hviscosity(i,j)+                            &
     &                      SmagorCoef*omn(i,j)*DefRate
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply boundary conditions
!-----------------------------------------------------------------------
!
!  Periodic boundary conditions.
!
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        visc3d_r)
      RETURN
      END SUBROUTINE hmixing_tile
      END MODULE hmixing_mod
