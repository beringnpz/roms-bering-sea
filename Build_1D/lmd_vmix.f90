      MODULE lmd_vmix_mod
!
!svn $Id: lmd_vmix.F 895 2009-01-12 21:06:20Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine computes the vertical mixing coefficients for       !
!  momentum and tracers  at the ocean surface boundary layer and       !
!  interior using the Large, McWilliams and Doney  (1994) mixing       !
!  scheme.                                                             !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review      !
!      and model with a nonlocal boundary layer parameterization,      !
!      Reviews of Geophysics, 32,363-403.                              !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: lmd_vmix
      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_vmix (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
      USE lmd_skpp_mod, ONLY : lmd_skpp
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
      CALL wclock_on (ng, iNLM, 18)
      CALL lmd_vmix_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng),                                     &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    OCEAN(ng) % rho,                              &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    MIXING(ng) % bvf,                             &
     &                    MIXING(ng) % Akt,                             &
     &                    MIXING(ng) % Akv)
      CALL lmd_skpp (ng, tile)
      CALL lmd_finish (ng, tile)
      CALL wclock_off (ng, iNLM, 18)
      RETURN
      END SUBROUTINE lmd_vmix
!
!***********************************************************************
      SUBROUTINE lmd_vmix_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
     &                          Hz,                                     &
     &                          z_r,                                    &
     &                          rho, u, v,                              &
     &                          bvf, Akt, Akv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      Integer, intent(in) :: nstp
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8) :: cff, lmd_iwm, lmd_iws, nu_sx, nu_sxc, shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Rig
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV
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
! Compute gradient Richardson number.
!-----------------------------------------------------------------------
!
!  Compute gradient Richardson number at horizontal RHO-points and
!  vertical W-points.  If zero or very small velocity shear, bound
!  computation by a large negative value.
!
      DO k=1,N(ng)-1
        DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear2=(cff*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+             &
     &                   u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)))**2+        &
     &             (cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+             &
     &                   v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))**2
            Rig(i,j,k)=bvf(i,j,k)/(shear2+eps)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute "interior" viscosities and diffusivities everywhere as
!  the superposition of three processes: local Richardson number
!  instability due to resolved vertical shear, internal wave
!  breaking, and double diffusion.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
!
!  Compute interior diffusivity due to shear instability mixing.
!
            cff=MIN(1.0_r8,MAX(0.0_r8,Rig(i,j,k))/lmd_Ri0)
            nu_sx=1.0_r8-cff*cff
            nu_sx=nu_sx*nu_sx*nu_sx
!
!  The shear mixing should be also a function of the actual magnitude
!  of the shear, see Polzin (1996, JPO, 1409-1425).
!
            shear2=bvf(i,j,k)/(Rig(i,j,k)+eps)
            cff=shear2*shear2/(shear2*shear2+16.0E-10_r8)
            nu_sx=cff*nu_sx
!
!  Compute interior diffusivity due to wave breaking (Gargett and
!  Holloway.
!
            cff=1.0_r8/SQRT(MAX(bvf(i,j,k),1.0E-7_r8))
            lmd_iwm=1.0E-6_r8*cff
            lmd_iws=1.0E-7_r8*cff
!           lmd_iwm=lmd_nuwm
!           lmd_iws=lmd_nuws
!
! Sum contributions due to internal wave breaking, shear instability
! and convective diffusivity due to shear instability.
!
            Akv(i,j,k)=lmd_iwm+lmd_nu0m*nu_sx
            Akt(i,j,k,itemp)=lmd_iws+lmd_nu0s*nu_sx
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE lmd_vmix_tile
!
!***********************************************************************
      SUBROUTINE lmd_finish (ng, tile)
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
      CALL lmd_finish_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      nstp(ng),                                   &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      MIXING(ng) % bvf,                           &
     &                      MIXING(ng) % Akt,                           &
     &                      MIXING(ng) % Akv)
      RETURN
      END SUBROUTINE lmd_finish
!
!***********************************************************************
      SUBROUTINE lmd_finish_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            nstp,                                 &
     &                            Hz,                                   &
     &                            z_r,                                  &
     &                            bvf, Akt, Akv)
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
      Integer, intent(in) :: nstp
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8) :: cff, lmd_iwm, lmd_iws, nu_sx, nu_sxc, shear2
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
!  Compute "interior" viscosities and diffusivities everywhere as
!  the superposition of three processes: local Richardson number
!  instability due to resolved vertical shear, internal wave
!  breaking, and double diffusion.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
!
!  Compute interior convective diffusivity due to static instability
!  mixing.
!
            cff=MAX(bvf(i,j,k),lmd_bvfcon)
            cff=MIN(1.0_r8,(lmd_bvfcon-cff)/lmd_bvfcon)
            nu_sxc=1.0_r8-cff*cff
            nu_sxc=nu_sxc*nu_sxc*nu_sxc
!
! Sum contributions due to internal wave breaking, shear instability
! and convective diffusivity due to shear instability.
!
            Akv(i,j,k)=Akv(i,j,k)+lmd_nu0c*nu_sxc
            Akt(i,j,k,itemp)=Akt(i,j,k,itemp)+lmd_nu0c*nu_sxc
            Akt(i,j,k,isalt)=Akt(i,j,k,isalt)+lmd_nu0c*nu_sxc
          END DO
        END DO
      END DO
!
!  Apply boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  Akv)
      DO itrc=1,NAT
        CALL bc_w3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    Akt(:,:,:,itrc))
      END DO
      RETURN
      END SUBROUTINE lmd_finish_tile
      END MODULE lmd_vmix_mod
