      MODULE uv3dmix_mod
!
!svn $Id: uv3dmix.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes horizontal viscosity of momentum.             !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC uv3dmix2
      CONTAINS
      SUBROUTINE uv3dmix2 (ng, tile)
!
!svn $Id: uv3dmix2_s.h 895 2009-01-12 21:06:20Z kate $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes harmonic mixing of momentum, along constant   !
!  S-surfaces,  from the horizontal divergence of the stress tensor.   !
!  A transverse  isotropy  is  assumed so the stress tensor is split   !
!  into vertical and horizontal subtensors.                            !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Wajsowicz, R.C, 1993: A consistent formulation of the           !
!         anisotropic stress tensor for use in models of the           !
!         large-scale ocean circulation, JCP, 105, 333-338.            !
!                                                                      !
!      Sadourny, R. and K. Maynard, 1997: Formulations of              !
!         lateral diffusion in geophysical fluid dynamics              !
!         models, In Numerical Methods of Atmospheric and              !
!         Oceanic Modelling. Lin, Laprise, and Ritchie,                !
!         Eds., NRC Research Press, 547-556.                           !
!                                                                      !
!      Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic              !
!         friction with a Smagorinsky-like viscosity for               !
!         use in large-scale eddy-permitting ocean models,             !
!         Monthly Weather Rev., 128, 8, 2935-2946.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
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
      CALL wclock_on (ng, iNLM, 30)
      CALL uv3dmix2_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nnew(ng),                           &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % om_p,                              &
     &                    GRID(ng) % om_r,                              &
     &                    GRID(ng) % on_p,                              &
     &                    GRID(ng) % on_r,                              &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pmon_p,                            &
     &                    GRID(ng) % pmon_r,                            &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % pnom_p,                            &
     &                    GRID(ng) % pnom_r,                            &
     &                    MIXING(ng) % visc3d_r,                        &
     &                    COUPLING(ng) % rufrc,                         &
     &                    COUPLING(ng) % rvfrc,                         &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v)
      CALL wclock_off (ng, iNLM, 30)
      RETURN
      END SUBROUTINE uv3dmix2
!
!***********************************************************************
      SUBROUTINE uv3dmix2_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nnew,                             &
     &                          Hz,                                     &
     &                          om_p, om_r, on_p, on_r,                 &
     &                          pm, pmon_p, pmon_r,                     &
     &                          pn, pnom_p, pnom_r,                     &
     &                          visc3d_r,                               &
     &                          rufrc, rvfrc, u, v)
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
      integer, intent(in) :: nrhs, nnew
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: visc3d_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, cff1, cff2
      real(r8) :: visc_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
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
!  Compute horizontal harmonic viscosity along constant S-surfaces.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=Hz(i,j,k)*0.5_r8*                                       &
     &          (pmon_r(i,j)*                                           &
     &           ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-                &
     &            (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-               &
     &           pnom_r(i,j)*                                           &
     &           ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-                &
     &            (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))
            UFx(i,j)=on_r(i,j)*on_r(i,j)*visc3d_r(i,j,k)*cff
            VFe(i,j)=om_r(i,j)*om_r(i,j)*visc3d_r(i,j,k)*cff
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            cff=0.125_r8*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+                    &
     &                    Hz(i-1,j-1,k)+Hz(i,j-1,k))*                   &
     &          (pmon_p(i,j)*                                           &
     &           ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-              &
     &            (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+             &
     &           pnom_p(i,j)*                                           &
     &           ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-              &
     &            (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))
            visc_p=0.25_r8*(visc3d_r(i-1,j-1,k)+visc3d_r(i-1,j,k)+      &
     &                      visc3d_r(i  ,j-1,k)+visc3d_r(i  ,j,k))
            UFe(i,j)=om_p(i,j)*om_p(i,j)*visc_p*cff
            VFx(i,j)=on_p(i,j)*on_p(i,j)*visc_p*cff
          END DO
        END DO
!
! Time-step harmonic, S-surfaces viscosity term. Notice that momentum
! at this stage is HzU and HzV and has m2/s units. Add contribution for
! barotropic forcing terms.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
            cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
     &                   (UFx(i,j  )-UFx(i-1,j))+                       &
     &                   (pm(i-1,j)+pm(i,j))*                           &
     &                   (UFe(i,j+1)-UFe(i  ,j)))
            cff2=dt(ng)*cff*cff1
            rufrc(i,j)=rufrc(i,j)+cff1
            u(i,j,k,nnew)=u(i,j,k,nnew)+cff2
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                           &
     &                   (VFx(i+1,j)-VFx(i,j  ))-                       &
     &                   (pm(i,j-1)+pm(i,j))*                           &
     &                   (VFe(i  ,j)-VFe(i,j-1)))
            cff2=dt(ng)*cff*cff1
            rvfrc(i,j)=rvfrc(i,j)+cff1
            v(i,j,k,nnew)=v(i,j,k,nnew)+cff2
          END DO
        END DO
      END DO K_LOOP
      RETURN
      END SUBROUTINE uv3dmix2_tile
      END MODULE uv3dmix_mod
