#include "cppdefs.h"
      MODULE ice_frazil_mod
#if defined ICE_MODEL && defined ICE_THERMO
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the frazil ice growth in the water when the
!  water temperature gets below freezing. It adjusts the water
!  temperature and salinity accordingly.
!
!  Reference: Steele et al. (1989). JPO, 19, 139-147.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_frazil

      CONTAINS

      SUBROUTINE ice_frazil (ng, tile)

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_stepping

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 44)
# endif
!
      CALL ice_frazil_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nnew(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % wfr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 44)
# endif
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      subroutine ice_frazil_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nnew,                                 &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            Hz, z_r, rho, t, wfr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: wfr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(out) :: wfr(LBi:UBi,LBj:UBj)
# endif
!
! Local variable definitions
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
!    real(r8), allocatable           :: buffer(:)
!    character (len=3), allocatable  :: op_handle(:)
# endif
      integer :: i, j, k, itrc

      real(r8), parameter :: Lhat = 79.2_r8
      real(r8), parameter :: r = 0.5_r8

      real(r8) :: t_freeze
      real(r8) :: s1
      real(r8) :: z1
      real(r8) :: gamma_k
      real(r8) :: t_fr

!ajh new kate update
      real(r8) :: sice
      real(r8) :: ice_dens
      real(r8) :: delta_wfr
      real(r8), parameter :: depressT = -0.0543_r8
      real(r8), parameter :: sice_ref = 3.2_r8
!

!  Inline functions
!  Freezing temperature (Gill, 1982)
!     t_freeze(s1,z1) = -0.0575*s1 + 1.710523d-3*sqrt(s1)**3
!    &       - 2.154996d-4*s1*s1 + 0.000753*z1
!  Freezing temperature (Steele et al. 1989)
      t_freeze(s1,z1) = -0.0543*s1 + 0.000759*z1

# include "set_bounds.h"

!ajh new kate update
      ice_dens = rhoice(ng)
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wfr(i,j) = 0.0_r8
        ENDDO
      ENDDO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
# ifdef MASKING
            IF (rmask(i,j) .ne. 0.0_r8) THEN
# endif
# ifdef WET_DRY
            IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
# endif
!ajh new kate update
              t_fr = t_freeze(t(i,j,k,nnew,isalt),z_r(i,j,k))
              IF (t(i,j,k,nnew,itemp) .lt. t_fr) THEN
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                     (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r)     &
     &                         - depressT * t(i,j,k,nnew,isalt))
                IF (gamma_k .lt. 0.0_r8) THEN
                  print *, 'trouble in ice_frazil', i, j, k,            &
     &             t(i,j,k,nnew,itemp), t(i,j,k,nnew,isalt),            &
     &             t_fr, wfr(i,j), gamma_k, Hz(i,j,k)
                  exit_flag = 1
                END IF
                wfr(i,j) = wfr(i,j) + gamma_k * Hz(i,j,k) *             &
     &                    (rho0 + rho(i,j,k) ) / ice_dens
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) + gamma_k *   &
     &                 (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                sice = MIN(sice_ref,t(i,j,k,nnew,isalt))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &                  (1.0_r8 + gamma_k) - sice*gamma_k
              ELSE IF (wfr(i,j) > 0 .and.                               &
     &                 t(i,j,k,nnew,itemp) .gt. t_fr) THEN
! Use heat at this level to melt some ice from below.
! gamma_k becomes negative here.
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                     (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r)     &
     &                         - depressT * t(i,j,k,nnew,isalt))
                delta_wfr = gamma_k * Hz(i,j,k) *                       &
     &                    (rho0 + rho(i,j,k) ) / rhoice(ng)
                IF ((wfr(i,j) + delta_wfr) > 0) THEN
                  wfr(i,j) = wfr(i,j) + delta_wfr
                ELSE
                  gamma_k = -wfr(i,j) * rhoice(ng) /                    &
     &                     (Hz(i,j,k)*(rho0+rho(i,j,k)))
                  wfr(i,j) = 0.0_r8
                ENDIF
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) + gamma_k *   &
     &                 (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &                                  (1.0_r8 + gamma_k)
!end new kate update
              END IF
# ifdef WET_DRY
            END IF
# endif
# ifdef MASKING
            END IF
# endif
          END DO
          wfr(i,j) = wfr(i,j)/dt(ng)
          IF (wfr(i,j) .lt. 0.0_r8) THEN
            print *, 'trouble in ice_frazil', i, j,                     &
     &         t(i,j,N(ng),nnew,itemp), t(i,j,N(ng),nnew,isalt),        &
     &         wfr(i,j), gamma_k, Hz(i,j,N(ng)),                        &
     &         t_fr, Lhat,z_r(i,j,N(ng))
          END IF
        END DO
      END DO
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wfr)
#ifdef DISTRIBUTE
!       buffer(1) = exit_flag
!       op_handle(1) = 'MAX'
!       CALL mp_reduce (ng, iNLM, 1, buffer, op_handle)
!       exit_flag = int(buffer(1))
      
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    wfr)

#endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!
!  Apply periodic boundary conditions.
!
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,itrc))
      END DO
# endif
# ifdef DISTRIBUTE
      DO itrc=1,NT(ng)
        CALL mp_exchange3d (ng, tile, iNLM, 1,                          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,itrc))
      END DO
# endif
      RETURN
      END SUBROUTINE ice_frazil_tile

#endif
      END MODULE ice_frazil_mod

