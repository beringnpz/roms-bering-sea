#include "cppdefs.h"
      MODULE tl_variability_mod

#if defined TANGENT && defined FOUR_DVAR
# ifdef EW_PERIODIC
#  define IR_RANGE Istr-1,Iend+1
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JR_RANGE Jstr-1,Jend+1
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: tl_variability.F 1015 2009-07-07 21:49:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine converts tangent error correlations, C, to background  !
!  error covariances,  B,  by multiplying by the standard deviations,  !
!  S.                                                                  !
!                                                                      !
!                                                                      !
!  The background error covariance is defined as:                      !
!                                                                      !
!    B = S C S                                                         !
!                                                                      !
!    C = C^(1/2) C^(T/2)                                               !
!                                                                      !
!    C^(1/2) = G L^(1/2) W^(-1/2)                               TLM    !
!    C^(T/2) = W^(-1/2) L^(T/2) G                               ADM    !
!                                                                      !
!  where                                                               !
!                                                                      !
!    B : background-error covariance matrix                            !
!    S : diagonal matrix of background-error standard deviations       !
!    C : symmetric matrix of background-error correlations             !
!    G : normalization coefficients matrix (convert to correlations)   !
!    L : tangent linear and adjoint diffusion operators                !
!    W : diagonal matrix of local area or volume metrics               !
!                                                                      !
!  Here, T/2 denote the transpose of a squared-root factor.            !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteo.    !
!      Soc, 127, 1815-1846.                                            !
!                                                                      !
!======================================================================!
!
      USE mod_kinds

      implicit none

      PRIVATE
      PUBLIC :: tl_variability

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_variability (ng, tile, Linp, Lweak)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL tl_variability_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Linp, Lweak,                            &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % e_t_obc,                 &
     &                          BOUNDARY(ng) % e_u_obc,                 &
     &                          BOUNDARY(ng) % e_v_obc,                 &
#  endif
     &                          BOUNDARY(ng) % e_ubar_obc,              &
     &                          BOUNDARY(ng) % e_vbar_obc,              &
     &                          BOUNDARY(ng) % e_zeta_obc,              &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % e_sustr,                   &
     &                          FORCES(ng) % e_svstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % e_stflx,                   &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % e_t,                        &
     &                          OCEAN(ng) % e_u,                        &
     &                          OCEAN(ng) % e_v,                        &
# else
     &                          OCEAN(ng) % e_ubar,                     &
     &                          OCEAN(ng) % e_vbar,                     &
# endif
     &                          OCEAN(ng) % e_zeta,                     &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % tl_t_obc,                &
     &                          BOUNDARY(ng) % tl_u_obc,                &
     &                          BOUNDARY(ng) % tl_v_obc,                &
#  endif
     &                          BOUNDARY(ng) % tl_ubar_obc,             &
     &                          BOUNDARY(ng) % tl_vbar_obc,             &
     &                          BOUNDARY(ng) % tl_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % tl_ustr,                   &
     &                          FORCES(ng) % tl_vstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % tl_tflux,                  &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
# else
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
# endif
     &                          OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_variability
!
!***********************************************************************
      SUBROUTINE tl_variability_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, LBij, UBij,   &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Linp, Lweak,                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                t_obc_std, u_obc_std, v_obc_std,  &
#  endif
     &                                ubar_obc_std, vbar_obc_std,       &
     &                                zeta_obc_std,                     &
# endif
# ifdef ADJUST_WSTRESS
     &                                sustr_std, svstr_std,             &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                stflx_std,                        &
# endif
# ifdef SOLVE3D
     &                                t_std, u_std, v_std,              &
# else
     &                                ubar_std, vbar_std,               &
# endif
     &                                zeta_std,                         &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                tl_t_obc, tl_u_obc, tl_v_obc,     &
#  endif
     &                                tl_ubar_obc, tl_vbar_obc,         &
     &                                tl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                                tl_ustr, tl_vstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                tl_tflux,                         &
# endif
# ifdef SOLVE3D
     &                                tl_t, tl_u, tl_v,                 &
# else
     &                                tl_ubar, tl_vbar,                 &
# endif
     &                                tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:,:,:,:)
      real(r8), intent(in) :: u_obc_std(LBij:,:,:)
      real(r8), intent(in) :: v_obc_std(LBij:,:,:)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:,:)
      real(r8), intent(in) :: vbar_obc_std(LBij:,:)
      real(r8), intent(in) :: zeta_obc_std(LBij:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u_std(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v_std(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ubar_std(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

# else

#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(in) :: u_obc_std(LBij:UBij,N(ng),4)
      real(r8), intent(in) :: v_obc_std(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: vbar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: zeta_obc_std(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: u_std(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(in) :: v_std(LBi:UBi,LBj:UBj,N(ng),NSA)
#  else
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj,NSA)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:UBi,LBj:UBj,NSA)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
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
# endif
      integer :: i, ib, ir, j, rec
# ifdef SOLVE3D
      integer :: it, k
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine error covariance standard deviation factors to use.
!-----------------------------------------------------------------------
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint
      ELSE
        rec=1                        ! strong constraint
      END IF
!
!-----------------------------------------------------------------------
!  Initial conditions and model error covariance: Multiply tangent
!  linear state by its corresponding standard deviation.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)*zeta_std(i,j,rec)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_zeta(:,:,Linp))
# endif
# ifndef SOLVE3D
!
!  2D U-momentum.
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          tl_ubar(i,j,Linp)=tl_ubar(i,j,Linp)*ubar_std(i,j,rec)
        END DO                                                          
      END DO                                                            
      DO j=JV_RANGE                                                     
        DO i=IV_RANGE                                                   
          tl_vbar(i,j,Linp)=tl_vbar(i,j,Linp)*vbar_std(i,j,rec)
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_ubar(:,:,Linp),                            &
     &                    tl_vbar(:,:,Linp))
#  endif
# endif
# ifdef SOLVE3D
!
!  3D U-momentum. 
!
      DO k=1,N(ng)
        DO j=JU_RANGE
          DO i=IU_RANGE
            tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)*u_std(i,j,k,rec)
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)*v_std(i,j,k,rec)
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_u(:,:,:,Linp),                             &
     &                    tl_v(:,:,:,Linp))
#  endif
!
!  Tangent linear tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              tl_t(i,j,k,Linp,it)=tl_t(i,j,k,Linp,it)*                  &
     &                            t_std(i,j,k,rec,it)
            END DO
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_t(:,:,:,Linp,:))
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Boundary conditions and model error covariance: Multiply tangent
!  linear state by its corresponding standard deviation.
!-----------------------------------------------------------------------
!
!  Free-surface open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isFsur,ng)))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isFsur,ng)).and.WESTERN_EDGE) THEN
            ib=iwest
            DO j=Jstr,Jend
              tl_zeta_obc(j,ib,ir,Linp)=tl_zeta_obc(j,ib,ir,Linp)*      &
     &                                  zeta_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO j=Jstr,Jend
              tl_zeta_obc(j,ib,ir,Linp)=tl_zeta_obc(j,ib,ir,Linp)*      &
     &                                  zeta_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO i=Istr,Iend
              tl_zeta_obc(i,ib,ir,Linp)=tl_zeta_obc(i,ib,ir,Linp)*      &
     &                                  zeta_obc_std(i,ib)
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO i=Istr,Iend
              tl_zeta_obc(i,ib,ir,Linp)=tl_zeta_obc(i,ib,ir,Linp)*      &
     &                                  zeta_obc_std(i,ib)
            END DO
          END IF
#  ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isFsur,ng)) THEN
              CALL mp_exchange2d_bry (ng, tile, iTLM, 1, ib,            &
     &                                LBij, UBij,                       &
     &                                NghostPoints,                     &
     &                                EWperiodic, NSperiodic,           &
     &                                tl_zeta_obc(:,ib,ir,Linp))
            END IF
          END DO
#  endif
        END DO
      END IF
!
!  2D U-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isUbar,ng)))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUbar,ng)).and.WESTERN_EDGE) THEN
            ib=iwest
            DO j=Jstr,Jend
              tl_ubar_obc(j,ib,ir,Linp)=tl_ubar_obc(j,ib,ir,Linp)*      &
     &                                  ubar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO j=Jstr,Jend
              tl_ubar_obc(j,ib,ir,Linp)=tl_ubar_obc(j,ib,ir,Linp)*      &
     &                                  ubar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO i=IstrU,Iend
              tl_ubar_obc(i,ib,ir,Linp)=tl_ubar_obc(i,ib,ir,Linp)*      &
     &                                  ubar_obc_std(i,ib)
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO i=IstrU,Iend
              tl_ubar_obc(i,ib,ir,Linp)=tl_ubar_obc(i,ib,ir,Linp)*      &
     &                                  ubar_obc_std(i,ib)
            END DO
          END IF
#  ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isUbar,ng)) THEN
              CALL mp_exchange2d_bry (ng, tile, iTLM, 1, ib,            &
     &                                LBij, UBij,                       &
     &                                NghostPoints,                     &
     &                                EWperiodic, NSperiodic,           &
     &                                tl_ubar_obc(:,ib,ir,Linp))
            END IF
          END DO
#  endif
        END DO
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isVbar,ng)))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVbar,ng)).and.WESTERN_EDGE) THEN
            ib=iwest
            DO j=JstrV,Jend
              tl_vbar_obc(j,ib,ir,Linp)=tl_vbar_obc(j,ib,ir,Linp)*      &
     &                                  vbar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO j=JstrV,Jend
              tl_vbar_obc(j,ib,ir,Linp)=tl_vbar_obc(j,ib,ir,Linp)*      &
     &                                  vbar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO i=Istr,Iend
              tl_vbar_obc(i,ib,ir,Linp)=tl_vbar_obc(i,ib,ir,Linp)*      &
     &                                  vbar_obc_std(i,ib)
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO i=Istr,Iend
              tl_vbar_obc(i,ib,ir,Linp)=tl_vbar_obc(i,ib,ir,Linp)*      &
     &                                  vbar_obc_std(i,ib)
            END DO
          END IF
#  ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isVbar,ng)) THEN
              CALL mp_exchange2d_bry (ng, tile, iTLM, 1, ib,            &
     &                                LBij, UBij,                       &
     &                                NghostPoints,                     &
     &                                EWperiodic, NSperiodic,           &
     &                                tl_vbar_obc(:,ib,ir,Linp))
            END IF
          END DO
#  endif
        END DO
      END IF

#  ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isUvel,ng)))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isUvel,ng)).and.WESTERN_EDGE) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
                tl_u_obc(j,k,ib,ir,Linp)=tl_u_obc(j,k,ib,ir,Linp)*      &
     &                                   u_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                tl_u_obc(j,k,ib,ir,Linp)=tl_u_obc(j,k,ib,ir,Linp)*      &
     &                                   u_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                tl_u_obc(i,k,ib,ir,Linp)=tl_u_obc(i,k,ib,ir,Linp)*      &
     &                                   u_obc_std(i,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                tl_u_obc(i,k,ib,ir,Linp)=tl_u_obc(i,k,ib,ir,Linp)*      &
     &                                   u_obc_std(i,k,ib)
              END DO
            END DO
          END IF
#   ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isUvel,ng)) THEN
              CALL mp_exchange3d_bry (ng, tile, iTLM, 1, ib,            &
     &                                LBij, UBij, 1, N(ng),             &
     &                                NghostPoints,                     &
     &                                EWperiodic, NSperiodic,           &
     &                                tl_u_obc(:,:,ib,ir,Linp))
            END IF
          END DO
#   endif
        END DO
      END IF
!
!  3D V-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isVvel,ng)))) THEN
        DO ir=1,Nbrec(ng)
          IF ((Lobc(iwest,isVvel,ng)).and.WESTERN_EDGE) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=JstrV,Jend
                tl_v_obc(j,k,ib,ir,Linp)=tl_v_obc(j,k,ib,ir,Linp)*      &
     &                                   v_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.EASTERN_EDGE) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                tl_v_obc(j,k,ib,ir,Linp)=tl_v_obc(j,k,ib,ir,Linp)*      &
     &                                   v_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.SOUTHERN_EDGE) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_v_obc(i,k,ib,ir,Linp)=tl_v_obc(i,k,ib,ir,Linp)*      &
     &                                   v_obc_std(i,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.NORTHERN_EDGE) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                tl_v_obc(i,k,ib,ir,Linp)=tl_v_obc(i,k,ib,ir,Linp)*      &
     &                                   v_obc_std(i,k,ib)
              END DO
            END DO
          END IF
#   ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isVvel,ng)) THEN
              CALL mp_exchange3d_bry (ng, tile, iADM, 1, ib,            &
     &                                LBij, UBij, 1, N(ng),             &
     &                                NghostPoints,                     &
     &                                EWperiodic, NSperiodic,           &
     &                                tl_v_obc(:,:,ib,ir,Linp))
            END IF
          END DO
#   endif
        END DO
      END IF
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (.not.Lweak.and.(ANY(Lobc(:,isTvar(it),ng)))) THEN
          DO ir=1,Nbrec(ng)
            IF ((Lobc(iwest,isTvar(it),ng)).and.WESTERN_EDGE) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_t_obc(j,k,ib,ir,Linp,it)=                          &
     &                                    tl_t_obc(j,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(j,k,ib,it)
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.EASTERN_EDGE) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  tl_t_obc(j,k,ib,ir,Linp,it)=                          &
     &                                    tl_t_obc(j,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(j,k,ib,it)
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.SOUTHERN_EDGE) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_t_obc(i,k,ib,ir,Linp,it)=                          &
     &                                    tl_t_obc(i,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(i,k,ib,it)
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.NORTHERN_EDGE) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  tl_t_obc(i,k,ib,ir,Linp,it)=                          &
     &                                    tl_t_obc(i,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(i,k,ib,it)
                END DO
              END DO
            END IF
#   ifdef DISTRIBUTE
            DO ib=1,4
              IF (Lobc(ib,isTvar(it),ng)) THEN
                CALL mp_exchange3d_bry (ng, tile, iADM, 1, ib,          &
     &                                  LBij, UBij, 1, N(ng),           &
     &                                  NghostPoints,                   &
     &                                  EWperiodic, NSperiodic,         &
     &                                  tl_t_obc(:,:,ib,ir,Linp,it))
              END IF
            END DO
#   endif
          END DO
        END IF
      END DO
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Initial conditions and model error covariance: Multiply tangent
!  linear state by its corresponding standard deviation.
!-----------------------------------------------------------------------

#  ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      IF (.not.Lweak) THEN
        DO ir=1,Nfrec(ng)
          DO j=JU_RANGE
            DO i=IU_RANGE
              tl_ustr(i,j,ir,Linp)=tl_ustr(i,j,ir,Linp)*sustr_std(i,j)
            END DO
          END DO
          DO j=JV_RANGE
            DO i=IV_RANGE
              tl_vstr(i,j,ir,Linp)=tl_vstr(i,j,ir,Linp)*svstr_std(i,j)
            END DO
          END DO
        END DO
#   ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      tl_ustr(:,:,:,Linp),                        &
     &                      tl_vstr(:,:,:,Linp))
#   endif
      END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Tangent linear surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO it=1,NT(ng)
          IF (Lstflux(it,ng)) THEN
            DO ir=1,Nfrec(ng)
              DO j=JR_RANGE
                DO i=IR_RANGE
                  tl_tflux(i,j,ir,Linp,it)=tl_tflux(i,j,ir,Linp,it)*    &
     &                                     stflx_std(i,j,it)
                END DO
              END DO
            END DO
#   ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, Nfrec(ng),       &
     &                          NghostPoints, EWperiodic, NSperiodic,   &
     &                          tl_tflux(:,:,:,Linp,it))
#   endif
          END IF
        END DO
      END IF
#  endif
# endif

# undef IR_RANGE
# undef IU_RANGE
# undef IV_RANGE
# undef JR_RANGE
# undef JU_RANGE
# undef JV_RANGE

      RETURN
      END SUBROUTINE tl_variability_tile
#endif
      END MODULE tl_variability_mod
