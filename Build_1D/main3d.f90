      SUBROUTINE main3d (ng)
!
!svn $Id: main3d.F 1003 2009-06-19 23:52:29Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for a single time step.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE biology_mod, ONLY : biology
      USE ccsm_flux_mod, ONLY : ccsm_flux
      USE diag_mod, ONLY : diag
      USE hmixing_mod, ONLY : hmixing
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
      USE lmd_vmix_mod, ONLY : lmd_vmix
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE rhs3d_mod, ONLY : rhs3d
      USE set_avg_mod, ONLY : set_avg
      USE set_avg2_mod, ONLY : set_avg2
      USE set_depth_mod, ONLY : set_depth
      USE set_massflux_mod, ONLY : set_massflux
      USE set_vbc_mod, ONLY : set_vbc
      USE set_zeta_mod, ONLY : set_zeta
      USE step2d_mod, ONLY : step2d
      USE step3d_t_mod, ONLY : step3d_t
      USE step3d_uv_mod, ONLY : step3d_uv
      USE wvelocity_mod, ONLY : wvelocity
!
      USE mod_forces
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: my_iif, next_indx1, subs, tile, thread
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations.
!=======================================================================
!
!  Set time indices and time clock.
!
      nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
      nnew(ng)=3-nstp(ng)
      nrhs(ng)=nstp(ng)
      CALL time_string (time(ng), time_code(ng))
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL set_data (ng, tile)
        END DO
      END DO
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ini_zeta (ng, tile, iNLM)
          END DO
        END DO
!
!  Initialize other state variables.
!
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL ini_fields (ng, tile, iNLM)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL set_massflux (ng, tile)
          CALL rho_eos (ng, tile)
          CALL diag (ng, tile)
        END DO
      END DO
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ccsm_flux (ng, tile)
          CALL set_vbc (ng, tile)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL lmd_vmix (ng, tile)
          CALL hmixing (ng, tile)
          CALL omega (ng, tile)
          CALL wvelocity (ng, tile, nstp(ng))
!
!--------------------------------------------------------------------
!  Compute deformation-dependent horizontal eddy viscosities
!  using Smagorinsky(1969) formulation
!--------------------------------------------------------------------
!
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1     ! irreversible loop
          CALL set_zeta (ng, tile)
          CALL set_avg (ng, tile)
          CALL set_avg2 (ng, tile)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL rhs3d (ng, tile)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
      DO my_iif=1,nfast(ng)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
        next_indx1=3-indx1(ng)
        IF (.not.PREDICTOR_2D_STEP(ng)) THEN
          PREDICTOR_2D_STEP(ng)=.TRUE.
          iif(ng)=my_iif
          IF (iif(ng).eq.1) THEN
            kstp(ng)=indx1(ng)
          ELSE
            kstp(ng)=3-indx1(ng)
          END IF
          knew(ng)=3
          krhs(ng)=indx1(ng)
        END IF
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.  No actual time-stepping is
!  performed during the auxiliary (nfast+1) time-step. It is needed
!  to finalize the fast-time averaging of 2D fields, if any, and
!  compute the new time-evolving depths.
!
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL step2d (ng, tile)
          END DO
        END DO
!
!  Set time indices for corrector step.
!
        IF (PREDICTOR_2D_STEP(ng)) THEN
          PREDICTOR_2D_STEP(ng)=.FALSE.
          knew(ng)=next_indx1
          kstp(ng)=3-knew(ng)
          krhs(ng)=3
          IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
        END IF
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
        IF (iif(ng).lt.(nfast(ng)+1)) THEN
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL step2d (ng, tile)
            END DO
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.  This call was moved from "step2d" to here.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL set_depth (ng, tile)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL step3d_uv (ng, tile)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL omega (ng, tile)
         CALL biology (ng, tile)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL step3d_t (ng, tile)
        END DO
      END DO
      RETURN
      END SUBROUTINE main3d
