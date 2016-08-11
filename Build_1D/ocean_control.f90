      MODULE ocean_control_mod
!
!svn $Id: nl_ocean.h 1060 2009-09-12 00:25:38Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Nonlinear model:                                          !
!                                                                      !
!  This driver executes ROMS/TOMS standard nonlinear model.  It        !
!  controls the initialization, time-stepping, and finalization        !
!  of the nonlinear model execution following ESMF conventions:        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize
      CONTAINS
      SUBROUTINE ROMS_initialize (first, MyCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first
      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.
      integer :: ng, thread,ic
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
        END IF
        DO ng=1,Ngrids
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
        END DO
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize model state variables for all nested/composed grids.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL initial (ng)
	IF (exit_flag.ne.NoError) RETURN
        CALL get_data (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Initialize run or ensemble counter.
!
      Nrun=1
!
!  Substract a time-step to model time after initialization because the
!  main time-stepping driver always add a single time-step.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,20) ng, ntstart(ng), ntend(ng)
        END IF
        time(ng)=time(ng)-dt(ng)
      END DO
 10   FORMAT (' Process Information:',/)
 20   FORMAT ('NL ROMS/TOMS: started time-stepping:',                   &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')')
      RETURN
      END SUBROUTINE ROMS_initialize
      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine runs ROMS/TOMS nonlinear model from specified starting !
!  (Tstr) to ending (Tend) time-steps.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, dimension(Ngrids) :: Tstr   ! starting time-step
      integer, dimension(Ngrids) :: Tend   ! ending   time-step
!
!  Local variable declarations.
!
      integer :: ng, my_iic
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------
!
      NEST_LOOP : DO ng=1,Ngrids
        NL_LOOP : DO my_iic=Tstr(ng),Tend(ng)
          iic(ng)=my_iic
          time(ng)=time(ng)+dt(ng)
          tdays(ng)=time(ng)*sec2day
!          print*,'iic in nl_ocean=',iic(ng)
!          print*,'time in nl_ocean=',time
!          print*,'tdays=',tdays
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
          CALL get_data (ng)
          IF (exit_flag.ne.NoError) RETURN
          CALL main3d (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO NL_LOOP
      END DO NEST_LOOP
      RETURN
      END SUBROUTINE ROMS_run
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear model execution.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      DO ng=1,Ngrids
        IF (LwrtRST(ng).and.(exit_flag.eq.1)) THEN
          IF (Master) WRITE (stdout,10)
 10       FORMAT (/,' Blowing-up: Saving latest model state into ',     & 
     &              ' RESTART file',/)
          IF (LcycleRST(ng).and.(NrecRST(ng).ge.2)) THEN
            tRSTindx(ng)=2
            LcycleRST(ng)=.FALSE.
          END IF
          blowup=exit_flag
          exit_flag=NoError
          CALL wrt_rst (ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF
      DO ng=1,Ngrids
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iNLM, 0)
        END DO
      END DO
!
!  Close IO files.
!
      CALL close_io
      RETURN
      END SUBROUTINE ROMS_finalize
      END MODULE ocean_control_mod
