      MODULE ocean_control_mod
!
!svn $Id: fte_ocean.h 999 2009-06-09 23:48:31Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Finite Time Eigenmodes (FTE) Driver:                      !
!                                                                      !
!  This driver computes the finite time eigenmodes of the propagator   !
!  R(0,t)  linearized about a time  evolving  circulation.  They are   !
!  often referred to as the finite time normal modes and are used to   !
!  measure the asymptotic stability of the circulation.                !
!                                                                      !
!  These  routines  control the  initialization,  time-stepping,  and  !
!  finalization of  ROMS/TOMS  model following ESMF conventions:       !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Moore, A.M. et al., 2004: A comprehensive ocean prediction and    !
!      analysis system based on the tangent linear and adjoint of a    !
!      regional ocean model, Ocean Modelling, 7, 227-258.              !
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
#ifdef AIR_OCEAN 
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

      integer :: ng, thread

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MyCOMM)) THEN
        OCN_COMM_WORLD=MyCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
#endif
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
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
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

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine computes the eigenvectors of the propagator R(0,t)     !
!  for either autonomous  or  non-autonomous ocean circulations. A     !
!  single integration of an arbitrary perturbation state vector "u"    !
!  forward in time over the interval  [0,t]  by the tangent linear     !
!  model: R(0,t)*u.  The eigenspectrum of  R(0,t) is computed with     !
!  the Arnoldi algorithm using ARPACK library:                         !
!                                                                      !
!  Lehoucq, R.B., D.C. Sorensen, and C. Yang, 1997:  ARPACK user's     !
!    guide:  solution  of  large  scale  eigenvalue  problems with     !
!    implicit restarted Arnoldi Methods, Rice University, 140p.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_storage
!
      USE propagator_mod
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: ITERATE, Lcomplex, LwrtGST

      integer :: NconvRitz, i, iter, ng

#ifdef DISTRIBUTE
      real(r8), external :: dlapy2, pdnorm2
#else
      real(r8), external :: dlapy2, dnrm2
#endif

      character (len=55) :: string
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
!  Initialize tangent linear for all grids first in order to compute
!  the size of the state vector, Nstate.  This size is computed in
!  routine "wpoints".
!
      DO ng=1,Ngrids
        CALL tl_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Currently, only non-nested applications are considered.  Otherwise,
!  a different structure for mod_storage is needed. 
!
      NEST_LOOP : DO ng=1,Ngrids

        IF (ng.eq.1) THEN
          CALL allocate_storage (ng)
        END IF
!
!  Initialize various parameters.
!
        Nrun=0

        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LwrtPER(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        nTLM(ng)=ntimes(ng)
!
!-----------------------------------------------------------------------
!  Implicit Restarted Arnoldi Method (IRAM) for the computation of
!  optimal perturbation Ritz eigenfunctions.
!-----------------------------------------------------------------------
!
        ITERATE=.TRUE.

        Lrvec=.TRUE.              ! Compute Ritz vectors
        ido=0                     ! reverse communication flag
        bmat='I'                  ! standard eigenvalue problem
        which='LM'                ! compute NEV largest eigenvalues
        howmany='A'               ! compute NEV Ritz vectors
        info=0                    ! random initial residual vector
        iparam(1)=1               ! exact shifts
        iparam(3)=MaxIterGST      ! maximum number of Arnoldi iterations
        iparam(4)=1               ! block size in the recurrence
        iparam(7)=1               ! type of eigenproblem being solved
!
!  ARPACK debugging parameters.
!
        logfil=stdout             ! output logical unit
        ndigit=-3                 ! number of decimal digits
        msaupd=1                  ! iterations, timings, Ritz
        msaup2=1                  ! norms, Ritz values
        msaitr=0
        mseigt=0
        msapps=0
        msgets=0
        mseupd=0
!
!  Determine size of the eigenproblem (Nsize) and size of work space
!  array SworkL (LworkL).
!
        Nsize=Nend(ng)-Nstr(ng)+1
        Nconv=0

#ifdef CHECKPOINTING
!
!  If restart, read in check pointing data GST restart NetCDF file.
!  Otherwise, create check pointing restart NetCDF file.
!
        IF (LrstGST) THEN
          CALL get_gst (ng, iTLM)
          ido=-2
        ELSE
          CALL def_gst (ng, iTLM)
        END IF
        IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Iterate until either convergence or maximum iterations has been
!  exceeded.
!
        iter=0
!
        ITER_LOOP : DO WHILE (ITERATE)
          iter=iter+1
!
!  Reverse communication interface.
!
#ifdef PROFILE
          CALL wclock_on (ng, iTLM, 38)
#endif
#ifdef DISTRIBUTE
          CALL pdnaupd (OCN_COMM_WORLD,                                 &
     &                  ido, bmat, Nsize, which, NEV, Ritz_tol,         &
     &                  resid(Nstr(ng):), NCV, Bvec(Nstr(ng):,1),       &
     &                  Nsize, iparam, ipntr,                           &
     &                  SworkD, SworkL, LworkL, info)
#else
          CALL dnaupd (ido, bmat, Nsize, which, NEV, Ritz_tol,          &
     &                 resid, NCV, Bvec, Nsize, iparam, ipntr,          &
     &                 SworkD, SworkL, LworkL, info)
#endif
#ifdef PROFILE
          CALL wclock_off (ng, iTLM, 38)
#endif
#ifdef CHECKPOINTING2
!
!  If appropriate, write out check point data into GST restart NetCDF
!  file. Notice that the restart data is always saved if MaxIterGST
!  is reached without convergence. It is also saved when convergence
!  is achieved (ido=99).
!
          IF ((MOD(iter,nGST).eq.0).or.(iter.ge.MaxIterGST).or.         &
              (ido.eq.99)) THEN
            CALL wrt_gst (ng, iTLM)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#endif
!
!  Terminate computations if maximum number of iterations is reached.
!  This will faciliate splitting the analysis in several computational
!  cycles using the restart option.
!
          IF ((iter.ge.MaxIterGST).and.(ido.ne.99)) THEN
            ITERATE=.FALSE.
            EXIT ITER_LOOP
          END IF
!
!  Perform matrix-vector operation:  R`(t,0)u
!
          IF (ABS(ido).eq.1) THEN
            NrecTLM(ng)=0
            tTLMindx(ng)=0
            Nconv=iaup2(4)
            CALL propagator (ng, Nstr(ng), Nend(ng),                   &
     &                       SworkD(ipntr(1):), SworkD(ipntr(2):))
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            IF (info.ne.0) THEN
              IF (Master) THEN
                CALL IRAM_error (info, 1, string)
                WRITE (stdout,10) 'DNAUPD', TRIM(string),               &
     &                            ', info = ', info
              END IF
              RETURN
            ELSE
!
!  Compute Ritz vectors. (The only choice left is IDO=99).
!
              IF (Master) THEN
                WRITE (stdout,20) 'Number of converged Ritz values:',   &
     &                            iparam(5)
                WRITE (stdout,20) 'Number of Arnoldi iterations taken:',&
     &                            iparam(3)
              END IF
#ifdef PROFILE
              CALL wclock_on (ng, iTLM, 38)
#endif
#ifdef DISTRIBUTE
              CALL pdneupd (OCN_COMM_WORLD,                             &
     &                      Lrvec, howmany, select,                     &
     &                      RvalueR, RvalueI, Rvector(Nstr(ng):,1),     &
     &                      Nsize, sigmaR, sigmaI, SworkEV, bmat,       &
     &                      Nsize, which, NEV, Ritz_tol,                &
     &                      resid(Nstr(ng):), NCV, Bvec(Nstr(ng):,1),   &
     &                      Nsize, iparam, ipntr,                       &
     &                      SworkD, SworkL, LworkL, info)
#else
              CALL dneupd (Lrvec, howmany, select,                      &
     &                     RvalueR, RvalueI, Rvector, Nsize,            &
     &                     sigmaR, sigmaI, SworkEV, bmat, Nsize,        &
     &                     which, NEV, Ritz_tol,                        &
     &                     resid, NCV, Bvec,                            &
     &                     Nsize, iparam, ipntr,                        &
     &                     SworkD, SworkL, LworkL, info)
#endif
#ifdef PROFILE
              CALL wclock_off (ng, iTLM, 38)
#endif
              IF (info.ne.0) THEN
                IF (Master) THEN
                  CALL IRAM_error (info, 2, string)
                  WRITE (stdout,10) 'DNEUPD', TRIM(string),             &
     &                              ', info = ', info
                END IF
                RETURN
              ELSE
!
!  Check residuals (Euclidean norm) and Ritz values.  Activate writing
!  of the initial and final perturbation for each eigenvector into
!  tangent history NetCDF file.
!
                Nrun=0
                NrecTLM(ng)=0
                tTLMindx(ng)=0
!!              ndefTLM(ng)=ntimes(ng)     ! one file per eigenvector
                NconvRitz=iparam(5)
                Lcomplex=.TRUE.
!
                DO i=1,NconvRitz
!
!  Ritz value is real.
!
                  IF (RvalueI(i).eq.0.0_r8) THEN
                    CALL propagator (ng, Nstr(ng), Nend(ng),            &
     &                               Rvector(Nstr(ng):,i), SworkD)
                    IF (exit_flag.ne.NoError) RETURN
                    CALL daxpy (Nsize, -RvalueR(i),                     &
     &                          Rvector(Nstr(ng):,i), 1, SworkD, 1)
#ifdef DISTRIBUTE
                    norm(i)=pdnorm2(OCN_COMM_WORLD, Nsize, SworkD, 1)
#else
                    norm(i)=dnrm2(Nstate(ng), SworkD, 1)
#endif
!
!  Ritz value is complex.
!
                  ELSE IF (Lcomplex) THEN
                    CALL propagator (ng, Nstr(ng), Nend(ng),            &
     &                               Rvector(Nstr(ng):,i  ), SworkD)
                    IF (exit_flag.ne.NoError) RETURN
                    CALL daxpy (Nsize, -RvalueR(i),                     &
     &                          Rvector(Nstr(ng):,i  ), 1, SworkD, 1)
                    CALL daxpy (Nsize,  RvalueI(i),                     &
     &                          Rvector(Nstr(ng):,i+1), 1, SworkD, 1)
#ifdef DISTRIBUTE
                    norm(i)=pdnorm2(OCN_COMM_WORLD, Nsize, SworkD, 1)
#else
                    norm(i)=dnrm2(Nstate(ng), SworkD, 1)
#endif
                    CALL propagator (ng, Nstr(ng), Nend(ng),            &
     &                               Rvector(Nstr(ng):,i+1), SworkD)
                    IF (exit_flag.ne.NoError) RETURN
                    CALL daxpy (Nsize, -RvalueI(i),                     &
     &                          Rvector(Nstr(ng):,i  ), 1, SworkD, 1)
                    CALL daxpy (Nsize, -RvalueR(i),                     &
     &                          Rvector(Nstr(ng):,i+1), 1, SworkD, 1)
#ifdef DISTRIBUTE
                    norm(i)=dlapy2(norm(i),                             &
     &                             pdnorm2(OCN_COMM_WORLD,              &
     &                                     Nsize, SworkD, 1))
#else
                    norm(i)=dlapy2(norm(i),                             &
     &                             dnrm2(Nstate(ng), SworkD, 1))
#endif
                    norm(i+1)=norm(i)
                    Lcomplex=.FALSE.
                  ELSE
                    Lcomplex=.TRUE.
                  END IF
                  IF (Master) THEN
                    WRITE (stdout,30) i, norm(i), RvalueR(i), RvalueI(i)
                  END IF
!
!  Write out Ritz eigenvalues and Ritz eigenvector Euclidean norm to
!  NetCDF file(s).
!
                  IF (LwrtTLM(ng)) THEN
                    CALL netcdf_put_fvar (ng, iTLM, TLMname(ng),        &
     &                                    'Ritz_rvalue', RvalueR(i:),   &
     &                                    start = (/i/),                &
     &                                    total = (/1/),                &
     &                                    ncid = ncTLMid(ng))
                    IF (exit_flag.ne. NoError) RETURN

                    CALL netcdf_put_fvar (ng, iTLM, TLMname(ng),        &
     &                                    'Ritz_ivalue', RvalueI(i:),   &
     &                                    start = (/i/),                &
     &                                    total = (/1/),                &
     &                                    ncid = ncTLMid(ng))
                    IF (exit_flag.ne. NoError) RETURN

                    CALL netcdf_put_fvar (ng, iTLM, TLMname(ng),        &
     &                                    'Ritz_norm', norm(i:),        &
     &                                    start = (/i/),                &
     &                                    total = (/1/),                &
     &                                    ncid = ncTLMid(ng))
                    IF (exit_flag.ne. NoError) RETURN
                  END IF
                END DO
              END IF
            END IF
            ITERATE=.FALSE.
          END IF

#ifdef CHECKPOINTING
!
!  If appropriate, write out check point data into GST restart NetCDF
!  file. Notice that the restart data is always saved if MaxIterGST
!  is reached without convergence. It is also saved when convergence
!  is achieved (ido=99).
!
          IF ((MOD(iter,nGST).eq.0).or.(iter.ge.MaxIterGST).or.         &
              ((ido.eq.99).and.LwrtGST)) THEN
            CALL wrt_gst (ng, iTLM)
            IF (exit_flag.ne.NoError) RETURN
            IF (ido.eq.99) LwrtGST=.FALSE.
          END IF
#endif

        END DO ITER_LOOP

      END DO NEST_LOOP
!
 10   FORMAT (/,1x,'Error in ',a,1x,a,a,1x,i5,/)
 20   FORMAT (/,a,1x,i2,/)
 30   FORMAT (4x,i4.4,'-th residual ',1p,e14.6,0p,                      &
     &        '  Ritz values ',1pe14.6,0p,2x,1pe14.6)

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear and adjoint models      !
!  execution.                                                          !
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
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE ROMS_finalize

      SUBROUTINE IRAM_error (info, icall, string)
!
!=======================================================================
!                                                                      !
!  This routine decodes internal error messages from the Implicit      !
!  Restarted Arnoldi Method (IRAM) for the computation of optimal      !
!  perturbation Ritz eigenfunctions.                                   !
!                                                                      !
!=======================================================================
!
!
!  imported variable declarations.
!
      integer, intent(in) :: info, icall

      character (len=*), intent(out) :: string
!
!-----------------------------------------------------------------------
!  Decode error message from IRAM.
!-----------------------------------------------------------------------
!
      IF (info.eq.0)  THEN
        string='Normal exit                                            '
      ELSE IF (info.eq.1) THEN
        IF (icall.eq.1) THEN
          string='Maximum number of iterations taken                   '
        ELSE
          string='Could not reorder Schur vectors                      '
        END IF
      ELSE IF (info.eq.3) THEN
        string='No shifts could be applied during an IRAM cycle        '
      ELSE IF (info.eq.-1) THEN
        string='Nstate must be positive                                '
      ELSE IF (info.eq.-2) THEN
        string='NEV must be positive                                   '
      ELSE IF (info.eq.-3) THEN
        string='NCV must be greater NEV and less than or equal Nstate  '
      ELSE IF (info.eq.-4) THEN
        string='Maximum number of iterations must be greater than zero '
      ELSE IF (info.eq.-5) THEN
        string='WHICH must be one of LM or SM                          '
      ELSE IF (info.eq.-6) THEN
        string='BMAT must be one of I or G                             '
      ELSE IF (info.eq.-7) THEN
        string='Length of private work array SworkL is not sufficient  '
      ELSE IF (info.eq.-8) THEN
        IF (icall.eq.1) THEN
          string='Error return from LAPACK eigenvalue calculation      '
        ELSE
          string='Error in DLAHQR in the Shurn vectors calculation     '
        END IF
      ELSE IF (info.eq.-9) THEN
        IF (icall.eq.1) THEN
          string='Starting vector is zero'
        ELSE
          string='Error in DTREVC in the eigenvectors calculation      '
        END IF
      ELSE IF (info.eq.-10) THEN
        string='IPARAM(7) must be 1, 2, 3, 4, 5                        '
      ELSE IF (info.eq.-11) THEN
        string='IPARAM(7) = 1 and BMAT = G are incompatable            '
      ELSE IF (info.eq.-12) THEN
        IF (icall.eq.1) THEN
          string='IPARAM(1) must be equal to 0 or 1                    '
        ELSE
          string='HOWMANY = S not yet implemented                      '
        END IF
      ELSE IF (info.eq.-13) THEN
        string='HOWMANY must be one of A or P if Lrvec = .TRUE.        '
      ELSE IF (info.eq.-14) THEN
        string='Did not find any eigenvalues to sufficient accuaracy   '
      ELSE IF (info.eq.-15) THEN
        string='Different count of converge Ritz values in DNEUPD      '
      ELSE IF (info.eq.-9999) THEN
        string='Could not build and Arnoldi factorization              '
      END IF

      RETURN
      END SUBROUTINE IRAM_error

      END MODULE ocean_control_mod
