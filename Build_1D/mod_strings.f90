      MODULE mod_strings
!
!svn $Id: mod_strings.F 933 2009-02-24 19:25:01Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  cdt         F90/F95 compiler used.                                  !
!  fflags      F90/F95 compiler flags.                                 !
!  title       Title of model run.                                     !
!  Coptions    Activated C-preprocessing options.                      !
!  StateMsg    Model state processing messages.                        !
!  Pregion     Model regions identifiers used for time profiling.      !
!                                                                      !
!=======================================================================
!
        implicit none
        character (len=80)   :: title
        character (len=2048) :: Coptions
        integer, parameter :: Nregion = 54
        character (len=55), dimension(11) :: StateMsg =                 &
     &    (/'Read state initial conditions,               ',            &
     &      'Read previous state initial conditions,      ',            &
     &      'Read previous adjoint state solution,        ',            &
     &      'Read latest adjoint state solution,          ',            &
     &      'Read initial/model normalization factors,    ',            &
     &      'Read correlation standard deviation,         ',            &
     &      'Read impulse forcing,                        ',            &
     &      'Read v-space increments,                     ',            &
     &      'Read background state,                       ',            &
     &      'Read boundary normalization factors,         ',            &
     &      'Read forcing normalization factors,          '/)
        character (len=50), dimension(Nregion) :: Pregion =             &
     &    (/'Initialization ...................................',       &
     &      'OI data assimilation .............................',       &
     &      'Reading of input data ............................',       &
     &      'Processing of input data .........................',       &
     &      'Processing of output time averaged data ..........',       &
     &      'Computation of vertical boundary conditions ......',       &
     &      'Computation of global information integrals ......',       &
     &      'Writing of output data ...........................',       &
     &      'Model 2D kernel ..................................',       &
     &      'Lagrangian floats trajectories ...................',       &
     &      'Tidal forcing ....................................',       &
     &      '2D/3D coupling, vertical metrics .................',       &
     &      'Omega vertical velocity ..........................',       &
     &      'Equation of state for seawater ...................',       &
     &      'Biological module, source/sink terms .............',       &
     &      'Sediment tranport module, source/sink terms ......',       &
     &      'Atmosphere-Ocean bulk flux parameterization ......',       &
     &      'KPP vertical mixing parameterization .............',       &
     &      'GLS vertical mixing parameterization .............',       &
     &      'My2.5 vertical mixing parameterization ...........',       &
     &      '3D equations right-side terms ....................',       &
     &      '3D equations predictor step ......................',       &
     &      'Pressure gradient ................................',       &
     &      'Harmonic mixing of tracers, S-surfaces ...........',       &
     &      'Harmonic mixing of tracers, geopotentials ........',       &
     &      'Harmonic mixing of tracers, isopycnals ...........',       &
     &      'Biharmonic mixing of tracers, S-surfaces .........',       &
     &      'Biharmonic mixing of tracers, geopotentials ......',       &
     &      'Biharmonic mixing of tracers, isopycnals .........',       &
     &      'Harmonic stress tensor, S-surfaces ...............',       &
     &      'Harmonic stress tensor, geopotentials ............',       &
     &      'Biharmonic stress tensor, S-surfaces .............',       &
     &      'Biharmonic stress tensor, geopotentials ..........',       &
     &      'Corrector time-step for 3D momentum ..............',       &
     &      'Corrector time-step for tracers ..................',       &
     &      'Two-way Atmosphere-Ocean models coupling .........',       &
     &      'Bottom boundary layer module .....................',       &
     &      'GST Analysis eigenproblem solution ...............',       &
     &      'Message Passage: 2D halo exchanges ...............',       &
     &      'Message Passage: 3D halo exchanges ...............',       &
     &      'Message Passage: 4D halo exchanges ...............',       &
     &      'Message Passage: data broadcast ..................',       &
     &      'Message Passage: data reduction ..................',       &
     &      'Message Passage: data gathering ..................',       &
     &      'Message Passage: data scattering..................',       &
     &      'Message Passage: boundary data gathering .........',       &
     &      'Message Passage: point data gathering ............',       &
     &      'Message Passage: multi-model coupling ............',       &
     &      'Ice thermodynamics................................',       &
     &      'Ice rheology coefficients.........................',       &
     &      'Generate coefficients for ice dynamics solver.....',       &
     &      'Generate RHS for ice dynamics solver..............',       &
     &      'Iterative solver of ice dynamics..................',       &
     &      'Advection of ice tracers..........................'/)
        character (len=80) :: my_os = "Linux"
        character (len=80) :: my_cpu = "x86_64"
        character (len=80) :: my_fort = "pgi"
        character (len=80) :: my_fc = "/opt/pgi/linux86-64/latest/bin/pgf90"
        character (len=160) :: my_fflags = " -O3 -Mfree"
      END MODULE mod_strings
