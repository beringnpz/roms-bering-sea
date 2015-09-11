#!/bin/sh
#
# This script provides a template for ROMS .in file.  It should be called 
# as follows:
# 
#   source bering10kTemplate.sh > myfile.in
#
# The following variables should be defined at the time this script is 
# run:
#
#   RBSdir:   Path to roms-bering-sea clone
#
#   ii:       i-direction partition for ROMS grid (NtileI)
#
#   jj:       j-direction partition for ROMS grid (NtileJ)
#
#   NTIMES:   number of time steps (DT = 600 currently)
#
#   NRREC:    restart switch
#
#   DSTART:   timestamp for initial time
#
#   INIFILE:  Input initial conditions file (can be a restart file)
#
#   HISFILE:  Output history file
#
#   AVGNAME:  Output averages file
#
#   sebsbio:  Input biological parameters file name

# The following is the contents of the .in file (with substitutions)
cat <<EOF
!  ROMS/TOMS Standard Input parameters.
!
!==============================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                              ===
!========================================================= Hernan G. Arango ===
!                                                                             !
! Input parameters can be entered in ANY order, provided that the parameter   !
! KEYWORD (usually, upper case) is typed correctly followed by "="  or "=="   !
! symbols. Any comment lines are allowed and must begin with an exclamation   !
! mark (!) in column one.  Comments may  appear to the right of a parameter   !
! specification to improve documentation.  All comments will ignored during   !
! reading.  Blank lines are also allowed and ignored. Continuation lines in   !
! a parameter specification are allowed and must be preceded by a backslash   !
! (\).  In some instances, more than one value is required for a parameter.   !
! If fewer values are provided, the  last value  is assigned for the entire   !
! parameter array.  The multiplication symbol (*),  without blank spaces in   !
! between, is allowed for a parameter specification.  For example, in a two   !
! grids nested application:                                                   !
!                                                                             !
!    AKT_BAK == 2*1.0d-6  2*5.0d-6              ! m2/s                        !
!                                                                             !
! indicates that the first two entries of array AKT_BAK,  in fortran column-  !
! major order, will have the same value of "1.0d-6" for grid 1,  whereas the  !
! next two entries will have the same value of "5.0d-6" for grid 2.           !
!                                                                             !
! In multiple levels of nesting and/or multiple connected domains  step-ups,  !
! "Ngrids" entries are expected for some of these parameters.  In such case,  !
! the order of the entries for a parameter is extremely important.  It  must  !
! follow the same order (1:Ngrids) as in the state variable declaration. The  !
! USER may follow the above guidelines for specifying his/her values.  These  !
! parameters are marked by "==" plural symbol after the KEYWORD.              !
!                                                                             !
!==============================================================================
!

! Application title.

       TITLE = Bering Sea 10km Grid

! C-preprocessing Flag.

    MyAppCPP = NEP5

! Input variable information file name.  This file needs to be processed
! first so all information arrays can be initialized properly.

     VARNAME = ${RBSdir}/Apps/NEP/GK_varinfo.dat

! Grid dimension parameters. See notes below in the Glossary for how to set
! these parameters correctly.

         Lm == 180           ! Number of I-direction INTERIOR RHO-points
         Mm == 256           ! Number of J-direction INTERIOR RHO-points
          N == 10            ! Number of vertical levels

       Nbed =  0             ! Number of sediment bed layers

        NAT =  2             ! Number of active tracers (usually, 2)
        NPT =  190           ! Number of inactive passive tracers
        NCS =  0             ! Number of cohesive (mud) sediment tracers
        NNS =  0             ! Number of non-cohesive (sand) sediment tracers

! Domain decomposition parameters for serial, distributed-memory or
! shared-memory configurations used to determine tile horizontal range
! indices (Istr,Iend) and (Jstr,Jend), [1:Ngrids].

      NtileI == ${ii}            ! I-direction partition 
      NtileJ == ${jj}            ! J-direction partition

! Time-Stepping parameters.

       NTIMES == ${NTIMES}
           DT == 600.
      NDTFAST == 40      

! Model iteration loops parameters.

       ERstr =  1
       ERend =  1
      Nouter =  1
      Ninner =  1
  Nintervals =  1

! Number of eigenvalues (NEV) and eigenvectors (NCV) to compute for the
! Lanczos/Arnoldi problem in the Generalized Stability Theory (GST)
! analysis. NCV must be greater than NEV (see documentation below).

         NEV =  2                               ! Number of eigenvalues
         NCV =  10                              ! Number of eigenvectors

! Input/Output parameters.

       NRREC =  ${NRREC}   !-1 for restart !0 for new
   LcycleRST == T
        NRST == 144 !432 !216  !192
        NSTA == 144
        NFLT == 144 !432 !216  !192
       NINFO == 1

! Output history, average, diagnostic files parameters.

     LDEFOUT == T
        NHIS == 144 
     NDEFHIS == 1440 
      NTSAVG == 1
        NAVG == 1008  
     NDEFAVG == 10080 
      NTSDIA == 1
        NDIA == 10
     NDEFDIA == 0

! Output tangent linear and adjoint models parameters.

   LcycleTLM == F
        NTLM == 72
     NDEFTLM == 0
   LcycleADJ == F
        NADJ == 72
     NDEFADJ == 0

! Output check pointing GST restart parameters.

     LrstGST =  F                               ! GST restart switch
  MaxIterGST =  500                             ! maximun number of iterations
        NGST =  10                              ! check pointing interval

! Relative accuracy of the Ritz values computed in the GST analysis.

    Ritz_tol =  1.0d-15

! Harmonic/biharmonic horizontal diffusion of tracer: [1:NAT+NPT,Ngrids].

        TNU2 == 250*25.0d0                          ! m2/s
        TNU4 == 250*0.0d0                         ! m4/s

! Harmononic/biharmonic, horizontal viscosity coefficient: [Ngrids].

       VISC2 == 25.0d0                          ! m2/s
       VISC4 == 0.0d0                           ! m4/s

! Vertical mixing coefficients for active tracers: [1:NAT+NPT,Ngrids]

     AKT_BAK == 250*1.0d-6                   ! m2/s

! Vertical mixing coefficient for momentum: [Ngrids].

     AKV_BAK == 1.0d-5                          ! m2/s

! Turbulent closure parameters.

     AKK_BAK == 5.0d-6                          ! m2/s
     AKP_BAK == 5.0d-6                          ! m2/s
      TKENU2 == 0.0d0                           ! m2/s
      TKENU4 == 0.0d0                           ! m4/s

! Generic length-scale turbulence closure parameters.

!      GLS_P == 3.0d0                           ! K-epsilon
!      GLS_M == 1.5d0
!      GLS_N == -1.0d0
!   GLS_Kmin == 7.6d-6
!   GLS_Pmin == 1.0d-12

!   GLS_CMU0 == 0.5477d0
!     GLS_C1 == 1.44d0
!     GLS_C2 == 1.92d0
!    GLS_C3M == -0.4d0
!    GLS_C3P == 1.0d0
!   GLS_SIGK == 1.0d0
!   GLS_SIGP == 1.30d0

       GLS_P == -1.0d0                           ! K-omega
       GLS_M == 0.5d0
       GLS_N == -1.0d0
    GLS_Kmin == 7.6d-6
    GLS_Pmin == 1.0d-12

    GLS_CMU0 == 0.5477d0
      GLS_C1 == 0.555d0
      GLS_C2 == 0.833d0
     GLS_C3M == -0.6d0
     GLS_C3P == 1.0d0
    GLS_SIGK == 2.0d0
    GLS_SIGP == 2.0d0

! Constants used in momentum stress computation.
                                     
        RDRG == 3.0d-04                    ! m/s
       RDRG2 == 3.0d-03                    ! nondimensional
         Zob == 0.02d0                     ! m
         Zos == 0.02d0                     ! m

! Height (m) of atmospheric measurements for Bulk fluxes parameterization.

      BLK_ZQ == 10.0d0                     ! air humidity
      BLK_ZT == 10.0d0                     ! air temperature
      BLK_ZW == 10.0d0                     ! winds

! Minimum depth for wetting and drying.

       DCRIT == 0.50d0                     ! m

! Various parameters.

       WTYPE == 1
     LEVSFRC == 15
     LEVBFRC == 1

! Vertical S-coordinates parameters, [1:Ngrids].

     THETA_S == 5.0d0                      ! 0 < THETA_S < 20
     THETA_B == 0.4d0                      ! 0 < THETA_B < 1
      TCLINE == 10.0d0                     ! m

! Mean Density and Brunt-Vaisala frequency.

        RHO0 =  1025.0d0                   ! kg/m3
     BVF_BAK =  1.0d-4                     ! 1/s2

! Time-stamp assigned for model initialization, reference time
! origin for tidal forcing, and model reference time for output
! NetCDF units attribute.

      DSTART =  ${DSTART}    
  TIDE_START =  -693962.0d0                ! days
    TIME_REF =  19000101.0d0               ! yyyymmdd.dd

! Nudging/relaxation time scales, inverse scales will be computed
! internally, [1:Ngrids].

       TNUDG == 360.0d0                   ! days
       ZNUDG == 360.0d0                      ! days
      M2NUDG == 360.0d0                      ! days
      M3NUDG == 360.0d0                      ! days

! Factor between passive (outflow) and active (inflow) open boundary
! conditions, [1:Ngrids]. If OBCFAC > 1, nudging on inflow is stronger
! than on outflow (recommended).

      OBCFAC == 120.0d0                      ! nondimensional

! Linear equation of State parameters:

          R0 == 1027.0d0                   ! kg/m3
          T0 == 10.0d0                     ! Celsius
          S0 == 35.0d0                     ! PSU
       TCOEF == 1.7d-4                     ! 1/Celsius
       SCOEF == 7.6d-4                     ! 1/PSU

! Slipperiness parameter: 1.0 (free slip) or -1.0 (no slip)

      GAMMA2 == 1.0d0

! Starting (DstrS) and ending (DendS) day for adjoint sensitivity forcing.
! DstrS must be less or equal to DendS. If both values are zero, their
! values are reset internally to the full range of the adjoint integration.

       DstrS == 0.0d0                      ! starting day
       DendS == 0.0d0                      ! ending day

! Starting and ending vertical levels of the 3D adjoint state variables
! whose sensitivity is required.

       KstrS == 1                          ! starting level
       KendS == 1                          ! ending level

! Logical switches (TRUE/FALSE) to specify the adjoint state variables
! whose sensitivity is required.

Lstate(isFsur) == F                        ! free-surface
Lstate(isUbar) == F                        ! 2D U-momentum
Lstate(isVbar) == F                        ! 2D V-momentum
Lstate(isUvel) == F                        ! 3D U-momentum
Lstate(isVvel) == F                        ! 3D V-momentum

! Logical switches (TRUE/FALSE) to specify the adjoint state tracer
! variables whose sensitivity is required (NT values are expected).

Lstate(isTvar) == F F                      ! tracers
! Stochastic optimals time decorrelation scale (days) assumed for
! red noise processes.

    SO_decay == 2.0d0                      ! days

! Logical switches (TRUE/FALSE) to specify the state surface forcing
! variable whose stochastic optimals is required.

SOstate(isUstr) == T                       ! surface u-stress
SOstate(isVstr) == T                       ! surface v-stress

! Logical switches (TRUE/FALSE) to specify the surface tracer forcing
! variable whose stochastic optimals is required (NT values are expected).

SOstate(isTsur) == F F                     ! surface tracer flux

! Stochastic optimals surface forcing standard deviation for
! dimensionalization.

SO_sdev(isUstr) == 1.0d0                   ! surface u-stress
SO_sdev(isVstr) == 1.0d0                   ! surface v-stress
SO_sdev(isTsur) == 1.0d0 1.0d0             ! NT surface tracer flux

! Logical switches (TRUE/FALSE) to activate writing of fields into
! HISTORY output file.

Hout(idUvel) == T                          ! 3D U-velocity
Hout(idVvel) == T                          ! 3D V-velocity
Hout(idWvel) == T                          ! 3D W-velocity
Hout(idOvel) == T                          ! omega vertical velocity
Hout(idUbar) == T                          ! 2D U-velocity
Hout(idVbar) == T                          ! 2D V-velocity
Hout(idFsur) == T                          ! free-surface

Hout(idTvar) == T T                        ! temperature and salinity

Hout(idUair) == F                          ! surface U-wind
Hout(idVair) == F                          ! surface V-wind
Hout(idUsms) == T                          ! surface U-stress
Hout(idVsms) == T                          ! surface V-stress
Hout(idUbms) == T                          ! bottom U-stress
Hout(idVbms) == T                          ! bottom V-stress

Hout(idUbrs) == F                          ! bottom U-current stress
Hout(idVbrs) == F                          ! bottom V-current stress
Hout(idUbws) == F                          ! bottom U-wave stress
Hout(idVbws) == F                          ! bottom V-wave stress
Hout(idUbcs) == F                          ! bottom max wave-current U-stress
Hout(idVbcs) == F                          ! bottom max wave-current V-stress
Hout(idUbot) == F                          ! bed wave orbital U-velocity
Hout(idVbot) == F                          ! bed wave orbital V-velocity
Hout(idUbur) == F                          ! bottom U-velocity above bed
Hout(idVbvr) == F                          ! bottom V-velocity above bed

Hout(idTsur) == T T                        ! surface net heat and salt flux
Hout(idLhea) == T                          ! latent heat flux
Hout(idShea) == T                          ! sensible heat flux
Hout(idLrad) == T                          ! longwave radiation flux
Hout(idSrad) == T                          ! shortwave radiation flux
Hout(idevap) == F                          ! evaporation rate
Hout(idrain) == F                          ! precipitation rate

Hout(idDano) == F                          ! density anomaly
Hout(idVvis) == T                          ! vertical viscosity
Hout(idTdif) == T                          ! vertical T-diffusion
Hout(idSdif) == F                          ! vertical Salinity diffusion
Hout(idHsbl) == T                          ! depth of surface boundary layer
Hout(idHbbl) == T                          ! depth of bottom boundary layer
Hout(idMtke) == F                          ! turbulent kinetic energy
Hout(idMtls) == F                          ! turbulent length scale

! Logical switches (TRUE/FALSE) to activate writing of ice prognostic
! variables into HISTORY output file.
  Hout(idUice) == T
 Hout(idVice)  == T
 Hout(idAice)  == T
 Hout(idHice)  == T
 Hout(idTice)  == T
 Hout(idHsno)  == T
 Hout(idTimid) == T
 Hout(idSfwat) == T
 Hout(idTauiw) == T
 Hout(idChuiw) == T
Hout(idAgeice) == T
 Hout(idSig11) == T
 Hout(idSig12) == T
 Hout(idSig22) == T
  Hout(idS0mk) == T
  Hout(idT0mk) == T

! Logical switches (TRUE/FALSE) to activate writing of extra inert passive
! tracers other than biological and sediment tracers. An inert passive tracer
! is one that it is only advected and diffused. Other processes are ignored.
! These tracers include, for example, dyes, pollutants, oil spills, etc.
! NPT values are expected. However, these switches can be activated using
! compact parameter specification.

Hout(inert) == T                          ! inert passive tracers

! Logical switches (TRUE/FALSE) to activate writing of exposed sediment
! layer properties into HISTORY output file.  Currently, MBOTP properties
! are expected for the bottom boundary layer and/or sediment models:
!
!   Hout(idBott(isd50)),  isd50 = 1        ! mean grain diameter
!   Hout(idBott(idens)),  idens = 2        ! mean grain density
!   Hout(idBott(iwsed)),  iwsed = 3        ! mean settling velocity
!   Hout(idBott(itauc)),  itauc = 4        ! critical erosion stress
!   Hout(idBott(irlen)),  irlen = 5        ! ripple length
!   Hout(idBott(irhgt)),  irhgt = 6        ! ripple height
!   Hout(idBott(ibwav)),  ibwav = 7        ! wave excursion amplitude
!   Hout(idBott(izdef)),  izdef = 8        ! default bottom roughness
!   Hout(idBott(izapp)),  izapp = 9        ! apparent bottom roughness
!   Hout(idBott(izNik)),  izNik = 10       ! Nikuradse bottom roughness
!   Hout(idBott(izbio)),  izbio = 11       ! biological bottom roughness
!   Hout(idBott(izbfm)),  izbfm = 12       ! bed form bottom roughness
!   Hout(idBott(izbld)),  izbld = 13       ! bed load bottom roughness
!   Hout(idBott(izwbl)),  izwbl = 14       ! wave bottom roughness
!   Hout(idBott(iactv)),  iactv = 15       ! active layer thickness
!   Hout(idBott(ishgt)),  ishgt = 16       ! saltation height
!
!                                 1 1 1 1 1 1 1
!               1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6

Hout(idBott) == T T T T T T T T T F F F F F F F

! Generic User parameters, [1:NUSER].

       NUSER =  0
        USER =  0.d0

! Input NetCDF file names, [1:Ngrids].


     GRDNAME == /home/aydink/ncfiles/nc_in/Bering_grid_withFeast.nc
     ININAME == ${INIFILE}
     IRPNAME == ocean_irp.nc
     IADNAME == ocean_iad.nc
!    CLMNAME == /wrkdir/blamey/Bering_10k/Bering10k_FeClm.nc
     BRYNAME == /home/aydink/ncfiles/nc_in/run42allyrs_bering10_als_biobry_tracers_new_iron_10layer.nc
!    FWDNAME == /ptmp/enrique/ROMS/NEP4/Inputs/ocean_fwd.nc
     ADSNAME == ocean_ads.nc

! Input forcing NetCDF file name(s).  The USER has the option to enter
! several files names per each nested grid.  For example, the USER may
! have a different files for wind products, heat fluxes, rivers, tides,
! etc.  The model will scan the file list and will read the needed data
! from the first file in the list containing the forcing field. Therefore,
! the order of the file names is very important.  If multiple forcing
! files per grid, enter first all the file names for grid 1, then grid 2,
! and so on.  Use a single line per entry with a continuation (\) symbol
! at the each entry, except the last one.


    NFFILES == 12                          ! number of forcing files

         FRCNAME == /home/aydink/ncfiles/CORE2/Pair.recent.Bering.nc   
                    /home/aydink/ncfiles/CORE2/lwrad.recent.Bering.nc  
                    /home/aydink/ncfiles/CORE2/Qair.recent.Bering.nc  
                    /home/aydink/ncfiles/CORE2/swrad.recent.Bering.nc  
                    /home/aydink/ncfiles/CORE2/tair_all.Bering.nc  
                    /home/aydink/ncfiles/CORE2/Uwind.recent.Bering.nc  
                    /home/aydink/ncfiles/CORE2/Vwind.recent.Bering.nc   
                    /home/aydink/ncfiles/CORE2/runoff.daitren.iaf.10FEB2011.als.nc  
                    /home/aydink/ncfiles/CORE2/rain.1948-2006.Bering.nc  
                    /home/aydink/ncfiles/CORE2/Bering_tides_otps.nc 
                    /home/aydink/ncfiles/CORE2/sss_fill_2004.nc
                    ./catch_weekly_link.nc 
                 
!        FRCNAME == /wrkdir/blamey/Bering_10k/Files/CCSM_1999_Pair.nc   \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_1999_lwrad_down.nc  \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_1999_Qair.nc   \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_1999_swrad.nc  \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_1999_Tair.nc   \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_1999_Uwind.nc  \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_1999_Vwind.nc  \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_Runoff.nc   \
!                   /wrkdir/blamey/Bering_10k/Files/CCSM_rain.nc  \
!                   /wrkdir/blamey/Bering_10k/Files2_26_09/Bering_tides_otps.nc  

! Output NetCDF file names, [1:Ngrids].

!        GSTNAME == ocean_gst.nc
         RSTNAME == /home/aydink/ncfiles/weekly/out_rst.nc
         HISNAME == ${HISFILE}
!        TLMNAME == nep5_tlm.nc
!        TLFNAME == ocean_tlf.nc
!        ADJNAME == nep5_tlm.nc
         AVGNAME == ${AVGFILE}
!        DIANAME == /wrkdir/blamey/BSmulti2_dia.nc
         STANAME == /home/aydink/ncfiles/weekly/out_sta.nc
         FLTNAME == /home/aydink/ncfiles/weekly/out_flt.nc
         DIANAME == /home/aydink/ncfiles/weekly/out_dia.nc

         GSTNAME == /tmp/aydtmp/dumm1_gst.nc
!        RSTNAME == /tmp/aydtmp/dumm2_rst.nc
!        HISNAME == /tmp/aydtmp/dumm3_his.nc
         TLMNAME == /tmp/aydtmp/dumm4_tlm.nc
         TLFNAME == /tmp/aydtmp/dumm5_tlf.nc
         ADJNAME == /tmp/aydtmp/dumm6_tlm.nc
!        AVGNAME == /tmp/aydtmp/dumm7_avg.nc
!        DIANAME == /tmp/aydtmp/dumm8_dia.nc
!        STANAME == /tmp/aydtmp/dumm9_sta.nc
!        FLTNAME == /tmp/aydtmp/dumm10_flt.nc
         
! Input ASCII parameter filenames.

         APARNAM =  External/assimilation.in
         
         SPOSNAM =  Apps/GK_FEAST/GK_stations_bering_10k.in
         IPARNAM =  Apps/Bering_10k/ice.in
         BPARNAM =  ${sebsbio}
         
!        SPOSNAM =  Apps/FEAST/stations_bering_10k.in
!        FPOSNAM =  Apps/FEAST/floatBS_Exp.in
!        IPARNAM =  Apps/FEAST/ice.in
!        BPARNAM =  Apps/FEAST/MERGE_sebs_bio.in
         
         SPARNAM =  External/sediment.in
         USRNAME =  External/MyFile.dat
!
!  GLOSSARY:
!  =========
!
!------------------------------------------------------------------------------
! Application title (string with a maximum of eighty characters) and
! C-preprocessing flag.
!------------------------------------------------------------------------------
!
!  TITLE       Application title.
!
!  MyAppCPP    Application C-preprocession option.
!
!------------------------------------------------------------------------------
! Variable information file name (string with a maximum of eighty characters).
!------------------------------------------------------------------------------
!
!  VARNAME     Input/Output variable information file name.  This file need to
!              be processed first so all information arrays and indices can be
!              initialized properly in "mod_ncparam.F".
!
!------------------------------------------------------------------------------
! Grid dimension parameters.
!------------------------------------------------------------------------------+!
! These parameters are very important since it determine the grid of the
! application to solve. They need to be read first in order to dynamically
! allocate all model variables.
!
! WARNING: It is trivial and posible to change these dimension parameters in
! -------  idealized applications via analytical expressions. However, in
! realistic applications any change to these parameters requires redoing all
! input NetCDF files.
!
!  Lm          Number of INTERIOR grid RHO-points in the XI-direction for
!                each nested grid, [1:Ngrids]. If using NetCDF files as
!                input, Lm=xi_rho-2 where "xi_rho" is the NetCDF file
!                dimension of RHO-points. Recall that all RHO-point
!                variables have a computational I-range of [0:Lm+1].
!
!  Mm          Number of INTERIOR grid RHO-points in the ETA-direction for
!                each nested grid, [1:Ngrids]. If using NetCDF files as
!                input, Mm=eta_rho-2 where "eta_rho" is the NetCDF file
!                dimension of RHO-points. Recall that all RHO-point
!                variables have a computational J-range of [0:Mm+1].
!
!  N           Number of vertical terrain-following levels at RHO-points,
!                [1:Ngrids].
!
!  Nbed        Number of sediment bed layers, [1:Ngrids]. This parameter
!                is only relevant if CPP option SEDIMENT is activated.
!
!                Mm+1  ___________________                _______  Kw = N
!                     |                   |              |       |
!                  Mm |   _____________   |              |       | Kr = N
!                     |  |             |  |              |_______|
!                     |  |             |  |              |       |
!                  Jr |  |             |  |              |       |
!                     |  |             |  |              |_______|
!                     |  |             |  |              |       |
!                   1 |  |_____________|  |              |       |
!                     |                   |              |_______|
!                   0 |___________________|              |       |
!                              Ir                        |       | 1
!                     0  1            Lm  Lm+1    h(i,j) |_______|
!                                                        ::::::::: 0
!                                                        :::::::::
!                                                        ::::::::: Nbed-1
!                                                        ::::::::: Nbed
!
!  NAT         Number of active tracer type variables. Usually, NAT=2 for
!                potential temperature and salinity.
!
!  NPT         Number of inert (dyes, age, etc) passive tracer type variables
!                to advect and diffuse only. This parameter is only relevant
!                if CPP option T_PASSIVE is activated.
!
!  NCS         Number of cohesive (mud) sediment tracer type variables. This
!                parameter is only relevant if CPP option SEDIMENT is
!                activated.
!
!  NNS         Number of non-cohesive (sand) sediment tracer type variables.
!                This parameter is only relevant if CPP option SEDIMENT is
!                activated.
!
!              The total of sediment tracers is NST=NCS+NNS. Notice that
!              NST must be greater than zero (NST>0).
!
!------------------------------------------------------------------------------
! Domain tile partition parameters.
!------------------------------------------------------------------------------
!
! Model tile decomposition parameters for serial and parallel configurations
! which are used to determine tile horizontal range indices (Istr,Iend) and
! (Jstr,Jend). In some computers, it is advantageous to have tile partitions
! in serial applications.
!
!  NtileI      Number of domain partitions in the I-direction (XI-coordinate).
!              It must be equal or greater than one.
!
!  NtileJ      Number of domain partitions in the J-direction (ETA-coordinate).
!              It must be equal or greater than one.
!
!  WARNING:    In shared-memory (OpenMP), the product of NtileI and NtileJ must
!              be a MULTIPLE of the number of parallel threads specified with
!              the OpenMP environmental variable OMP_NUM_THREADS.
!
!              In distributed-memory (MPI), the product of NtileI and NtileJ
!              must be EQUAL to the number of parallel nodes specified during
!              execution with the "mprun" or "mpirun" command.
!
!------------------------------------------------------------------------------
! Time-Stepping parameters.
!------------------------------------------------------------------------------
!
!  NTIMES      Total number time-steps in current run.  If 3D configuration,
!              NTIMES is the total of baroclinic time-steps.  If only 2D
!              configuration, NTIMES is the total of barotropic time-steps.
!
!  DT          Time-Step size in seconds.  If 3D configuration, DT is the
!              size of baroclinic time-step.  If only 2D configuration, DT
!              is the size of the barotropic time-step.
!
!  NDTFAST     Number of barotropic time-steps between each baroclinic time
!              step. If only 2D configuration, NDTFAST should be unity since
!              there is not need to splitting time-stepping.
!
!------------------------------------------------------------------------------
! Model iteration loops parameters.
!------------------------------------------------------------------------------
!
!  ERstr       Starting ensemble run (perturbation or iteration) number.
!
!  ERend       Ending   ensemble run (perturbation or iteration) number.
!
!  Nouter      Maximum number of 4DVAR outer loop iterations.
!
!  Ninner      Maximum number of 4DVAR inner loop iterations.
!
!  Nintervals  Number of time interval divisions for stochastic optimals
!              computations. It must be a multiple of NTIMES. The tangent
!              linear model (TLM) and the adjoint model (ADM) are integrated
!              forward and backward in different intervals.  For example,
!              if Nintervals=3,
!
!              1               NTIMES/3         2*NTIMES/3           NTIMES
!              +..................+..................+..................+
!              <========================================================> (1)
!                                 <=====================================> (2)
!                                                    <==================> (3)
!
!              In the first iteration (1), the TLM is integrated forward from
!              1 to NTIMES and the ADM is integrated backward from NTIMES to 1.
!              In the second iteration (2), the TLM is integrated forward from
!              NTIMES/3 to NTIMES and the ADM is integrated backward from
!              NTIMES to NTIMES/3. And so on.
!
!------------------------------------------------------------------------------
!  Eigenproblem parameters.
!------------------------------------------------------------------------------
!
!  NEV         Number of eigenvalues to compute for the Lanczos/Arnoldi
!              problem.  Notice that the model memory requirement increases
!              substantially as NEV increases.  The GST requires NEV+1
!              copies of the model state vector.  The memory requirements
!              are decreased in distributed-memory applications.
!
!  NCV         Number of eigenvectors to compute for the Lanczos/Arnoldi
!              problem. NCV must be greater than NEV.
!
!  At present, there is no a-priori analysis to guide the selection of NCV
!  relative to NEV.  The only formal requirement is that NCV > NEV. However
!  in optimal perturbations, it is recommended to have NCV greater than or
!  equal to 2*NEV. In Finite Time Eigenmodes (FTE) and Adjoint Finite Time
!  Eigenmodes (AFTE) the requirement is to have NCV greater than or equal to
!  2*NEV+1.
!
!  The efficiency of calculations depends critically on the combination of
!  NEV and NCV.  If NEV is large (greater than 10 say), you can use NCV=2*NEV+1
!  but for NEV small (less than 6) it will be inefficient to use NCV=2*NEV+1.
!  In complicated applications, you can start with NEV=2 and NCV=10. Otherwise,
!  it will iterate for very long time.
!
!------------------------------------------------------------------------------
! Input/Output parameters.
!------------------------------------------------------------------------------
!
!  NRREC       Switch to indicate re-start from a previous solution.  Use
!              NRREC=0 for new solutions. In a re-start solution, NRREC
!              is the time index of the re-start NetCDF file assigned for
!              initialization.  If NRREC is negative (said NRREC=-1), the
!              model will re-start from the most recent time record. That
!              is, the initialization record is assigned internally.
!              Notice that it is also possible to re-start from a history
!              or time-averaged NetCDF files.  If a history file is used
!              for re-start, it must contains all the necessary primitive
!              variables at all levels.
!
!  LcycleRST   Logical switch (T/F) used to recycle time records in output
!              re-start file.  If TRUE,  only the latest two re-start time
!              records are maintained.  If FALSE, all re-start fields are
!              saved every NRST time-steps without recycling.  The re-start
!              fields are written at all levels in double precision.
!
!  NRST        Number of time-steps between writing of re-start fields.
!
!  NSTA        Number of time-steps between writing data into stations file.
!              Station data is written at all levels.
!
!  NFLT        Number of time-steps between writing data into floats file.
!
!  NINFO       Number of time-steps between print of single line information
!              to standard output.  If also determines the interval between
!              computation of global energy diagnostics.
!
!------------------------------------------------------------------------------
!  Output history and average files parameters.
!------------------------------------------------------------------------------
!
!  LDEFOUT     Logical switch (T/F) used to create new output files when
!              initializing from a re-start file, abs(NRREC) > 0.  If TRUE
!              and applicable, a new history, average, diagnostic and
!              station files are created during the initialization stage.
!              If FALSE and applicable, data is appended to an existing
!              history, average, diagnostic and station files.  See also
!              parameters NDEFHIS, NDEFAVG and NDEFDIA below.
!
!  NHIS        Number of time-steps between writing fields into history file.
!
!  NDEFHIS     Number of time-steps between the creation of new history file.
!              If NDEFHIS=0, the model will only process one history file.
!              This feature is useful for long simulations when history files
!              get too large; it creates a new file every NDEFHIS time-steps.
!
!  NTSAVG      Starting time-step for the accumulation of output time-averaged
!              data.
!
!  NAVG        Number of time-steps between writing time-averaged data
!              into averages file.  Averaged date is written for all fields.
!
!  NDEFAVG     Number of time-steps between the creation of new average
!              file.  If NDEFAVG=0, the model will only process one average
!              file.  This feature is useful for long simulations when
!              average files get too large; it creates a new file every
!              NDEFAVG time-steps.
!
!  NTSDIA      Starting time-step for the accumulation of output time-averaged
!              diagnostics data.
!
!  NDIA        Number of time-steps between writing time-averaged diagnostics
!              data into diagnostics file.  Averaged date is written for all
!              fields.
!
!  NDEFDIA     Number of time-steps between the creation of new time-averaged
!              diagnostics file.  If NDEFDIA=0, the model will only process one
!              diagnostics file.  This feature is useful for long simulations
!              when diagnostics files get too large; it creates a new file
!              every NDEFDIA time-steps.
!
!------------------------------------------------------------------------------
!  Output tangent linear and adjoint model parameters.
!------------------------------------------------------------------------------
!
!  LcycleTLM   Logical switch (T/F) used to recycle time records in output
!              tangent linear file.  If TRUE, only the latest two time
!              records are maintained.  If FALSE, all tangent linear fields
!              are saved every NTLM time-steps without recycling.
!
!  NTLM        Number of time-steps between writing fields into tangent linear
!              model file.
!
!  NDEFTLM     Number of time-steps between the creation of new tangent linear
!              file. If NDEFTLM=0, the model will only process one tangent
!              linear file. This feature is useful for long simulations when
!              output NetCDF files get too large; it creates a new file every
!              NDEFTLM time-steps.
!
!  LcycleADJ   Logical switch (T/F) used to recycle time records in output
!              adjoint file.  If TRUE, only the latest two time records are
!              maintained.  If FALSE, all tangent linear fields re saved
!              every NADJ time-steps without recycling.
!
!  NADJ        Number of time-steps between writing fields into adjoint model
!              file.
!
!  NDEFADJ     Number of time-steps between the creation of new adjoint file.
!              If NDEFADJ=0, the model will only process one adjoint file.
!              This feature is useful for long simulations when output NetCDF
!              files get too large; it creates a new file every NDEFADJ
!              time-steps.
!
!------------------------------------------------------------------------------
!  Generalized Stability Theory (GST) analysis parameters.
!------------------------------------------------------------------------------
!
!  LrstGST     Logical switch (TRUE/FALSE) to restart GST analysis. If TRUE,
!              the check pointing data is read in from the GST restart NetCDF
!              file.  If FALSE and applicable, the check pointing GST data is
!              saved and overwritten every NGST iterations of the algorithm.
!
!  MaxIterGST  Maximum number of GST algorithm iterations.
!
!  NGST        Number of GST iterations between storing of check pointing
!              data into NetCDF file. The restart data is always saved if
!              MaxIterGST is reached without convergence. It is also saved
!              when convergence is achieved. It is always a good idea to
!              save the check pointing data at regular intervals so there
!              is a mechanism to recover from an unexpected interruption
!              in this very expensive computation. The check pointing data
!              can be used also to recompute the Ritz vectors by changing
!              some of the parameters, like convergence criteria (Ritz_tol)
!              and number of Arnoldi iterations (iparam(3)).
!
!  Ritz_tol    Relative accuracy of the Ritz values computed in the GST
!              analysis.
!
!------------------------------------------------------------------------------
! Harmonic/Biharmonic horizontal diffusion for active tracers.
!------------------------------------------------------------------------------
!
!  TNU2        Lateral, harmonic, constant, mixing coefficient (m2/s) for
!              active (NAT) and inert (NPT) tracer variables.  If variable
!              horizontal diffusion is activated, TNU2 is the mixing
!              coefficient for the largest grid-cell in the domain.
!
!  TNU4        Lateral, biharmonic, constant, mixing coefficient (m4/s) for
!              active (NAT) and inert (NPT) tracer variables.  If variable
!              horizontal diffusion is activated, TNU4 is the mixing
!              coefficient for the largest grid-cell in the domain.
!
!------------------------------------------------------------------------------
! Harmonic/biharmonic horizontal viscosity coefficients.
!------------------------------------------------------------------------------
!
!  VISC2       Lateral, harmonic, constant, mixing coefficient (m2/s) for
!              momentum.  If variable horizontal viscosity is activated, UVNU2
!              is the mixing coefficient for the largest grid-cell in the
!              domain.
!
!  VISC4       Lateral, biharmonic, constant mixing coefficient (m4/s) for
!              momentum. If variable horizontal viscosity is activated, UVNU4
!              is the mixing coefficient for the largest grid-cell in the
!              domain.
!
!------------------------------------------------------------------------------
! Vertical mixing coefficients for active tracers.
!------------------------------------------------------------------------------
!
!  AKT_BAK     Background vertical mixing coefficient (m2/s) for active
!              (NAT) and inert (NPT) tracer variables.
!
!------------------------------------------------------------------------------
! Vertical mixing coefficient for momentum.
!------------------------------------------------------------------------------
!
!  AKV_BAK     Background vertical mixing coefficient (m2/s) for momentum.
!
!------------------------------------------------------------------------------
! Turbulent closure parameters.
!------------------------------------------------------------------------------
!
!  AKK_BAK     Background vertical mixing coefficient (m2/s) for turbulent
!              kinetic energy.
!
!  AKP_BAK     Background vertical mixing coefficient (m2/s) for turbulent
!              generic statistical field, "psi".
!
!  TKENU2      Lateral, harmonic, constant, mixing coefficient (m2/s) for
!              turbulent closure variables.
!
!  TKENU4      Lateral, biharmonic, constant mixing coefficient (m4/s) for
!              turbulent closure variables.
!
!------------------------------------------------------------------------------
! Generic length-scale turbulence closure parameters.
!------------------------------------------------------------------------------
!
!  GLS_P       Stability exponent (non-dimensional).
!
!  GLS_M       Turbulent kinetic energy exponent (non-dimensional).
!
!  GLS_N       Turbulent length scale exponent (non-dimensional).
!
!  GLS_Kmin    Minimum value of specific turbulent kinetic energy
!
!  GLS_Pmin    Minimum Value of dissipation.
!
! Closure independent constraint parameters (non-dimensional):
!
!  GLS_CMU0    Stability coefficient.
!
!  GLS_C1      Shear production coefficient.
!
!  GLS_C2      Dissipation coefficient.
!
!  GLS_C3M     Buoyancy production coefficient (minus).
!
!  GLS_C3P     Buoyancy production coefficient (plus).
!
!  GLS_SIGK    Constant Schmidt number (non-dimensional) for turbulent
!              kinetic energy diffusivity.
!
!  GLS_SIGP    Constant Schmidt number (non-dimensional) for turbulent
!              generic statistical field, "psi".
!
! Suggested values for various parameterizations:
!
!              MY2.5         K-epsilon    K-omega      K-omega      K-tao
!
!      GLS_P = 0.d0          3.0d0       -1.0d0       -1.0d0       -3.0d0
!      GLS_M = 1.d0          1.5d0        0.5d0        0.5d0        0.5d0
!      GLS_N = 1.d0         -1.0d0       -1.0d0       -1.0d0        1.0d0
!   GLS_Kmin = 5.0d-6        7.6d-6       7.6d-6       7.6d-6       7.6d-6
!   GLS_Pmin = 5.0d-6        1.0d-12      1.0d-12      1.0d-12      1.0d-12
!
!   GLS_CMU0 = 0.5544d0      0.5477d0     0.5477d0     0.5477d0     0.5477d0
!     GLS_C1 = 0.9d0         1.44d0       0.555d0      0.52d0       0.173d0
!     GLS_C2 = 0.5d0         1.92d0       0.833d0      0.8d0        0.225d0
!    GLS_C3M = 0.9d0        -0.4d0       -0.6d0       -0.6d0        0.0d0
!    GLS_C3P = 0.9d0         1.0d0        1.0d0        1.0d0        0.0d0
!   GLS_SIGK = 1.96d0        1.0d0        2.0d0        2.0d0        1.46d0
!   GLS_SIGP = 1.96d0        1.30d0       2.0d0        2.0d0       10.8d0
!
!------------------------------------------------------------------------------
! Constants used in the computation of momentum stress.
!------------------------------------------------------------------------------
!
!  RDRG        Linear bottom drag coefficient (m/s).
!
!  RDRG2       Quadratic bottom drag coefficient.
!
!  Zob         Bottom roughness (m).
!
!  Zos         Surface roughness (m).
!
!------------------------------------------------------------------------------
! Height of atmospheric measurements for bulk fluxes parameterization.
!------------------------------------------------------------------------------
!
!  BLK_ZQ      Height (m) of surface air humidity measurement.  Usually,
!                recorded at 10 m.
!
!  BLK_ZT      Height (m) of surface air temperature measurement.  Usually,
!                recorded at 2 or 10 m.
!
!  BLK_ZW      Height (m) of surface winds measurement. Usually, recorded
!                at 10 m.
!
!------------------------------------------------------------------------------
! Wetting and drying parameters.
!------------------------------------------------------------------------------
!
!  DCRIT       Minimum depth (m) for wetting and drying.
!
!------------------------------------------------------------------------------
! Jerlow Water type.
!------------------------------------------------------------------------------
!
!  WTYPE       Jerlov water type: an integer value from 1 to 5.
!
!------------------------------------------------------------------------------
! Body-force parameters. Used when CPP option BODYFORCE is activated.
!------------------------------------------------------------------------------
!
!  LEVSFRC     Deepest level to apply surface momentum stress as a body-force.
!
!  LEVBFRC     Shallowest level to apply bottom momentum stress as a body-force.
!
!------------------------------------------------------------------------------
! Vertical S-coordinates parameters.
!------------------------------------------------------------------------------
!
!  THETA_S     S-coordinate surface control parameter, [0 < theta_s < 20].
!
!  THETA_B     S-coordinate bottom  control parameter, [0 < theta_b < 1].
!
!  TCLINE      Width (m) of surface or bottom boundary layer in which
!              higher vertical resolution is required during stretching.
!
!              WARNING:  Users need to experiment with these parameters. We
!                        have found out that the model goes unstable with
!                        high values of THETA_S.  In steep and very tall
!                        topography, it is recommended to use THETA_S < 3.0.
!
!------------------------------------------------------------------------------
! Mean Density and background Brunt-Vaisala frequency.
!------------------------------------------------------------------------------
!
!  RHO0        Mean density (Kg/m3) used when the Boussinesq approximation
!              is inferred.
!
!  BVF_BAK     Background Brunt-Vaisala frequency squared (1/s2). Typical
!              values for the ocean range (as a function of depth) from
!              1.0E-4 to 1.0E-6.
!
!------------------------------------------------------------------------------
! Time Stamps.
!------------------------------------------------------------------------------
!
!  DSTART      Time stamp assigned to model initialization (days).  Usually
!              a Calendar linear coordinate, like modified Julian Day.  For
!              Example:
!
!                       Julian Day = 1  for  Nov 25, 0:0:0 4713 BCE
!              modified Julian Day = 1  for  May 24, 0:0:0 1968  CE GMT
!
!              It is called truncated or modified Julian day because an offset
!              of 2440000 needs to be added.
!
!  TIDE_START  Reference time origin for tidal forcing (days). This is the
!              time used when processing input tidal model data. It is needed
!              in routine "set_tides" to compute the correct phase lag with
!              respect ROMS/TOMS initialization time.
!
!  TIME_REF    Reference time (yyyymmdd.f) used to compute relative time:
!              elapsed time interval since reference-time.  The "units"
!              attribute takes the form "time-unit since reference-time".
!              This parameter also provides information about the calendar
!              used:
!
!              If TIME_REF = -2, model time and DSTART are in modified Julian
!              days units.  The "units" attribute is:
!
!                      'time-units since 1968-05-23 00:00:00 GMT'
!
!              If TIME_REF = -1, model time and DSTART are in a calendar
!              with 360 days in every year (30 days each month).  The "units"
!              attribute is:
!
!                      'time-units since 0000-01-01 00:00:00'
!
!              If TIME_REF = 0, model time and DSTART are in a common year
!              calendar with 365.25 days.  The "units" attribute is:
!
!                      'time-units since 0000-01-01 00:00:00'
!
!              If TIME_REF > 0, model time and DSTART are the elapsed time
!              units since specified reference time.  For example,
!              TIME_REF=20020115.5 will yield the following attribute:
!
!                      'time-units since 2002-01-15 12:00:00'
!
!------------------------------------------------------------------------------
! Nudging/relaxation time scales, inverse scales will be computed internally.
!------------------------------------------------------------------------------
!
! When passive/active open boundary conditions are activated, these nudging
! values correspond to the passive (outflow) nudging time scales.
!
!  TNUDG       Nudging time scale (days) for active tracer variables.
!              (1:NAT+NPT,1:Ngrids) values are expected.
!
!  ZNUDG       Nudging time scale (days) for free-surface.
!
!  M2NUDG      Nudging time scale (days) for 2D momentum.
!
!  M3NUDG      Nudging time scale (days) for 3D momentum.
!
!  OBCFAC      Factor between passive (outflow) and active (inflow) open
!              boundary conditions.  The nudging time scales for the
!              active (inflow) conditions are obtained by multiplying
!              the passive values by OBCFAC. If OBCFAC > 1, nudging on
!              inflow is stronger than on outflow (recommended).
!
!------------------------------------------------------------------------------
! Linear equation of State parameters.
!------------------------------------------------------------------------------
!
! Ignoring pressure, the linear equation of state is:
!
!              rho(:,:,:) = R0 - R0 * TCOEF * (t(:,:,:,:,itemp) - T0)
!                              + R0 * SCOEF * (t(:,:,:,:,isalt) - S0)
!
!              Typical values:     R0 = 1027.0  kg/m3               
!                                  T0 = 10.0    Celsius
!                                  S0 = 35.0    PSU
!                               TCOEF = 1.7d-4  1/Celsius
!                               SCOEF = 7.6d-4  1/PSU
!
!  R0          Background density value (Kg/m3) used in Linear Equation of
!              State.
!
!  T0          Background potential temperature (Celsius) constant.
!
!  S0          Background salinity (PSU) constant.
!
!  TCOEF       Thermal expansion coefficient in Linear Equation of State.
!
!  SCOEF       Saline contraction coefficient in Linear Equation of State.
!
!------------------------------------------------------------------------------
! Slipperiness parameter.
!------------------------------------------------------------------------------
!
!  GAMMA2      Slipperiness variable, either 1.0 (free slip) or -1.0 (no slip).
!
!------------------------------------------------------------------------------
!  Adjoint sensitivity parameters.
!------------------------------------------------------------------------------
!
!  DstrS       Starting day for adjoint sensitivity forcing.
!
!  DendS       Ending   day for adjoint sensitivity forcing.
!
!              The adjoint forcing is applied at every time step according to
!              desired state functional stored in the adjoint sensitivity
!              NetCDF file. DstrS must be less or equal to DendS. If both
!              values are zero, their values are reset internally to the full
!              range of the adjoint integration.
!
!  KstrS       Starting vertical level of the 3D adjoint state variables whose
!                sensitivity is required.
!  KendS       Ending   vertical level of the 3D adjoint state variables whose
!                sensitivity is required.
!
!  Lstate      Logical switches (TRUE/FALSE) to specify the adjoint state
!                variables whose sensitivity is required.
!
!                Lstate(isFsur):   Free-surface
!                Lstate(isUbar):   2D U-momentum
!                Lstate(isVbar):   2D V-momentum
!                Lstate(isUvel):   3D U-momentum
!                Lstate(isVvel):   3D V-momentum
!                Lstate(isTvar):   Traces (NT values expected)
!
!------------------------------------------------------------------------------
!  Stochastic optimals parameters.
!------------------------------------------------------------------------------
!
!  SO_decay    Stochastic optimals time decorrelation scale (days) assumed
!                for red noise processes.
!
!  SOstate     Logical switches (TRUE/FALSE) to specify the state surface
!                forcing variable whose stochastic optimals is required.
!
!                SOstate(isustr):  surface u-stress
!                SOstate(isvstr):  surface v-stress
!                SOstate(isTsur):  surface tracer flux (NT values expected)
!
!  SO_sdev     Stochastic optimals surface forcing standard deviation for
!                dimensionalization.
!
!                SO_sdev(isustr):  surface u-stress
!                SO_sdev(isvstr):  surface v-stress
!                SO_sdev(isTsur):  surface tracer flux (NT values expected)
!
!------------------------------------------------------------------------------
! Logical switches (T/F) to activate writing of fields into HISTORY file.
!------------------------------------------------------------------------------
!
!  Hout(idUvel)  Write out 3D U-velocity component.
!  Hout(idVvel)  Write out 3D V-velocity component.
!  Hout(idWvel)  Write out 3D W-velocity component.
!  Hout(idOvel)  Write out 3D omega vertical velocity.
!  Hout(idUbar)  Write out 2D U-velocity component.
!  Hout(idVbar)  Write out 2D V-velocity component.
!  Hout(idFsur)  Write out free-surface.
!
!  Hout(idTvar)  Write out active (NAT) tracers: temperature and salinity.
!
!  Hout(idUsms)  Write out surface U-momentum stress.
!  Hout(idVsms)  Write out surface V-momentum stress.
!  Hout(idUbms)  Write out bottom  U-momentum stress.
!  Hout(idVbms)  Write out bottom  V-momentum stress.
!
!  Hout(idUbrs)  Write out current-induced, U-momentum stress.
!  Hout(idVbrs)  Write out current-induced, V-momentum stress.
!  Hout(idUbws)  Write out wind-induced, bottom U-wave stress.
!  Hout(idVbws)  Write out wind-induced, bottom V-wave stress.
!  Hout(idUbcs)  Write out bottom maximum wave and current U-stress.
!  Hout(idVbcs)  Write out bottom maximum wave and current V-stress.
!
!  Hout(idUbot)  Write out wind-induced, bed wave orbital U-velocity.
!  Hout(idVbot)  Write out wind-induced, bed wave orbital V-velocity.
!  Hout(idUbur)  Write out bottom U-velocity above bed.
!  Hout(idVbvr)  Write out bottom V-velocity above bed.
!
!  Hout(idTsur)  Write out surface net heat and salt flux
!  Hout(idLhea)  Write out latent heat flux.
!  Hout(idShea)  Write out sensible heat flux.
!  Hout(idLrad)  Write out long-wave radiation flux.
!  Hout(idSrad)  Write out short-wave radiation flux.
!  Hout(idevap)  Write out evaporation rate.
!  Hout(idrain)  Write out precipitation rate.
!
!  Hout(idDano)  Write out density anomaly.
!  Hout(idVvis)  Write out vertical viscosity coefficient.
!  Hout(idTdif)  Write out vertical diffusion coefficient of temperature.
!  Hout(idSdif)  Write out vertical diffusion coefficient of salinity.
!  Hout(idHsbl)  Write out depth of oceanic surface boundary layer.
!  Hout(idHbbl)  Write out depth of oceanic bottom boundary layer.
!  Hout(idMtke)  Write out turbulent kinetic energy.
!  Hout(idMtls)  Write out turbulent kinetic energy times length scale.
!
!  Hout(inert)   Write out extra inert passive tracers.
!
!  Hout(idBott)  Write out exposed sediment layer properties, 1:MBOTP.
!
!------------------------------------------------------------------------------
! Generic User parameters.
!------------------------------------------------------------------------------
!
!  NUSER       Number of User parameters to consider (integer).
!  USER        Vector containing user parameters (real array). This array
!                is used with the SANITY_CHECK to test the correctness of
!                the tangent linear adjoint models.  It contains information
!                of the model variable and grid point to perturb:
!
!                INT(user(1)):  tangent state variable to perturb
!                INT(user(2)):  adjoint state variable to perturb
!                               [isFsur=1] free-surface 
!                               [isUbar=2] 2D U-momentum
!                               [isVbar=3] 2D V-momentum
!                               [isUvel=4] 3D U-momentum
!                               [isVvel=5] 3D V-momentum
!                               [isTvar=6] Firt tracer (temperature)
!                               [   ...  ]
!                               [isTvar=?] Last tracer
!
!                INT(user(3)):  I-index of tangent variable to perturb
!                INT(user(4)):  I-index of adjoint variable to perturb
!                INT(user(5)):  J-index of tangent variable to perturb
!                INT(user(6)):  J-index of adjoint variable to perturb
!                INT(user(7)):  K-index of tangent variable to perturb, if 3D
!                INT(user(8)):  K-index of adjoint variable to perturb, if 3D
!
!                Set tangent and adjoint parameters to the same values
!                if perturbing and reporting the same variable.
!
!------------------------------------------------------------------------------
! Input/output NetCDF file names (string with a maximum of eighty characters).
!------------------------------------------------------------------------------
!
!  GRDNAME     Input grid file name.
!  ININAME     Input nonlinear initial conditions file name. It can be a
!                re-start file.
!  IRPNAME     Input representer model initial conditions file name.
!  ITLNAME     Input tangent linear model initial conditions file name.
!  IADNAME     Input adjoint model initial conditions file name.
!  FRCNAME     Input forcing fields file name.
!  CLMNAME     Input climatology fields file name.
!  BRYNAME     Input open boundary data file name.
!  FWDNAME     Input forward solution fields file name.
!  ADSNAME     Input adjoint sensitivity functional file name.
!
!  GSTNAME     Output GST analysis re-start file name.
!  RSTNAME     Output re-start file name.
!  HISNAME     Output history file name.
!  TLFNAME     Output impulse forcing for tangent linear (TLM and RPM) models.
!  TLMNAME     Output tangent linear file name.
!  ADJNAME     Output adjoint file name.
!  AVGNAME     Output averages file name.
!  DIANAME     Output diagnostics file name.
!  STANAME     Output stations file name.
!  FLTNAME     Output floats file name.
!
!------------------------------------------------------------------------------
! Input ASCII parameters file names.
!------------------------------------------------------------------------------
!
!  APARNAM     Input assimilation parameters file name.
!  SPOSNAM     Input stations positions file name.
!  FPOSNAM     Input initial drifters positions file name.
!  BPARNAM     Input biological parameters file name.
!  SPARNAM     Input sediment transport parameters file name.
!  USRNAME     USER's input generic file name.
!
EOF