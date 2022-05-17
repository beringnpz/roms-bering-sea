/*
** svn $Id: basin.h 8 2007-02-06 19:00:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for Northeast Pacific (NEP5) simulation
**  By default this will run a physics-only variant.  To add the BESTNPZ 
**  biogeochmistry, add "-DBESTNPZ" to MY_CPP_FLAGS; to add the full 
**  FEAST food web, add "-DBESTNPZ -DFEAST"
*/


#define  BERING_10K

#undef NETCDF4
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

/* general */

#define CURVGRID
#define MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#ifdef SOLVE3D
# undef SPLINES
#endif
#undef FLOATS
#define STATIONS
#undef WET_DRY

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_PASSIVE
#endif
 
/* salinity nudging */

#define SCORRECTION

/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL
# ifdef ICE_MODEL
#  define  ICE_THERMO
#  define  ICE_MK
#  undef   ICE_ALB_EC92
#  undef   ICE_SMOOTH
#  define  ICE_MOMENTUM
#  define  ICE_MOM_BULK
#  define  ICE_EVP
#  define  ICE_ADVECT
#  define  ICE_SMOLAR
#  define  ICE_UPWIND
#  define  ICE_BULK_FLUXES
#  define  ANA_AIOBC
#  define  ANA_HIOBC
#  define  ANA_HSNOBC
# endif
#endif

/* output stuff */
 
#define NO_WRITE_GRID
#undef OUT_DOUBLE
#define RST_SINGLE
#define AVERAGES
#undef AVERAGES2
#ifdef SOLVE3D
# undef AVERAGES_DETIDE
# define AVERAGES_AKT
# define AVERAGES_AKS
# define AVERAGES_AKV
# define AVERAGES_FLUXES
# undef AVERAGES_QUADRATIC
# undef DIAGNOSTICS_TS
#endif
#undef DIAGNOSTICS_UV
 
/* advection, dissipation, pressure grad, etc. */
 
#ifdef SOLVE3D
# define DJ_GRADPS
#endif
 
#define UV_ADV
#define UV_COR
#define UV_SADVECTION
 
#ifdef SOLVE3D
# define TS_C4HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif
 
#define UV_VIS2
#define UV_SMAGORINSKY
#define VISC_3DCOEF
#define MIX_S_UV
#define VISC_GRID
#define SPONGE

#ifdef SOLVE3D
# define TS_DIF2
# define MIX_GEO_TS
# define DIFF_GRID
#endif
 
 
/* vertical mixing */
 
#ifdef SOLVE3D
# define SOLAR_SOURCE
 
# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  undef LMD_BKPP
#  define LMD_NONLOCAL
#  define LMD_SHAPIRO
#  undef LMD_DDMIX
# endif
 
# undef GLS_MIXING
# undef MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
#endif
 
/* surface forcing */
 
#ifdef SOLVE3D
# define CORE_FORCING
# define BULK_FLUXES
# define CCSM_FLUXES
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define  DIURNAL_SRFLUX 
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO
#  define ALBEDO_CURVE  
#  if defined NETLONGWAVE /* A B10K-specific option, to keep LONGWAVE_OUT the default */
#   undef LONGWAVE_OUT
#   undef LONGWAVE
#  else
#   define LONGWAVE_OUT
#   undef LONGWAVE
#  endif
# endif
#endif
 
/* surface and side corrections */
 
#ifdef SOLVE3D
# undef SRELAXATION
# undef QCORRECTION
#endif
 
#ifdef SOLVE3D
# undef TCLIMATOLOGY
# undef TCLM_NUDGING
#endif
 
/* point sources (rivers, line sources) */
/* Using Runoff instead now             */

#ifdef SOLVE3D
# define RUNOFF
# define UV_PSOURCE
# define ANA_PSOURCE
# undef TS_PSOURCE
#endif
 
/* tides */
 
#define LTIDES
#ifdef LTIDES
# undef FILTERED   /*need on eventually*/
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES

# define UV_LDRAG
# define RDRG_GRID
# define DRAG_LIMITER
# undef UV_QDRAG
#else
# define UV_QDRAG
#endif
 
/* Boundary conditions...careful with grid orientation */
 
#define EASTERN_WALL
#define NORTHERN_WALL
#undef WESTERN_WALL
#undef SOUTHERN_WALL
 
#define RADIATION_2D
 
#ifndef NORTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# ifdef SOLVE3D
#  define NORTH_M3RADIATION
#  define NORTH_M3NUDGING
#  define NORTH_TRADIATION
#  define NORTH_TNUDGING
#  define NORTH_MIGRADIENT
# endif
#endif
 
#ifndef WESTERN_WALL
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# ifdef SOLVE3D
#  define WEST_M3RADIATION
#  define WEST_M3NUDGING
#  define WEST_TRADIATION
#  define WEST_TNUDGING
#  define WEST_MIGRADIENT
# endif
#endif
 
#ifndef SOUTHERN_WALL
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# ifdef SOLVE3D
#  define SOUTH_M3RADIATION
#  define SOUTH_M3NUDGING
#  define SOUTH_TRADIATION
#  define SOUTH_TNUDGING
#  define SOUTH_MIGRADIENT
# endif
#endif
 
#ifndef EASTERN_WALL
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# ifdef SOLVE3D
#  define EAST_M3RADIATION
#  define EAST_M3NUDGING
#  define EAST_TRADIATION
#  define EAST_TNUDGING
#  define EAST_MIGRADIENT
# endif
#endif
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif

/*
**  Biological model options.
*/

/* In order to allow use of this header file with different biology, I'm 
   leaving out any explicit define/undef options for BESTNPZ and FEAST.  
   These will be defined externally via MY_CPP_FLAGS.
*/ 

#undef NEMURO
#undef BIO_GOANPZ        /* Sarah Hinckley's 11 box model */

#if defined BEST_NPZ || defined BIO_GOANPZ || defined PASSIVE_TRACERS
# undef  BIOFLUX           /* sum Nitrogen fluxes between boxes */
# define ANA_BIOLOGY       /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define DIAPAUSE          /* Enable Neocalanus seasonal vertical migration */
# undef FLOAT_VWALK
#endif

#if defined NEMURO
# undef ANA_BIOLOGY        /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define IRON_LIMIT        /* Add iron as passive 11th tracer */
# define IRON_RELAX
# undef  IRON_RSIN
# define BIO_SEDIMENT
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
# undef  ANA_BIOSWRAD
# undef  DIAGNOSTICS_BIO
# undef  BIO_SEDIMENT
#endif

#ifdef BEST_NPZ

# define BIO_NUDGE_ONLY

# undef LIMIT_BIO_AKT
# if !defined(NEWSHADE) && !defined(NEWSHADESHALLOW) && !defined(COKELET)
#   define COKELET         /* updated light attenuation, based on  */
# endif
# undef KODIAK_IRAD        /* Generate irradiance with curve matching Kodiak data 
                              Else use Sarah Hinckly originl code   */
# define JELLY
# define IRON_LIMIT        /* Add iron  */
# define BENTHIC           /*FENNEL or BENTHIC or TRAP*/
# define ICE_BIO
# undef CLIM_ICE_1D
# define TCLM_NUDGING      /* Nudging of tracer climatology for iron */
# define ANA_TCLIMA        /* analytical tracers climatology for iron */
# define TCLIMATOLOGY      /* Processing of tracer climatology for iron */
# if defined CARBON
#  define CARBON_FLUX      /* For river fluxes of DIC,TA */
#  define OXYGEN           /* For oxygen cycling */
#  define STATIONARY2
# endif
# define STATIONARY
# undef PROD3 
# undef PROD2
# undef SINKVAR            /* for variable sinking rate*/
# undef DENMAN
# undef CORRECT_TEMP_BIAS  /* corrects ROMS temp for biology only */
#endif

/*
#undef  OFFLINE_BIOLOGY     define if offline simulation of bio tracers 
#if defined OFFLINE_BIOLOGY
# define AKSCLIMATOLOGY     Processing of AKS climatology 
# undef ANA_AKSCLIMA        Processing of AKS climatology 
#endif
*/

/* Here's the new stuff for feast*/

#ifdef FEAST
# ifdef PROD3
#  define FEAST_DAT        /* only use this if using 60-layer model */
# endif
# undef FEAST_DAT          /* so undef in the 10-layer model */
# define FEAST_FORCING
# define T_PASSIVE
# define ANA_PASSIVE
# undef ANA_PASSIVE
# undef FLOATS
# define ANA_BPFLUX
# define ANA_SPFLUX
# undef FEAST_NOEXCHANGE
#endif



