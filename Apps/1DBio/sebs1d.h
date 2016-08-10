
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
*/

#define BIO_1D 

#undef NETCDF4
#undef PARALLEL_IO
#undef OFFLINE_FLOATS

/* general */

#undef CURVGRID
#undef MASKING
#define NONLIN_EOS
#define SOLVE3D
#define SALINITY
#ifdef SOLVE3D
# undef SPLINES
#endif
#undef FLOATS
#undef STATIONS
#undef WET_DRY

#undef T_PASSIVE
#ifdef T_PASSIVE
# define ANA_PASSIVE
#endif
 
/* ice */

#ifdef SOLVE3D
# undef ICE_MODEL
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
#define AVERAGES2
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
#define UV_COR   /*!!!!!!!!!!!!!!!!!!!*/
#define UV_SADVECTION
 
#ifdef SOLVE3D
# define TS_C4HADVECTION
# define TS_C4VADVECTION
# undef TS_MPDATA
#endif
 
#define UV_VIS2  /*!!!!!!!!!!!!!!!!!!!*/
#define UV_SMAGORINSKY /*!!!!!!!!!!!!!!!!!!!*/
#define VISC_3DCOEF/*!!!!!!!!!!!!!!!!!!!*/
#define MIX_S_UV
#define VISC_GRID
#define SPONGE

#ifdef SOLVE3D
# define TS_DIF2
# undef MIX_GEO_TS /*define*/
# define MIX_ISO_TS  /*undef*/
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
#  define LMD_SHAPIRO /*!!!!!!!!!!!!!!!!!!!*/
#  undef LMD_DDMIX /*!!!!!!!!!!!!!!!!!!!*/
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
#  define LONGWAVE_OUT
#  define DIURNAL_SRFLUX  /*was undefined */
#  define EMINUSP
#  undef ANA_SRFLUX
#  undef ALBEDO
#  define ALBEDO_CURVE
#  undef LONGWAVE
# endif
#endif
 
 
/*
**  Select analytical fields configuration: define if using any of the
**  following options.  Set the appropriate analytical expression in
**  file "analytical.F".
*/


#undef  ANA_BPFLUX      /* analytical bottom passive tracers fluxes */
#undef  ANA_BSEDIM      /* analytical bottom sediment grain size & density */
#define ANA_BSFLUX      /* analytical bottom salinity flux */
#define ANA_BTFLUX      /* analytical bottom temperature flux */
#undef  ANA_CLOUD       /* analytical cloud fraction */
#undef  ANA_DIAG        /* Customized diagnostics */
#undef  ANA_FSOBC       /* analytical free-surface boundary conditions */
#undef  ANA_GRID        /* analytical model grid set-up */
#undef  ANA_HUMIDITY    /* analytical surface air humidity */
#undef  ANA_INITIAL     /* analytical initial conditions */
#define  ANA_M2CLIMA     /* analytical 2D momentum climatology */
#undef  ANA_M2OBC       /* analytical 2D momentum boundary conditions */
#undef  ANA_M3CLIMA     /* analytical 3D momentum climatology */
#undef  ANA_M3OBC       /* analytical 3D momentum boundary conditions */
#undef  ANA_PAIR        /* analytical surface air pressure */
#undef  ANA_PSOURCE     /* analytical point Sources/Sinks */
#undef  ANA_RAIN        /* analytical rain fall rate */
#undef  ANA_SEDIMENT    /* analytical sediment initial conditions */
#undef  ANA_SMFLUX      /* analytical surface momentum stress */
#undef  ANA_SPFLUX      /* analytical surface passive tracers fluxes */
#undef  ANA_SPINNING    /* analytical time-varying rotation force */
#undef ANA_SRFLUX      /* analytical surface shortwave radiation flux */
#undef  ANA_SSFLUX      /* analytical surface salinity flux */
#undef  ANA_SSH         /* analytical sea surface height */
#undef  ANA_SSS         /* analytical sea surface salinity */
#undef  ANA_SST         /* analytical SST and dQdSST */
#undef  ANA_STFLUX      /* analytical surface temperature flux */
#undef  ANA_TAIR        /* analytical surface air temperature */
#undef  ANA_TCLIMA      /* analytical tracers climatology */
#undef  ANA_TOBC        /* analytical tracers boundary conditions */
#undef  ANA_VMIX        /* analytical vertical mixing coefficients */
#undef  ANA_WINDS       /* analytical surface winds */
#undef  ANA_WWAVE       /* analytical wind induced waves */


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
 
/* Using Runoff instead now */
#ifdef SOLVE3D
# undef RUNOFF
# undef UV_PSOURCE
# undef ANA_PSOURCE
# undef TS_PSOURCE
#endif
 
/* tides */
 
#undef LTIDES   /* define */
#ifdef LTIDES
# define FILTERED
# define SSH_TIDES
# define UV_TIDES
# define ADD_FSOBC
# define ADD_M2OBC
# undef RAMP_TIDES
# define TIDES_ASTRO
# define POT_TIDES
# undef UV_LDRAG
# undef RDRG_GRID
# define DRAG_LIMITER
# define UV_QDRAG
#else
# define UV_QDRAG
#endif
 
/* Boundary conditions...careful with grid orientation */
 

 
#define RADIATION_2D  /* !!!!!!! */

#define EW_PERIODIC       /* East-West periodic boundaries */
#define NS_PERIODIC       /* North-South periodic boundaries */ 
 
#undef  EASTERN_WALL      /* Eastern edge, closed wall condition */
#undef  WESTERN_WALL      /* Western edge, closed wall condition */
#undef  NORTHERN_WALL     /* Northern edge, closed wall condition */
#undef  SOUTHERN_WALL     /* Southern edge, closed wall condition */


 
 
/*                                                                 
**  Turn ON or OFF options for reading and processing of climatological                                                               
**  fields.  The nudging of climatology data is primarily used in sponge
**  areas.
*/ 
#define M2CLIMATOLOGY   /* Processing of 2D momentum climatology */
#undef  M3CLIMATOLOGY   /* Processing of 3D momentum climatology   define */
#undef  OCLIMATOLOGY    /* Processing of omega climatology    define*/      
#define TCLIMATOLOGY    /* Processing of tracer climatology   define*/
#undef  ZCLIMATOLOGY    /* Processing of SSH climatology      define*/        
#define  M2CLM_NUDGING   /* Nudging of 2D momentum climatology */
#undef  M3CLM_NUDGING   /* Nudging of 3D momentum climatology  define   */
#define TCLM_NUDGING    /* Nudging of tracer climatology       define */
#undef  ZCLM_NUDGING    /* Nudging of SSH climatology          define   */
 
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

#undef  NPZD1              /* Craig Lewis's 4 box model */
#undef  BIO_GOANPZ         /* Sarah Hinckley's 10 box model */
#define BEST_NPZ         /* Georgina Gibsons BEST NPZ model  */
#undef  PASSIVE_TRACERS    /* add 5 tracer boxes that are passive */

#if defined BEST_NPZ || defined BIO_GOANPZ || defined PASSIVE_TRACERS
# undef  BIOFLUX           /* sum Nitrogen fluxes between boxes */
# define  ANA_BIOLOGY       /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define DIAPAUSE          /* Enable Neocalanus seasonal vertical migration */
#endif



#ifdef BEST_NPZ
#  undef LIMIT_BIO_AKT
#  undef NEWSHADE   /*formulation for self shading in PAR calc basrd on Morel*/
#  define COKELET /*formulation for PAR based on Ned Cokelet obs in Bering*/
#  undef KODIAK_IRAD /* Generate irradiance with curve matching Kodiak data
                       Else use Sarah Hinckly originl code   */
#  define JELLY
#  undef IRON_LIMIT        /* Add iron  */
#  define BENTHIC /*FENNEL or BENTHIC or TRAP*/
#  undef ICE_BIO
#  undef CLIM_ICE_1D
#  define  TCLM_NUDGING    /* Nudging of tracer climatology does all inc pys */
#  undef ANA_TCLIMA     /* analytical tracers climatology  */
#  define TCLIMATOLOGY   /* Processing of tracer climatology  */
#  define FULL_NUDGE
#  undef  BIO_ONLY_NUDGE /*switch for nunge of bio tracers to climatology NO3, NH4, Fe*/
#  undef INI_NO3_WOA
#  undef STATIONARY
#  undef STATIONARY2
#  undef PROD3
#  undef PROD2
#  undef SINKVAR      /* for variable sinking rate*/
#  undef DENMAN
#  undef CORRECT_TEMP_BIAS /* corrects ROMS temp for biology only */
#  undef EUP_VM  
# endif

# undef  OFFLINE_BIOLOGY   /* define if offline simulation of bio tracers */
#   if defined OFFLINE_BIOLOGY
#    define AKSCLIMATOLOGY   /* Processing of AKS climatology */
#    undef ANA_AKSCLIMA      /* Processing of AKS climatology */
#   endif


