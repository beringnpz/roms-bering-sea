/*
** svn $Id: test_head.h 895 2009-01-12 21:06:20Z kate $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Sediment Test Headland Case.
**
** Application flag:   TEST_HEAD
** Input scripts:      ocean_test_head.in
**                     coupling_test_head.in
**                     sediment_test_head.in
*/

#undef  CURVGRID
#define WET_DRY
#define NEARSHORE_MELLOR
#define SWAN_COUPLING
#define MASKING
#define OUT_DOUBLE
#define UV_ADV
#undef  UV_COR
#define TS_MPDATA
#define DJ_GRADPS
#undef  SALINITY
#define SOLVE3D
#define SPLINES
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EAST_M3RADIATION
#define WEST_M3RADIATION
#define EAST_TGRADIENT
#define WEST_TGRADIENT

#define SSH_TIDES
#define WEST_FSCLAMPED
#define EAST_FSCLAMPED

#define ANA_FSOBC
#define ANA_M2OBC

#define UV_TIDES
#define FSOBC_REDUCED
#define WEST_M2REDUCED
#define EAST_M2REDUCED

#define UV_LOGDRAG
#undef  UV_QDRAG
#undef  MB_BBL
#undef  SG_BBL
#undef  SSW_BBL

#ifdef SWAN_COUPLING
# define MCT_LIB
#endif

#ifdef SG_BBL
# define SG_CALC_ZNOT
# undef  SG_LOGINT
#endif
#ifdef MB_BBL
# define MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_LOGINT
#endif
#if defined MB_BBL || defined SG_BBL || defined SSW_BBL || defined NEARSHORE_MELLOR
# undef ANA_WWAVE
# ifndef ANA_WWAVE
# endif
#endif

#ifdef SOLVE3D
# define SEDIMENT
# ifdef SEDIMENT
#  define SED_MORPH
#  undef  SED_DENS
#  define SUSPLOAD
#  define BEDLOAD_MPM
#  undef  BEDLOAD_SOULSBY
#  define AVERAGES_BEDLOAD
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  undef ANA_SEDIMENT
# endif
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX

# undef  ANA_VMIX
# define MY25_MIXING
# undef  GLS_MIXING
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  undef  CRAIG_BANNER
#  undef  CHARNOK
#  undef  ZOS_HSIG
#  undef  TKE_WAVEDISS
# endif
#endif

#undef  ANA_GRID
#undef  ANA_INITIAL
#define ANA_SMFLUX

