/*
** svn $Id: lmd_test.h 895 2009-01-12 21:06:20Z kate $
*******************************************************************************
** Copyright (c) 2002-2009 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for K-Profile Parameterization Test.
**
** Application flag:   LMD_TEST
** Input script:       ocean_lmd_test.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define NONLIN_EOS
#define SALINITY
#define AVERAGES
#define AVERAGES_AKS
#define AVERAGES_AKT
#define AVERAGES_AKV
#define STATIONS
#define SOLVE3D
#define SPLINES
#define EW_PERIODIC
#define NS_PERIODIC
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_DDMIX
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
