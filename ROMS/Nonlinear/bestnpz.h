#include "cppdefs.h"
      SUBROUTINE biology (ng,tile)

!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Georgina Gibsons BESTNPZ Code July 2010
!
!  Modified from Sarah Hinckleys GOANPZ code                           !
!  Implemented by Craig Lewis (CVL)                                    !
!  Modified by Liz Dobbins and Sarah Hinckley                          !
!                                                                      !
!=======================================================================

      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_mixing

#if defined BERING_10K
#if defined ICE_BIO
      USE mod_ice
#endif
#if defined FEAST
      USE mod_feast
#endif
#endif
# if defined CLIM_ICE_1D
      USE mod_clima
#endif

      integer, intent(in) :: ng, tile

#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN

#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   N(ng), NT(ng),                                 &
     &                   nnew(ng), nstp(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % t                                  &
#if defined BENTHIC
     &                   ,OCEAN(ng) % bt                                &
#endif
#if defined FEAST
     &                   ,OCEAN(ng) % u                                 &
     &                   ,OCEAN(ng) % v                                 &
     &                   ,GFEAST(ng)                                    &
#endif
#if defined ICE_BIO
# ifdef CLIM_ICE_1D
     &                   ,OCEAN(ng) % it                                &
     &                   ,OCEAN(ng) % itL                               &
     &                   ,CLIMA(ng) % tclm                              &
# elif defined BERING_10K
     &                   ,ICE(ng) % IcePhL                              &
     &                   ,ICE(ng) % IceNO3                              &
     &                   ,ICE(ng) % IceNH4                              &
     &                   ,ICE(ng) % IceLog                              &
     &                   ,ICE(ng) % ti                                  &
     &                   ,ICE(ng) % hi                                  &
     &                   ,ICE(ng) % ai                                  &
     &                   ,ICE(ng) % ageice                              &
     &                   ,ICE(ng) % ui                                  &
     &                   ,ICE(ng) % vi                                  &
# endif
#endif
#ifdef STATIONARY
     &                   ,OCEAN(ng) % st                                &
     &                   ,NTS(ng)                                       &
#endif
#ifdef STATIONARY2
     &                   ,OCEAN(ng) % st2                               &
     &                   ,NTS2(ng)                                      &
#endif
#ifdef PROD3
     &                   ,OCEAN(ng) % pt3                               &
     &                   ,NPT3(ng)                                      &
#endif
#ifdef PROD2
     &                   ,OCEAN(ng) % pt2                               &
     &                   ,NPT2(ng)                                      &
#endif
#ifdef BIOFLUX
     &                  ,OCEAN(ng) % bflx                               &
#endif
     &                  ,MIXING(ng) % Akt                               &
     &                              )

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology

!************************************************************************
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         UBk, UBt,                                &
     &                         nnew, nstp,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w, srflx, t                   &
#if defined BENTHIC
     &                            ,bt                                   &
#endif
#if defined FEAST
     &                           ,u,v                                   &
     &                           ,GF                                &
#endif
#if defined ICE_BIO
# ifdef CLIM_ICE_1D
     &                            ,it                                   &
     &                            ,itL                                  &
     &                            ,tclm                                 &
# elif defined BERING_10K
     &                            , IcePhL                              &
     &                            , IceNO3                              &
     &                            , IceNH4                              &
     &                            , IceLog                              &
     &                            ,ti                                   &
     &                            ,hi                                   &
     &                            ,ai                                   &
     &                            ,ageice                               &
     &                            ,ui,vi                                &
# endif
#endif
#ifdef STATIONARY
     &                          ,st                                     &
     &                          ,UBst                                   &
#endif
#ifdef STATIONARY2
     &                          ,st2                                    &
     &                          ,UBst2                                  &
#endif
#ifdef PROD3
     &                          ,pt3                                    &
     &                          ,UBpt3                                  &
#endif
#ifdef PROD2
     &                          ,pt2                                    &
     &                          ,UBpt2                                  &
#endif
#ifdef BIOFLUX
     &                          ,bflx                                   &
#endif
     &                          ,Akt                                    &

     &                                  )

     !==================================================================
     !  VARIABLE DECLARATIONS
     !==================================================================

      USE mod_param
      USE mod_biology
      USE mod_scalars
      USE mod_ocean
      USE mod_grid
#if defined FEAST
      USE mod_feast
#endif
#if defined BERING_10K

# if defined ICE_BIO
#  ifdef BERING_10K
      USE mod_ice
      USE IcePhLbc_mod, ONLY : IcePhLbc_tile
      USE IceNO3bc_mod, ONLY : IceNO3bc_tile
      USE IceNH4bc_mod, ONLY : IceNH4bc_tile
#  endif
# endif

# ifdef DISTRIBUTE
      USE mp_exchange_mod

# endif

#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod
#endif

#if defined CLIM_ICE_1D
      USE mod_clima
#endif

      implicit none

      !  Imported variable declarations.

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: UBk, UBt
      integer, intent(in) :: nnew, nstp
!
#ifdef STATIONARY
      integer, intent(in) :: UBst
#endif
#if defined STATIONARY2
      integer, intent(in) :: UBst2
#endif
#if defined PROD3
      integer, intent(in) :: UBpt3
#endif
#if defined PROD2
      integer, intent(in) :: UBpt2
#endif

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)

# ifdef STATIONARY
      real(r8), intent(inout) :: st(LBi:,LBj:,:,:,:)
# endif
# ifdef STATIONARY2
      real(r8), intent(inout) :: st2(LBi:,LBj:,:,:)
# endif
# ifdef PROD3
      real(r8), intent(inout) :: pt3(LBi:,LBj:,:,:,:)
# endif
# ifdef PROD2
      real(r8), intent(inout) :: pt2(LBi:,LBj:,:,:)
# endif
# if defined BENTHIC
      real(r8), intent(inout) :: bt(LBi:,LBj:,:,:,:)
# endif
# if defined FEAST
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,UBk,3),v(LBi:UBi,LBj:UBj,UBk,3)
      TYPE (T_FEAST) :: GF
      !real(r8), intent(inout) :: ghfal(LBi:UBi,LBj:UBj,nrates,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
      !real(r8), intent(inout) :: ghfl(LBi:UBi,LBj:UBj,nrates,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS)
      !real(r8), intent(inout) :: ghfsp(LBi:UBi,LBj:UBj,nrates,NUM_SIMPLE_SPECIES)
      !real(r8), intent(inout) :: ozm(LBi:UBi,LBj:UBj,UBk,NUM_PLANKTON)
      !real(r8), intent(inout) :: ofdat(LBi:UBi,LBj:UBj,UBk,NFDAT)
      !real(r8), intent(inout) :: incatch (LBi:UBi,LBj:UBj,NUM_GEARS,TOT_FEAST)
# endif
# if defined ICE_BIO
#  ifdef CLIM_ICE_1D
      real(r8), intent(inout) :: it(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: itL(LBi:,LBj:,:,:)

#  elif defined BERING_10K
      real(r8), intent(in) :: ti(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: ageice(LBi:,LBj:,:)
      real(r8), intent(inout) :: IcePhL(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNO3(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNH4(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceLog(LBi:,LBj:,:)
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
#  endif

#  ifdef BIOFLUX
      real(r8), intent(inout) :: bflx(:,:)
#  endif
# endif

#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# ifdef STATIONARY
      real(r8), intent(inout) :: st(LBi:UBi,LBj:UBj,UBk,3,NTS(ng))
# endif
# ifdef STATIONARY2
      real(r8), intent(inout) :: st2(LBi:UBi,LBj:UBj,3,NTS2(ng))
# endif
# ifdef PROD3
      real(r8), intent(inout) :: pt3(LBi:UBi,LBj:UBj,UBk,3,NPT3(ng))
# endif
# ifdef PROD2
      real(r8), intent(inout) :: pt2(LBi:UBi,LBj:UBj,3,NPT2(ng))
# endif
# if defined BENTHIC
      real(r8), intent(inout) :: bt(LBi:UBi,LBj:UBj,UBk,3,1)
# endif
# if defined FEAST
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,UBk,3),v(LBi:UBi,LBj:UBj,UBk,3)
# endif
# if defined ICE_BIO
#  ifdef CLIM_ICE_1D
      real(r8), intent(inout) :: it(LBi:UBi,LBj:UBj,3,1)
      real(r8), intent(inout) :: itL(LBi:UBi,LBj:UBj,3,1)

      real(r8), intent(inout) ::tclmG(LBi:UBi,LBj:UBj,UBk,3,NH(ng)+2)
      real(r8), intent(inout) ::tclm(LBi:UBi,LBj:UBj,UBk,NT(ng)+2)
#  elif defined BERING_10K
      real(r8), intent(in) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ageice(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IcePhL(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNO3(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNH4(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceLog(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
#  endif
# endif
# ifdef BIOFLUX
      real(r8), intent(inout) :: bflx(LBi:UBi,LBj:UBj)
# endif
#endif
#if defined FEAST

      real(r8) :: ROMS_depth(LBi:UBi,LBj:UBj,UBk)
      real(r8) :: ROMS_edges(LBi:UBi,LBj:UBj,UBk+1)
      real(r8) :: ROMS_temp(LBi:UBi,LBj:UBj,UBk)
      real(r8) :: ROMS_zoop(LBi:UBi,LBj:UBj,UBk,NUM_PLANKTON)
      !real :: ZoopFishDeath(LBi:UBi,LBj:UBj,UBk,NUM_PLANKTON)
      !real :: out_zoop_mort(LBi:UBi,LBj:UBj,UBk,NUM_PLANKTON)
      !real :: tvar(LBi:UBi,LBj:UBj,NUM_TVAR)

      real(r8) :: fal(LBi:UBi,LBj:UBj,nfvaral,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
      real(r8) :: fl(LBi:UBi,LBj:UBj,nfvarl,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS )
      real(r8) :: fsp(LBi:UBi,LBj:UBj,nfvar,NUM_SIMPLE_SPECIES)
      real(r8) :: ratepar(nrates)
      real(r8) :: t_new(LBi:Ubi,LBj:UBj)
      !real(r8) :: t_old(LBi:Ubi,LBj:UBj)
      real(r8) :: base_speed, diff_coeff
      real(r8) :: delN, Ftot, Fmax, BB, CC, delCF, delCAL
!     real(r8) :: predSumCop, predSumNCaS, predSumEupS, predSumNCaO, predSumEupO

      real(r8) :: diff_xiup(LBi:Ubi,LBj:UBj)
      real(r8) :: diff_xidn(LBi:Ubi,LBj:UBj)
      real(r8) :: diff_etup(LBi:Ubi,LBj:UBj)
      real(r8) :: diff_etdn(LBi:Ubi,LBj:UBj)
      real(r8) :: u_move(LBi:Ubi,LBj:UBj)
      real(r8) :: v_move(LBi:Ubi,LBj:UBj)
      real(r8) :: u_up(LBi:Ubi,LBj:UBj)
      real(r8) :: u_dn(LBi:Ubi,LBj:UBj)
      real(r8) :: v_up(LBi:Ubi,LBj:UBj)
      real(r8) :: v_dn(LBi:Ubi,LBj:UBj)

      real(r8), dimension(LBi:Ubi,LBj:UBj) :: t_left, CFmat, CAmat, CA_left, CF_left, CFnew, CAnew
      real(r8), dimension(LBi:Ubi,LBj:UBj) :: gT_u_up, gT_u_dn, gT_v_up, gT_v_dn
      real(r8), dimension(LBi:Ubi,LBj:UBj) :: gCF_u_up, gCF_u_dn, gCF_v_up, gCF_v_dn
      real(r8), dimension(LBi:Ubi,LBj:UBj) :: gCA_u_up, gCA_u_dn, gCA_v_up, gCA_v_dn
      integer :: nfeast,nng, feast_calls
      !integer :: dilo,dihi,djlo,djhi
      !integer :: itrczoop(10),ip
      integer :: fvaral,fvarl,fvar,spal,spl,sp,spy,lc,ac,gr,isp,izoop,klev
      integer :: ictr,jctr
      real(r8) :: ftstp
      real(r8) :: NNstart, CFstart, CAstart, WW, Navail, Nloss, Ngloss, rec
      real(r8), dimension(TOT_LINKS):: NNbase,CFbase,CAbase,Npromote
      real(r8), dimension(TOT_FEAST):: eggs
      !real :: damper,dampex
      !real :: cff
      !real :: FXF(IminS:ImaxS,JminS:JmaxS),FEF(IminS:ImaxS,JminS:JmaxS)
      !real :: happy(IminS:ImaxS,JminS:JmaxS)
      !real :: happy2(IminS:ImaxS,JminS:JmaxS)
      !real :: happy3(IminS:ImaxS,JminS:JmaxS)
      !real :: hape(IminS:ImaxS,JminS:JmaxS),hapx(IminS:ImaxS,JminS:JmaxS)
      !real :: hapcente(IminS:ImaxS,JminS:JmaxS),hapcentx(IminS:ImaxS,JminS:JmaxS)
      !real :: v_swim(IminS:ImaxS,JminS:JmaxS),u_swim(IminS:ImaxS,JminS:JmaxS)
      !real :: hapgrad,maxswim,actswim,zoofac,crowdfac,locfac,tempfac
      !real :: mid_var1,mid_var2,mid_var3,mid_var4,var1,var2,var3,var4

      real :: CurD

#endif
      real(r8) :: predSumCop, predSumNCaS, predSumEupS, predSumNCaO, predSumEupO

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)

      !  Local variable declarations.

      integer :: i, j, k, ibio, ibio2,itr, itrmx, itrc, itrc2
      real(r8) :: cff5,cff6,cff6b,cff7,cff8,cff9,cff10,cff11
#ifdef FEAST
      integer :: iv, spn,lcp,elder,younger
#endif
#if defined BENTHIC
      integer :: ibioB
      real(r8) :: bf,fbase,TSS,Ifs,atss,btss,SF
      real(r8) :: avgD,avgDF,avgPS,avgPL,dw,wcPS,wcPL,wcD,wcDF,PSsum
      real(r8) :: totD, totDF, totPS, totPL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: frac1, frac2
      real(r8) ::sumD,sumDF,sumPL
#endif
#ifdef ICE_BIO
      integer :: ibioBI
#endif
      integer :: Iter,is
      integer :: iday, month, year

      real(r8) :: cff0,cff1, cff1b,cff2,cff3,cff4,dz
      real(r8) :: TFMZS,TFMZL,TFCop,TFNCa,TFEup,TFJel
      real(r8) :: Drate, Pmax,NOup, NHup,offset
      real(r8) :: dtdays,Ra,Rf
      real(r8) :: LightLim,NOLim,NHLim,IronLim
      real(r8) :: hour,yday,lat,k_phy,Dl,Par1,k_extV,k_chlV
      real(r8) :: Sal1,Temp1
      real(r8) :: ParMax
      real(r8) :: BasalMetMZL, BasalMetCop, BasalMetNC, BasalMetCM, BasalMetEup
      real(r8) :: Iron1,kfePh,respPh,BasalMet,BasalMetJel
      real(r8) :: PON,Dep1,Nitrif,NH4R
      real(r8) :: NitrifMax,DLNitrif

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: DBio
!       real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak
#if defined PROD3
      real(r8), dimension(IminS:ImaxS,N(ng),NBPT3) :: Prod
#endif
#if defined PROD2
      real(r8), dimension(IminS:ImaxS,NBPT2) :: Prod2
#endif
#if defined STATIONARY
      real(r8), dimension(IminS:ImaxS,N(ng),NBTS) :: Stat3
#endif
#if defined STATIONARY2
      real(r8), dimension(IminS:ImaxS,NBTS2) :: Stat2
#endif

#if defined BENTHIC
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: BioB
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: DBioB
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: Bio_bakB
#endif
#if defined ICE_BIO
      real(r8), dimension(IminS:ImaxS,NIceT(ng)) :: BioBI
      real(r8), dimension(IminS:ImaxS,NIceT(ng)) :: DBioBI
      real(r8), dimension(IminS:ImaxS,NIceT(ng)) :: Bio_bakBI
#endif
#if defined BIOFLUX
      real(r8), dimension(NT(ng),NT(ng)) :: BioFlx
#endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS) :: PARs

      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dens
      real(r8), dimension(N(ng)) :: DensV
      real(r8), dimension(N(ng)) :: ZW_V
      real(r8), dimension(IminS:ImaxS) :: StabParam
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TestVal
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncCop
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncNeo
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncEup

#ifdef JELLY
      real(r8),dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng))::TempFuncJel
#endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: HzL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: z_wL

      real(r8), dimension(IminS:ImaxS,N(ng)) :: sinkIN,sinkOUT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: riseIN,riseOUT
#ifdef ICE_BIO
      real(r8) :: aiceIfrac,aiceNfrac,dhicedt,trs,cwi,twi
      real(r8) ::grow1, GROWAice,reN,fNO3,RAi0,RgAi
      real(r8) :: sb, gesi
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng),3):: aib
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ice_thick, ice_status
#endif
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif

#ifdef DIAPAUSE
      logical :: downwardNC = .false., upwardNC = .false.
      logical :: downwardCM = .false., upwardCM = .false.
#endif

      real(r8), parameter :: eps  = 1.0E-20_r8
      real(r8), parameter :: minv = 0.0E-20_r8

      real(r8) :: Alpha
      real(r8) :: ALPHA_N,ALPHA_P, kN, kP
      real(r8) ::respNC, respCM, eCM, eNC

      ! Vertical movement

      real(r8), dimension(1,N(ng)) :: dBtmp
      real(r8) :: flxtmp

      real(r8) :: RSNC, RENC, SSNC, SENC, RSCM, RECM, SSCM, SECM

      ! Bio tracer setup

      integer  :: iiNO3,    iiNH4,    iiPhS,  iiPhL,  iiMZS, iiMZL, iiCop
      integer  :: iiNCaS,   iiNCaO,   iiEupS, iiEupO, iiDet, iiDetF
      integer  :: iiJel,    iiFe,     iiBen,  iiBenDet
      integer  :: iiIcePhL, iiIceNO3, iiIceNH4
      real(r8), dimension(IminS:ImaxS,N(ng),20) :: Bio3d, Bio2d, Bio_bak ! TODO: add DBio
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Temp, Salt

      ! Intermediate fluxes

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gpp_NO3_PhS, Gpp_NO3_PhL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gpp_NH4_PhS, Gpp_NH4_PhL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_PhS_MZL, Gra_PhL_MZL, Ege_MZL_Det
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_PhS_Cop, Gra_PhL_Cop, Gra_MZL_Cop, Gra_IPhL_Cop, Ege_Cop_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_PhS_NCaS, Gra_PhL_NCaS, Gra_MZL_NCaS, Gra_IPhL_NCaS, Ege_NCaS_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_PhS_NCaO, Gra_PhL_NCaO, Gra_MZL_NCaO, Gra_IPhL_NCaO, Ege_NCaO_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_PhS_EupS, Gra_PhL_EupS, Gra_MZL_EupS, Gra_Cop_EupS, Gra_IPhL_EupS, Gra_Det_EupS, Gra_DetF_EupS, Ege_EupS_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_PhS_EupO, Gra_PhL_EupO, Gra_MZL_EupO, Gra_Cop_EupO, Gra_IPhL_EupO, Gra_Det_EupO, Gra_DetF_EupO, Ege_EupO_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_Cop_Jel, Gra_EupS_Jel, Gra_EupO_Jel, Gra_NCaS_Jel, Gra_NCaO_Jel, Ege_Jel_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Mor_PhS_Det, Mor_PhL_Det, Mor_MZL_Det, Mor_Cop_DetF, Mor_NCaS_DetF, Mor_EupS_DetF, Mor_NCaO_DetF, Mor_EupO_DetF, Mor_Jel_DetF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Res_PhS_NH4, Res_PhL_NH4,Res_MZL_NH4, Res_Cop_NH4, Res_NCaS_NH4, Res_NCaO_NH4, Res_EupS_NH4, Res_EupO_NH4, Res_Jel_NH4
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dec_Det_NH4, Dec_DetF_NH4, Dec_NH4_NO3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_Det_Ben,Gra_DetF_Ben, Gra_PhS_Ben, Gra_PhL_Ben, Gra_BenDet_Ben, Exc_Ben_NH4, Exc_Ben_BenDet, Res_Ben_NH4, Mor_Ben_BenDet, Pre_Ben_BenDet, Dec_BenDet_NH4
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gpp_INO3_IPhL, Gpp_INH4_IPhL, Res_IPhL_INH4, Mor_IPhL_INH4, Dec_INH4_INO3, Twi_IPhL_PhL, Twi_INO3_NO3, Twi_INH4_NH4

      ! Biological source/sinks

      real(r8) :: LightLimS, NOLimS, NHLimS, IronLimS
      real(r8) :: LightLimL, NOLimL, NHLimL, IronLimL
      real(r8) :: alphaPhSv, alphaPhLv, DrateS, DrateL, PmaxS, PmaxL, PmaxsS, PmaxsL
      real(r8) :: IcePhlAvail
      real(r8), dimension(IminS:ImaxS,N(ng)) :: BasMetMZL, BasMetCop, BasMetNC, BasMetCM, BasMetEup
      real(r8) :: ParW

      !==================================================================
      !  SOME SETUP APPLICABLE TO ALL GRID CELLS
      !==================================================================

#include "set_bounds.h"

      ! Extract year, yday, month, day, and hour

      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)

      ! A few parameters...

      dtdays = dt(ng)*sec2day/REAL(BioIter(ng),r8)  ! time step, in days

#ifdef DIAPAUSE

      ! Copepod diapause is determined by time of year, based on sinking/
      ! rising day-of-year input parameters.  Set movement
      ! direction flags for on- and offshore large copepods here, and
      ! lower respiration rates if they're in the diapause (downward)
      ! phase.
      ! NCaS = CM = mostly C. marshallae, on-shelf
      ! NCaO = NC = mostly Neocalanus, off-shelf

      RSNC = MOD(RiseStart, 366.0_r8)
      RENC = MOD(RiseEnd,   366.0_r8)
      SSNC = MOD(SinkStart, 366.0_r8)
      SENC = MOD(SinkEnd,   366.0_r8)

      if ((RiseStartCM .eq. 0.0_r8) .and. (RiseEndCM .eq. 0.0_r8) .and. &
        (SinkStartCM .eq. 0.0_r8) .and. (SinkEndCM .eq. 0.0_r8)) then

        ! All 0 is the shortcut for lagging the onshelf group movement
        ! 1 month behind the offshelf group (this was the original
        ! hard-coded behavior, and I wanted to maintain
        ! back-compatibility with an input parameter file that doesn't
        ! include the newer Rise/SinkCM parameters)

        RSCM = MOD(RiseStart + 30, 366.0_r8)
        RECM = MOD(RiseEnd   + 30, 366.0_r8)
        SSCM = MOD(SinkStart + 30, 366.0_r8)
        SECM = MOD(SinkEnd   + 30, 366.0_r8)

      else

        RSCM = MOD(RiseStartCM, 366.0_r8)
        RECM = MOD(RiseEndCM,   366.0_r8)
        SSCM = MOD(SinkStartCM, 366.0_r8)
        SECM = MOD(SinkEndCM,   366.0_r8)

      endif

      upwardNC =   ((RSNC.lt.RENC) .and.                                &
     &              (yday.ge.RSNC .and. yday.le.RENC))                  &
     &             .or.                                                 &
     &             ((RSNC.gt.RENC) .and.                                &
     &              (yday.ge.RSNC .or.  yday.le.RENC))

      upwardCM =   ((RSCM.lt.RECM) .and.                                &
     &              (yday.ge.RSCM .and. yday.le.RECM))                  &
     &             .or.                                                 &
     &             ((RSCM.gt.RECM) .and.                                &
     &              (yday.ge.RSCM .or.  yday.le.RECM))

      downwardNC = ((SSNC.lt.SENC) .and.                                &
     &              (yday.ge.SSNC .and. yday.le.SENC))                  &
     &             .or.                                                 &
     &             ((SSNC.gt.SENC) .and.                                &
     &              (yday.ge.SSNC .or.  yday.le.SENC))

      downwardCM = ((SSCM.lt.SECM) .and.                                &
     &              (yday.ge.SSCM .and. yday.le.SECM))                  &
     &             .or.                                                 &
     &             ((SSCM.gt.SECM) .and.                                &
     &              (yday.ge.SSCM .or.  yday.le.SECM))

      if (downwardNC) then
        respNC = respNCa * 0.1_r8
        eNC = 0
      else
        respNC = respNCa
        eNC = eNCa
      end if

      if (downwardCM) then
        respCM = respNCa * 0.1_r8
        eCM = 0
      else
        respCM = respNCa
        eCM = eNCa
      end if

#endif


      !==================================================================
      !  JLOOP: BEGIN HORIZONTAL LOOP
      !==================================================================

      J_LOOP : DO j=Jstr,Jend

        !---------------------------------
        ! Biological tracer variable setup
        !---------------------------------

        ! The various biological state variables are passed into this
        ! function in a few different arrays, depending on whether the
        ! variable was part of the original GOANPZ model, added with the
        ! benthic submodel, or added with the ice model (and for ice,
        ! whether we're running ROMS with a full ice model or with 1D
        ! climatological ice).
        !
        ! To make long-term maintenance of this code easier, we'll place
        ! all these variables into a single i x k x var array, where
        ! i = horizontal grid cell looping dimension, k = depth (counting
        ! from bottom to top), and var is the biological state variable
        ! index.

        ! First, some handy indices, so I don't have to switch around
        ! between the 3 different sets used for pelagic, benthic, and ice
        ! variables in the input arrays.

        iiNO3    = 1
        iiNH4    = 2
        iiPhS    = 3
        iiPhL    = 4
        iiMZS    = 5
        iiMZL    = 6
        iiCop    = 7
        iiNCaS   = 8
        iiEupS   = 9
        iiNCaO   = 10
        iiEupO   = 11
        iiDet    = 12
        iiDetF   = 13
        iiJel    = 14
        iiFe     = 15
        iiBen    = 16
        iiBenDet = 17
        iiIcePhL = 18
        iiIceNO3 = 19
        iiIceNH4 = 20

        ! All state variables will be saved in two different versions:
        ! per-volume (Bio3d) and per-area (Bio2d).  This redundancy makes
        ! for clearer (for human readers) code.

        Bio3d = 0 ! Initialize to 0
        Bio2d = 0

        ! Pelagic variables: These are originally stored in per-volume
        ! concentrations in each water column layer. NO3 and NH4 are in
        ! mmol N m^-3, Fe is in umol Fe m^-3, and the rest are in mg C
        ! m^-3.

        DO itrc=1,NBT  ! Pelagic variables
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio3d(i,k,itrc)=max(t(i,j,k,nstp,idbio(itrc)),0.0_r8)
              Bio2d(i,k,itrc)=Bio3d(i,k,itrc)*Hz(i,j,k)
            END DO
          END DO
        END DO

        ! Benthic variables: These are originally stored in per-area
        ! concentrations in each benthic layer, in mg C m^-2.  The
        ! benthic layers are of unknown thickness.
        !
        ! At the moment, BESTNPZ hard-codes the number of benthic layers
        ! (NBL) to 1; for bookkeeping purposes, we're going to store
        ! benthic biomass in the bottom layer of our Bio3d/2d arrays.  If
        ! we ever change the number of benthic layers, we may need to
        ! rethink this schema.

        DO itrc=1,NBEN
          ibioB=idben(itrc) ! Note: idben(i) = i
          DO k=1,NBL(ng) ! Note: For BESTNPZ, NBL = 1 is hard-coded in mod_param.F
            DO i=Istr,Iend
              Bio2d(i,k,itrc+NBT)=bt(i,j,k,nstp,ibioB) ! TODO: not restricted to >0?
              Bio3d(i,k,itrc+NBT)=Bio2d(i,k,itrc+NBT)/Hz(i,j,k)
            END DO
          END DO
        END DO

#ifdef ICE_BIO
        ! Before we get to the ice variables, we'll collect some info
        ! about the ice itself: ice thickness, and status (i.e. whether
        ! ice has appeared or disappeared between this time step and the
        ! last).
        !   ice_status =  2  ice present at this time step and previous
        !   ice_status =  1  ice appeared at this step
        !   ice_status =  0  no ice at either
        !   ice_status = -1  ice disappeared at this step

        DO i=Istr,Iend
# if defined CLIM_ICE_1D

          ! Ice thickness

          ice_thick(i,j) = MAX(0.0_r8,tclm(3,3,N(ng),i1CI))

          ! Ice status

          if (ice_thick(i,j).gt.aidz) THEN
            itL(i,j,nstp,iIceLog) =1.0_r8
          else
            itL(i,j,nstp,iIceLog) =-1.0_r8
          endif

          if     (itL(i,j,nstp,iIceLog).gt.0 .and. itL(i,j,nnew,iIceLog).le.0) THEN
            ice_status(i,j) = 1.0
          elseif (itL(i,j,nstp,iIceLog).gt.0 .and. itL(i,j,nnew,iIceLog).gt.0) THEN
            ice_status(i,j) = 2.0
          elseif (itL(i,j,nstp,iIceLog).le.0 .and. itL(i,j,nnew,iIceLog).gt.0) THEN
            ice_status(i,j) = -1.0
          else
            ice_status(i,j) = 0.0
          endif

# elif defined BERING_10K

          ! Ice thickness

          if (hi(i,j,nstp).gt.0.0_r8)THEN
            ice_thick(i,j) = hi(i,j,nstp)
          else
            ice_thick(i,j)=0.0_r8
          end if

          ! Ice status

          cff1=IceLog(i,j,nnew)
          cff2=IceLog(i,j,nstp)

          IceLog(i,j,nnew)=cff2
          IceLog(i,j,nstp)=cff1

          if     (IceLog(i,j,nstp).gt.0 .and. IceLog(i,j,nnew).le.0) THEN
            ice_status(i,j) = 1.0
          elseif (IceLog(i,j,nstp).gt.0 .and. IceLog(i,j,nnew).gt.0) THEN
            ice_status(i,j) = 2.0
          elseif (IceLog(i,j,nstp).le.0 .and. IceLog(i,j,nnew).gt.0) THEN
            ice_status(i,j) = -1.0
          else
            ice_status(i,j) = 0.0
          endif

# endif
        END DO
#endif

        ! Ice variables: Ice variables are passed into this function in
        ! different arrays depending on whether the model is running with
        ! climatological one-dimensional ice, or coupled to a full ice
        ! model. In either case, they're originally stored in per-volume
        ! concentrations in a single ice skeletal layer of prescribed
        ! thickness (aidz).  For bookkeeping, we'll put this in the top
        ! layer of the Bio3d/2d arrays, but remember that the conversion
        ! factor assumes a different layer thickness than the top water
        ! layer.

        ! If there ice is present in both this step and the last, start
        ! by extracting ice biomass from the main ice bio tracer array.
        ! Otherwise, no biomass to start (we'll deal with changing ice in
        ! a moment).

#ifdef ICE_BIO
        DO i=Istr,Iend
          if (ice_status(i,j) .eq. 2.0) then

# ifdef CLIM_ICE_1D
            Bio3d(i,N(ng),iiIcePhL) = max(it(i,j,nstp,idice(1)), 0.0_r8)
            Bio3d(i,N(ng),iiIceNO3) = max(it(i,j,nstp,idice(2)), 0.0_r8)
            Bio3d(i,N(ng),iiIceNH4) = max(it(i,j,nstp,idice(3)), 0.0_r8)
# elif defined BERING_10K
            Bio3d(i,N(ng),iiIcePhL) = max(0.0_r8, IcePhL(i,j,nstp))
            Bio3d(i,N(ng),iiIceNO3) = max(0.0_r8, IceNO3(i,j,nstp))
            Bio3d(i,N(ng),iiIceNH4) = max(0.0_r8, IceNH4(i,j,nstp))
# endif
          endif

          Bio2d(i,N(ng),iiIcePhL) = Bio3d(i,N(ng),iiIcePhL)*aidz
          Bio2d(i,N(ng),iiIceNO3) = Bio3d(i,N(ng),iiIceNO3)*aidz
          Bio2d(i,N(ng),iiIceNH4) = Bio3d(i,N(ng),iiIceNH4)*aidz

        END DO
#endif

        ! Temperature and salinity, for easier reference

        Temp = t(Istr:Iend,j,1:N(ng),nstp,itemp)
        Salt = t(Istr:Iend,j,1:N(ng),nstp,isalt)

#ifdef CORRECT_TEMP_BIAS
        Temp = Temp - 1.94_r8 ! bias correction for bio only, not fed back
#endif

        ! Initialize the rate of change, dB/dt, to 0 for all elements.
        ! Same for all intermediate flux arrays.  Note that these fluxes
        ! will hold the 2D equivalent of all the fluxes (i.e. per area,
        ! rather than per volume); this makes it easier to keep track of
        ! things that are moving between different-sized layers (e.g. ice
        ! to surface layer, or benthos to water column)

        DBio = 0 ! Initializes entire array to 0

        Gpp_NO3_PhS    = 0
        Gpp_NO3_PhL    = 0
        Gpp_NH4_PhS    = 0
        Gpp_NH4_PhL    = 0
        Gra_PhS_MZL    = 0
        Gra_PhL_MZL    = 0
        Ege_MZL_Det    = 0
        Gra_PhS_Cop    = 0
        Gra_PhL_Cop    = 0
        Gra_MZL_Cop    = 0
        Gra_IPhL_Cop   = 0
        Ege_Cop_DetF   = 0
        Gra_PhS_NCaS   = 0
        Gra_PhL_NCaS   = 0
        Gra_MZL_NCaS   = 0
        Gra_IPhL_NCaS  = 0
        Ege_NCaS_DetF  = 0
        Gra_PhS_NCaO   = 0
        Gra_PhL_NCaO   = 0
        Gra_MZL_NCaO   = 0
        Gra_IPhL_NCaO  = 0
        Ege_NCaO_DetF  = 0
        Gra_PhS_EupS   = 0
        Gra_PhL_EupS   = 0
        Gra_MZL_EupS   = 0
        Gra_Cop_EupS   = 0
        Gra_IPhL_EupS  = 0
        Gra_Det_EupS   = 0
        Gra_DetF_EupS  = 0
        Ege_EupS_DetF  = 0
        Gra_PhS_EupO   = 0
        Gra_PhL_EupO   = 0
        Gra_MZL_EupO   = 0
        Gra_Cop_EupO   = 0
        Gra_IPhL_EupO  = 0
        Gra_Det_EupO   = 0
        Gra_DetF_EupO  = 0
        Ege_EupO_DetF  = 0
        Gra_Cop_Jel    = 0
        Gra_EupS_Jel   = 0
        Gra_EupO_Jel   = 0
        Gra_NCaS_Jel   = 0
        Gra_NCaO_Jel   = 0
        Ege_Jel_DetF   = 0
        Mor_PhS_Det    = 0
        Mor_PhL_Det    = 0
        Mor_MZL_Det    = 0
        Mor_Cop_DetF   = 0
        Mor_NCaS_DetF  = 0
        Mor_EupS_DetF  = 0
        Mor_NCaO_DetF  = 0
        Mor_EupO_DetF  = 0
        Mor_Jel_DetF   = 0
        Res_PhS_NH4    = 0
        Res_PhL_NH4    = 0
        Res_MZL_NH4    = 0
        Res_Cop_NH4    = 0
        Res_NCaS_NH4   = 0
        Res_NCaO_NH4   = 0
        Res_EupS_NH4   = 0
        Res_EupO_NH4   = 0
        Res_Jel_NH4    = 0
        Dec_Det_NH4    = 0
        Dec_DetF_NH4   = 0
        Dec_NH4_NO3    = 0
        Gra_Det_Ben    = 0
        Gra_DetF_Ben   = 0
        Gra_PhS_Ben    = 0
        Gra_PhL_Ben    = 0
        Gra_BenDet_Ben = 0
        Exc_Ben_NH4    = 0
        Exc_Ben_BenDet = 0
        Res_Ben_NH4    = 0
        Mor_Ben_BenDet = 0
        Pre_Ben_BenDet = 0
        Dec_BenDet_NH4 = 0
        Gpp_INO3_IPhL  = 0
        Gpp_INH4_IPhL  = 0
        Res_IPhL_INH4  = 0
        Mor_IPhL_INH4  = 0
        Dec_INH4_INO3  = 0
        Twi_IPhL_PhL   = 0
        Twi_INO3_NO3   = 0
        Twi_INH4_NH4   = 0




        ! Save a copy of the original biomass

        Bio_bak = Bio2d

#ifdef ICE_BIO
        ! Move tracers between surface water layer and ice skeletal layer
        ! if ice appeared or disappeared

        DO i=Istr,Iend
          if (ice_status(i,j) .eq. 1.0) then

            ! If new ice appeared, assume the biomass from the surface
            ! layer in the previous step is now spread evenly across the
            ! surface water column and the ice skeletal layer.

            Bio3d(i,N(ng),iiPhL) = Bio2d(i,N(ng),iiPhL)/(Hz(i,j,N(ng))+aidz)
            Bio3d(i,N(ng),iiNO3) = Bio2d(i,N(ng),iiNO3)/(Hz(i,j,N(ng))+aidz)
            Bio3d(i,N(ng),iiNH4) = Bio2d(i,N(ng),iiNH4)/(Hz(i,j,N(ng))+aidz)

            Bio3d(i,N(ng),iiIcePhL) = Bio3d(i,N(ng),iiPhl)
            Bio3d(i,N(ng),iiIceNO3) = Bio3d(i,N(ng),iiNO3)
            Bio3d(i,N(ng),iiIceNH4) = Bio3d(i,N(ng),iiNH4)

            Bio2d(i,N(ng),iiPhL)    = Bio3d(i,N(ng),iiPhl) * Hz(i,j,N(ng))
            Bio2d(i,N(ng),iiNO3)    = Bio3d(i,N(ng),iiNO3) * Hz(i,j,N(ng))
            Bio2d(i,N(ng),iiNH4)    = Bio3d(i,N(ng),iiNH4) * Hz(i,j,N(ng))

            Bio2d(i,N(ng),iiIcePhL) = Bio3d(i,N(ng),iiIcePhl) * aidz
            Bio2d(i,N(ng),iiIceNO3) = Bio3d(i,N(ng),iiIceNO3) * aidz
            Bio2d(i,N(ng),iiIceNH4) = Bio3d(i,N(ng),iiIceNH4) * aidz

          elseif (ice_status(i,j) .eq. -1.0) then

            ! If ice disappeared, biomass that was in the ice gets dumped
            ! into the water surface layer

            Bio2d(i,N(ng),iiPhL) = Bio2d(i,N(ng),iiPhL) + Bio2d(i,N(ng),iiIcePhL)
            Bio2d(i,N(ng),iiNO3) = Bio2d(i,N(ng),iiNO3) + Bio2d(i,N(ng),iiIceNO3)
            Bio2d(i,N(ng),iiNH4) = Bio2d(i,N(ng),iiNH4) + Bio2d(i,N(ng),iiIceNH4)

            Bio2d(i,N(ng),iiIcePhL) = 0.0_r8
            Bio2d(i,N(ng),iiIceNO3) = 0.0_r8
            Bio2d(i,N(ng),iiIceNH4) = 0.0_r8

            Bio3d(i,N(ng),iiPhL) = Bio2d(i,N(ng),iiPhL)/Hz(i,j,N(ng))
            Bio3d(i,N(ng),iiNO3) = Bio2d(i,N(ng),iiNO3)/Hz(i,j,N(ng))
            Bio3d(i,N(ng),iiNH4) = Bio2d(i,N(ng),iiNH4)/Hz(i,j,N(ng))

            Bio3d(i,N(ng),iiIcePhL) = 0.0_r8
            Bio3d(i,N(ng),iiIceNO3) = 0.0_r8
            Bio3d(i,N(ng),iiIceNH4) = 0.0_r8

          endif
        END DO
#endif

        ! Calculate inverse layer thickness

        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO

        !-------------------------------------
        ! Calculate Day Length and Surface PAR
        !-------------------------------------

        DO i=Istr,Iend

#if defined DIURNAL_SRFLUX
          ! Calculate Day Length:  Day Length is already accounted for
          ! in ANA_SWRAD so disable correction

          Dl = 24.0_r8
#else
          ! Day Length calculation (orig from R. Davis) from latitude and
          ! declination.

          lat = GRID(ng) % latr(i,j)
          cff1 = 2.0_r8 * pi * ( yday-1.0_r8 ) / 365.0_r8
          cff2 = 0.006918_r8 - 0.399912_r8*cos(cff1)                    &
     &           + 0.070257_r8*sin(cff1) - 0.006758_r8*cos(2*cff1)      &
     &           + 0.000907_r8*sin(2*cff1) - 0.002697_r8*cos(3*cff1)    &
     &           + 0.00148_r8*sin(3*cff1) ! Solar declination from Oberhuber (1988) (COADS documentation)
          cff3 = lat * pi /180.0_r8
          IF ( abs( -tan(cff3)*tan(cff2) ) .le. 1.0_r8 ) THEN
            cff1 = acos( -tan(cff3)*tan(cff2) ) * 180.0_r8 / pi
            Dl = 2.0_r8 / 15.0_r8 * cff1
          ELSE
            IF ( yday.gt.90.0_r8 .and. yday.lt.270.0_r8 ) THEN
              Dl = 24.0_r8
            ELSE
              Dl = 0.0_r8
            END IF
          END IF
#endif

#ifdef KODIAK_IRAD

          ! Calculate PAR at the surface
          ! Eyeball fit of data from Hinckley''s ezeroday.dat (E d-1 m-2)
          ! TODO: this option seems to be missing the actual setting of PARs(i)??

          cff2 = 41.0_r8 - 35.0_r8                                      &
     &           * COS( ( 12.0_r8 + yday) * 2.0_r8 * pi / 365.0_r8 )
#else

          ! Calculate PAR at the surface: use Shortwave radiation
          ! ( = surface solar irradiance) converted from deg C m/s to
          ! E/m2/day   rho0=1025 Cp=3985

          PARs(i) =PARfrac(ng) * srflx(i,j) * rho0 * Cp * 0.394848_r8
#endif
        END DO

        !------------------------------------------
        ! Calculate light decay in the water column
        !------------------------------------------

#ifdef NEWSHADE
        ! Georgina Gibsons version after Morel 1988 (in Loukos 1997)

        DO i=Istr,Iend

          cff10=grid(ng) % h(i,j)
          k_extV= k_ext+2.00_r8*exp(-cff10*.05)

          cff0=PARs(i)

          DO k=N(ng),1,-1

            dz=0.5_r8*(z_w(i,j,k)-z_w(i,j,k-1))
            cff5=(Bio3d(i,k,iiPhS)/ccr)+ (Bio3d(i,k,iiPhL)/ccrPhL)
            cff2 = (k_chl*(cff5)**(0.428_r8))
            !cff2 = min(0.05_r8,max(0.0067_r8,(k_chl*(cff5)**(-0.428_r8))))
            PAR(i,k) = cff0 * EXP(-(k_extV+cff2)*dz)
            cff0=cff0 * EXP(-(k_extV+cff2)*dz*2.0_r8)

          END DO
        END DO

#elif defined COKELET
        ! Version from Ned Cokelet

        DO i=Istr,Iend

          cff10=grid(ng) % h(i,j)
          !  k_extV= k_ext+k_extZ*exp(-cff10*.05)
          k_extV= k_ext
          cff0=PARs(i)

          DO k=N(ng),1,-1

            dz=0.5_r8*(z_w(i,j,k)-z_w(i,j,k-1))
            cff5=(Bio3d(i,k,iiPhS)/ccr)+ (Bio3d(i,k,iiPhL)/ccrPhL)
            cff2 = (k_chlA*(cff5)**(k_chlB))

            PAR(i,k) = cff0 * EXP(-(k_extV+cff2)*dz)
            cff0=cff0 * EXP(-(k_extV+cff2)*dz*2.0_r8)

          END DO
        END DO
#else
        ! Version from Sarah Hinckley old C code
        DO k=N(ng),1,-1
          DO i=Istr,Iend
            cff3 = z_r(i,j,k)+2.5_r8
            IF ( cff3 .gt. -71.0_r8 ) THEN
              cff1 = k_ext + k_chl *                                 &
     &                  ( Bio3d(i,k,iiPhS) + Bio3d(i,k,iiPhL) ) / ccr
            ELSE
                cff1 = 0.077_r8
            END IF
            PAR(i,k) = PARfrac(ng) * cff2 * exp( cff1 * cff3 )


          END DO
        END DO
#endif


        !================================================================
        !  Begin time loop (if BioIter > 1, this divides the main time
        !  step into smaller steps for biological calculations)
        !================================================================

        ITER_LOOP: DO Iter=1,BioIter(ng)

          !==============================================================
          !  Biological Source/Sink terms.
          !==============================================================

          !------------------------------
          ! Phytoplankton production
          !------------------------------

          LightLimS = 1.0_r8
          NOLimS    = 1.0_r8
          NHLimS    = 1.0_r8
          IronLimS  = 1.0_r8
          LightLimL = 1.0_r8
          NOLimL    = 1.0_r8
          NHLimL    = 1.0_r8
          IronLimL  = 1.0_r8

          DO k=1,N(ng)
            DO i=Istr,Iend

              ! Slope of P-I curve

              if (PARs(i).lt.30.0) then
                alphaPhSv = 18
                alphaPhLv = 10
              elseif (PARs(i).gt.40.0) then
                alphaPhSv = 5.6
                alphaPhLv = 2.2
              else
                alphaPhSv = 18.0-((18.0-5.6)/(40.0-30.0))*(PARs(i)-30.0)
                alphaPhLv = 10.0-((10.0-2.2)/(40.0-30.0))*(PARs(i)-30.0)
              end if

              ! Maximum uptake rate

              DrateS = DiS * 10.0_r8 ** (DpS * Temp(i,k))
              DrateL = DiL * 10.0_r8 ** (DpL * Temp(i,k))

              PmaxS = (2.0_r8 ** DrateS - 1.0_r8 )   ! maximum daily mass specific growth rate FROST (1987)
              PmaxL = (2.0_r8 ** DrateL - 1.0_r8 )

              PmaxsS=PmaxS*ccr                       ! max chla specific growth rate from FROST (1987)
              PmaxsL=PmaxL*ccrPhL

#ifdef DENMAN
              ! N03 limitation following Denman

              NOLimS = (Bio3d(i,k,iiNO3) + Bio3d(i,k,iiNH4))/(k1PhS + Bio3d(i,k,iiNO3) + Bio3d(i,k,iiNH4))
              NOLimL = (Bio3d(i,k,iiNO3) + Bio3d(i,k,iiNH4))/(k1PhL + Bio3d(i,k,iiNO3) + Bio3d(i,k,iiNH4))
#else
              ! NO3 limitation following Lomas (Marine Biology 1999)

              NOLimS = (Bio3d(i,k,iiNO3)/(k1PhS + Bio3d(i,k,iiNO3))) * (1-(0.8_r8*Bio3d(i,k,iiNH4)/(k2PhS + Bio3d(i,k,iiNH4))))
              NOLimL = (Bio3d(i,k,iiNO3)/(k1PhL + Bio3d(i,k,iiNO3))) * (1-(0.8_r8*Bio3d(i,k,iiNH4)/(k2PhL + Bio3d(i,k,iiNH4))))
#endif
#ifdef IRON_LIMIT

              ! Iron limitation

              IronLimS = eps + Bio3d(i,k,iiFe)/(kfePhS + Bio3d(i,k,iiFe))*(kfePhS + FeCritPS)/FeCritPS
              IronLimL = eps + Bio3d(i,k,iiFe)/(kfePhL + Bio3d(i,k,iiFe))*(kfePhL + FeCritPL)/FeCritPL
#endif

              ! Light limitation

#ifdef DENMAN
              LightLimS = TANH(alphaPhSv * PAR(i,k) / PmaxsS)
              LightLimL = TANH(alphaPhLv * PAR(i,k) / PmaxsL) ! TODO: alphaPhLv was Alpha in orig code, but the code to set its value was commented out.
#else
              OffSet = 0.0_r8
              LightLimS = TANH(alphaPhSv * MAX(PAR(i,k) - OffSet,0.0_r8)/PmaxsS)
              LightLimL = TANH(alphaPhLv * MAX(PAR(i,k) - OffSet,0.0_r8)/PmaxsL)
#endif
#ifdef DENMAN
              ! Uptake of NO3 and NH4

              cff1 = 0.2/(0.2+Bio3d(i,k,iiNH4))
              cff2 = cff1 * Bio3d(i,k,iiNO3)/(Bio3d(i,k,iiNO3)+Bio3d(i,k,iiNH4))

              cff3 = PmaxS*LightLimS*NOLimS*IronLimS
              Gpp_NO3_PhS(i,k) = Bio3d(i,k,iiPhS) * cff3 * cff2       ! mg C m^-3 d^-1
              Gpp_NH4_PhS(i,k) = Bio3d(i,k,iiPhS) * cff3 * (1 - cff2) ! mg C m^-3 d^-1

              cff3 = PmaxL*LightLimL*NOLimL*IronLimL
              Gpp_NO3_PhL(i,k) = Bio3d(i,k,iiPhL) * cff3 * cff2       ! mg C m^-3 d^-1
              Gpp_NH4_PhL(i,k) = Bio3d(i,k,iiPhL) * cff3 * (1 - cff2) ! mg C m^-3 d^-1
#else
              ! Ammonium limitation

              NHLimS = Bio3d(i,k,iiNH4) / ( k2PhS + Bio3d(i,k,iiNH4) )
              NHLimL = Bio3d(i,k,iiNH4) / ( k2PhL + Bio3d(i,k,iiNH4) )
              if((NOLimS+NHLimS).gt.1.0_r8) then
                NHLimS= 1.0_r8-NOLimS
              endif
              if((NOLimL+NHLimL).gt.1.0_r8) then
                NHLimL= 1.0_r8-NOLimL
              endif

              ! NO3 uptake

              Gpp_NO3_PhS(i,k) = MAX(0.0_r8,(Bio3d(i,k,iiPhS)/ccr)   *PmaxsS*MIN(LightLimS, NOLimS, IronLimS))  ! mg C m^-3 d^-1
              Gpp_NO3_PhL(i,k) = MAX(0.0_r8,(Bio3d(i,k,iiPhL)/ccrPhL)*PmaxsL*MIN(LightLimL, NOLimL, IronLimL))

              ! NH4 uptake

              Gpp_NH4_PhS(i,k) = MAX(0.0_r8,(Bio3d(i,k,iiPhS)/ccr)    * PmaxsS * MIN(LightLimS, NHLimS)) ! mg C m^-3 d^-1
              Gpp_NH4_PhL(i,k) = MAX(0.0_r8,(Bio3d(i,k,iiPhL)/ccrPhL) * PmaxsL * MIN(LightLimL ,NHLimL))

#endif
              ! Convert intermediate fluxes from volumetric to per area

              Gpp_NO3_PhS(i,k) = Gpp_NO3_PhS(i,k) * Hz(i,j,k) ! mg C m^-2 d^-1
              Gpp_NH4_PhS(i,k) = Gpp_NH4_PhS(i,k) * Hz(i,j,k) ! mg C m^-2 d^-1
              Gpp_NO3_PhL(i,k) = Gpp_NO3_PhL(i,k) * Hz(i,j,k) ! mg C m^-2 d^-1
              Gpp_NH4_PhL(i,k) = Gpp_NH4_PhL(i,k) * Hz(i,j,k) ! mg C m^-2 d^-1

              ! If doubling rate is 0, ignore the above (would be more
              ! efficient to check before calculating, but then I would
              ! have to split up the large/small calcs; doing it this way
              ! for the sake of maintenance and debugging)

              if (DiS.le.0_r8) THEN
                Gpp_NO3_PhS(i,k) = 0; ! mg C m^-2 d^-1
                Gpp_NH4_PhS(i,k) = 0; ! mg C m^-2 d^-1
              endif
              if (DiL.le.0_r8) THEN
                Gpp_NO3_PhL(i,k) = 0; ! mg C m^-2 d^-1
                Gpp_NH4_PhL(i,k) = 0; ! mg C m^-2 d^-1
              endif

            END DO
          END DO

          !------------------------------
          ! Grazing and predation
          !------------------------------

          DO k=1,N(ng)
            DO i=Istr,Iend

#ifdef ICE_BIO
              ! Amount of ice algae available, assuming that all algae in
              ! the bottom 2 cm (aidz=0.02) of ice is relocated to the
              ! top water layer.  Will be 0 for all but k = N(ng)

              IcePhlAvail = Bio2d(i,k,iiIcePhl)/Hz(i,j,k) ! mg C m^-3
#else
              ! Placeholder to simplify code in no-ice runs

              IcePhlAvail = 0.0_r8
#endif

              ! Microzooplankton

              cff1 = fpPhSMZL * Bio3d(i,k,iiPhS)**2 +                   &
     &               fpPhLMZL * Bio3d(i,k,iiPhL)**2
              cff2 = eMZL * Bio3d(i,k,iiMZL) / (fMZL + cff1)
              cff3 = Q10MZL**((Temp(i,k)-Q10MZLT)/10.0_r8)

              Gra_PhS_MZL(i,k) = fpPhSMZL * (Bio3d(i,k,iiPhS)**2) * cff2 * cff3 ! mg C m^-3
              Gra_PhL_MZL(i,k) = fpPhLMZL * (Bio3d(i,k,iiPhL)**2) * cff2 * cff3

              Ege_MZL_Det(i,k) = (1.0_r8 - gammaMZL) * cff1 * cff2 * cff3 ! mg C m^-2

              ! Copepods

              cff1 = fpPhSCop * Bio3d(i,k,iiPhS)**2                     &
     &             + fpPhLCop * Bio3d(i,k,iiPhL)**2                     &
     &             + fpMZLCop * Bio3d(i,k,iiMZL)**2                     &
     &             + fpPhLCop * (IcePhlAvail)**2

              cff2 = eCop * Bio3d(i,k,iiCop) / (fCop + cff1)
              cff3 = Q10Cop**((Temp(i,k)-Q10CopT)/10.0_r8)

              if (cff1.lt.0.01_r8) THEN ! Starvation response, used in Res below
                BasMetCop(i,k) = respCop*cff1/0.01_r8
              else
                BasMetCop(i,k) = respCop
              endif

              Gra_PhS_Cop(i,k)  = fpPhSCop * (Bio3d(i,k,iiPhS)**2) * cff2 * cff3
              Gra_PhL_Cop(i,k)  = fpPhLCop * (Bio3d(i,k,iiPhL)**2) * cff2 * cff3
              Gra_MZL_Cop(i,k)  = fpMZLCop * (Bio3d(i,k,iiMZL)**2) * cff2 * cff3
              Gra_IPhL_Cop(i,k) = fpPhLCop * (IcePhlAvail)**2  * cff2 * cff3

              Ege_Cop_DetF(i,k) = (1.0_r8 - gammaCop) * cff1 * cff2 * cff3

              ! On-shelf Neocalanus

              cff1 = fpPhSNCa * Bio3d(i,k,iiPhS)**2                     &
     &             + fpPhLNCa * Bio3d(i,k,iiPhL)**2                     &
     &             + fpMZLNCa * Bio3d(i,k,iiMZL)**2                     &
     &             + fpPhLNCa * (IcePhlAvail)**2

              cff2 = eNCa * Bio3d(i,k,iiNCaS) / (fNCa + cff1)
              cff3 = Q10NCa ** ((Temp(i,k)-Q10NCaT)/10.0_r8)

              if (cff1.lt.0.01_r8) THEN ! Starvation response, used in Res below
                BasMetNC(i,k) = respNC*cff1/0.01_r8
                BasMetCM(i,k) = respCM*cff1/0.01_r8
              else
                BasMetNC(i,k) = respNC
                BasMetCM(i,k) = respCM
              endif

              Gra_PhS_NCaS(i,k)  = fpPhSNCa * Bio3d(i,k,iiPhS)**2 * cff2 * cff3
              Gra_PhL_NCaS(i,k)  = fpPhLNCa * Bio3d(i,k,iiPhL)**2 * cff2 * cff3
              Gra_MZL_NCaS(i,k)  = fpMZLNCa * Bio3d(i,k,iiMZL)**2 * cff2 * cff3
              Gra_IPhL_NCaS(i,k) = fpPhLNCa * (IcePhlAvail)**2 * cff2 * cff3

              Ege_NCaS_DetF(i,k) = (1.0_r8 - gammaNCa) * cff1 * cff2 * cff3

              ! Off-shelf Neocalanus

              cff2 = eNCa * Bio3d(i,k,iiNCaO) / (fNCa + cff1)

              Gra_PhS_NCaO(i,k)  = fpPhSNCa * Bio3d(i,k,iiPhS)**2 * cff2 * cff3
              Gra_PhL_NCaO(i,k)  = fpPhLNCa * Bio3d(i,k,iiPhL)**2 * cff2 * cff3
              Gra_MZL_NCaO(i,k)  = fpMZLNCa * Bio3d(i,k,iiMZL)**2 * cff2 * cff3
              Gra_IPhL_NCaO(i,k) = fpPhLNCa * (IcePhlAvail)**2 * cff2 * cff3

              Ege_NCaO_DetF(i,k) = (1.0_r8 - gammaNCa) * cff1 * cff2 * cff3

              ! On-shelf euphausiids

              cff1 = fpPhSEup * Bio3d(i,k,iiPhS)**2                      &
     &             + fpPhLEup * Bio3d(i,k,iiPhL)**2                      &
     &             + fpMZLEup * Bio3d(i,k,iiMZL)**2                      &
     &             + fpCopEup * Bio3d(i,k,iiCop)**2                      &
     &             + fpPhLEup * (IcePhlAvail)**2    ! live food

              cff0 = fpDetEup * Bio3d(i,k,iiDet)**2                      &
     &             + fpDetEup * Bio3d(i,k,iiDetF)**2 ! detrital food

              cff2 = eEup * Bio3d(i,k,iEupS) / (fEup + cff1 + cff0)
              cff3 = Q10Eup ** ((Temp(i,k)-Q10EupT) / 10.0_r8)

              cff4 = 1.0_r8 - (0.5_r8 + 0.5_r8*tanh((grid(ng)%h(i,j) - 200_r8)/.3_r8)) ! depth-limiter, stops growth if they move deeper than 200m TODO

              if (cff1 .lt. 0.01_r8)THEN
                BasMetEup = respEup*cff1/0.01_r8 ! TODO: doesn't consider detrital food?
              else
                BasMetEup = respEup
              endif

              Gra_PhS_EupS(i,k)  = fpPhSEup * Bio3d(i,k,iiPhS)**2  * cff2 * cff3 * cff4
              Gra_PhL_EupS(i,k)  = fpPhLEup * Bio3d(i,k,iiPhL)**2  * cff2 * cff3 * cff4
              Gra_MZL_EupS(i,k)  = fpMZLEup * Bio3d(i,k,iiMZL)**2  * cff2 * cff3 * cff4
              Gra_Cop_EupS(i,k)  = fpCopEup * Bio3d(i,k,iiCop)**2  * cff2 * cff3 * cff4
              Gra_IPhL_EupS(i,k) = fpPhLEup * (IcePhlAvail)**2     * cff2 * cff3 * cff4
              Gra_Det_EupS(i,k)  = fpDetEup * Bio3d(i,k,iiDet)**2  * cff2 * cff3 * cff4
              Gra_DetF_EupS(i,k) = fpDetEup * Bio3d(i,k,iiDetF)**2 * cff2 * cff3 * cff4

              Ege_EupS_DetF(i,k) = ((1.0_r8 - gammaEup) * cff1 +        &
      &                             (1.0_r8 - 0.3_r8)   * cff0) *       &
      &                            cff2 * cff3 * cff4

              ! Off-shelf euphausiids

              cff2 = eEup * Bio3d(i,k,iEupO) / (fEup + cff1 + cff0)

              cff4 = 1.0_r8 - (0.5_r8 + 0.5_r8*tanh((200_r8 - grid(ng)%h(i,j))/.3_r8)) ! depth-limiter, stops growth if they move shallower than 200m TODO

              Gra_PhS_EupO(i,k)  = fpPhSEup * Bio3d(i,k,iiPhS)**2  * cff2 * cff3 * cff4
              Gra_PhL_EupO(i,k)  = fpPhLEup * Bio3d(i,k,iiPhL)**2  * cff2 * cff3 * cff4
              Gra_MZL_EupO(i,k)  = fpMZLEup * Bio3d(i,k,iiMZL)**2  * cff2 * cff3 * cff4
              Gra_Cop_EupO(i,k)  = fpCopEup * Bio3d(i,k,iiCop)**2  * cff2 * cff3 * cff4
              Gra_IPhL_EupO(i,k) = fpPhLEup * (IcePhlAvail)**2     * cff2 * cff3 * cff4
              Gra_Det_EupO(i,k)  = fpDetEup * Bio3d(i,k,iiDet)**2  * cff2 * cff3 * cff4
              Gra_DetF_EupO(i,k) = fpDetEup * Bio3d(i,k,iiDetF)**2 * cff2 * cff3 * cff4

              Ege_EupO_DetF(i,k) = ((1.0_r8 - gammaEup) * cff1 +        &
      &                             (1.0_r8 - 0.3_r8)   * cff0) *       &
      &                            cff2 * cff3 * cff4

              ! Jellyfish

              cff1 = fpCopJel * Bio3d(i,k,iiCop)**2 +                   &
    &                fpNCaJel * Bio3d(i,k,iiNCaS)**2 +                  &
    &                fpNCaJel * Bio3d(i,k,iiNCaO)**2 +                  &
    &                fpEupJel * Bio3d(i,k,iiEupS)**2 +                  &
    &                fpEupJel * Bio3d(i,k,iiEupO)**2

              cff2 = eJel * Bio3d(i,k,iiJel) / (fJel + cff1)
              cff3= Q10Jele ** ((Temp(i,k)-Q10JelTe) / 10.0_r8)

              Gra_Cop_Jel(i,k)  = fpCopJel * Bio3d(i,k,iiCop)**2  * cff2 * cff3
              Gra_NCaS_Jel(i,k) = fpNCaJel * Bio3d(i,k,iiNCaS)**2 * cff2 * cff3
              Gra_NCaO_Jel(i,k) = fpNCaJel * Bio3d(i,k,iiNCaO)**2 * cff2 * cff3
              Gra_EupS_Jel(i,k) = fpEupJel * Bio3d(i,k,iiEupS)**2 * cff2 * cff3
              Gra_EupO_Jel(i,k) = fpEupJel * Bio3d(i,k,iiEupO)**2 * cff2 * cff3

              ! Note: mentioned in docs that gammaJel can be >1 to allow
              ! for an outside food source.  However, GG's code doesn't
              ! spefify how the flux to detritus might change in that
              ! case (as written currently, that extra would come out of
              ! the DetF biomass via a negative egestion flux) TODO

              Ege_Jel_DetF(i,k) = (1.0_r8 - gammaJel) * cff1 * cff2 * cff3

            END DO
          END DO

          ! Convert all grazing and egestion fluxes from volumetric to
          ! integrated over layer

          DO k=1,N(ng)
            DO i=Istr,Iend
              Gra_PhS_MZL  (i,k)  = Gra_PhS_MZL  (i,k) * Hz(i,j,k)
              Gra_PhL_MZL  (i,k)  = Gra_PhL_MZL  (i,k) * Hz(i,j,k)
              Ege_MZL_Det  (i,k)  = Ege_MZL_Det  (i,k) * Hz(i,j,k)
              Gra_PhS_Cop  (i,k)  = Gra_PhS_Cop  (i,k) * Hz(i,j,k)
              Gra_PhL_Cop  (i,k)  = Gra_PhL_Cop  (i,k) * Hz(i,j,k)
              Gra_MZL_Cop  (i,k)  = Gra_MZL_Cop  (i,k) * Hz(i,j,k)
              Gra_IPhL_Cop (i,k)  = Gra_IPhL_Cop (i,k) * Hz(i,j,k)
              Ege_Cop_DetF (i,k)  = Ege_Cop_DetF (i,k) * Hz(i,j,k)
              Gra_PhS_NCaS (i,k)  = Gra_PhS_NCaS (i,k) * Hz(i,j,k)
              Gra_PhL_NCaS (i,k)  = Gra_PhL_NCaS (i,k) * Hz(i,j,k)
              Gra_MZL_NCaS (i,k)  = Gra_MZL_NCaS (i,k) * Hz(i,j,k)
              Gra_IPhL_NCaS(i,k)  = Gra_IPhL_NCaS(i,k) * Hz(i,j,k)
              Ege_NCaS_DetF(i,k)  = Ege_NCaS_DetF(i,k) * Hz(i,j,k)
              Gra_PhS_NCaO (i,k)  = Gra_PhS_NCaO (i,k) * Hz(i,j,k)
              Gra_PhL_NCaO (i,k)  = Gra_PhL_NCaO (i,k) * Hz(i,j,k)
              Gra_MZL_NCaO (i,k)  = Gra_MZL_NCaO (i,k) * Hz(i,j,k)
              Gra_IPhL_NCaO(i,k)  = Gra_IPhL_NCaO(i,k) * Hz(i,j,k)
              Ege_NCaO_DetF(i,k)  = Ege_NCaO_DetF(i,k) * Hz(i,j,k)
              Gra_PhS_EupS (i,k)  = Gra_PhS_EupS (i,k) * Hz(i,j,k)
              Gra_PhL_EupS (i,k)  = Gra_PhL_EupS (i,k) * Hz(i,j,k)
              Gra_MZL_EupS (i,k)  = Gra_MZL_EupS (i,k) * Hz(i,j,k)
              Gra_Cop_EupS (i,k)  = Gra_Cop_EupS (i,k) * Hz(i,j,k)
              Gra_IPhL_EupS(i,k)  = Gra_IPhL_EupS(i,k) * Hz(i,j,k)
              Gra_Det_EupS (i,k)  = Gra_Det_EupS (i,k) * Hz(i,j,k)
              Gra_DetF_EupS(i,k)  = Gra_DetF_EupS(i,k) * Hz(i,j,k)
              Ege_EupS_DetF(i,k)  = Ege_EupS_DetF(i,k) * Hz(i,j,k)
              Gra_PhS_EupO (i,k)  = Gra_PhS_EupO (i,k) * Hz(i,j,k)
              Gra_PhL_EupO (i,k)  = Gra_PhL_EupO (i,k) * Hz(i,j,k)
              Gra_MZL_EupO (i,k)  = Gra_MZL_EupO (i,k) * Hz(i,j,k)
              Gra_Cop_EupO (i,k)  = Gra_Cop_EupO (i,k) * Hz(i,j,k)
              Gra_IPhL_EupO(i,k)  = Gra_IPhL_EupO(i,k) * Hz(i,j,k)
              Gra_Det_EupO (i,k)  = Gra_Det_EupO (i,k) * Hz(i,j,k)
              Gra_DetF_EupO(i,k)  = Gra_DetF_EupO(i,k) * Hz(i,j,k)
              Ege_EupO_DetF(i,k)  = Ege_EupO_DetF(i,k) * Hz(i,j,k)
              Gra_Cop_Jel  (i,k)  = Gra_Cop_Jel  (i,k) * Hz(i,j,k)
              Gra_EupS_Jel (i,k)  = Gra_EupS_Jel (i,k) * Hz(i,j,k)
              Gra_EupO_Jel (i,k)  = Gra_EupO_Jel (i,k) * Hz(i,j,k)
              Gra_NCaS_Jel (i,k)  = Gra_NCaS_Jel (i,k) * Hz(i,j,k)
              Gra_NCaO_Jel (i,k)  = Gra_NCaO_Jel (i,k) * Hz(i,j,k)
              Ege_Jel_DetF (i,k)  = Ege_Jel_DetF (i,k) * Hz(i,j,k)
            END DO
          END DO

          !------------------------------
          ! Mortality and senescence
          !------------------------------

          ! TODO: might not need the loops here, but the GF%zoop_force
          ! dimensions complicate things so I'm keeping it

          ! TODO: change so exponent is set by user?  Would simplify the
          ! mXXX (linear) vs mpredXXX (quadratic) coefficient choice, but
          ! then the user would need to be careful that they changed the
          ! parameter and exponent in tandem

          DO k=1,N(ng)
            DO i=Istr,Iend

              ! Phytoplankton (linear senescence)

              Mor_PhS_Det(i,k) = mPhS * Bio3d(i,k,iiPhS)
              Mor_PhL_Det(i,k) = mPhL * Bio3d(i,k,iiPhL)

              ! Microzooplankton (quadratic mortality, with hard-coded
              ! option for linear)

!             Mor_MZL_Det(i,k) = mMZL*Bio3d(i,k,iiMZL)          ! linear
              Mor_MZL_Det(i,k) = mpredMZL*Bio3d(i,k,iiMZL)**2   ! quadratic

#ifdef fixedPRED
              ! TODO: original DBio(i,k,iXXX) = DBio(i,k,iXXX) - 0.5*Hz(i,j,k)/dtdays
              ! Implies coefficient units of mg C * day * m^-4???  Typo?
              ! Supposed to be constant rate, or maybe constant fraction
              ! of biomass?  Assuming the former for now.
              Mor_Cop_DetF(i,k)  = 0.5
              Mor_NCaS_DetF(i,k) = 0.5
              Mor_EupS_DetF(i,k) = 1
              Mor_NCaO_DetF(i,k) = 0.5
              Mor_EupO_DetF(i,k) = 1
#else
              TFEup = Q10Eup ** ((Temp(i,k)-Q10EupT) / 10.0_r8)
# ifdef FEAST
              ! Mesozooplankton (based on predation by fish)

              Mor_Cop_DetF(i,k)  = TFEup*(mpredCop + fpredCop  * GF%zoop_force(1,1,i,j,1))*Bio3d(i,k,iiCop)**2
              Mor_NCaS_DetF(i,k) = TFEup*(mpredNca + fpredNcaS * GF%zoop_force(1,2,i,j,1))*Bio3d(i,k,iiNCaS)**2
              Mor_EupS_DetF(i,k) = TFEup*(mpredEup + fpredEupS * GF%zoop_force(1,4,i,j,1))*Bio3d(i,k,iiEupS)**2
              Mor_NCaO_DetF(i,k) = TFEup*(mpredNca + fpredNcaO * GF%zoop_force(1,3,i,j,1))*Bio3d(i,k,iiNCaO)**2
              Mor_EupO_DetF(i,k) = TFEup*(mpredEup + fpredEupO * GF%zoop_force(1,5,i,j,1))*Bio3d(i,k,iiEupO)**2
# else
              ! Mesozooplankton (quadratic predation closure)
              Mor_Cop_DetF(i,k)  = TFEup*(mpredCop)*Bio3d(i,k,iiCop)**2
              Mor_NCaS_DetF(i,k) = TFEup*(mpredNca)*Bio3d(i,k,iiNCaS)**2
              Mor_EupS_DetF(i,k) = TFEup*(mpredEup)*Bio3d(i,k,iiEupS)**2
              Mor_NCaO_DetF(i,k) = TFEup*(mpredNca)*Bio3d(i,k,iiNCaO)**2
              Mor_EupO_DetF(i,k) = TFEup*(mpredEup)*Bio3d(i,k,iiEupO)**2
# endif
#endif

              ! Jellyfish (quadratic predation closure)

              Mor_Jel_DetF(i,k) = mpredJel*Bio3d(i,k,iiJel)**2

            END DO
          END DO

          ! Convert mortality fluxes from volumetric to integrated over
          ! layer

          DO k=1,N(ng)
            DO i=Istr,Iend
              Mor_PhS_Det(i,k)   = Mor_PhS_Det(i,k)   * Hz(i,j,k)
              Mor_PhL_Det(i,k)   = Mor_PhL_Det(i,k)   * Hz(i,j,k)
              Mor_MZL_Det(i,k)   = Mor_MZL_Det(i,k)   * Hz(i,j,k)
              Mor_Cop_DetF(i,k)  = Mor_Cop_DetF(i,k)  * Hz(i,j,k)
              Mor_NCaS_DetF(i,k) = Mor_NCaS_DetF(i,k) * Hz(i,j,k)
              Mor_EupS_DetF(i,k) = Mor_EupS_DetF(i,k) * Hz(i,j,k)
              Mor_NCaO_DetF(i,k) = Mor_NCaO_DetF(i,k) * Hz(i,j,k)
              Mor_EupO_DetF(i,k) = Mor_EupO_DetF(i,k) * Hz(i,j,k)
              Mor_Jel_DetF(i,k)  = Mor_Jel_DetF(i,k)  * Hz(i,j,k)
            END DO
          END DO

          !------------------------------
          ! Respiration
          !------------------------------

          ! Phytoplankton

          Res_PhS_NH4 = exp(KtBm_PhS * (Temp - TmaxPhS)) * respPhS * Bio3d(:,:,iiPhS)
          Res_PhL_NH4 = exp(KtBm_PhL * (Temp - TmaxPhL)) * respPhL * Bio3d(:,:,iiPhL)

          ! Microzooplankton

          Res_MZL_NH4 = exp(KtBm_MZL * (Temp - TmaxMZL)) * respMZL * Bio3d(:,:,iiMZL)

          ! Mesozooplankton (BasMetXXX is respXXX w/ starvation response)

          Res_Cop_NH4  = exp(ktbmC * (Temp - TrefC)) * BasMetCop * Bio3d(:,:,iiCop)
          Res_NCaS_NH4 = exp(ktbmN * (Temp - TrefN)) * BasMetCM  * Bio3d(:,:,iiNCaS)
          Res_NCaO_NH4 = exp(ktbmN * (Temp - TrefN)) * BasMetNC  * Bio3d(:,:,iiNCaO)
          Res_EupS_NH4 = exp(ktbmE * (Temp - TrefE)) * BasMetEup * Bio3d(:,:,iiEupS)
          Res_EupO_NH4 = exp(ktbmE * (Temp - TrefE)) * BasMetEup * Bio3d(:,:,iiEupO)

          ! Jellyfish

          Res_Jel_NH4 = Q10Jelr ** ((Temp-Q10JelTr)/10.0_r8) * respJel * Bio3d(:,:,iiJel)

          ! Convert respiration fluxes from volumetric to integrated over
          ! layer

          DO k=1,N(ng)
            DO i=Istr,Iend
              Res_PhS_NH4(i,k)  = Res_PhS_NH4(i,k)  * Hz(i,j,k)
              Res_PhL_NH4(i,k)  = Res_PhL_NH4(i,k)  * Hz(i,j,k)
              Res_MZL_NH4(i,k)  = Res_MZL_NH4(i,k)  * Hz(i,j,k)
              Res_Cop_NH4(i,k)  = Res_Cop_NH4(i,k)  * Hz(i,j,k)
              Res_NCaS_NH4(i,k) = Res_NCaS_NH4(i,k) * Hz(i,j,k)
              Res_NCaO_NH4(i,k) = Res_NCaO_NH4(i,k) * Hz(i,j,k)
              Res_EupS_NH4(i,k) = Res_EupS_NH4(i,k) * Hz(i,j,k)
              Res_EupO_NH4(i,k) = Res_EupO_NH4(i,k) * Hz(i,j,k)
              Res_Jel_NH4(i,k)  = Res_Jel_NH4(i,k)  * Hz(i,j,k)
            END DO
          END DO

          !------------------------------
          ! Decomposition, nitrification,
          ! and remineralization
          !------------------------------

          DO k=1,N(ng)
            DO i=Istr,Iend

              ! Detrital remineralization

              PON = Bio3d(i,k,iiDet)*xi  ! Particulate organic nitrogen in Det
              Dec_Det_NH4(i,k) = (Pv0 * exp(PvT*Temp(i,k)) * PON) ! mmol N m^-3 d^-1

              PON = Bio3d(i,k,iiDetF)*xi  ! Particulate organic nitrogen in DetF
              Dec_DetF_NH4(i,k) = (Pv0 * exp(PvT*Temp(i,k)) * PON) ! mmol N m^-3 d^-1

              ! Nitrification

              ParW = PAR(i,k)/0.394848_r8 ! convert to W TODO: is this supposed to be PAR(i,k) or PARs(i)?
              NitrifMax = Nitr0 * exp(-ktntr*(Temp(i,k) - ToptNtr)**2)     ! Arhonditsis 2005 temperature dependence
              DLNitrif = (1 - MAX(0.0_r8, (ParW - tI0)/(KI + ParW - tI0))) ! Fennel light dependence
              cff1 = Bio3d(i,k,iiNH4)/(KNH4Nit +Bio3d(i,k,iiNH4))          ! Arhonditsis saturation

              Dec_NH4_NO3(i,k) = NitrifMax * Bio3d(i,k,iiNH4) * DLNitrif * cff1 !  mmol N m^-3 d^-1

            END DO
          END DO

          ! Convert decomp fluxes from volumetric to integrated over
          ! layer, and from N to C for consistency with other fluxes

          DO k=1,N(ng)
            DO i=Istr,Iend
              Dec_Det_NH4(i,k)  = Dec_Det_NH4(i,k)  * Hz(i,j,k)/xi
              Dec_DetF_NH4(i,k) = Dec_DetF_NH4(i,k) * Hz(i,j,k)/xi
              Dec_NH4_NO3(i,k)  = Dec_NH4_NO3(i,k)  * Hz(i,j,k)/xi
            END DO
          END DO


#ifdef BENTHIC
          !-----------------
          !Benthic Sub Model
          !-----------------

          DO i=Istr,Iend

            ! Pelagic food accessible to benthic infauna

            dw = 1.0_r8 ! assume bottom 1 m is accessible

            totD  = 0.0_r8
            totDF = 0.0_r8
            totPS = 0.0_r8
            totPL = 0.0_r8

            cff2 = 0.0_r8 ! accounted-for height above bottom
            DO k = 1,N(ng)

              ! Fraction of this layer contributing to benthic feeding

              cff1 = max(min(dw, cff2 + Hz(i,j,k)) - cff2, 0.0_r8) ! m
              frac1(i,k) = cff1/Hz(i,j,k)

              ! Fraction of benthic feeding coming from this level

              frac2(i,k) = cff1/dw

              ! Food available to benthos

              totD  = totD  + Bio2d(i,k,iiDet) *frac1(i,k)
              totDF = totDF + Bio2d(i,k,iiDetF)*frac1(i,k)
              totPS = totPS + Bio2d(i,k,iiPhS) *frac1(i,k)
              totPL = totPL + Bio2d(i,k,iiPhL) *frac1(i,k)

              cff2 = cff2 + Hz(i,j,k)
            END DO

            ! Potential food available from water column

            cff1=(prefD *totD /((prefD *totD )+LupP))*prefD *totD
            cff2=(prefD *totDF/((prefD *totDF)+LupP))*prefD *totDF
            cff3=(prefPS*totPS/((prefPS*totPS)+LupP))*prefPS*totPS
            cff4=(prefPL*totPL/((prefPL*totPL)+LupP))*prefPL*totPL

            cff6 = cff1+cff2+cff3+cff4 ! Total pelagic food

            ! Potential food available from  sea floor

            cff5 = (prefD * Bio2d(i,1,iiBenDet) /                       &
     &             (prefD * Bio2d(i,1,iiBenDet) + LupD)) *              &
                    prefD * Bio2d(i,1,iiBenDet)

            ! Temperature mediation (for feeding and mortality)

            cff0 = q10r**((Temp(i,1)-T0benr)/10.0_r8)

            ! Total uptake of each food category

            cff7  = min(cff1,(cff0*cff1*BioB(i,k,iBen)*Rup/(cff6+KupP))) ! D
            cff8  = min(cff2,(cff0*cff2*BioB(i,k,iBen)*Rup/(cff6+KupP))) ! DF
            cff9  = min(cff3,(cff0*cff3*BioB(i,k,iBen)*Rup/(cff6+KupP))) ! PS
            cff10 = min(cff4,(cff0*cff4*BioB(i,k,iBen)*Rup/(cff6+KupP))) ! PL
            cff11 = min(cff5,(cff0*cff5*BioB(i,k,iBen)*Rup/(cff5+KupD))) ! BenDet

            ! Distribute pelagic feeding losses to appropriate water
            ! column layers
            ! TODO: will need to account for multi-layers-to-bottom-layer
            ! for input to Ben when calculating total Ben DBio

            DO k = 1,N(ng)

              Gra_Det_Ben(i,k)  = cff7  * frac2(i,k) ! mg C m^-2 d^-1
              Gra_DetF_Ben(i,k) = cff8  * frac2(i,k)
              Gra_PhS_Ben(i,k)  = cff9  * frac2(i,k)
              Gra_PhL_Ben(i,k)  = cff10 * frac2(i,k)

            END DO

            ! Benthic feeding takes place in bottom layer for bookkeeping
            ! purposes

            Gra_BenDet_Ben(i,1) = cff11 ! mg C m^-2 d^-1

            ! Assume all excretion occurs in the bottom layer too.  Half
            ! goes to NH4 and half to BenDet

            Exc_Ben_BenDet(i,1) = (eexD * (cff7 + cff8 + cff11) +       &
     &                             eex  * (cff9 + cff10)) * 0.5_r8
            Exc_Ben_NH4(i,1) = Exc_Ben_BenDet(i,1)

            ! Respiration (also takes place in bottom layer)

            cff3 = cff0 * Bio2d(i,1,iiBen) * Rres
            cff4 = ((1_r8 - eexD) * (cff7 + cff8 + cff11) +             &
     &              (1_r8 - eex)  * (cff9 + cff10)) * Qres

            Res_Ben_NH4(i,1) = cff3 + cff4 ! mg C m^-2 d^-1

            ! Mortality

            Mor_Ben_BenDet(i,1) = rmort*Bio2d(i,1,iiBen)*cff0 ! mg C m^-2 d^-1

            ! Additional predation TODO: is it right that this goes immediately to BenDet?
            ! This additional loss is due to undefined predation on
            ! benthic infauna, which is assumed to make its way back to
            ! the benthic detritus pool

            Pre_Ben_BenDet(i,1) = Bio2d(i,1,iiBen)**2 * cff0 * BenPred ! mg C m^-2 d^-1

            ! Benthic remineralization: assumes only the top 25% is
            ! available to remineralize to NH4 (in bottom layer)

            PON = Bio3d(i,k,iiBenDet)*0.25*xi  ! Benthic Particulate organic nitrogen
            cff1 = Pv0*exp(PvT*Temp(i,1))*PON  ! Kawamiya 2000, mmol N m^-3

            Dec_BenDet_NH4(i,1) = cff1*Hz(i,j,k)/xi ! mg C m^-2

          END DO
#endif

#ifdef ICE_BIO
          !-----------------
          ! Ice Sub Model
          !-----------------

          DO i=Istr,Iend
            if (ice_status(i,j) .gt. 0.0_r8) then

              ! Ice algae production limitation terms

              Temp1 = Temp(i,N(ng)) ! Assume temperature of top layer = temp ice skeletal layer
              Par1  = PARs(i)       ! surface light

              aiceIfrac = (1-exp(-alphaIb*Par1))*exp(-betaI*Par1) ! light limitation

              cff1 = Bio3d(i,N(ng),iiIceNO3)/(ksnut1 + Bio3d(i,N(ng),iiIceNO3)) ! NO3 limitation
              cff2 = Bio3d(i,N(ng),iiIceNH4)/(ksnut2 + Bio3d(i,N(ng),iiIceNH4)) ! NH4 limitation
              aiceNfrac = cff1*exp(-inhib*Bio3d(i,N(ng),iiIceNH4)) + cff2       ! N limitation
              fNO3      = cff1*exp(-inhib*Bio3d(i,N(ng),iiIceNH4))/aiceNfrac    ! f-ratio

# ifdef BERING_10K
              ! Ice algae growth is also limited by suboptimal brine
              ! salinity in the ice.  This value isn't tracked explicitly
              ! by the ice model, so instead we use the brine salinty vs
              ! ice temperature polynomial fit from Arrigo 1993 Appendix
              ! A to estimate brine salinity.

              if (ti(i,j,nstp) .ge. -22.9_r8) THEN
                cff1=-3.9921
                cff2=-22.7
                cff3=-1.0015
                cff4=-0.019956
              else if (ti(i,j,nstp) .gt. -44.0_r8  .AND. ti(i,j,nstp) .lt. -22.9_r8) THEN
                cff1=206.24
                cff2=-1.8907
                cff3=-0.060868
                cff4=-0.0010247
              else
                cff1=-4442.1
                cff2=-277.86
                cff3=-5.501
                cff4=-0.03669
              endif

              sb = cff1 + cff2*ti(i,j,nstp) + cff3*ti(i,j,nstp)**2 +    &
     &             cff4*ti(i,j,nstp)**3 ! brine salinity

              ! Salinity impact on ice algal growth (gesi) is determined
              ! by a polynomial fit to brine salinty (Arrigo 1993
              ! Appendix B)

              gesi = max(0.0_r8, (1.1e-2+3.012e-2*sb                    &
     &             +1.0342e-3*sb**2                                     &
     &             -4.6033e-5*sb**3                                     &
     &             +4.926e-7*sb**4                                      &
     &             -1.659e-9*sb**5              ))

# else
              ! When running without an accompanying ice model, assume no
              ! salinity limitation on ice algae growth

              gesi=1.0_r8
# endif

              ! Ice algae production

              grow1 = mu0*exp(0.0633*Temp1)
              GROWAice=grow1 * min(aiceNfrac,aiceIfrac) * gesi

              Gpp_INO3_IPhL(i,N(ng)) = (GrowAice * Bio3d(i,N(ng),iiIcePhL) *    fNO3 )*aidz ! mg C m^-2 d^-1
              Gpp_INH4_IPhL(i,N(ng)) = (GrowAice * Bio3d(i,N(ng),iiIcePhL) * (1-fNO3))*aidz ! mg C m^-2 d^-1

              ! Ice algae respiration

              RAi0 = R0i*mu0*exp(0.0633*Temp1)

              Res_IPhL_INH4(i,N(ng)) = (RAi0 * Bio3d(i,N(ng),iiIcePhL))*aidz ! mg C m^-2 d^-1

              ! Ice algae mortality

              RgAi = rg0*exp(rg*Temp1)

              Mor_IPhL_INH4(i,N(ng)) = (Bio3d(i,N(ng),iiIcePhL)*RgAi)*aidz ! mg C m^-2 d^-1

              ! Nitrification

              Dec_INH4_INO3(i,N(ng)) = (annit*Bio3d(i,N(ng),iiIceNH4)/xi)*aidz ! mg C m^-2 d^-1

              ! Ice/water convective exchange covers transfer of algae,
              ! NO3, and NH4 between the ice and surface water based on
              ! water exchange between the two layers, following Jin et
              ! al. 2006.  The water-ice interface transport (twi) rate
              ! is determined based on a polynomial fit with rate of
              ! change of ice thickness.

# if defined CLIM_ICE_1D
              dhicedt=it(i,j,nnew,iIceZ)-it(i,j,nstp,iIceZ) ! change in ice thickness over this time step (m)
# elif defined BERING_10K
              dhicedt=hi(i,j,nstp)-hi(i,j,nnew) ! change in ice thickness over this time step (m)
# endif
              dhicedt=dhicedt*sec2day/dtdays ! convert to m/s for polynomial fit

              IF (dhicedt.lt.0) THEN ! Ice is melting

                trs=4.49e-6*ABS(dhicedt)-1.39e-5*ABS(dhicedt)**2
                trs=trs*86400   ! convert back to m/d
                twi=720*trs

              ELSE                   ! Ice is growing
                trs=9.667e-11+4.49e-6*dhicedt-1.39e-5*dhicedt**2
                trs=trs*86400   ! convert back to m/d
                twi=72*trs
              ENDIF

              ! IcePhL can get washed out of ice, but not in, so assume
              ! [PhL] = 0 for this exchange

              Twi_IPhL_PhL(i,N(ng)) = twi * Bio3d(i,N(ng),iiIcePhL) ! mg C m^-2 d^-1

              ! NO3 and NH4 have two-way exchange, based on gradient
              ! across the ice/water interface.  Note that I'm using
              ! ice-to-water as the naming convention for this flux, but
              ! it may be negative, implying the reverse direction.

              Twi_INO3_NO3(i,N(ng)) = twi * (Bio3d(i,N(ng),iiIceNO3) - Bio3d(i,N(ng),iiNO3))/xi ! mg C m^-2 d^-1
              Twi_INH4_NH4(i,N(ng)) = twi * (Bio3d(i,N(ng),iiIceNH4) - Bio3d(i,N(ng),iiNH4))/xi ! mg C m^-2 d^-1

            endif
          END DO
#endif

          !------------------------------
          ! Combine bio source/sinks
          !------------------------------

          DBio(:,:,iiNO3   ) = (Dec_NH4_NO3                             &
     &                       +  Twi_INO3_NO3                            &
     &                       -  Gpp_NO3_PhS                             &
     &                       -  Gpp_NO3_PhL)*xi*dtdays ! NO3: mmolN m^-2

          DBio(:,:,iiNH4   ) = (Res_PhS_NH4                             &
     &                       +  Res_PhL_NH4                             &
     &                       +  Res_MZL_NH4                             &
     &                       +  Res_Cop_NH4                             &
     &                       +  Res_NCaS_NH4                            &
     &                       +  Res_NCaO_NH4                            &
     &                       +  Res_EupS_NH4                            &
     &                       +  Res_EupO_NH4                            &
     &                       +  Res_Jel_NH4                             &
     &                       +  Dec_Det_NH4                             &
     &                       +  Dec_DetF_NH4                            &
     &                       +  Exc_Ben_NH4                             &
     &                       +  Res_Ben_NH4                             &
     &                       +  Dec_BenDet_NH4                          &
     &                       +  Twi_INH4_NH4                            &
     &                       -  Gpp_NH4_PhS                             &
     &                       -  Gpp_NH4_PhL                             &
     &                       -  Dec_NH4_NO3)*xi*dtdays ! NH4: mmol N m^-2

          DBio(:,:,iiPhS   ) = (Gpp_NO3_PhS                             &
     &                       +  Gpp_NH4_PhS                             &
     &                       -  Gra_PhS_MZL                             &
     &                       -  Gra_PhS_Cop                             &
     &                       -  Gra_PhS_NCaS                            &
     &                       -  Gra_PhS_NCaO                            &
     &                       -  Gra_PhS_EupS                            &
     &                       -  Gra_PhS_EupO                            &
     &                       -  Mor_PhS_Det                             &
     &                       -  Res_PhS_NH4                             &
     &                       -  Gra_PhS_Ben)*dtdays ! PhS: mg C m^-2

          DBio(:,:,iiPhL   ) = (Gpp_NO3_PhL                             &
     &                       +  Gpp_NH4_PhL                             &
     &                       +  Twi_IPhL_PhL                            &
     &                       -  Gra_PhL_MZL                             &
     &                       -  Gra_PhL_Cop                             &
     &                       -  Gra_PhL_NCaS                            &
     &                       -  Gra_PhL_NCaO                            &
     &                       -  Gra_PhL_EupS                            &
     &                       -  Gra_PhL_EupO                            &
     &                       -  Mor_PhL_Det                             &
     &                       -  Res_PhL_NH4                             &
     &                       -  Gra_PhL_Ben)*dtdays ! PhL: mg C m^-2

          DBio(:,:,iiMZL   ) = (Gra_PhS_MZL                             &
     &                       +  Gra_PhL_MZL                             &
     &                       -  Ege_MZL_Det                             &
     &                       -  Gra_MZL_Cop                             &
     &                       -  Gra_MZL_NCaS                            &
     &                       -  Gra_MZL_NCaO                            &
     &                       -  Gra_MZL_EupS                            &
     &                       -  Gra_MZL_EupO                            &
     &                       -  Mor_MZL_Det                             &
     &                       -  Res_MZL_NH4)*dtdays ! MZL: mg C m^-2

          DBio(:,:,iiCop   ) = (Gra_PhS_Cop                             &
     &                       +  Gra_PhL_Cop                             &
     &                       +  Gra_MZL_Cop                             &
     &                       +  Gra_IPhL_Cop                            &
     &                       -  Ege_Cop_DetF                            &
     &                       -  Gra_Cop_EupS                            &
     &                       -  Gra_Cop_EupO                            &
     &                       -  Gra_Cop_Jel                             &
     &                       -  Mor_Cop_DetF                            &
     &                       -  Res_Cop_NH4)*dtdays ! Cop: mg C m^-2

          DBio(:,:,iiNCaS  ) = (Gra_PhS_NCaS                            &
     &                       +  Gra_PhL_NCaS                            &
     &                       +  Gra_MZL_NCaS                            &
     &                       +  Gra_IPhL_NCaS                           &
     &                       -  Ege_NCaS_DetF                           &
     &                       -  Gra_NCaS_Jel                            &
     &                       -  Mor_NCaS_DetF                           &
     &                       -  Res_NCaS_NH4)*dtdays ! NCaS: mg C m^-2

          DBio(:,:,iiEupS  ) = (Gra_PhS_EupS                            &
     &                       +  Gra_PhL_EupS                            &
     &                       +  Gra_MZL_EupS                            &
     &                       +  Gra_Cop_EupS                            &
     &                       +  Gra_IPhL_EupS                           &
     &                       +  Gra_Det_EupS                            &
     &                       +  Gra_DetF_EupS                           &
     &                       -  Ege_EupS_DetF                           &
     &                       -  Gra_EupS_Jel                            &
     &                       -  Mor_EupS_DetF                           &
     &                       -  Res_EupS_NH4)*dtdays ! EupS: mg C m^-2

          DBio(:,:,iiNCaO  ) = (Gra_PhS_NCaO                            &
     &                       +  Gra_PhL_NCaO                            &
     &                       +  Gra_MZL_NCaO                            &
     &                       +  Gra_IPhL_NCaO                           &
     &                       -  Ege_NCaO_DetF                           &
     &                       -  Gra_NCaO_Jel                            &
     &                       -  Mor_NCaO_DetF                           &
     &                       -  Res_NCaO_NH4)*dtdays ! NCaO: mg C m^-2

          DBio(:,:,iiEupO  ) = (Gra_PhS_EupO                            &
     &                        +  Gra_PhL_EupO                           &
     &                        +  Gra_MZL_EupO                           &
     &                        +  Gra_Cop_EupO                           &
     &                        +  Gra_IPhL_EupO                          &
     &                        +  Gra_Det_EupO                           &
     &                        +  Gra_DetF_EupO                          &
     &                        -  Ege_EupO_DetF                          &
     &                        -  Gra_EupO_Jel                           &
     &                        -  Mor_EupO_DetF                          &
     &                        -  Res_EupO_NH4)*dtdays ! EupO: mg C m^-2

          DBio(:,:,iiDet   )  = (Ege_MZL_Det                            &
     &                        +  Mor_PhS_Det                            &
     &                        +  Mor_PhL_Det                            &
     &                        +  Mor_MZL_Det                            &
     &                        -  Gra_Det_EupS                           &
     &                        -  Gra_Det_EupO                           &
     &                        -  Dec_Det_NH4                            &
     &                        -  Gra_Det_Ben)*dtdays ! Det: mg C m^-2

          DBio(:,:,iiDetF  )  = (Ege_Cop_DetF                           &
     &                        +  Ege_NCaS_DetF                          &
     &                        +  Ege_NCaO_DetF                          &
     &                        +  Ege_EupS_DetF                          &
     &                        +  Ege_EupO_DetF                          &
     &                        +  Ege_Jel_DetF                           &
     &                        +  Mor_Cop_DetF                           &
     &                        +  Mor_NCaS_DetF                          &
     &                        +  Mor_EupS_DetF                          &
     &                        +  Mor_NCaO_DetF                          &
     &                        +  Mor_EupO_DetF                          &
     &                        +  Mor_Jel_DetF                           &
     &                        -  Gra_DetF_EupS                          &
     &                        -  Gra_DetF_EupO                          &
     &                        -  Dec_DetF_NH4                           &
     &                        -  Gra_DetF_Ben)*dtdays ! DetF: mg C m^-2

          DBio(:,:,iiJel   )  = (Gra_Cop_Jel                            &
     &                        +  Gra_EupS_Jel                           &
     &                        +  Gra_EupO_Jel                           &
     &                        +  Gra_NCaS_Jel                           &
     &                        +  Gra_NCaO_Jel                           &
     &                        -  Ege_Jel_DetF                           &
     &                        -  Mor_Jel_DetF                           &
     &                        -  Res_Jel_NH4)*dtdays ! Jel: mg C m^-2

          DBio(:,:,iiFe    )  = (                                       &
     &                        -  Gpp_NO3_PhS                            &
     &                        -  Gpp_NO3_PhL)*FeC*dtdays ! Fe: umol Fe m^-2

          DBio(:,:,iiBen   )  = (Gra_Det_Ben                            &
     &                        +  Gra_DetF_Ben                           &
     &                        +  Gra_PhS_Ben                            &
     &                        +  Gra_PhL_Ben                            &
     &                        +  Gra_BenDet_Ben                         &
     &                        -  Exc_Ben_NH4                            &
     &                        -  Exc_Ben_BenDet                         &
     &                        -  Res_Ben_NH4                            &
     &                        -  Mor_Ben_BenDet                         &
     &                        -  Pre_Ben_BenDet)*dtdays ! Ben: mg C m^-2

          DBio(:,:,iiBenDet)  = (Exc_Ben_BenDet                         &
     &                        +  Mor_Ben_BenDet                         &
     &                        +  Pre_Ben_BenDet                         &
     &                        -  Gra_BenDet_Ben                         &
     &                        -  Dec_BenDet_NH4)*dtdays ! BenDet: mg C m^-2

          DBio(:,:,iiIcePhL)  = (Gpp_INO3_IPhL                          &
     &                        +  Gpp_INH4_IPhL                          &
     &                        -  Gra_IPhL_Cop                           &
     &                        -  Gra_IPhL_NCaS                          &
     &                        -  Gra_IPhL_NCaO                          &
     &                        -  Gra_IPhL_EupS                          &
     &                        -  Gra_IPhL_EupO                          &
     &                        -  Res_IPhL_INH4                          &
     &                        -  Mor_IPhL_INH4                          &
     &                        -  Twi_IPhL_PhL)*dtdays ! IcePhL: mg C m^-2

          DBio(:,:,iiIceNO3)  = (Dec_INH4_INO3                          &
     &                        -  Gpp_INO3_IPhL                          &
     &                        -  Twi_INO3_NO3)*xi*dtdays ! IceNO3: mmol N m^-2

          DBio(:,:,iiIceNH4)  = (Res_IPhL_INH4                          &
     &                        +  Mor_IPhL_INH4                          &
     &                        -  Gpp_INH4_IPhL                          &
     &                        -  Dec_INH4_INO3                          &
     &                        -  Twi_INH4_NH4)*xi*dtdays ! IceNH4: mmol N m^-2


          ! TODO: Collect net production, bflux values

          ! Add DBio terms to existing biomass

          Bio2d = Bio2d + DBio

          ! Infauna (Ben) group can receive flux from water column layers.
          ! Move these additions to the bottom layer now, consistent with
          ! the initial setup of the Bio2d and Bio3d arrays.

          Bio2d(:,1,iiBen) = sum(Bio2d(:,:,iiBen), DIM=2)
          Bio2d(:,2:N(ng),iiBen) = 0

          ! TODO: Eliminate negatives?  Hopefully processes are
          ! formulated to prevent any, but very fast overturning might
          ! result in numerical issues.  Brute force zero traps will
          ! eliminate conservation of mass, so I'd prefer to look into
          ! increasing BioIter if this is a problem


          ! Sync volumetric version to the updated per-area values

          DO i=Istr,Iend
            DO k = 1,N(ng)
              DO itrc = 1,NBT+NBEN ! Pelagic (and benthic, for bookkeeping)
                Bio3d(i,k,itrc) = Bio2d(i,k,itrc)/Hz(i,j,k)
              END DO
              DO itrc = 18,20 ! Ice
                Bio3d(i,k,itrc) = Bio2d(i,k,itrc)/aidz
              END DO
            END DO
          END DO

          !==============================================================
          ! Vertical Movement
          !==============================================================

          ! This section includes all vertical movement of state
          ! variables, including sinking of phytoplankton and particulate
          ! detritus, large copepod seasonal diapause, and (evenuntually)
          ! euphausiid diel vertical migration.


          ! Initialize temporary arrays to 0

          dBtmp = 0
          flxtmp = 0

          ! Small phytoplankton: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wPhS, Bio3d(i,:,iiPhS), dBtmp,         &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiPhS) = Bio3d(i,1:N(ng),iiPhS) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiBenDet) = Bio2d(i,1,iiBenDet) + flxtmp*0.79_r8

          END DO


          ! Large phytoplankton: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wPhL, Bio3d(i,:,iiPhL), dBtmp,         &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiPhL) = Bio3d(i,1:N(ng),iiPhL) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiBenDet) = Bio2d(i,1,iiBenDet) + flxtmp*0.79_r8

          END DO

          ! Slow-sinking detritus: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wDet, Bio3d(i,:,iiDet), dBtmp,         &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiDet) = Bio3d(i,1:N(ng),iiDet) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiBenDet) = Bio2d(i,1,iiBenDet) + flxtmp*0.79_r8

          END DO

          ! Fast-sinking detritus: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wDetF, Bio3d(i,:,iiDetF), dBtmp,       &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiDetF) = Bio3d(i,1:N(ng),iiDetF) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiBenDet) = Bio2d(i,1,iiBenDet) + flxtmp*0.79_r8

          END DO

          ! On-shelf large copepods (NCaS i.e. CM): Move up and down
          ! based on dates set in input file.  Downward movement is
          ! stopped at 200 m or halfway through the bottom layer,
          ! whichever is shallower; upward movement is stopped halfway
          ! through the top layer. No biomass should cross the bottom or
          ! surface boundary.

          DO i=Istr,Iend

            if (downwardCM) then

              call BioVert(N(ng), -wNCsink, Bio3d(i,:,iiNCaS), dBtmp,   &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     max((z_w(i,j,0)+z_w(i,j,1))/2, -200.0_r8), flxtmp)
              Bio3d(i,1:N(ng),iiNCaS) = Bio3d(i,1:N(ng),iiNCaS) + dBtmp(1,1:N(ng))

            else if (upwardCM) then

              call BioVert(N(ng), wNCrise, Bio3d(i,:,iiNCaS), dBtmp,    &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     (z_w(i,j,N(ng)-1)+z_w(i,j,N(ng)))/2, flxtmp)
              Bio3d(i,1:N(ng),iiNCaS) = Bio3d(i,1:N(ng),iiNCaS) + dBtmp(1,1:N(ng))

            end if

          END DO

          ! Off-shelf large copepods (NCaO i.e. NC): Move up and down
          ! based on dates set in input file.  Downward movement is
          ! stopped at 400 m or halfway through the bottom layer,
          ! whichever is shallower; upward movement is stopped halfway
          ! through the top layer.  If the water depth is less than 400 m,
          ! copepods are assumed to die when they hit the bottom, and
          ! biomass is transferred to benthic detritus.

          DO i=Istr,Iend

            if (downwardNC) then

              call BioVert(N(ng), -wNCsink, Bio3d(i,:,iiNCaO), dBtmp,   &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     -400.0_r8, flxtmp)
              Bio3d(i,1:N(ng),iiNCaO) = Bio3d(i,1:N(ng),iiNCaO) + dBtmp(1,1:N(ng))
              Bio2d(i,1,iBenDet) = Bio2d(i,1,iBenDet) + flxtmp

            else if (upwardNC) then

              call BioVert(N(ng), wNCrise, Bio3d(i,:,iiNCaO), dBtmp,    &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     (z_w(i,j,N(ng)-1)+z_w(i,j,N(ng)))/2, flxtmp)
              Bio3d(i,1:N(ng),iiNCaO) = Bio3d(i,1:N(ng),iiNCaO) + dBtmp(1,1:N(ng))

            end if

          END DO

          ! TODO: overhauled to here

          !==============================================================
          ! Update Bio array
          !==============================================================
          DO i=Istr,Iend

            DO itrc=1,NBT
              ibio=idbio(itrc)
              DO k=1,N(ng)

                Bio(i,k,ibio)=Bio(i,k,ibio)+DBio(i,k,ibio)

              END DO
            END DO
          END DO
#ifdef BENTHIC
          DO itrc=1,NBEN
            ibioB=idben(itrc)
            DO k=1,NBL(ng)
              DO i=Istr,Iend

                BioB(i,k,ibioB)=BioB(i,k,ibioB)+DBioB(i,k,ibioB)

              END DO
            END DO
          END DO

#endif
#ifdef ICE_BIO
# if defined CLIM_ICE_1D

          DO i=Istr,Iend

            if(itL(i,j,nstp,1).gt.0_r8)THEN
              BioBI(i,iIceNO3)=BioBI(i,iIceNO3)+DBioBI(i,iIceNO3)
              BioBI(i,iIceNH4)=BioBI(i,iIceNH4)+DBioBI(i,iIceNH4)
              BioBI(i,iIcePhL)=BioBI(i,iIcePhL)+DBioBI(i,iIcePhL)

            else
              BioBI(i,iIceNO3)=0_r8
              BioBI(i,iIceNH4)=0_r8
              BioBI(i,iIcePhL)=0_r8
            endif

          END DO

# elif defined BERING_10K

          DO i=Istr,Iend

            if (IceLog(i,j,nstp).gt.0 )THEN

              IcePhL(i,j,nnew) = IcePhL(i,j,nstp) + DBioBI(i,iIcePhL)
              IceNO3(i,j,nnew) = IceNO3(i,j,nstp) + DBioBI(i,iIceNO3)
              IceNH4(i,j,nnew) = IceNH4(i,j,nstp) + DBioBI(i,iIceNH4)
            else
              IcePhL(i,j,nnew) = 0_r8
              IceNO3(i,j,nnew) = 0_r8
              IceNH4(i,j,nnew) = 0_r8
            endif

#  ifdef STATIONARY2


!           Stat2(i,5)=IceNH4(i,j,nstp)
!           Stat2(i,6)=IceNH4(i,j,nnew)
#  endif


          END DO

# endif
#endif


        END DO ITER_LOOP

        !=============================================
        !  Update global tracer variables (m Tunits).
        !=============================================

        DO i=Istr,Iend
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)

              t(i,j,k,nnew,ibio)=MAX(t(i,j,k,nnew,ibio)+                &
      &                               (Bio(i,k,ibio)-Bio_bak(i,k,ibio)) &
      &                               *Hz(i,j,k)                        &
      &                               ,0.0_r8)

#ifdef TS_MPDATA
              t(i,j,k,3,ibio)=t(i,j,k,nnew,ibio)*Hz_inv(i,k)

#endif

              t(i,j,k,nnew,iMZS)=0.0_r8

            END DO
          END DO
        END DO

#ifdef BENTHIC

        DO itrc=1,NBEN
          ibioB=idben(itrc)
          DO k=1,NBL(ng)
            DO i=Istr,Iend

              !  check indexing here. think it ok.
              bt(i,j,k,nnew,ibioB)=MAX(bt(i,j,k,nstp,ibioB)+            &
     &                           (BioB(i,k,ibioB)-Bio_bakB(i,k,ibioB))  &
     &                           ,0.0_r8)

# ifdef MASKING
              bt(i,j,k,nnew,ibioB)=  bt(i,j,k,nnew,ibioB)*rmask(i,j)
# endif
            END DO
          END DO
        END DO

#endif
#ifdef ICE_BIO
# if defined CLIM_ICE_1D


        DO i=Istr,Iend

          if(itL(i,j,nstp,1).gt.0_r8)THEN

            it(i,j,nnew,iIceNO3)=MAX(it(i,j,nstp,iIceNO3)+              &
     &                           BioBI(i,iIceNO3)-Bio_bakBI(i,iIceNO3)  &
     &                               ,0.0_r8)
            it(i,j,nnew,iIceNH4)=MAX(it(i,j,nstp,iIceNH4)+              &
     &                           BioBI(i,iIceNH4)-Bio_bakBI(i,iIceNH4)  &
     &                               ,0.0_r8)
            it(i,j,nnew,iIcePhL)=MAX(it(i,j,nstp,iIcePhL)+              &
     &                           BioBI(i,iIcePhL)-Bio_bakBI(i,iIcePhL)  &
     &                               ,0.0_r8)

          else

            it(i,j,nnew,iIcePhL)=0_r8
            it(i,j,nnew,iIceNH4)=0_r8
            it(i,j,nnew,iIceNO3)=0_r8
          endif


          itL(i,j,nnew,iIceLog)=itL(i,j,nstp,iIceLog)


#  ifdef MASKING
          it(i,j,nnew,iIcePhL)  = it(i,j,nnew,iIcePhL)*rmask(i,j)
          it(i,j,nnew,iIceNH4)  = it(i,j,nnew,iIceNH4)*rmask(i,j)
          it(i,j,nnew,iIceNO3)  = it(i,j,nnew,iIceNO3)*rmask(i,j)
          itL(i,j,nnew,iIceLog) = itL(i,j,nnew,iIceLog)*rmask(i,j)
#  endif

        END DO

# elif defined BERING_10K
        DO i=Istr,Iend
          if (IceLog(i,j,nstp).ge.0_r8 )THEN

            !  ajh added zero trap on these

            IcePhL(i,j,nnew) = max(0.,IcePhL(i,j,nnew))
            IceNO3(i,j,nnew) =  max(0.,IceNO3(i,j,nnew))
            IceNH4(i,j,nnew) =  max(0.,IceNH4(i,j,nnew))
          else
            IcePhL(i,j,nnew) =0_r8
            IceNO3(i,j,nnew) =  0_r8
            IceNH4(i,j,nnew) =  0_r8
          endif

#  ifdef MASKING
          IcePhL(i,j,nnew) = IcePhL(i,j,nnew)*rmask(i,j)
          IceNO3(i,j,nnew) =  IceNO3(i,j,nnew)*rmask(i,j)
          IceNH4(i,j,nnew) =  IceNH4(i,j,nnew)*rmask(i,j)
#  endif
            IceLog(i,j,nnew)=IceLog(i,j,nstp)

        END DO

# endif


#endif


#ifdef STATIONARY
        DO k=1,N(ng)

          DO i=Istr,Iend

             st(i,j,k,nstp,1)  =    Stat3(i,k,1)
             st(i,j,k,nstp,2)  =    Stat3(i,k,2)
             st(i,j,k,nstp,3)  =    Stat3(i,k,3)
             st(i,j,k,nstp,4)  =    Stat3(i,k,4)
             st(i,j,k,nstp,5)  =    Stat3(i,k,5)
             st(i,j,k,nstp,6)  =    Stat3(i,k,6)
             st(i,j,k,nstp,7)  =    Stat3(i,k,7)
             st(i,j,k,nstp,8)  =    Stat3(i,k,8)
             st(i,j,k,nstp,9)  =    Stat3(i,k,9)
             st(i,j,k,nstp,10) =    Stat3(i,k,10)
             st(i,j,k,nstp,11) =    Stat3(i,k,11)
             st(i,j,k,nstp,12) =    Stat3(i,k,12)
             st(i,j,k,nstp,13) =    Stat3(i,k,13)
             st(i,j,k,nstp,14) =    Stat3(i,k,14)
             st(i,j,k,nstp,15) =    Stat3(i,k,15)
             st(i,j,k,nstp,16) =    Stat3(i,k,16)


          END DO
        END DO
#endif
#ifdef STATIONARY2

        DO i=Istr,Iend
          st2(i,j,nstp,1) =   Stat2(i,1)
          st2(i,j,nstp,2) =   Stat2(i,2)
          st2(i,j,nstp,3) =   Stat2(i,3)
          st2(i,j,nstp,4) =   Stat2(i,4)
          st2(i,j,nstp,5) =   Stat2(i,5)
          st2(i,j,nstp,6) =   Stat2(i,6)
          st2(i,j,nstp,7) =   Stat2(i,7)
          st2(i,j,nstp,8) =   Stat2(i,8)
# ifdef MASKING
          st2(i,j,nstp,1) =  st2(i,j,nstp,1)*rmask(i,j)
          st2(i,j,nstp,2) =  st2(i,j,nstp,2)*rmask(i,j)
          st2(i,j,nstp,3) =  st2(i,j,nstp,3)*rmask(i,j)
          st2(i,j,nstp,4) =  st2(i,j,nstp,4)*rmask(i,j)
          st2(i,j,nstp,5) =  st2(i,j,nstp,5)*rmask(i,j)
          st2(i,j,nstp,6) =  st2(i,j,nstp,6)*rmask(i,j)
          st2(i,j,nstp,7) =  st2(i,j,nstp,7)*rmask(i,j)
          st2(i,j,nstp,8) =  st2(i,j,nstp,8)*rmask(i,j)

# endif
        END DO
#endif
#ifdef PROD3
        DO k=1,N(ng)
          DO i=Istr,Iend
            pt3(i,j,k,nnew,iPhSprd) = pt3(i,j,k,nstp,iPhSprd) +         &
     &                             Prod(i,k,iPhSPrd)

            pt3(i,j,k,nnew,iPhLprd) = pt3(i,j,k,nstp,iPhLprd) +         &
     &                             Prod(i,k,iPhLPrd)
            pt3(i,j,k,nnew,iMZSprd) = pt3(i,j,k,nstp,iMZSprd) +         &
     &                             Prod(i,k,iMZSPrd)
            pt3(i,j,k,nnew,iMZLprd) = pt3(i,j,k,nstp,iMZLprd) +         &
     &                             Prod(i,k,iMZLPrd)
            pt3(i,j,k,nnew,iCopPrd) = pt3(i,j,k,nstp,iCopPrd) +         &
     &                             Prod(i,k,iCopPrd)
            pt3(i,j,k,nnew,iNCaPrd) = pt3(i,j,k,nstp,iNCaPrd) +         &
     &                             Prod(i,k,iNCaPrd)
            pt3(i,j,k,nnew,iEupPrd) = pt3(i,j,k,nstp,iEupPrd) +         &
     &                             Prod(i,k,iEupPrd)
# ifdef JELLY

            pt3(i,j,k,nnew,iJelPrd) = pt3(i,j,k,nstp,iJelPrd) +         &
     &                             Prod(i,k,iJelPrd)
# endif
# ifdef FEAST

            ! Note that fish zoopmort are lagged by one timestep due to state update

!           pt3(i,j,k,nnew,iFCopMort)   = ZoopFishDeath(i,j,k,1)
!           pt3(i,j,k,nnew,iFNCaSMort)  = ZoopFishDeath(i,j,k,2)
!           pt3(i,j,k,nnew,iFNCaOMort)  = ZoopFishDeath(i,j,k,3)
!           pt3(i,j,k,nnew,iFEupSMort)  = ZoopFishDeath(i,j,k,4)
!           pt3(i,j,k,nnew,iFEupOMort)  = ZoopFishDeath(i,j,k,5)
            pt3(i,j,k,nnew,iFishOne)   = GF%ofdat(i,j,k,1)
            pt3(i,j,k,nnew,iFishTwo)   = GF%ofdat(i,j,k,2)
            pt3(i,j,k,nnew,iFishThree) = GF%ofdat(i,j,k,3)
            pt3(i,j,k,nnew,iFishFour)  = GF%ofdat(i,j,k,4)
            pt3(i,j,k,nnew,iFishFive)  = GF%ofdat(i,j,k,5)
            pt3(i,j,k,nnew,iFishSix)   = GF%ofdat(i,j,k,6)
            pt3(i,j,k,nnew,iFishSeven) = GF%ofdat(i,j,k,7)
            pt3(i,j,k,nnew,iFishEight) = GF%ofdat(i,j,k,8)
# endif
          END DO
        END DO
#endif


#ifdef PROD2

        DO i=Istr,Iend
# ifdef BENTHIC

          pt2(i,j,nnew,iBenPrd) = pt2(i,j,nstp,iBenPrd)                 &
     &                           +Prod2(i,iBenPrd)
# endif

# ifdef ICE_BIO
          pt2(i,j,nnew,iIAPrd) = pt2(i,j,nstp,iIAPrd)+Prod2(i,iIAPrd)

# endif
          pt2(i,j,nnew,iXPrd) = pt2(i,j,nstp,iXPrd)+Prod2(i,iXPrd)
        END DO
#endif


      END DO J_LOOP

#if defined EW_PERIODIC || defined NS_PERIODIC

      ! Apply periodic boundary conditions.

      DO itrc=1,NBT
        ibio=idbio(itrc)

        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,ibio))
      END DO
#endif

#ifdef DISTRIBUTE

      ! Exchange boundary data.


      !ajh
      !added block on this for passives when feast is present
      !NOTE will need to exchange when passives/fish are changed elsewhere
# ifdef FEAST_NOEXCHANGE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NAT,         &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,:))

      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NAT+NPT+1, NT(ng),                    &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,:))

# else
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,:))
# endif
!
# ifdef STATIONARY
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NBTS,        &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    st(:,:,:,nnew,:))
# endif

#endif


#if defined ICE_BIO
# ifdef BERING_10K

      CALL IcePhLbc_tile (ng, tile,                                    &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp, nnew,                                     &
     &                ui, vi, IcePhL)
      CALL IceNO3bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp, nnew,                                     &
     &                ui, vi, IceNO3)
      CALL IceNH4bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp, nnew,                                     &
     &                ui, vi, IceNH4)

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    IcePhL(:,:,nnew), IceNO3(:,:,nnew),         &
     &                    IceNH4(:,:,nnew))

#  endif






# endif
#endif

#if defined FEAST
# include "feast_step.h"
#endif

      RETURN
      END SUBROUTINE biology_tile

!=====================================================================
!  BioVert: biological tracer vertical movement
!  K. Kearney, 03/2017
!
!  This routine handles all vertical movement of tracers, regardless
!  of direction.  It provides a wrapper around the BIOSINK routine.
!  I've chosen to use this wrapper rather than separate sinking and
!  rising functions to make code maintenance easier (any changes to
!  the vertical movement scheme can be isolated to the BIOSINK
!  subroutine)
!
! Input (and *output) variables:
!
!   nn:      1 x 1,    number of vertical layers
!   wBio:    1 x 1,    vertical velocity (negative down, positive up,
!                      m/d)
!   Bio:     1 x nn,   concentration in each layer (mass/m^3)
!  *dBioOut: 1 x nn,   change in concentration per layer (mass/m^3)
!   HzL:     1 x nn,   thickness of each layer (m)
!   dtdays:  1 x 1,    time step (days)
!   z_wL:    1 x nn+1, depth of layer edges (m, negative down)
!   zlimit:  1 x 1,    depth limit for movement (m, negative down).  If
!                      value is outside the range of z_wL, no limit is
!                      imposed and material can pass through the bottom
!                      edge.
!  *flx:     1 x 1,    flux lost through the bottom edge for sinking, or
!                      top edge for rising (mass/m^2)
!
!=====================================================================
      subroutine BioVert(nn,wBio,Bio,dBioOut,HzL,dtdays,z_wL,zlimit,flx)

      USE mod_kinds
      implicit none

      integer,  intent(in) :: nn
      real(r8), intent(in) :: wBio    ! negative down, positive up
      real(r8), intent(in) :: zlimit
      real(r8), intent(in) :: dtdays
      real(r8), intent(in) :: HzL(1,1,nn)
      real(r8), intent(in) :: z_wL(1,1,0:nn)
      real(r8), intent(in) :: Bio(1,nn)

      real(r8) :: Bflip(1,nn)
      real(r8) :: zwflip(1,0:nn)
      real(r8) :: Hzflip(1,nn)
      real(r8) :: dBflip(1,nn)
      real(r8) :: zlimflip
      integer k

      real(r8), intent(out) :: dBioOut(1,nn)
      real(r8), intent(out) :: flx

      if (wBio .le. 0) then ! Sinking

        call BIOSINK(nn,-wBio,Bio,dBioOut,HzL,dtdays,z_wL,zlimit,flx)

      else ! Rising

        ! Flip the water column upside down.  This new water column has its
        ! surface at z = 0 and bottom at -(water depth + free surface height).
        ! Layer thicknesses and biomass concentrations are the same as
        ! before, but in reverse order.

        DO k = 1,nn
          Bflip(1,k) = Bio(1,nn+1-k) ! flip
          Hzflip(1,k) = HzL(1,1,nn+1-k) ! flip
        END DO
        DO k = 0,nn
          zwflip(1,k) = z_wL(1,1,0) - z_wL(1,1,nn-k) ! make surface the bottom
        END DO

        zlimflip = z_wL(1,1,0) - zlimit ! relocate in flipped layers

        dBflip = 0

        call BIOSINK(nn,wBio,Bflip,dBflip,Hzflip,dtdays,zwflip,zlimflip,flx)

        ! Flip db back

        DO k = 1,nn
          dBioOut(1,k) = dBflip(1,nn+1-k) ! flip back
        END DO

      endif

      ! Calculate how much mass was lost

!       flx = 0
!       do k = 1,nn
!         flx = flx - (HzL(1,1,k)*dBioOut(1,k))
!       end do

      END SUBROUTINE BioVert

!=====================================================================
! BIOSINK  particle sinking subroutine After J. Warner sed sink code
! G. Gibson July 2008
!
! K.Kearney 2017: Small modifications to input/output parameters (was
! experimenting with this code vs the GOANPZ zlimit formulation, so I
! wanted to match the inputs between the two... plus added the flx
! output).  This version now allows for explicit limiting of sinking
! depth.  I also eliminated the i-loop in this code, moving that to the
! main calling routine; this was done to make testing easier outside of
! the ROMS tiling scheme (probably at the expense of some speed, but I'm
! okay with that).
!
! Input (and *output) variables:
!
!   nn:      1 x 1,    number of vertical layers
!   wBio:    1 x 1,    sinking rate (positive value, m/d)
!   Bio:     1 x nn,   concentration in each layer (mass/m^3)
!  *dBioOut: 1 x nn,   change in concentration per layer (mass/m^3)
!   HzL:     1 x nn,   thickness of each layer (m)
!   dtdays:  1 x 1,    time step (days)
!   z_wL:    1 x nn+1, depth of layer edges (m, negative down)
!   zlimit:  1 x 1,    depth limit for sinking (m, negative down).  If
!                      value is outside the range of z_wL, no limit is
!                      imposed and material can pass through the bottom
!                      edge.
!  *flx:     1 x 1,    flux lost through the bottom edge (mass/m^2)
!
!=====================================================================
      subroutine BIOSINK(nn,wBio,Bio,dBioOut,HzL,dtdays,z_wL,zlimit,flx)
!
      USE mod_kinds
      implicit none
!
      integer, intent(in)  :: nn
      real(r8), intent(in) :: wBio
      real(r8), intent(in) :: zlimit
      real(r8), intent(in) :: z_wL(1,0:nn)
      real(r8), intent(in) :: Bio(1,nn)
      real(r8), intent(in) :: HzL(1,nn)
      real(r8), intent(in) :: dtdays

      real(r8), intent(out) :: dBioOut(1,nn)
      real(r8), intent(out) :: flx

      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff

      real(r8)::  wBiod(1,0:nn)
      real(r8) :: FC(1,0:nn)

      real(r8) :: Hz_inv(1,nn)
      real(r8) :: Hz_inv2(1,nn)
      real(r8) :: Hz_inv3(1,nn)

      integer :: ksource(1,nn)

      real(r8) :: qR(1,nn)
      real(r8) :: qL(1,nn)
      real(r8) :: WL(1,nn)
      real(r8) :: WR(1,nn)
      real(r8), dimension(1,nn) :: qc

!  Compute inverse thickness to avoid repeated divisions.
!

      DO k=1,nn
       DO i=1,1
         Hz_inv(i,k)=1.0_r8/HzL(i,k)
       END DO
      END DO
      DO k=1,nn-1
       DO i=1,1
         Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=2,nn-1
       DO i=1,1
         Hz_inv3(i,k)=1.0_r8/(HzL(i,k-1)+HzL(i,k)+HzL(i,k+1))
       END DO
      END DO

      DO k=1,nn
         DO i=1,1
           qc(i,k)=Bio(i,k)
         END DO
      END DO

!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!

      DO k=nn-1,1,-1
        DO i=1,1

          FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
        END DO
      END DO

!     DO k=nn-1,1,-1
!       DO i=1,1
!         print*,'LBi=',LBi,'UBi=',UBi
!         print*,'i=',i,'k=',k
!         if(i.le.UBi)THEN
!           FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
!         endif
!       END DO
!     END DO


      DO k=2,nn-1
        DO i=1,1
          dltR=HzL(i,k)*FC(i,k)
          dltL=HzL(i,k)*FC(i,k-1)
          cff=HzL(i,k-1)+2.0_r8*HzL(i,k)+HzL(i,k+1)
          cffR=cff*FC(i,k)
          cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
          cff=(dltR-dltL)*Hz_inv3(i,k)
          dltR=dltR-cff*HzL(i,k+1)
          dltL=dltL+cff*HzL(i,k-1)
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
          WR(i,k)=(2.0_r8*dltR-dltL)**2
          WL(i,k)=(dltR-2.0_r8*dltL)**2
        END DO
      END DO
      cff=1.0E-14_r8
      DO k=2,nn-2
         DO i=1,1
          dltL=MAX(cff,WL(i,k  ))
          dltR=MAX(cff,WR(i,k+1))
          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
          qL(i,k+1)=qR(i,k)
        END DO
      END DO

      DO i=1,1
        FC(i,nn)=0.0_r8              ! no-flux boundary condition
#if defined LINEAR_CONTINUATION
        qL(i,nn)=qR(i,nn-1)
        qR(i,nn)=2.0_r8*qc(i,nn)-qL(i,nn)
#elif defined NEUMANN
        qL(i,nn)=qR(i,nn-1)
        qR(i,nn)=1.5_r8*qc(i,nn)-0.5_r8*qL(i,nn)
#else
        qR(i,nn)=qc(i,nn)         ! default strictly monotonic
        qL(i,nn)=qc(i,nn)         ! conditions
        qR(i,nn-1)=qc(i,nn)
#endif
#if defined LINEAR_CONTINUATION
        qR(i,1)=qL(i,2)
        qL(i,1)=2.0_r8*qc(i,1)-qR(i,1)
#elif defined NEUMANN
        qR(i,1)=qL(i,2)
        qL(i,1)=1.5_r8*qc(i,1)-0.5_r8*qR(i,1)
#else
        qL(i,2)=qc(i,1)                 ! bottom grid boxes are
        qR(i,1)=qc(i,1)                 ! re-assumed to be
        qL(i,1)=qc(i,1)                 ! piecewise constant.
#endif
      END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
      DO k=1,nn
        DO i=1,1
          dltR=qR(i,k)-qc(i,k)
          dltL=qc(i,k)-qL(i,k)
          cffR=2.0_r8*dltR
          cffL=2.0_r8*dltL
          IF ((dltR*dltL).lt.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
        END DO
      END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by nn).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!

      DO k=1,nn
        DO i=1,1
          cff=dtdays*ABS(wBio)
          FC(i,k-1)=0.0_r8
          WL(i,k)=z_wL(i,k-1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO
      DO k=1,nn
        DO ks=k,nn-1
          DO i=1,1
            IF (WL(i,k).gt.z_wL(i,ks)) THEN
              ksource(i,k)=ks+1
              FC(i,k-1)=FC(i,k-1)+WR(i,ks)
            END IF
          END DO
        END DO
      END DO
!
!  Finalize computation of flux: add fractional part.
!
      DO k=1,nn
        DO i=1,1
          ks=ksource(i,k)
          cu=MIN(1.0_r8,(WL(i,k)-z_wL(i,ks-1))*Hz_inv(i,ks))
          FC(i,k-1)=FC(i,k-1)+                                          &
     &                  HzL(i,ks)*cu*                                   &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))

!
!  G.Gibson  - FC is the flux into the level
!            - should be 0 at the surface

          if (k.eq.nn) then
            FC(i,k)=0.0_r8
          endif
        END DO
      END DO

      ! Hard barrier on sinking

      IF ((zlimit.ge.minval(z_wL)) .and. (zlimit.le.maxval(z_wL))) then
        DO k=0,nn
          if (z_wL(1,k) .le. zlimit) then
            FC(1,k) = 0
          end if
        end do
      end if

      DO k=1,nn
        DO i=1,1
          dBioOut(i,k) = (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
        END DO
      END DO

      flx = FC(1,0)

      RETURN
      END SUBROUTINE BIOSINK




        Function ComputeDensity(Temp1,Sal1)
!----------------------------------------------------------------
! Computes the water column density from salinity and temperature
! Returns sigma-t
!----------------------------------------------------------------
        USE mod_kinds

        Real(r8) ComputeDensity
        Real(r8) Temp1, Sal1
        Real(r8) Sig
        Sig = 999.842594 + 0.06793952 * Temp1
        Sig = Sig - 0.00909529 * Temp1 ** 2 +                      &
     &          0.0001001685 * Temp1 ** 3
        Sig = Sig - 0.000001120083 * Temp1 ** 4 +                  &
     &          0.000000006536332 * Temp1 ** 5
        Sig = Sig + 0.824493 * Sal1 - 0.0040899 * Temp1 * Sal1
        Sig = Sig + 0.000076438 * Temp1 ** 2 * Sal1 -              &
     &          0.00000082467 * Temp1 ** 3 * Sal1
        Sig = Sig + 0.0000000053875 * Temp1 ** 4 * Sal1 -          &
     &          0.00572466 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00010227 * Temp1 * Sal1 ** (3 / 2) -         &
     &          0.0000016546 * Temp1 ** 2 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00048314 * Sal1 ** 2
         ComputeDensity = Sig - 1000
        End Function ComputeDensity
!===============================================================
        Function GetPhytoResp2(Temp1, Tref, KbmPh)
!------------------------------------------------------
! Computes the temperature correction for phytoplankton
! respiration according to Arhonditsis 2005.
!------------------------------------------------------
        USE mod_kinds

        Real(r8) GetPhytoResp2
        Real(r8) Temp1      !Temperature, passed
        Real(r8) Tref       !Reference temperature
        Real(r8) KbmPh      !Half saturation, temperature

        Real(r8) Resp       !Returned variable

        Resp = exp(KbmPh * (Temp1 - Tref))
        GetPhytoResp2 = Resp
        Return
        End Function GetPhytoResp2
!=====================================================================
      FUNCTION GetLightLimIronSml(alphaPh, PAR1, Pmax1,          &
     &     CrChlRatio1,IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIronSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIronSml = LightLim
      END FUNCTION GetLightLimIronSml
!=====================================================================
      FUNCTION GetLightLimSml(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLimSml = LightLim
      END FUNCTION GetLightLimSml
!=====================================================================
      FUNCTION GetLightLimIron(alphaPh, PAR1, Pmax1, CrChlRatio1,  &
     &     IronLim1, ParMax)
!------------------------------------------------------------------
! Light lim with varying alpha. Works with iron limitation. Alph is
! a function of the surface light intensity.
!------------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

       Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
       Real(r8) GetLightLimIron
       Real(r8) Alpha,LightLim,OffSet,ParMax

       !Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
       !if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = alphaPh
       end if
       LightLim = TANH( Alpha * Par1/Pmax1/CrChlRatio1/IronLim1)
       GetLightLimIron = LightLim

      END FUNCTION GetLightLimIron
!=======================================================================
      FUNCTION GetLightLim(alphaPh, PAR1, Pmax1, CrChlRatio1, ParMax)
!-----------------------------------------------------------------
! Generates a light lim with varying alphaPh without iron
!-----------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

       Real(r8) :: alphaPh,Pmax1,PAR1,CrChlRatio1
       Real(r8) GetLightLim
       Real(r8) Alpha,LightLim,OffSet,ParMax

!       Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
!       if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = alphaPh
       end if
       LightLim = TANH(alphaPh * PAR1/Pmax1/CrChlRatio1)

       GetLightLim = LightLim
      END FUNCTION GetLightLim
!==============================================================
        Function GetCopepodResp(Temp1,respVal,ktbm,Tref)
!--------------------------------------------------------------
! Computes copepod respiration according to Arhonditsis (2005).
!--------------------------------------------------------------
        USE mod_kinds
        USE mod_param

        Real(r8) GetCopepodResp
        real(r8) :: respVal
        Real(r8) Temp1       !Passed variable
!        Real(r8) :: bm  = 0.04   !Basal metabolic rate day**-1
        Real(r8) :: ktbm  != 0.05 !Temperature response degrees C**-1
        Real(r8) :: Tref  != 20   !Reference temperature degrees C
        Real(r8) Resp        !Returned variable

        Resp = respVal * exp(ktbm * (Temp1 - Tref))
        GetCopepodResp = Resp
        Return
        End Function GetCopepodResp
!==============================================================
        Function GetJelResp(Temp1,bmJ,ktbmJ,TrefJ)

        USE mod_param
        implicit none
!
        real(r8) :: GetJelResp
        real(r8) :: Temp1          !Passed variable
        real(r8) :: bmJ   != 0.04   !Basal metabolic rate day**-1
        real(r8) :: ktbmJ != 0.05 !Temperature response degrees C**-1
        real(r8) :: TrefJ != 20   !Reference temperature degrees C
        real(r8) :: Resp           !Returned variable

        Resp = bmJ * exp(ktbmJ * (Temp1 - TrefJ))
        GetJelResp = Resp
        Return
        End Function GetJelResp
!=================================================================
        Function GetBasalMetabolism(respPh,kfePh,Iron1)
!---------------------------------------------------------
! Computes an iron correction for the basal metabolism for
! the phytoplankton respiration calculation
!---------------------------------------------------------
        USE mod_kinds

        Real(r8) GetBasalMetabolism
        Real(r8) Iron1     !Iron concentration
        Real(r8) kfePh     !Half saturation for iron
        Real(r8) respPh    !Phytoplankton uncorrected basal metabolism
        Real(r8) BaseMet   !Phytoplankton basal metabolism

        BaseMet = Iron1/(kfePh + Iron1)
        BaseMet = BaseMet * ((kfePh +2)/2)
        GetBasalMetabolism = BaseMet * respPh

        Return
        End Function GetBasalMetabolism
!=========================================================
        Function GetNitrif(Temp1,Dep1,NH4R)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrif
                               !--------------------------------
        Real(r8) Temp1, Dep1, NH4R !Passed variables
        Real(r8) NH4conv  /14/     !mg N/mol N
        Real(r8) KNH4Nit  /0.08/   !Half Sat Con mg N/m3/day
        Real(r8) KTNitr   /0.002/  !Temperature responce dec C^2
        Real(r8) ToptNtr  /28/     !Optimum nitrification temp
        Real(r8) Zox      /20/     !50% nitrification depth
        Real(r8) Nexp     /6/      !Exponent to adjust profile shape
        Real(r8) NitrMax  /0.011/  !Maximum nitrification (mM/m3/d
                               !--------------------------------

        Real(r8) Nitr, DepCor

        NH4R = NH4R * NH4conv
        Nitr = NH4R/(KNH4Nit + NH4R)
        Nitr = Nitr * exp(-KTNitr*(Temp1 - ToptNtr)**2)
        DepCor = (Dep1**Nexp)/( (Zox**Nexp) + Dep1**Nexp)
        Nitr = (Nitr * DepCor) * NitrMax
        GetNitrif = Nitr
        Return
        End Function GetNitrif
!========================================================================
        Function GetNitrif2(Temp1,PAR1,NH4R)
        !---------------------------------------------------------
        !Computes nitrificaton from Kawamiya with light correction
        !from Fennel; Kawamiya (2000), Fennel (2006)
        !---------------------------------------------------------
        USE mod_kinds
        USE mod_param
         Real(r8) GetNitrif2
                               !--------------------------------
        Real(r8) :: Temp1, NH4R       !Passed variables
        Real(r8) :: I0 = 0.0095     !Threshold,light inhibition, W m-2
        Real(r8) :: KI = 4.0        !Half Saturation light intensity, W m-2
        Real(r8) :: KN0 = 0.03       !Nitrification at 0 deg C, day-1
        Real(r8) :: KNT = 0.0693     !Temperature coefficient
        Real(r8) :: ParW              !Par in watts
        Real(r8) :: NitrMax           !Maximum nitrification
        Real(r8) :: Nitr              !Nitrification
                               !---------------------------------
        Real(r8) :: cff1, PAR1

        !-----------------------------------
        !Temperature dependent nitrification
        !-----------------------------------
        KI = 1.5
        KN0 = 0.15
        !KNT = 0.07
        NitrMax = (KN0*Exp(KNT*Temp1))*NH4R
        !-----------------------------------
        !Convert PAR in E m-2 d-1 to W day-1
        !-----------------------------------
        ParW = PAR1/0.394848_r8
        !---------------------------------
        !Light correction of nitrification
        !---------------------------------
        cff1 = (ParW-I0)/(KI+ParW-I0)
        Nitr = NitrMax*(1-MAX(0.0_r8,cff1))
        GetNitrif2 = Nitr
        Return
        End Function GetNitrif2
!-----------------------------------------------------------

  Function GetNitrifLight(Par1,tI0,KI)
   USE mod_kinds

        Real(r8) GetNitrifLight
        Real(r8) :: Par1               !--------------------------------
        Real(r8) :: tI0                 !Threshold,light inhibition, W m-2
        Real(r8) :: KI                 !Half Saturation light intensity, W m-2
        Real(r8) :: ParW               !Par in watts
        Real(r8) :: cff1

      ParW = Par1/0.394848_r8 !convert PAR back to watts
            cff1 = (ParW-tI0)/(KI+ParW-tI0)

      GetNitrifLight=(1-MAX(0.0_r8,cff1))

        Return
        End Function GetNitrifLight
!-----------------------------------------------
   Function GetNitrifMaxK(Temp1,KnT,Nitr0)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrifMaxK

        Real(r8) :: Temp1         !--------------------------------
        Real(r8) :: KnT           !Passed variables
        Real(r8) :: Nitr0         !--------------------------------


  GetNitrifMaxK=Nitr0*exp(KnT*Temp1)

        Return
        End Function GetNitrifMaxK

!-----------------------------------------------------------

  Function GetNitrifMaxA(Nitr0, ktntr,Temp1,ToptNtr)
   USE mod_kinds

        Real(r8) GetNitrifMaxA

        Real(r8) :: Temp1
        Real(r8) :: Ktntr           !--------------------------------
        Real(r8) :: Nitr0           !Passed variables
        Real(r8) :: ToptNtr         !--------------------------------


  GetNitrifMaxA=Nitr0*exp(-ktntr*(Temp1 - ToptNtr)**2)


        Return
        End Function GetNitrifMaxA

!-----------------------------------------------------------
!=====================================================================
      FUNCTION GetLightLimIron2(alphaPh, PAR1, Pmax1, CrChlRatio1, &
     &     IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIron2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIron2 = LightLim
      END FUNCTION GetLightLimIron2
!=====================================================================
      FUNCTION GetLightLim2(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLim2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLim2 = LightLim
      END FUNCTION GetLightLim2


!===================================================================
     FUNCTION ComputeStability(ng,z_wL,Dens)
!--------------------------------------------------------
! Computes the stability parameter from density and depth
! Simpson stability parameter.
!--------------------------------------------------------
     USE mod_param
!
     implicit none
          integer, intent(in) ::  ng
     real(r8), intent(in) :: z_wL(0:N(ng))
     real(r8), intent(in) :: Dens(1:N(ng))
     Integer k,indx
     real ComputeStability
     real dZ,SumSigma,SumDep,RMean,Z1,Dep1,Dep2

     !SumSigma = Dens(30)
     !SumSigma = 0
     SumDep = 0
     indx = N(ng)     !30

     SumDep = z_wL(indx) * (-1.0_r8)
     SumSigma = Dens(indx) * SumDep
     indx = indx - 1
    DO k = 1, N(ng) - 1
       IF (z_wL(indx) .ge. -140) THEN
          dZ = z_wL(indx + 1) - z_wL(indx)
          !Print *, '**',Dens(k),z_wL(k),SumSigma,SumDep,dZ,k
          SumSigma = SumSigma + 0.5 * &
    &                (Dens(indx) + Dens(indx+1)) * dZ
          SumDep = SumDep + dZ
       END IF
       indx = indx - 1
     END DO
     !if (SumDep == 0) then
     !    print *,ng,z_wL(indx),Dens(indx)
     !end if
     RMean = SumSigma/SumDep
     !print *, '&&&&%%%%',RMean, SumDep, SumSigma,z_wL(1)
     !Print *, z_wL(30), Dens(1),Dens(30)
     indx = N(ng)-1          !29
     Z1 = 0.5 * (z_wL(indx))
     dZ = z_wL(indx)
     SumSigma = Z1 * (RMean - Dens(indx - 1)) * dZ
     !print *, 'xxx',SumSigma
     DO k = 1, N(ng) - 1
       IF (z_wL(indx) .ge. -140) THEN
          Dep1 = z_wL(indx) * (-1.0_r8)
          Dep2 = z_wL(indx -1) * (-1.0_r8)
          Z1 = 0.5 * (Dep1 + Dep2)
          dZ = (Dep2 - Dep1)
!         print *, '&&**',z_wL(indx),z_wL(indx + 1),dZ,Z1
!        print *, '&&**',SumSigma,z_wL(k),dens(k),Z1*(Rmean -   &
!     &   0.5*(Dens(k) + Dens(k+1))) * dZ, dZ,rmean
          SumSigma = SumSigma + Z1*(Rmean -                     &
    &         0.5*(Dens(indx) + Dens(indx + 1))) * dZ
       END IF
       indx = indx - 1
     END DO
     !Z1 = 0.5 * z_wL(N(ng))
     !dZ = z_wL(N(Ng))
     !SumSigma = SumSigma + Z1*(RMean - Dens(N(ng))) * dZ
     ComputeStability = -9.8 * SumSigma / SumDep

     END FUNCTION ComputeStability

