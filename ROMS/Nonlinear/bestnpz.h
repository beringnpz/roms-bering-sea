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

      real(r8), dimension(N(ng)) :: Btmp, Hztmp
      real(r8), dimension(0:N(ng)) :: zwtmp
      real(r8) :: sinkout2

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

      ! Biological source/sinks

      real(r8) :: LightLimS, NOLimS, NHLimS, IronLimS
      real(r8) :: LightLimL, NOLimL, NHLimL, IronLimL
      real(r8) :: alphaPhSv, alphaPhLv, DrateS, DrateL, PmaxS, PmaxL, PmaxsS, PmaxsL
      real(r8) :: IcePhlAvail
      real(r8), dimension(IminS:ImaxS,N(ng)) :: BasMetMZL, BasMetCop, BasMetNC, BasMetCM, BasMetEup

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

          ! TODO: overhauled to here (still need to compile check)

          !=====================================================
          ! Phytoplankton Linear Mortality and Senescence Terms
          !=====================================================
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!             cff1 = MAX( minmPhS , maxmPhS -                           &
!    &              ( maxmPhS - minmPhS ) * Bio(i,k,iNO3) / NcritPhS)
!             cff2 = MAX( minmPhL , maxmPhL -                           &
!    &              ( maxmPhL - minmPhL ) * Bio(i,k,iNO3) / NcritPhL)


              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &                         mPhS* Bio(i,k,iPhS) * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                     &
     &                         mPhL* Bio(i,k,iPhL) * dtdays

#ifdef STATIONARY
              if(i.eq.2) THEN
                if(j.eq.2) THEN
!                 Stat3(i,k,2)=  mPhS * dtdays * Bio(i,k,iPhL)
                endif
              endif
              if(i.eq.3) THEN
                if(j.eq.3) THEN
!                 Stat3(i,k,2)=  mPhL * dtdays * Bio(i,k,iPhL)
                endif
              endif
#endif
              !--------------------------------------------------
              !  Additions to detritus pool - phytoplankton mort
              !--------------------------------------------------

              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &                         ( mPhS * Bio(i,k,iPhS)) * dtdays
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &                         ( mPhL * Bio(i,k,iPhL)) * dtdays
#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iDet) = bflx(iPhS,iDet)                       &
     &            + mPhS*Bio(i,k,iPhS)* dtdays*xi
                bflx(iPhL,iDetF) = bflx(iPhL,iDetF)                     &
     &            + mPhL*Bio(i,k,iPhL)* dtdays*xi
              END IF
#endif
            END DO
          END DO

          !============================================================
          !  Microzooplankton Mortality - use only linear OR QUADRATIC
          !============================================================

          DO k=1,N(ng)
            DO i=Istr,Iend


              !  Linear

!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!     &                        mMZS * Bio(i,k,iMZS) * dtdays
!             DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
!     &                         mMZL * Bio(i,k,iMZL) * dtdays

              !  Quadratic

!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!    &                         mpredMZS*dtdays*Bio(i,k,iMZS)**2
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &                         mpredMZL*dtdays*Bio(i,k,iMZL)**2

#ifdef STATIONARY
              if(i.eq.2) THEN
                if(j.eq.2) THEN
!                 Stat3(i,k,7)=  mpredMZL*dtdays*Bio(i,k,iMZL)**2
                endif
              endif
#endif
              !----------------------------------------------------------
              !  Additions to detritus pool - natural microzoo mortality
              !----------------------------------------------------------

              ! if linear (George)

!             Bio(i,k,iDet) = DBio(i,k,iDet) +                          &
!    &                        mMZL * Bio(i,k,iMZL) * dtdays

              ! if quadratic (Ken)

              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
      &                        (mpredMZL * Bio(i,k,iMZL)**2 ) * dtdays

#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
!               bflx(iMZS,iDet)= bflx(iMZS,iDet)                        &
!     &           + mMZS * Bio(i,k,iMZS) * dtdays*xi

                bflx(iMZL,iDet)= bflx(iMZL,iDet)                        &
     &            + mpredMZL*( Bio(i,k,iMZL)**2) * dtdays*xi
              END IF
#endif
            END DO
          END DO

          !============================================
          !  Mesozooplankton Mortality (Closure terms)
          !============================================

          DO k=1,N(ng)
            DO i=Istr,Iend

              TFEup = Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8)

#ifdef fixedPRED
              DBio(i,k,iCop)  = DBio(i,k,iCop)  - 0.5*Hz(i,j,k)/dtdays
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) - 0.5*Hz(i,j,k)/dtdays
              DBio(i,k,iEupS) = DBio(i,k,iEupS) - 1  *Hz(i,j,k)/dtdays
              DBio(i,k,iNCaO) = DBio(i,k,iNCaS) - 0.5*Hz(i,j,k)/dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupS) - 1  *Hz(i,j,k)/dtdays

#else
# ifdef FEAST
              predSumCop  = TFEup*(mpredCop + fpredCop  * GF%zoop_force(1,1,i,j,1))*Bio(i,k,iCop)**2
              predSumNCaS = TFEup*(mpredNca + fpredNcaS * GF%zoop_force(1,2,i,j,1))*Bio(i,k,iNCaS)**2
              predSumEupS = TFEup*(mpredEup + fpredEupS * GF%zoop_force(1,4,i,j,1))*Bio(i,k,iEupS)**2
              predSumNCaO = TFEup*(mpredNca + fpredNcaO * GF%zoop_force(1,3,i,j,1))*Bio(i,k,iNCaO)**2
              predSumEupO = TFEup*(mpredEup + fpredEupO * GF%zoop_force(1,5,i,j,1))*Bio(i,k,iEupO)**2
# else
              predSumCop  = TFEup*(mpredCop)*Bio(i,k,iCop)**2
              predSumNCaS = TFEup*(mpredNca)*Bio(i,k,iNCaS)**2
              predSumEupS = TFEup*(mpredEup)*Bio(i,k,iEupS)**2
              predSumNCaO = TFEup*(mpredNca)*Bio(i,k,iNCaO)**2
              predSumEupO = TFEup*(mpredEup)*Bio(i,k,iEupO)**2
# endif

              DBio(i,k,iCop)  = DBio(i,k,iCop)  - predSumCop  * dtdays
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) - predSumNCaS * dtdays
              DBio(i,k,iEupS) = DBio(i,k,iEupS) - predSumEupS * dtdays
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) - predSumNCaO * dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupO) - predSumEupO * dtdays
#endif

#ifdef FEAST
# ifdef PROD3
              pt3(i,j,k,nnew,iQCopMort)  = predSumCop
              pt3(i,j,k,nnew,iQNCaSMort) = predSumNCaS
              pt3(i,j,k,nnew,iQEupSMort) = predSumEupS
              pt3(i,j,k,nnew,iQNCaOMort) = predSumNCaO
              pt3(i,j,k,nnew,iQEupOMort) = predSumEupO
# endif
#endif

#ifdef STATIONARY
!g            Stat3(i,k,8)=TFEup*mpredNCa*dtdays*Bio(i,k,iNCaS)**2
#endif

              !-----------------------------------
              !  Detritus from nonlinear mortality
              !-----------------------------------

              DBio(i,k,iDetF) = DBio(i,k,iDetF) + dtdays *              &
     &                          (predSumCop + predSumNCaS + predSumEupS &
     &                          + predSumNCaO + predSumEupO)

#ifdef STATIONARY

# ifdef fixedPRED
              Stat3(i,k,14) =                                           &
     &          + 0.5*Hz(i,j,k)/dtdays                                  &
     &          +0.5*Hz(i,j,k)/dtdays                                   &
     &          + 1*Hz(i,j,k)/dtdays                                    &
     &          0.5*Hz(i,j,k)/dtdays                                    &
     &          +1*Hz(i,j,k)/dtdays

# else

              Stat3(i,k,14) = dtdays *                                  &
     &                        (predSumCop + predSumNCaS + predSumEupS + &
     &                         predSumNCaO + predSumEupO)

# endif

#endif

#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iCop,iDetF)= bflx(iCop,iDetF)  + dtdays*xi*predSumCop
                bflx(iNcaS,iDetF)= bflx(iNCaS,iDetF)+ dtdays*xi*predSumNCaS
                bflx(iEupS,iDetF)= bflx(iEupS,iDetF)+ dtdays*xi*predSumEupS
                bflx(iNcaO,iDetF)= bflx(iNCaO,iDetF)+ dtdays*xi*predSumNCaO
                bflx(iEupO,iDetF)= bflx(iEupO,iDetF)+ dtdays*xi*predSumEupO
              END IF
#endif

#if defined JELLY
              DBio(i,k,iJel) = DBio(i,k,iJel) - mpredJel *              &
     &                         Bio(i,k,iJel)**2  * dtdays


# ifdef STATIONARY
!             Stat3(i,k,5) =mpredJel *  Bio(i,k,iJel)**2

# endif
              DBio(i,k,iDetF) = DBio(i,k,iDetF)                         &
     &                   + mpredJel * Bio(i,k,iJel)**2  * dtdays

# if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iJel,iDet)= bflx(iJel,iDet)+                       &
     &                           mpredJel*dtdays*Bio(i,k,iJel)*xi
              END IF
# endif

# ifdef STATIONARY
              IF (i.eq.3.and.j.eq.3) THEN
!g              Stat3(i,k,11)=mpredJel*Bio(i,k,iJel)**2  * dtdays

              END IF
# endif


#endif
            END DO
          END DO


          !=================================
          !Phytoplankton respiration losses
          !=================================

          DO k=1,N(ng)
            DO i=Istr,Iend

              BasalMet = respPhS

#ifdef IRON_LIMIT
              !----------------------------------------------
              !  Correct basal metabolism for iron limitation
              !----------------------------------------------

!             Iron1 = Bio(i,k,iFe)
!             respPh = respPhS
!             kfePh = kfePhS
!             BasalMet = GetBasalMetabolism(respPh,kfePh,Iron1)

#endif

              !-----------------------------------
              !  Arhonditsis temperature functions
              !-----------------------------------
!
              TempFuncPhS(i,k) = GetPhytoResp2(Temp1,TmaxPhS,KtBm_PhS)

              !-------------------------------------------------
              !  Change in concentration of Small Phytoplankton
              !-------------------------------------------------

              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &          TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)



              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &          xi * TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)


              !--------------------------------------------
              !  Primary production of Small phytoplankton
              !--------------------------------------------
#ifdef PROD3
              Prod(i,k,iPHSprd) = Prod(i,k,iPHSprd) -                   &
     &          TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)


#endif
#ifdef STATIONARY
              if(i.eq.2) THEN
                if(j.eq.2) THEN
    !             Stat3(i,k,1)=TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
                endif
              endif
#endif

              BasalMet = respPhL
#ifdef IRON_LIMIT

              !-----------------------------------------------
              !  Correct basal metabolism for iron limitation
              !-----------------------------------------------

!             respPh = respPhL
!             kfePh = kfePhL
!             BasalMet = GetBasalMetabolism(respPh,kfePh,Iron1)
#endif
              !------------------------------------
              !  Arhonditsis temperature functions
              !------------------------------------
              TempFuncPhL(i,k) = GetPhytoResp2(Temp1,TmaxPhL,KtBm_PhL)

              !-------------------------------------------------
              !  Change in concentration of Large Phytoplankton
              !-------------------------------------------------

              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &          TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
#ifdef STATIONARY

              if(i.eq.3) THEN
                if(j.eq.3) THEN
!                 Stat3(i,k,1)=TempFuncPhL(i,k)*BasalMet*dtdays!*Bio(i,k,iPhL)
                endif
              endif
#endif

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &           xi * TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)

              !--------------------------------------------
              !  Primary production of Large phytoplankton
              !--------------------------------------------

#ifdef PROD3
              Prod(i,k,iPHLprd) = Prod(i,k,iPHLprd) -                   &
     &          TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
#endif


#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iNH4)= bflx(iPhS,iNH4)+                       &
     &            TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)*xi

                bflx(iPhL,iNH4)= bflx(iPhL,iNH4)+                       &
     &            TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)*xi
              END IF
#endif

            END DO
          END DO

          !======================================
          !  Microzooplankton respiration losses
          !======================================
          DO k=1,N(ng)
            DO i=Istr,Iend

!             cff1 = fpPhSMZL * Bio(i,k,iPhS)+ fpPhLMZL * Bio(i,k,iPhL)

!             if(cff1.lt.1.0_r8)THEN
!               BasalMetMZL=0.0_r8
!             else
!               BasalMetMZL=respMZL
!             end if


              !-------------------------
              !  Small Microzooplankton
              !-------------------------

              !  Arhonditsis temperature functions

!             TempFuncMZS(i,k) = GetPhytoResp2(Temp1,TmaxMZS,           &
!    &                           KtBm_MZS)

!             BasalMet = respMZS
!             TFMZL = Q10MZS ** ( (Bio(i,k,itemp)-Q10MZST) / 10.0_r8 )
!
              !----------------------------------------------------
              !  Change in concentration of small microzooplankton
              !----------------------------------------------------

!             DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
!     &         TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)
!     &         TFMZS*respMZS*dtdays*Bio(i,k,iMZS)

              !------------------------------------
              !  Small Microzooplankton production
              !------------------------------------

#ifdef PROD3
!             Prod(i,k,iMZS) = Prod(i,k,iMZS) -                         &
!     &         TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)
!     &         TFMZS*respMZS*dtdays*Bio(i,k,iMZS)
#endif
              !-------------------------------------------------------------
              !  Add ammonium to correct for excretion related to metabolism
              !-------------------------------------------------------------

!             DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
!     &         xi*(TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS))
!     &         xi*(TFMZS*respMZS*dtdays*Bio(i,k,iMZS))

              !-------------------------
              !  Large Microzooplankton
              !-------------------------

              !  Arhonditsis temperature functions


              TFMZL = exp(KtBm_MZL * (Bio(i,k,itemp) - TmaxMZL))

              BasalMetMZL = respMZL

              !----------------------------------------------------
              !  Change in concentration of large microzooplankton
              !----------------------------------------------------
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
!    &          TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)
     &          TFMZL* BasalMetMZL*dtdays*Bio(i,k,iMZL)
#ifdef STATIONARY
!             if(i.eq.2) THEN
!               if(j.eq.2) THEN
!                 Stat3(i,k,9)=TFMZL* BasalMetMZL*dtdays
!               endif
!             endif
#endif

              !---------------------------------------
              !  Large Microzooplankton net production
              !---------------------------------------

#ifdef PROD3
              Prod(i,k,iMZLprd) = Prod(i,k,iMZLprd) -                   &
     &             TFMZL* BasalMetMZL*dtdays*Bio(i,k,iMZL)
#endif
              !-------------------------------------------------------------
              !  Add ammonium to correct for excretion related to metabolism
              !-------------------------------------------------------------

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
!    &          xi*(TempFuncMZL(i,k)*BasalMetMZL*dtdays*Bio(i,k,iMZL))
     &          xi*(TFMZL* BasalMetMZL*dtdays*Bio(i,k,iMZL))
#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
!               bflx(iMZS,iNH4)= bflx(iMZS,iNH4) + TempFuncMZS(i,k) *   &
!    &            BasalMetMZL*dtdays*Bio(i,k,iMZS)*xi
                bflx(iMZL,iNH4)= bflx(iMZL,iNH4) + xi * (TFMZL*         &
     &            BasalMetMZL*dtdays*Bio(i,k,iMZL))

              END IF
#endif

            END DO
          END DO

          !=====================================
          !  Mesozooplankton respiration losses
          !=====================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              cff1 = fpPhSCop * Bio(i,k,iPhS)+ fpPhLCop * Bio(i,k,iPhL) &
     &             + fpMZLCop * Bio(i,k,iMZL)



              if(cff1.lt.0.01_r8)THEN  !0.05
                BasalMetCop= respCop*cff1/0.01_r8
              else
                BasalMetCop= respCop
              endif

              !---------------------------------
              !  Copepod respiration correction
              !---------------------------------

              TFCop = exp(ktbmC * (Bio(i,k,itemp) - TrefC))

              !------------------------------------
              !  Neocalanus respiration correction
              !------------------------------------

              TFNCa = exp(ktbmN * (Bio(i,k,itemp) - TrefN))

              !-----------------------------------
              !  Euphausiid respiration correction
              !-----------------------------------

              TFEup = exp(ktbmE * (Bio(i,k,itemp) - TrefE))


              ! starvation response

              cff1 = fpPhSEup * Bio(i,k,iPhS)                           &
     &               + fpPhLEup * Bio(i,k,iPhL)                         &
!    &               + fpMZSEup * Bio(i,k,iMZS)                         &
     &               + fpMZLEup * Bio(i,k,iMZL)                         &
     &               + fpCopEup * Bio(i,k,iCop)

              cff2 = fpPhSNCa * Bio(i,k,iPhS)+ fpPhLNCa * Bio(i,k,iPhL) &
!    &               + fpMZSNCa * Bio(i,k,iMZS)                         &
     &               + fpMZLNCa * Bio(i,k,iMZL)




              if(cff1.lt.0.01_r8)THEN
                BasalMetEup=respEup*cff1/0.01_r8
              else
                BasalMetEup= respEup
              endif



              if(cff2.lt.0.01_r8)THEN
                BasalMetNC= respNC*cff2/0.01_r8
                BasalMetCM= respCM*cff2/0.01_r8
              else
                BasalMetNC=respNC
                BasalMetCM=respCM
              end if

              !---------------------------------------------------------
              !  Change in concentration from small copepod respiration
              !---------------------------------------------------------

              DBio(i,k,iCop) = DBio(i,k,iCop) -                         &
     &          BasalMetCop*TFCop*Bio(i,k,iCop)*dtdays


#ifdef PROD3
              Prod(i,k,iCopprd) = Prod(i,k,iCopprd) -                   &
     &          TFCop*BasalMetCop*Bio(i,k,iCop)*dtdays
#endif
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &          xi*(TFCop*BasalMetCop*dtdays*Bio(i,k,iCop))

              !----------------------------------------------------------
              !  Change in concentration from large copepods respiration
              !----------------------------------------------------------

              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) -                       &
     &            TFNCa*BasalMetCM*Bio(i,k,iNCaS)*dtdays

              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) -                       &
     &            TFNCa*BasalMetNC*Bio(i,k,iNCaO)*dtdays
#ifdef PROD3
              Prod(i,k,iNCaprd) = Prod(i,k,iNCaprd) -                   &
      &           TFNCa*BasalMetCM*Bio(i,k,iNCaS)*dtdays

              Prod(i,k,iNCaprd) = Prod(i,k,iNCaprd) -                   &
     &            TFNCa*BasalMetNC*Bio(i,k,iNCaO)*dtdays
#endif
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFNCa*BasalMetCM*dtdays*Bio(i,k,iNCaS))

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFNCa*BasalMetNC*dtdays*Bio(i,k,iNCaO))

              !------------------------------------------------------
              !  Change in concentration from euphausiid respiration
              !------------------------------------------------------

              DBio(i,k,iEupS) = DBio(i,k,iEupS) -                       &
     &            TFEup*BasalMetEup*Bio(i,k,iEupS)*dtdays

              DBio(i,k,iEupO) = DBio(i,k,iEupO) -                       &
     &            TFEup*BasalMetEup*Bio(i,k,iEupO)*dtdays
#ifdef PROD3
              Prod(i,k,iEupprd) = Prod(i,k,iEupprd) -                   &
     &            TFEup*BasalMetEup*Bio(i,k,iEupS)*dtdays
              Prod(i,k,iEupprd) = Prod(i,k,iEupprd) -                   &
     &            TFEup*BasalMetEup*Bio(i,k,iEupO)*dtdays
#endif
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFEup*BasalMetEup* dtdays*Bio(i,k,iEupS))

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFEup*BasalMetEup*dtdays*Bio(i,k,iEupO))


#ifdef STATIONARY
              Stat3(i,k,10)= TFNCa*BasalMetCM*Bio(i,k,iNCaS)*dtdays
#endif


#if defined JELLY

              BasalMetJel = respJel

              TFJel = Q10Jelr ** ( (Bio(i,k,itemp)-Q10JelTr)/10.0_r8)

!             TFJel =1.0_r8

              DBio(i,k,iJel) = DBio(i,k,iJel) -                         &
     &            TFJel *BasalMetJel* Bio(i,k,iJel)*dtdays

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            xi*(TFJel*BasalMetJel*dtdays*Bio(i,k,iJel))

# ifdef STATIONARY
!             Stat3(i,k,6) =  TFJel *BasalMetJel* Bio(i,k,iJel)

# endif

# ifdef STATIONARY
              IF (i.eq.3.and.j.eq.3) THEN
!               Stat3(i,k,12)= Ra !TFJel *BasalMetJel* Bio(i,k,iJel)*dtdays
              endif

# endif
# ifdef PROD3
              Prod(i,k,iJelprd) = Prod(i,k,iJelprd) -                   &
     &            TFJel*BasalMetJel*Bio(i,k,iJel)*dtdays
# endif


#endif

#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iCop,iNH4)  = bflx(iCop,iNH4)  +                   &
     &              TFCop*BasalMetCop*Bio(i,k,iCop) *dtdays*xi
                bflx(iNCaS,iNH4) = bflx(iNCaS,iNH4) +                   &
     &              TFNCa*BasalMetNCa*Bio(i,k,iNCaS)*dtdays*xi
                bflx(iNCaO,iNH4) = bflx(iNCaO,iNH4) +                   &
     &              TFNCa*BasalMetNCa*Bio(i,k,iNCaO)*dtdays*xi
                bflx(iEupS,iNH4) = bflx(iEupS,iNH4) +                   &
     &              TFEup*BasalMetEup*Bio(i,k,iEupS)*dtdays*xi
                bflx(iEupO,iNH4) = bflx(iEupO,iNH4) +                   &
     &              TFEup*BasalMetEup*Bio(i,k,iEupO)*dtdays*xi
                bflx(iJel,iNH4)  = bflx(iJel,iNH4)  +                   &
     &              TFJel*BasalMetJel*Bio(i,k,iJel) *dtdays*xi
              END IF
#endif
            END DO
          END DO

          !=============
          ! Molting:
          !=============

          !--------------------------------------------------------
          !  NOTE: It is unclear where molting equation came from.
          !  This is present only for euphausiids, not copepods
          !--------------------------------------------------------

!g        DO k=1,N(ng)
!g          DO i=Istr,Iend
!g            cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*     &
!g   &               Bio(i,k,iEupS)
!g            DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
!g            DBio(i,k,iEupS) = DBio(i,k,iEupS) - cff1 * dtdays
!g
!g            cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*     &
!g   &               Bio(i,k,iEupO)
!g            DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
!g            DBio(i,k,iEupO) = DBio(i,k,iEupO) - cff1 * dtdays
!g

!g#if defined BIOFLUX && defined BEST_NPZ
!g            IF (i.eq.3.and.j.eq.3) THEN
!g              bflx(iEupS,iDet)= bflx(iEupS,iDet) + cff1* dtdays*xi
!g              bflx(iEupO,iDet)= bflx(iEupO,iDet) + cff1* dtdays*xi
!g            END IF
!g#endif
!g          END DO
!g        END DO

          !==========================================
          ! Detrital Remineralization   (Det -> NH4)
          !==========================================

          DO k=1,N(ng)
            DO i=Istr,Iend

              Temp1 = Bio(i,k,itemp)
              !----------------------
              !  From Frost (1993).
              !----------------------
!             cff1 = regen * dgrad * Bio(i,k,iDet)
!             DBio(i,k,iNH4) = DBio(i,k,iNH4) + xi * cff1 * dtdays
!             DBio(i,k,iDet) = DBio(i,k,iDet) - cff1 * dtdays

              !-----------------------
              !  From Kawamiya(2000)
              !-----------------------
              PON = Bio(i,k,iDet)*xi  !Particulate organic nitrogen

              DBio(i,k,iDet) = DBio(i,k,iDet) -                         &
     &            ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            ((Pv0*exp(PvT*Temp1)*PON))*dtdays

#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iDet,iNH4)= bflx(iDet,iNH4) +                      &
     &            ((Pv0*exp(PvT*Temp1)*PON))*dtdays
              END IF
#endif

              PON = Bio(i,k,iDetF)*xi  !Particulate organic nitrogen

              DBio(i,k,iDetF) = DBio(i,k,iDetF) -                       &
     &            ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &            ((Pv0*exp(PvT*Temp1)*PON))*dtdays



#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iDetF,iNH4)= bflx(iDetF,iNH4) +                    &
     &            ((Pv0*exp(PvT*Temp1)*PON))*dtdays
              END IF
#endif
            END DO
          END DO
          !=============================
          !Nitrification  (NH4 -> NO3)
          !=============================
          DO k=1,N(ng)
            DO i=Istr,Iend


              !-----------------------
              !Temperature dependance
              !-----------------------

              ! Kawamiya 2000  NitrMax=Nitr0*exp(KnT*Temp1) - Ken

!             NitrifMax=GetNitrifMaxK(Temp1,KnT,Nitr0)

              ! Arhonditsis NitrMax=Nitr0*exp(-ktntr*(Temp1 - ToptNtr)**2)

              NitrifMax=GetNitrifMaxA(Nitr0, ktntr,Temp1,ToptNtr)

              !  No temperaure effects - NitrMax is constant

!             NitrifMax=Nitr0


              !-------------------------
              !  Light/Depth dependance
              !-------------------------

              !  Fennel

              DLNitrif=GetNitrifLight(Par1,tI0,KI)

              ! Denman

!             DLNitrif = (z_wL(k)**10_r8)/( (20_r8**10_r8) +            &
!    &                    z_wL(k)**10_r8)

              !  No Depth/Ligt effects

!               DLNitrif= 1.0_r8


              !-------------
              !  Saturation
              !-------------

              ! Arhonditsis

              cff1 = Bio(i,k,iNH4)/(KNH4Nit +Bio(i,k,iNH4))

              !  No saturation -ken

!               cff1 =1.0_r8


              DBio(i,k,iNH4) = DBio(i,k,iNH4)  - NitrifMax *            &
     &            Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays

              DBio(i,k,iNO3) = DBio(i,k,iNO3) + NitrifMax *             &
     &            Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays

#if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iNH4,iNO3)= bflx(iNH4,iNO3) +                      &
     &            NitrifMax  * Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays
              END IF
#endif
            END DO
          END DO

#ifdef BENTHIC
          !=================
          !Benthic Sub Model
          !=================

          DO i=Istr,Iend

            !----------------------------
            !  Growth of Benthic Infauna
            !----------------------------

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

              totD  = totD  + Bio(i,k,iDet)  * Hz(i,j,k) * frac1(i,k)
              totDF = totDF + Bio(i,k,iDetF) * Hz(i,j,k) * frac1(i,k)
              totPS = totPS + Bio(i,k,iPhS)  * Hz(i,j,k) * frac1(i,k)
              totPL = totPL + Bio(i,k,iPhL)  * Hz(i,j,k) * frac1(i,k)

              cff2 = cff2 + Hz(i,j,k)
            END DO

            ! Potential food available from water column

            cff1=(prefD *totD /((prefD *totD )+LupP))*prefD *totD
            cff2=(prefD *totDF/((prefD *totDF)+LupP))*prefD *totDF
            cff3=(prefPS*totPS/((prefPS*totPS)+LupP))*prefPS*totPS
            cff4=(prefPL*totPL/((prefPL*totPL)+LupP))*prefPL*totPL

            cff6 = cff1+cff2+cff3+cff4 ! Total pelagic food

            ! Potential food available from  sea floor

            cff5 = (prefD * BioB(i,k,iBenDet) / (prefD *                &
     &             BioB(i,k,iBenDet) + LupD)) * prefD * BioB(i,k,iBenDet)

            ! Uptake rates mediated by bottom layer temperature

             k = 1
             Temp1 = Bio(i,k,itemp)
             cff0 = q10r**((Temp1-T0benr)/10.0_r8)

            !  uptake of each food category

            cff7  = min(cff1,(cff0*cff1*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff8  = min(cff2,(cff0*cff2*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff9  = min(cff3,(cff0*cff3*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff10 = min(cff4,(cff0*cff4*BioB(i,k,iBen)*Rup/(cff6+KupP)))
            cff11 = min(cff5,(cff0*cff5*BioB(i,k,iBen)*Rup/(cff5+KupD)))

            ! Addition to benthos

            DBioB(i,1,iBen) = DBioB(i,1,iBen) + (cff7 +cff8+cff9+cff10+cff11) * dtdays

            ! Feeding losses from appropriate layers

            DBioB(i,k,iBenDet) = DBioB(i,k,iBenDet) - dtdays*cff11

            DO k = 1,N(ng)

              DBio(i,k,iDet)  = DBio(i,k,iDet)  - cff7  * frac2(i,k) * dtdays / Hz(i,j,k)
              DBio(i,k,iDetF) = DBio(i,k,iDetF) - cff8  * frac2(i,k) * dtdays / Hz(i,j,k)
              DBio(i,k,iPhS)  = DBio(i,k,iPhS)  - cff9  * frac2(i,k) * dtdays / Hz(i,j,k)
              DBio(i,k,iPhL)  = DBio(i,k,iPhL)  - cff10 * frac2(i,k) * dtdays / Hz(i,j,k)

            END DO

# if defined BIOFLUX && defined BEST_NPZ
            IF (i.eq.3.and.j.eq.3) THEN
              bflx(NT(ng)+2,NT(ng)+1)=bflx(NT(ng)+2,NT(ng)+1)           &
     &            +dtdays*cff9*xi
            ENDIF
# endif

            !--------------------
            ! Benthic Production
            !--------------------
!
# ifdef PROD2
            Prod2(i,iBenPrd)=Prod2(i,iBenPrd) + (cff11)*dtdays
            Prod2(i,iBenPrd)=Prod2(i,iBenPrd) + (cff7+cff8+cff9+cff10)*dtdays
# endif

# if defined BIOFLUX && defined BEST_NPZ
            IF (i.eq.3.and.j.eq.3) THEN

              bflx(iDet,NT(ng)+1) = bflx(iDet,NT(ng)+1) +               &
    &             dtdays*cff6/Hz(i,j,k)*xi
              bflx(iPhL,NT(ng)+1) = bflx(iPhL,NT(ng)+1) +               &
    &             dtdays*cff7/Hz(i,j,k)*xi
              bflx(iPhS,NT(ng)+1) = bflx(iPhS,NT(ng)+1) +               &
    &             dtdays*cff8/Hz(i,j,k)*xi

            END IF
# endif
            !------------
            !  Excretion
            !------------

            ! Assume all excretion occurs in the bottom layer

            k = 1 ! I'm trying to hard-code the bottom layer thing, but just in case I missed some...
            cff1=cff7 *eexD
            cff2=cff8 *eexD
            cff3=cff9 *eex
            cff4=cff10*eex
            cff5=cff11*eexD

            DBioB(i,1,iBen) = DBioB(i,1,iBen)                             &
     &                        -(cff1+cff2+cff3+cff4+cff5)*dtdays


            !  Material excreted is considered inorganic and not available
            !  for further secondary production

            DBio(i,1,iNH4)= DBio(i,1,iNH4)+ xi * dtdays * 0.5_r8        &
     &         *(cff1+cff2+cff3+cff4+cff5)/Hz(i,j,1)

            DBioB(i,1,iBenDet)= DBioB(i,1,iBenDet)+dtdays               &
     &         *(cff1+cff2+cff3+cff4+cff5)*0.5_r8


# if defined BIOFLUX && defined BEST_NPZ
            IF (i.eq.3.and.j.eq.3) THEN
              bflx(NT(ng)+1,NT(ng)+2)=bflx(NT(ng)+1,NT(ng)+2)           &
     &          + (cff6+cff7+cff8+cff9+cff10)*dtdays*xi
            ENDIF
# endif
            !--------------
            !  Respiration
            !--------------

            cff3=cff0*BioB(i,1,iBen)*Rres
            cff4=Qres*(((1_r8-eexD)*cff7) + ((1_r8-eexD)*cff8) +        &
    &            ((1_r8-eex)*cff9) + ((1_r8-eex)*cff10) +               &
    &            ((1_r8-eexD)*cff11))



            cff6=cff3+cff4

            DBioB(i,1,iBen)=DBioB(i,1,iBen) -cff6*dtdays
            DBio(i,1,iNH4)= DBio(i,1,iNH4)+ xi*dtdays*cff6/Hz(i,j,1)


# ifdef STATIONARY2
            Stat2(i,4)=cff3*dtdays
            Stat2(i,5)=cff4*dtdays

# endif


# if defined BIOFLUX && defined BEST_NPZ
            IF (i.eq.3.and.j.eq.3) THEN
              bflx(NT(ng)+1,iNH4)=bflx(NT(ng)+1,iNH4)                   &
     &          +xi*dtdays*cff6/Hz(i,j,1)
            ENDIF
# endif
            !---------------------
            !  Benthic Production
            !---------------------
# ifdef PROD2
            Prod2(i,iBenPrd) = Prod2(i,iBenPrd) - cff6*dtdays

# endif
            !------------
            !  Mortality
            !------------

            cff1= rmort*BioB(i,k,iBen)*cff0
            DBioB(i,k,iBen) = DBioB(i,k,iBen)- cff1 * dtdays
            DBioB(i,k,iBenDet)= DBioB(i,k,iBenDet)+cff1*dtdays
# ifdef STATIONARY2

            Stat2(i,6)= cff1 * dtdays
# endif
# if defined BIOFLUX && defined BEST_NPZ
            IF (i.eq.3.and.j.eq.3) THEN
              bflx(NT(ng)+1,NT(ng)+2)=bflx(NT(ng)+1,NT(ng)+2)           &
     &          + cff1*dtdays*xi
            ENDIF
# endif

            !------------
            !  Predation
            !------------

            DBioB(i,k,iBen) =DBioB(i,k,iBen)                            &
     &           -cff0*BenPred*dtdays*BioB(i,k,iBen)**2

            DBioB(i,k,iBenDet) =DBioB(i,k,iBenDet)                      &
     &           +cff0*BenPred*dtdays*BioB(i,k,iBen)**2


# ifdef STATIONARY2

            Stat2(i,7)=cff0*BenPred*dtdays*BioB(i,k,iBen)**2
# endif

            !-------------------------------------
            !  (Det -> NH4) temperature dependent
            !-------------------------------------

            ! Assume only the top 25% available -> NH4

            PON = BioB(i,k,iBenDet)*0.25*xi/Hz(i,j,k)  !Benthic Particulate organic nitrogen

!           cff5= q10**((Temp1-T0ben)/10.0_r8)
            cff1=Pv0*exp(PvT*Temp1)*PON       !-Kawamiya 2000
!           cff1=Pv0*cff5*PON

            !  Arhonditsis

            cff2 = Bio(i,k,iNH4)/(KNH4Nit +Bio(i,k,iNH4))

            DBioB(i,k,iBenDet) =DBioB(i,k,iBenDet)                      &
     &           - (cff1/xi)*dtdays*Hz(i,j,k)

            DBio(i,k,iNH4)= DBio(i,k,iNH4)                              &
     &           + (cff1)*dtdays


# if defined BIOFLUX && defined BEST_NPZ
            IF (i.eq.3.and.j.eq.3) THEN
              bflx(NT(ng)+2,iNH4)=bflx(NT(ng)+2,iNH4)                   &
     &            +  iremin*(cff1)*dtdays
            ENDIF
# endif

          END DO

#endif



          !
          !==============================================================
          ! Vertical Sinking of Particles(Phyto and Det)
          !==============================================================
          !
          ! Use Sinking Code adapted from J. Warner ROMS sediment code
          !   - this is similar to the approach taken in Fasham
          !     biology code but adapted to be used in subroutine
          !   - Also addapted to return particle flux in and out of each
          !     level - used in BENTHIC sub model
          !   - Incorporated Liz Dobbins zlimit to determine sinking rate
          !     zlimit =1 for constant sinking - use for Phyto and Det
          !     zlimit =-1*NcmaxZ for Neocalanus- sink rate is then attenuated
          !     as max sinking depth is approached
          !
          ! G.Gibson  July 2008
          ! K.Kearney March 2016: updated to correct mass-accumulation bug

          ! Note: all sinking and rising uses the subroutine TracerSink,
          ! which modifies Btmp and sinkout

          ! Initialize temporary arrays to 0

          Btmp = 0    ! tracer profile, (1:n)
          Hztmp = 0   ! layer thickness, (1:n)
          zwtmp = 0   ! layer edges, (0:n)
          sinkout2 = 0 ! loss out of bottom cell (1)

          ! Small phytoplankton: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call TracerSink(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wPhS, 1.0_r8, Bio(i,:,iPhS), Btmp, sinkout2)

            DO k = 1,N(ng)
              DBio(i,k,iPhS) = DBio(i,k,iPhS) + (Btmp(k) - Bio(i,k,iPhS))
            END DO
            DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8

          END DO


          ! Large phytoplankton: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call TracerSink(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wPhL, 1.0_r8, Bio(i,:,iPhL), Btmp, sinkout2)

            DO k = 1,N(ng)
              DBio(i,k,iPhL) = DBio(i,k,iPhL) + (Btmp(k) - Bio(i,k,iPhL))
            END DO
            DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8

          END DO

          ! Slow-sinking detritus: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call TracerSink(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wDet, 1.0_r8, Bio(i,:,iDet), Btmp, sinkout2)

            DO k = 1,N(ng)
              DBio(i,k,iDet) = DBio(i,k,iDet) + (Btmp(k) - Bio(i,k,iDet))
            END DO
            DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8

          END DO

          ! Fast-sinking detritus: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call TracerSink(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wDetF, 1.0_r8, Bio(i,:,iDetF), Btmp, sinkout2)

            DO k = 1,N(ng)
              DBio(i,k,iDetF) = DBio(i,k,iDetF) + (Btmp(k) - Bio(i,k,iDetF))
            END DO
            DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2*0.79_r8

          END DO

          ! On-shelf large copepods (NCaS i.e. CM): Move up and down
          ! based on dates set in input file.  Stop at either 200m or the
          ! water depth, whichever is shallower.  No biomass should cross
          ! the bottom or surface boundary.

          DO i=Istr,Iend

            if (downwardCM) then

              call TracerSink(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wNCsink, max((z_w(i,j,0)+z_w(i,j,1))/2, -200.0_r8), Bio(i,:,iNCaS), Btmp, sinkout2)

              DO k = 1,N(ng)
                DBio(i,k,iNCaS) = DBio(i,k,iNCaS) + (Btmp(k) - Bio(i,k,iNCaS))
              END DO

            else if (upwardCM) then

              call TracerRise(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wNCrise, (z_w(i,j,N(ng)-1) + z_w(i,j,N(ng)))/2, Bio(i,:,iNCaS), Btmp, sinkout2)

              DO k = 1,N(ng)
                DBio(i,k,iNCaS) = DBio(i,k,iNCaS) + (Btmp(k) - Bio(i,k,iNCaS))
              END DO

            end if

          END DO

          ! Off-shelf large copepods (NCaO i.e. NC): Move up and down
          ! based on dates set in input file.  Diapause to 400 m.
          ! If water is shallower than 400m, assume they die when they
          ! hit the bottom, and transfer biomass to benthic detritus.

          DO i=Istr,Iend

            if (downwardNC) then

              call TracerSink(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wNCsink, -400.0_r8, Bio(i,:,iNCaO), Btmp, sinkout2)

              DO k = 1,N(ng)
                DBio(i,k,iNCaO) = DBio(i,k,iNCaO) + (Btmp(k) - Bio(i,k,iNCaO))
              END DO

              DBio(i,1,iBenDet) = DBioB(i,1,iBenDet) + sinkout2

            else if (upwardNC) then

              call TracerRise(N(ng), Hz(i,j,:), z_w(i,j,:), dtdays, wNCrise, (z_w(i,j,N(ng)-1) + z_w(i,j,N(ng)))/2, Bio(i,:,iNCaO), Btmp, sinkout2)

              DO k = 1,N(ng)
                DBio(i,k,iNCaO) = DBio(i,k,iNCaO) + (Btmp(k) - Bio(i,k,iNCaO))
              END DO

            end if

          END DO



          !==============================================================
          !Ice Sub Model
          !==============================================================
!
#ifdef ICE_BIO

          DO i=Istr,Iend

# ifdef CLIM_ICE_1D

            IF (itL(i,j,nstp,iIceLog).gt.0.0) THEN   !0.02
              Temp1 = Bio(i,N(ng),itemp)

# elif defined BERING_10K

            IF (IceLog(i,j,nstp).gt.0.0_r8) THEN
              Temp1 = ti(i,j,nstp)
# endif

              Sal1 = Bio(i,N(ng),isalt)
              Par1 = PARs(i)

              !----------------------
              ! Growth of Ice  Algae
              !----------------------

              ! light limitation

              aiceIfrac=(1-exp(-alphaIb*Par1))*exp(-betaI*Par1)

! nutrient limitation

# ifdef CLIM_ICE_1D
              cff1=BioBI(i,iIceNO3)/(ksnut1+BioBI(i,iIceNO3))
              cff2=BioBI(i,iIceNH4)/(ksnut2+BioBI(i,iIceNH4))

              aiceNfrac=cff1*exp(-inhib*BioBI(i,iIceNH4))+cff2
              fNO3=(cff1*exp(-inhib*BioBI(i,iIceNH4)))/aiceNfrac

# elif defined BERING_10K
              cff1=IceNO3(i,j,nstp)/(ksnut1+IceNO3(i,j,nstp))
              cff2=IceNH4(i,j,nstp)/(ksnut2+IceNH4(i,j,nstp))

              aiceNfrac=cff1*exp(-inhib*IceNH4(i,j,nstp))+cff2
              fNO3=(cff1*exp(-inhib*IceNH4(i,j,nstp)))/aiceNfrac
# endif


              ! salinity impact (gesi) on ice algal growth

# ifdef BERING_10K
              ! use ice temperature to determine brine salinity (sb) -from Arrigo 1993

              if (ti(i,j,nstp)> -22.9_r8) THEN
                cff1=-3.9921
                cff2=-22.7
                cff3=-1.0015
                cff4=-0.019956
              else if (ti(i,j,nstp)> -44.0_r8  .AND. ti(i,j,nstp)< -22.9_r8) THEN
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

              sb=cff1+cff2* ti(i,j,nstp)+cff3*ti(i,j,nstp)**2           &
     &           +cff4* ti(i,j,nstp)**3



              gesi = max(0.0_r8,(1.1e-2+3.012e-2*sb                     &
     &             +1.0342e-3*sb**2                                     &
     &             -4.6033e-5*sb**3                                     &
     &             +4.926e-7*sb**4                                      &
     &             -1.659e-9*sb**5              ))

# else

              gesi=1.0_r8
# endif

              grow1=mu0*exp(0.0633*Temp1)

              ! Alternate temperature dependance code - need to look into more

              ! freezing point of seawater at sea level

!             cff1=-0.0575
!             cff2=1.710523E-3
!             cff3=-2.154996E-4
!             cff4=cff1*Sal1+cff2*Sal1**(3/2)+cff3*Sal1**2

!             grow1=mu0*exp(cff4-Temp1)


# ifdef STATIONARY2

!             Stat2(i,1)= Temp1
!             Stat2(i,2)= sb
!             Stat2(i,3)= gesi


# endif


              GROWAice=grow1*min(aiceNfrac,aiceIfrac)*gesi

# ifdef CLIM_ICE_1D
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &         + GROWAice*BioBI(i,iIcePhL)* dtdays

#  ifdef PROD2
              Prod2(i,iIAPrd) = Prod2(i,iIAPrd)                         &
     &         + GROWAice*BioBI(i,iIcePhL)* dtdays
#  endif
# elif defined BERING_10K
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &         + GROWAice*IcePhL(i,j,nstp)* dtdays

#  ifdef PROD2
              Prod2(i,iIAPrd) = Prod2(i,iIAPrd)                         &
     &         + GROWAice*IcePhL(i,j,nstp)* dtdays
#  endif
# endif

              !------------------------------
              !  Respiration of Ice Algae
              !------------------------------
!             RAi0=R0i*GROWAice
              RAi0=R0i*mu0*exp(0.0633*Temp1)
# ifdef CLIM_ICE_1D
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &         -BioBI(i,iIcePhL)*RAi0*dtdays

#  ifdef PROD2
              Prod2(i,iIAPrd) = Prod2(i,iIAPrd)                         &
     &         -BioBI(i,iIcePhL)*RAi0*dtdays
#  endif
#  ifdef STATIONARY2


#  endif
# elif defined BERING_10K
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &         -IcePhL(i,j,nstp)*RAi0*dtdays
#  ifdef PROD2
              Prod2(i,iIAPrd) = Prod2(i,iIAPrd)                         &
     &         -IcePhL(i,j,nstp)*RAi0*dtdays
#  endif
# endif

              ! ----------------------------------------
              !  mortality of Ice Algae
              ! ----------------------------------------
              RgAi=rg0*exp(rg*Temp1)

# ifdef CLIM_ICE_1D
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &         -BioBI(i,iIcePhL)*RgAi*dtdays

              DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                       &
     &         +BioBI(i,iIcePhL)*RgAi*dtdays*xi

# elif defined BERING_10K
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &          -IcePhL(i,j,nstp)*RgAi*dtdays

              DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                       &
     &          +IcePhL(i,j,nstp)*RgAi*dtdays*xi
# endif


              !--------------
              ! Nitrification
              !--------------


# ifdef CLIM_ICE_1D
              reN=annit*BioBI(i,iIceNH4)

# elif defined BERING_10K

              reN=annit*IceNH4(i,j,nstp)

# endif
              DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)+reN*dtdays
              DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)-reN*dtdays


# ifdef CLIM_ICE_1D
              cff2=BioBI(i,iIcePhL)*(GROWAice-RAi0)

# elif defined BERING_10K

              cff2=IcePhL(i,j,nstp)*(GROWAice-RAi0)
# endif

              DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)                       &
     &            -fNO3*cff2*xi*dtdays


              DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                       &
     &            -(1.0_r8-fNO3)*cff2*xi*dtdays

              !------------------------------
              ! Ice/water convective exchange
              !------------------------------

              ! Convective exchange rate is based on a polynomial fit to
              ! ice growth/melting rate

# if defined CLIM_ICE_1D
              dhicedt=it(i,j,nnew,iIceZ)-it(i,j,nstp,iIceZ)
# elif defined BERING_10K
              dhicedt=hi(i,j,nstp)-hi(i,j,nnew)
# endif
              dhicedt=dhicedt*sec2day/dtdays !convert to m/s

              trs=9.667e-11+4.49e-6*dhicedt-1.39e-5*dhicedt**2
              trs=trs*86400   !convert to m/d
              twi=72*trs

              IF (dhicedt.lt.0) THEN

                trs=4.49e-6*ABS(dhicedt)-1.39e-5*ABS(dhicedt)**2
                trs=trs*86400
                twi=720*trs

              ENDIF

              ! IcePhL can get washed out of ice, but not in, so assume
              ! [PhL] = 0 for this exchange

# if defined CLIM_ICE_1D
              DBioBI(i,iIcePhL)  = DBioBI(i,iIcePhL)  + (twi * -BioBI(i,iIcePhL)/aidz)*dtdays
              DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL) + (twi * BioBI(i,iIcePhL)/Hz(i,j,N(ng)))*dtdays
# elif defined BERING_10K
              DBioBI(i,iIcePhL)  = DBioBI(i,iIcePhL)  + (twi * -IcePhL(i,j,nstp)/aidz)*dtdays
              DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL) + (twi * IcePhL(i,j,nstp)/Hz(i,j,N(ng)))*dtdays
# endif

              ! NO3 and NH4 have two-way exchange, based on gradient
              ! across the ice/water interface

# if defined CLIM_ICE_1D
              DBioBI(i,iIceNO3)  = DBioBI(i,iIceNO3)  + (twi * (Bio(i,N(ng),iNO3) - BioBI(i,iIceNO3))/aidz)*dtdays
              DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3) + (twi * (BioBI(i,iIceNO3)) - Bio(i,N(ng),iNO3)/Hz(i,j,N(ng)))*dtdays

              DBioBI(i,iIceNH4)  = DBioBI(i,iIceNH4)  + (twi * (Bio(i,N(ng),iNH4) - BioBI(i,iIceNH4))/aidz)*dtdays
              DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4) + (twi * (BioBI(i,iIceNH4)) - Bio(i,N(ng),iNH4)/Hz(i,j,N(ng)))*dtdays

# elif defined BERING_10K
              DBioBI(i,iIceNO3)  = DBioBI(i,iIceNO3)  + (twi * (Bio(i,N(ng),iNO3) - IceNO3(i,j,nstp))/aidz)*dtdays
              DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3) + (twi * (IceNO3(i,j,nstp)) - Bio(i,N(ng),iNO3)/Hz(i,j,N(ng)))*dtdays

              DBioBI(i,iIceNH4)  = DBioBI(i,iIceNH4)  + (twi * (Bio(i,N(ng),iNH4) - IceNH4(i,j,nstp))/aidz)*dtdays
              DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4) + (twi * (IceNH4(i,j,nstp)) - Bio(i,N(ng),iNH4)/Hz(i,j,N(ng)))*dtdays

# endif

              ! TODO: add bflux associated with ice/water exchange

!               IF (twi.lt.0) THEN
!
! # if defined CLIM_ICE_1D
!                 DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                     &
!      &            + BioBI(i,iIcePhL)*(twi/aidz)*dtdays
!
! #  ifdef STATIONARY2
!
! !               Stat2(i,6)= BioBI(i,iIcePhL)*(twi/aidz)
! !               Stat2(i,7)= BioBI(i,iIcePhL)*twi/Hz(i,j,N(ng))
! #  endif
!
!                 DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL)                 &
!      &            -twi*BioBI(i,iIcePhL)*dtdays/Hz(i,j,N(ng))
!
! # elif defined BERING_10K
!
!                 DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                     &
!      &            + IcePhL(i,j,nstp)*(twi/aidz)*dtdays
!
!                 DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL)                 &
!      &            -twi*IcePhL(i,j,nstp)*dtdays/Hz(i,j,N(ng))
!
! # endif
!
!
!
!               ENDIF
! ! nutrient gradient between ice and water
! # if defined CLIM_ICE_1D
!               cff1=twi*(Bio(i,N(ng),iNO3)-BioBI(i,iIceNO3))*dtdays
!               cff2=twi*(Bio(i,N(ng),iNH4)-BioBI(i,iIceNH4))*dtdays
!
!               IF (twi.lt.0) THEN
!                 DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)                     &
!      &                      + BioBI(i,iIceNO3)*twi*dtdays
!
!
!                 DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                     &
!      &                      + BioBI(i,iIceNH4)*twi*dtdays
!
! # elif defined BERING_10K
!
!               cff1=twi*(Bio(i,N(ng),iNO3)-IceNO3(i,j,nstp))*dtdays
!               cff2=twi*(Bio(i,N(ng),iNH4)-IceNH4(i,j,nstp))*dtdays
!
!               IF (twi.lt.0) THEN
!                 DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)                     &
!      &                      + IceNO3(i,j,nstp)*twi*dtdays
!
!
!                 DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                     &
!      &                      + IceNH4(i,j,nstp)*twi*dtdays
! !               Stat2(i,1)= IceNO3(i,j,nstp)*twi*dtdays
! # endif
!
!                 DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)+cff1/Hz(i,j,N(ng))
!                 DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)+cff2/Hz(i,j,N(ng))
!
!
! # if defined CLIM_ICE_1D
! #  if defined BIOFLUX && defined BEST_NPZ
!
! #  endif
! # endif
!               ELSE IF (twi.gt.0) THEN
!                 if(cff1.gt.0_r8)THEN
! # if defined BERING_10K
!                   if(( IceNO3(i,j,nstp)+cff1/aidz).lt.Bio(i,N(ng),iNO3))THEN
!
! # elif defined CLIM_ICE_1D
!                   if((BioBI(i,iIceNO3)+cff1/aidz).lt.Bio(i,N(ng),iNO3))THEN
! # endif
!
!
!                     DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)+cff1/aidz
!                     DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)-cff1/Hz(i,j,N(ng))
!                   endif
!                 else if(cff1.lt.0_r8)THEN
! # if defined BERING_10K
!                   if((Bio(i,N(ng),iNO3)-cff1/Hz(i,j,N(ng))).lt.IceNO3(i,j,nstp))THEN
! # elif defined CLIM_ICE_1D
!                   if((Bio(i,N(ng),iNO3)-cff1/Hz(i,j,N(ng))).lt.BioBI(i,iIceNO3))THEN
! # endif
!                     DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)+cff1/aidz
!                     DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)-cff1/Hz(i,j,N(ng))
!                   endif
!                 endif
!
!                 if(cff2.gt.0_r8)THEN
! # if defined BERING_10K
!                   if((IceNH4(i,j,nstp)+cff2/aidz).lt.Bio(i,N(ng),iNH4))THEN
! # elif defined CLIM_ICE_1D
!                   if((BioBI(i,iIceNH4)+cff2/aidz).lt.Bio(i,N(ng),iNH4))THEN
! # endif
!                     DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)+cff2/aidz
!                     DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)-cff2/Hz(i,j,N(ng))
!                   endif
!                 else if(cff2.lt.0_r8)THEN
! # if defined BERING_10K
!                  if((Bio(i,N(ng),iNH4)-cff2/Hz(i,j,N(ng))).lt.IceNH4(i,j,nstp))THEN
! # elif defined CLIM_ICE_1D
!                  if((Bio(i,N(ng),iNH4)-cff2/Hz(i,j,N(ng))).lt.BioBI(i,iIceNH4))THEN
! # endif
!                     DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)+cff2/aidz
!                     DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNO3)-cff2/Hz(i,j,N(ng))
!                   endif
!                 endif
! !    endif
!
!
! # if defined CLIM_ICE_1D
! #  if defined BIOFLUX && defined BEST_NPZ
!                 IF (i.eq.3.and.j.eq.3) THEN
!                   bflx(iNO3,NT(ng)+4)=bflx(iNO3,NT(ng)+4)               & !NO3->NO3i
!      &                  + cff1/aidz
!
!                   bflx(iNH4,NT(ng)+5)=bflx(iNH4,NT(ng)+5)               & !NH4->NH4i
!      &                  + cff2/aidz
!                 ENDIF
! #  endif
! # endif
!               ENDIF
!
! !             DBioBI(i,iIcePhL)=10_r8
! !           ENDIF

# if defined CLIM_ICE_1D
            ELSE IF(itL(i,j,nstp,iIceLog).lt.0_r8.and.itL(i,j,nnew,iIceLog).gt.0_r8) THEN

              DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL)                   &
     &               + BioBI(i,iIcePhL)*aidz/Hz(i,j,N(ng))

              DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)                   &
     &               + BioBI(i,iIceNO3)*aidz/Hz(i,j,N(ng))
              DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)                   &
     &               + BioBI(i,iIceNH4)*aidz/Hz(i,j,N(ng))


              DBioBI(i,iIceNO3)=-BioBI(i,iIceNO3)
              DBioBI(i,iIceNH4)=-BioBI(i,iIceNH4)
              DBioBI(i,iIcePhL)=-BioBI(i,iIcePhL)
#  ifdef PROD2
              Prod2(i,iIAPrd) = 0_r8
#  endif



# elif defined BERING_10K

            ELSE IF (IceLog(i,j,nstp).lt.0_r8.and.IceLog(i,j,nnew).gt.0.0_r8) THEN

              DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL)                   &
     &               + IcePhL(i,j,nnew)*aidz/Hz(i,j,N(ng))

              if((Bio(i,N(ng),iNO3)+IceNO3(i,j,nnew)).le.IceNO3(i,j,nnew))THEN


                DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)                 &
     &               + IceNO3(i,j,nnew)*aidz/Hz(i,j,N(ng))
              endif

              if((Bio(i,N(ng),iNH4)+IceNH4(i,j,nnew)).le.IceNH4(i,j,nnew))THEN


                DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)                 &
     &               + IceNH4(i,j,nnew)*aidz/Hz(i,j,N(ng))
              endif

              DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)-IceNO3(i,j,nstp)
              DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)-IceNH4(i,j,nstp)
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)-IcePhL(i,j,nstp)
#  ifdef PROD2
              Prod2(i,iIAPrd) = 0_r8
#  endif

# endif
# if defined CLIM_ICE_1D
#  if defined BIOFLUX && defined BEST_NPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(NT(ng)+3,iPhL)=bflx(NT(ng)+3,iPhL)                 & !PLi->PhL
     &            + it(i,j,nnew,iIcePhL)*xi/aidz

                bflx(NT(ng)+4,iNO3)=bflx(NT(ng)+4,iNO3)                 & !NO3i->NO3
     &            + it(i,j,nnew,iIceNO3)/aidz

                bflx(NT(ng)+5,iNH4)=bflx(NT(ng)+5,iNH4)                 & !NH4i->NH4
     &            + it(i,j,nnew,iIceNH4)/aidz
              ENDIF
#  endif
# endif

            ELSE
              !no recent hitory of ice so icebio do nothing
# if defined CLIM_ICE_1D

              DO itrc=1,3 !NIB
                ibioBI=idice(itrc)
                DBioBI(i,ibioBI)=0_r8

              END DO
#  ifdef PROD2
              Prod2(i,iIAPrd) = 0_r8
#  endif
# elif defined BERING_10K
              DBioBI(i,iIceNO3)=0_r8
              DBioBI(i,iIceNH4)=0_r8
              DBioBI(i,iIcePhL)=0_r8
#  ifdef PROD2
              Prod2(i,iIAPrd) = 0_r8
#  endif
# endif

            ENDIF

            if(i.eq.3)THEN
!              Stat2(3,7)=  DBioBI(i,iIcePhL)
            endif

          END DO

#endif


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
#include "feast_step.h"
#endif

      RETURN
      END SUBROUTINE biology_tile

!=====================================================================
! BIOSINK_1  particle sinking subroutine After J. Warner sed sink code
! G. Gibson July 2008
!=====================================================================
      subroutine BIOSINK_1(ng,wBio,Bio,sinkIN,sinkOUT,HzL,dtdays,z_wL,  &
     &                       zlimit,LBi,UBi,IminS, ImaxS)

      USE mod_param
!
      implicit none
!

      integer, intent(in)     :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in)    :: wBio
      real(r8), intent(in)    :: zlimit

      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(inout) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))

      real(r8), intent(in)    :: dtdays


      real(r8), intent(out) :: sinkIN(IminS:ImaxS,N(ng))
      real(r8), intent(out) :: sinkOUT(IminS:ImaxS,N(ng))

      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff


      real(r8):: dBio(0:N(ng)), wBiod(LBi:UBi,0:N(ng))
      real(r8) :: FC(IminS:ImaxS,0:N(ng))

      real(r8) :: Hz_inv(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv2(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv3(IminS:ImaxS,N(ng))


      integer :: ksource(IminS:ImaxS,N(ng))

      real(r8) :: qR(IminS:ImaxS,N(ng))
      real(r8) :: qL(IminS:ImaxS,N(ng))
      real(r8) :: WL(IminS:ImaxS,N(ng))
      real(r8) :: WR(IminS:ImaxS,N(ng))
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc




      IF ( zlimit .lt. 0 ) THEN
       DO k=0,N(ng)
         DO i=LBi, UBi
          IF ( z_wL(i,k) .ge. zlimit ) THEN
            wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2/ &
     &        (zlimit/2)**2 )
          ELSE
            wBiod(i,k) = 0.0_r8
          END IF
        END DO
       END DO
      ELSE
        DO k=0,N(ng)
          DO i=LBi, UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
!
!
!  Compute inverse thickness to avoid repeated divisions.
!

      DO k=1,N(ng)
       DO i=LBi,UBi
         !if (HzL(i,k)==0) then
         !   print *,i,k, HzL(i,k)
         !end if
         Hz_inv(i,k)=1.0_r8/HzL(i,k)
       END DO
      END DO
      DO k=1,N(ng)-1
       DO i=LBi,UBi
         Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=2,N(ng)-1
       DO i=LBi,UBi
         Hz_inv3(i,k)=1.0_r8/(HzL(i,k-1)+HzL(i,k)+HzL(i,k+1))
       END DO
      END DO

      DO k=1,N(ng)
         DO i=LBi,UBi
           qc(i,k)=Bio(i,k)
         END DO
      END DO

!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!

      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi

          FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
        END DO
      END DO

!     DO k=N(ng)-1,1,-1
!       DO i=LBi,UBi
!         print*,'LBi=',LBi,'UBi=',UBi
!         print*,'i=',i,'k=',k
!         if(i.le.UBi)THEN
!           FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
!         endif
!       END DO
!     END DO


      DO k=2,N(ng)-1
        DO i=LBi,UBi
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
      DO k=2,N(ng)-2
         DO i=LBi,UBi
          dltL=MAX(cff,WL(i,k  ))
          dltR=MAX(cff,WR(i,k+1))
          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
          qL(i,k+1)=qR(i,k)
        END DO
      END DO

      DO i=LBi,UBi
        FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
#if defined LINEAR_CONTINUATION
        qL(i,N(ng))=qR(i,N(ng)-1)
        qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
#elif defined NEUMANN
        qL(i,N(ng))=qR(i,N(ng)-1)
        qR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
#else
        qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
        qL(i,N(ng))=qc(i,N(ng))         ! conditions
        qR(i,N(ng)-1)=qc(i,N(ng))
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
      DO k=1,N(ng)
        DO i=LBi,UBi
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
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!

      DO k=1,N(ng)
        DO i=LBi,UBi
          cff=dtdays*ABS(wBiod(i,k))
          FC(i,k-1)=0.0_r8
          WL(i,k)=z_wL(i,k-1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO
      DO k=1,N(ng)
        DO ks=k,N(ng)-1
          DO i=LBi,UBi
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
      DO k=1,N(ng)
        DO i=LBi,UBi
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

          if (k.eq.N(ng)) then
            FC(i,k)=0.0_r8
          endif
        END DO
      END DO

      DO k=1,N(ng)
        DO i=LBi,UBi
!
!  The Bio variables are now updated in the main subroutine
!
!         Bio(i,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
          sinkIN(i,k)=FC(i,k)
          sinkOUT(i,k)=FC(i,k-1)

          if(k.eq.N(ng))THEN
            sinkIN(i,k)=0.0_r8
          endif

        END DO
      END DO

      RETURN
      END SUBROUTINE BIOSINK_1
!=====================================================================
! BIORISE_1  rising particle subroutine a reversal of J. Warner sed sink code
!      plus and attenuation of rise rate based on closeness to max sink depth
! G.Gibson July 2008
!=====================================================================
     subroutine BIORISE_1(ng,wBio,Bio,rOUT,rIN,HzL,dtdays,z_wL,         &
     &                       LBi,UBi,zlimit,IminS, ImaxS)




      USE mod_param
!
      implicit none
!
!     real(r8), intent(in)    :: z_w(0:N(ng))
      integer, intent(in)     :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in)    :: wBio
      real(r8), intent(in)    :: zlimit
      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(in) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
!     real(r8), intent(in)    :: HzL(N(ng))
      real(r8), intent(in)    :: dtdays


      real(r8), intent(out) :: rIN(IminS:ImaxS,N(ng))
      real(r8), intent(out) :: rOUT(IminS:ImaxS,N(ng))

      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff


      real(r8):: dBio(0:N(ng)), wBiod(LBi:UBi,0:N(ng))
      real(r8) :: FC(IminS:ImaxS,0:N(ng))

      real(r8) :: Hz_inv(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv2(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv3(IminS:ImaxS,N(ng))


      integer :: ksource(IminS:ImaxS,N(ng))

      real(r8) :: qR(IminS:ImaxS,N(ng))
      real(r8) :: qL(IminS:ImaxS,N(ng))
      real(r8) :: WL(IminS:ImaxS,N(ng))
      real(r8) :: WR(IminS:ImaxS,N(ng))
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc


      IF ( zlimit .lt. 0 ) THEN
       DO k=0,N(ng)
         DO i=LBi, UBi
          IF ( z_wL(i,k) .ge. zlimit ) THEN
            wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2/        &
     &                   (zlimit/2)**2 )
          ELSE
            wBiod(i,k) = 0.0_r8
          END IF
        END DO
       END DO
      ELSE
        DO k=0,N(ng)
          DO i=LBi, UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
!
!
!  Compute inverse thickness to avoid repeated divisions.
!

      DO k=1,N(ng)
       DO i=LBi,UBi
         Hz_inv(i,k)=1.0_r8/HzL(i,k)
       END DO
      END DO
      DO k=1,N(ng)-1
       DO i=LBi,UBi
         Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k+1))
       END DO
      END DO
      DO k=2,N(ng)-1
       DO i=LBi,UBi
         Hz_inv3(i,k)=1.0_r8/(HzL(i,k-1)+HzL(i,k)+HzL(i,k+1))
       END DO
      END DO

      DO k=1,N(ng)
         DO i=LBi,UBi
           qc(i,k)=Bio(i,k)
         END DO
      END DO

! !
! !  Reconstruct vertical profile of suspended sediment "qc" in terms
! !  of a set of parabolic segments within each grid box. Then, compute
! !  semi-Lagrangian flux due to sinking.
! !


      DO k=2,N(ng),1
        DO i=LBi,UBi

          FC(i,k)=(qc(i,k-1)-qc(i,k))*Hz_inv2(i,k)

        END DO
      END DO

      DO k=N(ng)-1,2,-1
         DO i=LBi,UBi
          dltR=HzL(i,k)*FC(i,k)
          dltL=HzL(i,k)*FC(i,k+1)

          cff=HzL(i,k+1)+2.0_r8*HzL(i,k)+HzL(i,k-1)
          cffR=cff*FC(i,k)
          cffL=cff*FC(i,k+1)




! !
! !  Apply PPM monotonicity constraint to prevent oscillations within the
! !  grid box.
! !
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
! !
! !  Compute right and left side values (qR,qL) of parabolic segments
! !  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
! !
! !  NOTE: Although each parabolic segment is monotonic within its grid
! !        box, monotonicity of the whole profile is not guaranteed,
! !        because qL(k+1)-qR(k) may still have different sign than
! !        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
! !        are reconciled using WENO procedure.
! !
          cff=(dltR-dltL)*Hz_inv3(i,k)
          dltR=dltR-cff*HzL(i,k-1)
          dltL=dltL+cff*HzL(i,k+1)

          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
          WR(i,k)=(2.0_r8*dltR-dltL)**2
          WL(i,k)=(dltR-2.0_r8*dltL)**2
        END DO
      END DO

      cff=1.0E-14_r8

      DO k=N(ng),2,-1
        DO i=LBi,UBi
          dltL=MAX(cff,WL(i,k  ))

          dltR=MAX(cff,WR(i,k-1))

          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k-1))/(dltR+dltL)

          qL(i,k-1)=qR(i,k)
        END DO
      END DO

      DO i=LBi,UBi
             FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
#if defined LINEAR_CONTINUATION
             qL(i,N(ng))=qR(i,N(ng)-1)
             qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
#elif defined NEUMANN
             qL(i,N(ng))=qR(i,N(ng)-1)
             qR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
#else
             qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
             qL(i,N(ng))=qc(i,N(ng))         ! conditions
             qR(i,N(ng)-1)=qc(i,N(ng))
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
! !
! !  Apply monotonicity constraint again, since the reconciled interfacial
! !  values may cause a non-monotonic behavior of the parabolic segments
! !  inside the grid box.
! !
!

      DO k=N(ng),1,-1
        DO i=LBi,UBi
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
! !
! !  After this moment reconstruction is considered complete. The next
! !  stage is to compute vertical advective fluxes, FC. It is expected
! !  that sinking may occurs relatively fast, the algorithm is designed
! !  to be free of CFL criterion, which is achieved by allowing
! !  integration bounds for semi-Lagrangian advective flux to use as
! !  many grid boxes in upstream direction as necessary.
! !
! !  In the two code segments below, WL is the z-coordinate of the
! !  departure point for grid box interface z_w with the same indices;
! !  FC is the finite volume flux; ksource(:,k) is index of vertical
! !  grid box which contains the departure point (restricted by N(ng)).
! !  During the search: also add in content of whole grid boxes
! !  participating in FC.
! !

!
!

      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
          cff=dtdays*ABS(wBiod(i,k))
          FC(i,k)=0.0_r8

          WL(i,k)=z_wL(i,k+1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO

      DO k=N(ng)-1,1,-1
        DO ks=N(ng),k,1
          DO i=LBi,UBi
            IF (WL(i,k).gt.z_wL(i,ks)) THEN
              ksource(i,k)=ks-1

              FC(i,k+1)=FC(i,k+1)+WR(i,ks)

            END IF
          END DO
        END DO
      END DO

!
!
! !
! !  Finalize computation of flux: add fractional part.
! !
      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
          ks=ksource(i,k)

          cu=MIN(1.0_r8,(WL(i,k)-z_wL(i,ks+1))*Hz_inv(i,ks))
          FC(i,k+1)=FC(i,k+1)+                                          &
      &                  HzL(i,ks)*cu*                                  &
      &                  (qL(i,ks)+                                     &
      &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-               &
      &                       (1.5_r8-cu)*                              &
      &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
        END DO

     END DO


     DO k=1,N(ng),1
       DO i=LBi,UBi
         rIN(i,k)=MAX(0.0_r8,FC(i,k))
         if (k.eq.1) THEN
           rIN(i,1)=0.0_r8
         endif
       END DO
     END DO

     DO k=1,N(ng)-1,1
       DO i=LBi,UBi
         rOUT(i,k)=MAX(0.0_r8,FC(i,k+1))
         rOUT(i,N(ng))=0.0_r8
       END DO
     END DO



      RETURN

      END SUBROUTINE BIORISE_1








!=====================================================================
! TracerSink, TracerRise (and BioSink, called by both)
! New subroutine for vertical movement, replaces BIOSINK_1, BIORISE_1
!=====================================================================

!************************************************************************

      Subroutine TracerSink(nn, Hz, z_w, dtdays, wBio, zlimit, Bio, Bout, lost)
      ! Arguments:
      ! nn:     scalar, # layers
      ! Hz:     n x 1, layer thickness (m)
      ! z_w:    n+1 x 1, depths of layer edges (m, negative down)
      ! dtdays: time step (days)
      ! wBio:   sinking rate (m d^-1, positive down)
      ! zlimit: depth below which sinking stops (m, negative down, >=0 = no limit)
      ! Bio:    n x 1, biomass concentration input (g m^-3)
      ! Bout:   n x 1, biomass concentration output (g m^-3)
      ! lost:   scalar, amount of biomass lost across lower boundary (g m^-2)

      USE mod_param
      implicit none

      integer,  intent(in) :: nn
      real(r8), intent(in) :: wBio
      real(r8), intent(in) :: zlimit
      real(r8), intent(in) :: dtdays
      real(r8), intent(in)  :: Hz(1,1,nn)
      real(r8), intent(in)  :: z_w(1,1,0:nn)
      real(r8), intent(in)  :: Bio(1,nn)

      real(r8), intent(out) :: Bout(nn)
      real(r8), intent(out) :: lost

      integer  :: k
      real(r8) :: zwtmp(0:nn)
      real(r8) :: Hztmp(nn)

      ! Pack the array slices into appropriate 1D arrays

      DO k = 1,nn
        Bout(k) = Bio(1,k)
        Hztmp(k) = Hz(1,1,k)
      END DO
      DO k = 0,nn
        zwtmp(k) = z_w(1,1,k)
      END DO

      ! Call the sinking routine

      call BioSink(nn, Bout, wBio, Hztmp, dtdays, zwtmp, zlimit, lost)

      RETURN
      END SUBROUTINE TracerSink


!************************************************************************

      Subroutine TracerRise(nn, Hz, z_w, dtdays, wBio, zlimit, Bio, Bout, lost)
      ! Arguments:
      ! n:      scalar, # layers
      ! Hz:     1 x 1 x n, slice, layer thickness array (m)
      ! z_w:    1 x 1 x n+1 slice, depths of layer edges (m, negative down)
      ! dtdays: time step (days)
      ! wBio:   rising rate (m d^-1, positive up)
      ! zlimit: depth above which rising stops (m, negative down, >=0 = no limit)
      ! Bio:    1 x n x 1 slice, biomass concentration input (g m^-3)
      ! Bout:   n x 1, biomass concentration output (g m^-3)
      ! lost:   scalar, amount of biomass lost across upper boundary (g m^-2)

      USE mod_param
      implicit none

      integer,  intent(in)  :: nn
      real(r8), intent(in)  :: wBio
      real(r8), intent(in)  :: zlimit
      real(r8), intent(in)  :: dtdays
      real(r8), intent(in)  :: Hz(1,1,nn)
      real(r8), intent(in)  :: z_w(1,1,0:nn)
      real(r8), intent(in)  :: Bio(1,nn)

      real(r8), intent(out) :: Bout(nn)
      real(r8), intent(out) :: lost

      integer  :: k
      real(r8) :: zwtmp(0:nn)
      real(r8) :: Hztmp(nn)
      real(r8) :: Btmp(nn)
      real(r8) :: zlimflip

      ! Flip the water column upside down.  The bottom is now at z = 0
      ! and the surface at water depth+free surface height.

      DO k = 1,nn
        Btmp(k) = Bio(1,nn+1-k) ! flip
        Hztmp(k) = Hz(1,1,nn+1-k) ! flip
      END DO
      DO k = 0,nn
        zwtmp(k) = z_w(1,1,0) - z_w(1,1,nn-k) ! make surface the bottom
      END DO

      zlimflip = z_w(1,1,0) - zlimit

      ! Call the sinking routine

      call BioSink(nn, Btmp, wBio, Hztmp, dtdays, zwtmp, zlimflip, lost)

      ! Flip the water column back

      DO k = 1,nn
        Bout(k) = Btmp(nn+1-k) ! flip back
      END DO

      RETURN
      END SUBROUTINE TracerRise

!************************************************************************

      Subroutine BioSink(n, Btmp, wBio, HzL, dtdays, z_wL, zlimit, sinkout)

      !------------------------------------------------------------------
      ! Computes redistribution of a tracer in the water column due to
      ! sinking.  After J. Warner sed sink code.
      !
      !   n       = # layers
      !   Btmp    = n x 1 array, tracer concentration (bottom to top)
      !   wBio    = sinking rate (m d^-1)
      !   HzL     = n x 1 array, thickness of layers (m)
      !   dtdays  = time step (d)
      !   z_wL    = n+1 x 1 array, depth of layer edges (m, negative)
      !   zlimit  = maximum depth for sinking (m, negative)
      !   sinkout = amount lost out of bottom cell (concentration)
      !
      ! Modifies Btmp and sinkout in the calling program
      !------------------------------------------------------------------

      USE mod_param
      implicit none

      integer,  intent(in) :: n
      real(r8), intent(in) :: wBio
      real(r8), intent(in) :: zlimit
      real(r8), intent(in) :: dtdays
      real(r8), intent(in) :: z_wL(0:n)
      real(r8), intent(inout) :: Btmp(n)
      real(r8), intent(in) :: HzL(n)
      real(r8), intent(inout) :: sinkout

      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff

      real(r8):: dBio(0:n), wBiod(0:n)
      real(r8) :: FC(0:n)

      real(r8) :: Hz_inv(n)
      real(r8) :: Hz_inv2(n)
      real(r8) :: Hz_inv3(n)

      integer :: ksource(n)

      real(r8) :: qR(n)
      real(r8) :: qL(n)
      real(r8) :: WL(n)
      real(r8) :: WR(n)
      real(r8), dimension(n) :: qc


      ! Modify sinking rates as necessary

      IF ( zlimit .lt. 0 ) THEN
        DO k=0,n
          IF ( z_wL(k) .ge. zlimit ) THEN
            wBiod(k) = wBio*exp(-1*(z_wL(k)-(zlimit/2))**2/(zlimit/2)**2)
          ELSE
            wBiod(k) = 0.0_r8
          END IF
        END DO
      ELSE
        DO k=0,n
          wBiod(k) = wBio
        END DO
      END IF


      ! Compute inverse thickness to avoid repeated divisions.


      DO k=1,n
        Hz_inv(k)=1.0_r8/HzL(k)
      END DO
      DO k=1,n-1
        Hz_inv2(k)=1.0_r8/(HzL(k)+HzL(k+1))
      END DO
      DO k=2,n-1
        Hz_inv3(k)=1.0_r8/(HzL(k-1)+HzL(k)+HzL(k+1))
      END DO

      DO k=1,n
        qc(k)=Btmp(k)
      END DO


      ! Reconstruct vertical profile of suspended sediment "qc" in terms
      ! of a set of parabolic segments within each grid box. Then,
      ! compute semi-Lagrangian flux due to sinking.

      DO k=n-1,1,-1
          FC(k)=(qc(k+1)-qc(k))*Hz_inv2(k)
      END DO

      DO k=2,n-1
          dltR=HzL(k)*FC(k)
          dltL=HzL(k)*FC(k-1)
          cff=HzL(k-1)+2.0_r8*HzL(k)+HzL(k+1)
          cffR=cff*FC(k)
          cffL=cff*FC(k-1)

          ! Apply PPM monotonicity constraint to prevent oscillations
          ! within the grid box.

          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF

          ! Compute right and left side values (qR,qL) of parabolic
          ! segments within grid box Hz(k); (WR,WL) are measures of
          ! quadratic variations.
          !
          !  NOTE: Although each parabolic segment is monotonic within
          !        its grid box, monotonicity of the whole profile is not
          !        guaranteed, because qL(k+1)-qR(k) may still have
          !        different sign than qc(k+1)-qc(k).  This possibility
          !        is excluded, after qL and qR are reconciled using WENO
          !        procedure.

          cff=(dltR-dltL)*Hz_inv3(k)
          dltR=dltR-cff*HzL(k+1)
          dltL=dltL+cff*HzL(k-1)
          qR(k)=qc(k)+dltR
          qL(k)=qc(k)-dltL
          WR(k)=(2.0_r8*dltR-dltL)**2
          WL(k)=(dltR-2.0_r8*dltL)**2
      END DO

      cff=1.0E-14_r8

      DO k=2,n-2
        dltL=MAX(cff,WL(k  ))
        dltR=MAX(cff,WR(k+1))
        qR(k)=(dltR*qR(k)+dltL*qL(k+1))/(dltR+dltL)
        qL(k+1)=qR(k)
      END DO


      FC(n)=0.0_r8                  ! no-flux boundary condition
#if defined LINEAR_CONTINUATION
      qL(n)=qR(n-1)
      qR(n)=2.0_r8*qc(n)-qL(n)
#elif defined NEUMANN
      qL(n)=qR(n-1)
      qR(n)=1.5_r8*qc(n)-0.5_r8*qL(n)
#else
      qR(n)=qc(n)                   ! default strictly monotonic
      qL(n)=qc(n)                   ! conditions
      qR(n-1)=qc(n)
#endif
#if defined LINEAR_CONTINUATION
      qR(1)=qL(2)
      qL(1)=2.0_r8*qc(1)-qR(1)
#elif defined NEUMANN
      qR(1)=qL(2)
      qL(1)=1.5_r8*qc(1)-0.5_r8*qR(1)
#else
      qL(2)=qc(1)                   ! bottom grid boxes are
      qR(1)=qc(1)                   ! re-assumed to be
      qL(1)=qc(1)                   ! piecewise constant.
#endif


      ! Apply monotonicity constraint again, since the reconciled
      ! interfacial values may cause a non-monotonic behavior of the
      ! parabolic segments inside the grid box.

      DO k=1,n
          dltR=qR(k)-qc(k)
          dltL=qc(k)-qL(k)
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
          qR(k)=qc(k)+dltR
          qL(k)=qc(k)-dltL
      END DO

      ! After this moment reconstruction is considered complete. The next
      ! stage is to compute vertical advective fluxes, FC. It is expected
      ! that sinking may occurs relatively fast, the algorithm is designed
      ! to be free of CFL criterion, which is achieved by allowing
      ! integration bounds for semi-Lagrangian advective flux to use as
      ! many grid boxes in upstream direction as necessary.
      !
      ! In the two code segments below, WL is the z-coordinate of the
      ! departure point for grid box interface z_w with the same indices;
      ! FC is the finite volume flux; ksource(:,k) is index of vertical
      ! grid box which contains the departure point (restricted by N(ng)).
      ! During the search: also add in content of whole grid boxes
      ! participating in FC.

      DO k=1,n
        cff=dtdays*ABS(wBiod(k))
        FC(k-1)=0.0_r8
        WL(k)=z_wL(k-1)+cff
        WR(k)=HzL(k)*qc(k)
        ksource(k)=k
      END DO
      DO k=1,n
        DO ks=k,n-1
          IF (WL(k).gt.z_wL(ks)) THEN
            ksource(k)=ks+1
            FC(k-1)=FC(k-1)+WR(ks)
          END IF
        END DO
      END DO

      ! Finalize computation of flux: add fractional part.

      DO k=1,n

        ks=ksource(k)
        cu=MIN(1.0_r8,(WL(k)-z_wL(ks-1))*Hz_inv(ks))
        FC(k-1)=FC(k-1)+                                                &
     &                HzL(ks)*cu*                                       &
     &                (qL(ks)+                                          &
     &                 cu*(0.5_r8*(qR(ks)-qL(ks))-                      &
     &                     (1.5_r8-cu)*                                 &
     &                     (qR(ks)+qL(ks)-2.0_r8*qc(ks))))

      END DO

      ! New profile of tracer Btmp (mass per volume)

      DO k=1,n
        Btmp(k)=qc(k)+(FC(k)-FC(k-1))*Hz_inv(k)
      END DO

      ! Amount lost out of bottom cell (mass per area)

      sinkout = FC(0)

      RETURN
      END SUBROUTINE BioSink

!************************************************************************

!=======================================================================
!    DetSINK
!----------------
!====================================================================
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

