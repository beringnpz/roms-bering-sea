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

      ! TODO: is ICE_BIO bestnpz-specific?
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
     &                   nnew(ng), nstp(ng)                             &
#ifdef MASKING
     &                   ,GRID(ng) % rmask                              &
#endif
     &                   ,GRID(ng) % Hz                                 &
     &                   ,GRID(ng) % z_r                                &
     &                   ,GRID(ng) % z_w                                &
     &                   ,GRID(ng) % h                                  &
     &                   ,GRID(ng) % latr                               &
     &                   ,FORCES(ng) % srflx                            &
     &                   ,OCEAN(ng) % t                                 &
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
     &                         nnew, nstp                               &
#ifdef MASKING
     &                         ,rmask                                   &
#endif
     &                         ,Hz                                      &
     &                         ,z_r                                     &
     &                         ,z_w                                     &
     &                         ,h                                       &
     &                         ,lat                                     &
     &                         ,srflx                                   &
     &                         ,t                                       &
#if defined BENTHIC
     &                         ,bt                                      &
#endif
#if defined FEAST
     &                         ,u                                       &
     &                         ,v                                       &
     &                         ,GF                                      &
#endif
#if defined ICE_BIO
# ifdef CLIM_ICE_1D
     &                         ,it                                      &
     &                         ,itL                                     &
     &                         ,tclm                                    &
# elif defined BERING_10K
     &                         ,IcePhL                                  &
     &                         ,IceNO3                                  &
     &                         ,IceNH4                                  &
     &                         ,IceLog                                  &
     &                         ,ti                                      &
     &                         ,hi                                      &
     &                         ,ai                                      &
     &                         ,ageice                                  &
     &                         ,ui                                      &
     &                         ,vi                                      &
# endif
#endif
#ifdef STATIONARY
     &                         ,st                                      &
     &                         ,UBst                                    &
#endif
#ifdef STATIONARY2
     &                         ,st2                                     &
     &                         ,UBst2                                   &
#endif
     &                         ,Akt                                     &
     &                          )

     !==================================================================
     !  VARIABLE DECLARATIONS
     !==================================================================

     ! Modules

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
#   ifndef MATLABCOMPILE
      USE mod_ice
      USE IcePhLbc_mod, ONLY : IcePhLbc_tile
      USE IceNO3bc_mod, ONLY : IceNO3bc_tile
      USE IceNH4bc_mod, ONLY : IceNH4bc_tile
#   endif
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

      integer, intent(in)     :: ng, tile
      integer, intent(in)     :: LBi, UBi, LBj, UBj
      integer, intent(in)     :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in)     :: UBk, UBt
      integer, intent(in)     :: nnew, nstp
#ifdef STATIONARY
      integer, intent(in)     :: UBst
#endif
#if defined STATIONARY2
      integer, intent(in)     :: UBst2
#endif

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in)    :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in)    :: Hz(LBi:,LBj:,:)
      real(r8), intent(in)    :: z_r(LBi:,LBj:,:)
      real(r8), intent(in)    :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in)    :: h(LBi:,LBj:)
      real(r8), intent(in)    :: lat(LBi:,LBj:)
      real(r8), intent(in)    :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)

# ifdef STATIONARY
      real(r8), intent(inout) :: st(LBi:,LBj:,:,:,:)
# endif
# ifdef STATIONARY2
      real(r8), intent(inout) :: st2(LBi:,LBj:,:,:)
# endif
# if defined BENTHIC
      real(r8), intent(inout) :: bt(LBi:,LBj:,:,:,:)
# endif
# if defined FEAST
      real(r8), intent(in)    :: u(LBi:UBi,LBj:UBj,UBk,3),v(LBi:UBi,LBj:UBj,UBk,3)
      TYPE (T_FEAST)          :: GF
# endif
# if defined ICE_BIO
#  ifdef CLIM_ICE_1D
      real(r8), intent(inout) :: it(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: itL(LBi:,LBj:,:,:)

#  elif defined BERING_10K
      real(r8), intent(in)    :: ti(LBi:,LBj:,:)
      real(r8), intent(in)    :: hi(LBi:,LBj:,:)
      real(r8), intent(in)    :: ai(LBi:,LBj:,:)
      real(r8), intent(in)    :: ageice(LBi:,LBj:,:)
      real(r8), intent(inout) :: IcePhL(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNO3(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNH4(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceLog(LBi:,LBj:,:)
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
#  endif
# endif

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:) ! TODO why is this passed in?  Never used?

#else

# ifdef MASKING
      real(r8), intent(in)    :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in)    :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in)    :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in)    :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in)    :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: lat(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# ifdef STATIONARY
      real(r8), intent(inout) :: st(LBi:UBi,LBj:UBj,UBk,3,UBst)
# endif
# ifdef STATIONARY2
      real(r8), intent(inout) :: st2(LBi:UBi,LBj:UBj,3,UBst2)
# endif
# if defined BENTHIC
      real(r8), intent(inout) :: bt(LBi:UBi,LBj:UBj,NBL(ng),3,NBeT(ng))
# endif
# if defined FEAST
      real(r8), intent(in)    :: u(LBi:UBi,LBj:UBj,UBk,3),v(LBi:UBi,LBj:UBj,UBk,3)
# endif
# if defined ICE_BIO
#  ifdef CLIM_ICE_1D
      real(r8), intent(inout) :: it(LBi:UBi,LBj:UBj,3,1)
      real(r8), intent(inout) :: itL(LBi:UBi,LBj:UBj,3,1)

      real(r8), intent(inout) :: tclmG(LBi:UBi,LBj:UBj,UBk,3,NH(ng)+2)
      real(r8), intent(inout) :: tclm(LBi:UBi,LBj:UBj,UBk,NT(ng)+2)
#  elif defined BERING_10K
      real(r8), intent(in)    :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: ai(LBi:UBi,LBj:UBj,2)     ! TODO: never used?
      real(r8), intent(in)    :: ageice(LBi:UBi,LBj:UBj,2) ! TODO: never used?
      real(r8), intent(inout) :: IcePhL(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNO3(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNH4(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceLog(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
#  endif
# endif
      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#endif

      ! Local variable declarations

#if defined FEAST

      real(r8) :: ROMS_depth(LBi:UBi,LBj:UBj,UBk)
      real(r8) :: ROMS_edges(LBi:UBi,LBj:UBj,UBk+1)
      real(r8) :: ROMS_temp(LBi:UBi,LBj:UBj,UBk)
      real(r8) :: ROMS_zoop(LBi:UBi,LBj:UBj,UBk,NUM_PLANKTON)
      real(r8) :: fal(LBi:UBi,LBj:UBj,nfvaral,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
      real(r8) :: fl(LBi:UBi,LBj:UBj,nfvarl,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS )
      real(r8) :: fsp(LBi:UBi,LBj:UBj,nfvar,NUM_SIMPLE_SPECIES)
      real(r8) :: ratepar(nrates)
      real(r8) :: t_new(LBi:Ubi,LBj:UBj)
      real(r8) :: base_speed, diff_coeff
      real(r8) :: delN, Ftot, Fmax, BB, CC, delCF, delCAL
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
      integer :: fvaral,fvarl,fvar,spal,spl,sp,spy,lc,ac,gr,isp,izoop,klev
      integer :: ictr,jctr
      real(r8) :: ftstp
      real(r8) :: NNstart, CFstart, CAstart, WW, Navail, Nloss, Ngloss, rec
      real(r8), dimension(TOT_LINKS):: NNbase,CFbase,CAbase,Npromote
      real(r8), dimension(TOT_FEAST):: eggs
      real :: CurD

      integer :: iv, spn,lcp,elder,younger

#endif

      integer :: i, j, k, ibio, itrc
      real(r8) :: cff0,cff1, cff1b,cff2,cff3,cff4
      real(r8) :: cff5,cff6,cff6b,cff7,cff8,cff9,cff10,cff11


#if defined BENTHIC
      integer :: ibioB
      real(r8) :: dw
      real(r8) :: totD, totDF, totPS, totPL, totBD
      real(r8), dimension(N(ng)) :: frac1, mfromlayer
      real(r8), dimension(N(ng),4) :: frac2
#endif

      integer :: Iter,is
      integer :: iday, month, year
      real(r8) :: hour, yday

      real(r8) :: dz
      real(r8) :: TFEup
      real(r8) :: dtdays
      real(r8) :: Dl, Par1, k_extV, k_chlV
      real(r8) :: Temp1
      real(r8) :: PON
      real(r8) :: NitrifMax, DLNitrif

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS) :: PARs

      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR

#ifdef ICE_BIO
      real(r8) :: aiceIfrac, aiceNfrac, dhicedt, trs, twi
      real(r8) :: grow1, GROWAice, fNO3, RAi0, RgAi
      real(r8) :: sb, gesi
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick, ice_status
#endif

      ! Setup

      real(r8), dimension(1,N(ng)) :: dBtmp
      real(r8) :: flxtmp
      logical :: defaultCMdp

      real(r8) :: RSNC, RENC, SSNC, SENC, RSCM, RECM, SSCM, SECM
      real(r8) :: respNC, respCM, eCM, eNC
      real(r8) :: targetdepth

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: cffsw, Ifrac

      ! Bio tracer setup

      integer  :: iiNO3,    iiNH4,    iiPhS,  iiPhL,  iiMZS, iiMZL, iiCop
      integer  :: iiNCaS,   iiNCaO,   iiEupS, iiEupO, iiDet, iiDetF
      integer  :: iiJel,    iiFe,     iiBen,  iiDetBen
      integer  :: iiIcePhL, iiIceNO3, iiIceNH4
      real(r8), dimension(IminS:ImaxS,N(ng),20) :: Bio3d, Bio2d, Bio_bak, DBio
      real(r8), dimension(IminS:ImaxS,N(ng),20) :: extrabio
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
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Rem_Det_NH4, Rem_DetF_NH4, Nit_NH4_NO3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gra_Det_Ben,Gra_DetF_Ben, Gra_PhS_Ben, Gra_PhL_Ben, Gra_DetBen_Ben, Exc_Ben_NH4, Exc_Ben_DetBen, Res_Ben_NH4, Mor_Ben_DetBen, Pre_Ben_DetBen, Rem_DetBen_NH4
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Gpp_INO3_IPhL, Gpp_INH4_IPhL, Res_IPhL_INH4, Mor_IPhL_INH4, Nit_INH4_INO3, Twi_IPhL_PhL, Twi_INO3_NO3, Twi_INH4_NH4
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Ver_PhS_DetBen, Ver_PhS_Out, Ver_PhL_DetBen, Ver_PhL_Out, Ver_Det_DetBen, Ver_Det_Out, Ver_DetF_DetBen, Ver_DetF_Out, Ver_NCaO_DetBen, Ver_NCaS_DetF, Ver_NCaS_DetBen
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Frz_PhL_IPhL, Frz_NO3_INO3, Frz_NH4_INH4
      real(r8), dimension(IminS:ImaxS,N(ng)) :: prod_PhS, prod_PhL, prod_MZL, prod_Cop, prod_NCaS, prod_EupS, prod_NCaO, prod_EupO, prod_Jel, prod_Ben, prod_IcePhL

      ! Biological source/sinks

      real(r8) :: LightLimS, NOLimS, NHLimS, IronLimS
      real(r8) :: LightLimL, NOLimL, NHLimL, IronLimL
      real(r8) :: alphaPhSv, alphaPhLv, DrateS, DrateL, PmaxS, PmaxL, PmaxsS, PmaxsL
      real(r8) :: IcePhlAvail
      real(r8), dimension(IminS:ImaxS,N(ng)) :: BasMetMZL, BasMetCop, BasMetNC, BasMetCM, BasMetEup
      real(r8) :: ParW, OffSet
      real(r8) :: fracUsed

      ! Parameter default values

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

#include "set_bounds.h"

#ifdef MATLABCOMPILE
!       ! Check inputs (prints values of scalars, shapes of arrays)
!       print *, "ng     = ", ng
!       print *, "tile   = ", tile
!       print *, "LBi    = ", LBi
!       print *, "UBi    = ", UBi
!       print *, "LBj    = ", LBj
!       print *, "UBj    = ", UBj
!       print *, "IminS  = ", IminS
!       print *, "ImaxS  = ", ImaxS
!       print *, "JminS  = ", JminS
!       print *, "JmaxS  = ", JmaxS
!       print *, "UBk    = ", UBk
!       print *, "UBt    = ", UBt
!       print *, "nnew   = ", nnew
!       print *, "nstp   = ", nstp
! # ifdef MASKING
!       print *, "rmask  shape = ", shape(rmask )
! # endif
!       print *, "Hz     shape = ", shape(Hz    )
!       print *, "z_r    shape = ", shape(z_r   )
!       print *, "z_w    shape = ", shape(z_w   )
!       print *, "h      shape = ", shape(h     )
!       print *, "lat    shape = ", shape(lat   )
!       print *, "srflx  shape = ", shape(srflx )
!       print *, "t      shape = ", shape(t     )
! # if defined BENTHIC
!       print *, "bt     shape = ", shape(bt    )
! # endif
! # if defined FEAST
!       print *, "u      shape = ", shape(u     )
!       print *, "v      shape = ", shape(v     )
!       print *, "GF     shape = ", shape(GF    )
! # endif
! # if defined ICE_BIO
! #  ifdef CLIM_ICE_1D
!       print *, "it     shape = ", shape(it    )
!       print *, "itL    shape = ", shape(itL   )
!       print *, "tclm   shape = ", shape(tclm  )
! #  elif defined BERING_10K
!       print *, "IcePhL shape = ", shape(IcePhL)
!       print *, "IceNO3 shape = ", shape(IceNO3)
!       print *, "IceNH4 shape = ", shape(IceNH4)
!       print *, "IceLog shape = ", shape(IceLog)
!       print *, "ti     shape = ", shape(ti    )
!       print *, "hi     shape = ", shape(hi    )
!       print *, "ai     shape = ", shape(ai    )
!       print *, "ageice shape = ", shape(ageice)
!       print *, "ui     shape = ", shape(ui    )
!       print *, "vi     shape = ", shape(vi    )
! #  endif
! # endif
! # ifdef STATIONARY
!       print *, "st     shape = ", shape(st    )
!       print *, "UBst   = ", UBst
! # endif
! # ifdef STATIONARY2
!       print *, "st2    shape = ", shape(st2   )
!       print *, "UBst2  = ", UBst2
! # endif
!       print *, "Akt    shape = ", shape(Akt   )
#endif

      !==================================================================
      !  SOME SETUP APPLICABLE TO ALL GRID CELLS
      !==================================================================

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

      ! All 0 is the shortcut for lagging the onshelf group movement
      ! 1 month behind the offshelf group.  This maintains back-
      ! compatibility with input files that don't include the newer
      ! Rise/SinkCM parameters, since missing parameters are set to 0 by
      ! default.  The one-month lag was the hard-coded behavior before I
      ! added separate input parameters for the NCaS (i.e. CM) group.

      defaultCMdp = ((RiseStartCM .eq. 0.0_r8) .and.                    &
     &               (RiseEndCM   .eq. 0.0_r8) .and.                    &
     &               (SinkStartCM .eq. 0.0_r8) .and.                    &
     &               (SinkEndCM   .eq. 0.0_r8))

      if (defaultCMdp) then

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

      ! Light attenuation fraction

      CALL bestnpz_swfrac_tile(ng, tile,                                &
   &        LBi, UBi, LBj, UBj,                                         &
   &        IminS, ImaxS, JminS, JmaxS,                                 &
   &        h(IminS:ImaxS,JminS:JmaxS),                                 &
   &        z_w(IminS:ImaxS,JminS:JmaxS,0:N(ng)),                       &
   &        max(0.0_r8, t(IminS:ImaxS,JminS:JmaxS,1:N(ng),nstp,iPhS)),  &
   &        max(0.0_r8, t(IminS:ImaxS,JminS:JmaxS,1:N(ng),nstp,iPhL)),  &
   &        cffsw, Ifrac)

#endif

      !==================================================================
      !  JLOOP: BEGIN HORIZONTAL LOOP
      !==================================================================

      J_LOOP : DO j=Jstr,Jend

        !---------------------------------
        ! Biological tracer variable setup
        !---------------------------------

        ! The various biological state variables are passed into this
        ! function in a few different arrays, depending on whether it's a
        ! tracer variable in the water (subject to advection and
        ! diffusion), a benthic variable (no movement), or an ice
        ! variable (subject to movement within the ice submodel, if
        ! running with the ice submodel, or being set analytically in
        ! climatological 1D mode).
        !
        ! Further complicating things is the mix of units.  The main
        ! tracer array uses a mix of units: Tunits (i.e. tracer units, in
        ! this case mass m^-3) for t(:,:,:,nstp,:) but m Tunits for
        ! t(:,:,:,nnew,:) (b/c we don't yet know the layer thicknesses
        ! that will be used in the predictor step nnew).  The benthic
        ! Tunits are mass m^-2, and the ice Tunits are mass m^-3.
        !
        ! To make long-term maintenance of this code easier, within this
        ! J_LOOP section we'll work with i x k x var arrays, where
        ! i = horizontal grid cell looping dimension, k = depth (counting
        ! from bottom to top), and var is the biological state variable
        ! index.

        ! First, set up the var indices, so I don't have to switch around
        ! between the 3 different sets used for pelagic, benthic, and ice
        ! variables in the input arrays, or worry about which options
        ! (benthic, ice, jelly, etc) are turned on.

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
        iiDetBen = 17
        iiIcePhL = 18
        iiIceNO3 = 19
        iiIceNH4 = 20

        ! All state variables will be saved in two different versions:
        ! mass m^-3 (Bio3d) and mass m^-2 (Bio2d).  This redundancy
        ! makes for clearer (for human readers) code,

        Bio3d = 0.0_r8 ! Initialize to 0
        Bio2d = 0.0_r8

        ! Pelagic variables: These are originally stored in per-volume
        ! concentrations in each water column layer. NO3 and NH4 are in
        ! mmol N m^-3, Fe is in umol Fe m^-3, and the rest are in mg C
        ! m^-3.

        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio3d(i,k,iiNO3 ) = t(i,j,k,nstp,iNO3)
            Bio3d(i,k,iiNH4 ) = t(i,j,k,nstp,iNH4)
            Bio3d(i,k,iiPhS ) = t(i,j,k,nstp,iPhS)
            Bio3d(i,k,iiPhL ) = t(i,j,k,nstp,iPhL)
            Bio3d(i,k,iiMZS ) = t(i,j,k,nstp,iMZS)
            Bio3d(i,k,iiMZL ) = t(i,j,k,nstp,iMZL)
            Bio3d(i,k,iiCop ) = t(i,j,k,nstp,iCop)
            Bio3d(i,k,iiNCaS) = t(i,j,k,nstp,iNCaS)
            Bio3d(i,k,iiEupS) = t(i,j,k,nstp,iEupS)
            Bio3d(i,k,iiNCaO) = t(i,j,k,nstp,iNCaO)
            Bio3d(i,k,iiEupO) = t(i,j,k,nstp,iEupO)
            Bio3d(i,k,iiDet ) = t(i,j,k,nstp,iDet)
            Bio3d(i,k,iiDetF) = t(i,j,k,nstp,iDetF)
#ifdef JELLY
            Bio3d(i,k,iiJel ) = t(i,j,k,nstp,iJel)
#endif
#ifdef IRON_LIMIT
            Bio3d(i,k,iiFe  ) = t(i,j,k,nstp,iFe)
#endif
            DO itrc = iiNO3,iiFe
              Bio2d(i,k,itrc) = Bio3d(i,k,itrc)*Hz(i,j,k)
            END DO

#ifdef STATIONARY
            ! Rate of change due to processes outside this subroutine
            ! (Note: this calc depends on the step3d_t.F code where
            ! stationary diagnostic values are copied from nstp to nnew...
            ! if that changes, this will break.)

            st(i,j,k,nstp,131) = (Bio2d(i,k,iiNO3 ) - st(i,j,k,nnew,117))/dtdays ! mmol N m^-2 d^-1
            st(i,j,k,nstp,132) = (Bio2d(i,k,iiNH4 ) - st(i,j,k,nnew,118))/dtdays ! mmol N m^-2 d^-1
            st(i,j,k,nstp,133) = (Bio2d(i,k,iiPhS ) - st(i,j,k,nnew,119))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,134) = (Bio2d(i,k,iiPhL ) - st(i,j,k,nnew,120))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,135) = (Bio2d(i,k,iiMZL ) - st(i,j,k,nnew,121))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,136) = (Bio2d(i,k,iiCop ) - st(i,j,k,nnew,122))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,137) = (Bio2d(i,k,iiNCaS) - st(i,j,k,nnew,123))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,138) = (Bio2d(i,k,iiEupS) - st(i,j,k,nnew,124))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,139) = (Bio2d(i,k,iiNCaO) - st(i,j,k,nnew,125))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,140) = (Bio2d(i,k,iiEupO) - st(i,j,k,nnew,126))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,141) = (Bio2d(i,k,iiDet ) - st(i,j,k,nnew,127))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,142) = (Bio2d(i,k,iiDetF) - st(i,j,k,nnew,128))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,143) = (Bio2d(i,k,iiJel ) - st(i,j,k,nnew,129))/dtdays ! mg C m^-2 d^-1
            st(i,j,k,nstp,144) = (Bio2d(i,k,iiFe  ) - st(i,j,k,nnew,130))/dtdays ! umol Fe m^-2 d^-1
#endif
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

#ifdef BENTHIC
        DO k=1,NBL(ng) ! Note: For BESTNPZ, NBL = 1 is hard-coded in mod_param.F
          DO i=Istr,Iend
            Bio2d(i,k,iiBen   ) = bt(i,j,k,nstp,iBen)
            Bio2d(i,k,iiDetBen) = bt(i,j,k,nstp,iDetBen)

            DO itrc=iiBen,iiDetBen
              Bio3d(i,k,itrc) = Bio2d(i,k,itrc)/Hz(i,j,k)
            END DO
          END DO
        END DO
#endif

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
          ! (Note: Most arrays have data arranged such that time
          ! dimensions = [1 2 3]
          !            = [nstp nnew 3]
          !            = [current predictor intermediate]
          ! Ice needs to hold on to the previous step, so those arrays
          ! are [1 2] = [nprev nstp] = [previous current].  The flipflop
          ! of values below is simply done so we can keep all the nstp
          ! values at the same index location (1) and therefore access
          ! them with the nstp shorthand.  Side effect is the confusing
          ! shorthand that for ice, nnew refers to the previous step
          ! when reading, and the predictor step later when writing.

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

        ! Start by extracting ice biomass from the main ice bio tracer
        ! arrays (we'll deal with changes due to appearing or
        ! disappearing ice in a moment)

#ifdef ICE_BIO
        DO i=Istr,Iend

# ifdef CLIM_ICE_1D
            Bio3d(i,N(ng),iiIcePhL) = it(i,j,nstp,idice(1))
            Bio3d(i,N(ng),iiIceNO3) = it(i,j,nstp,idice(2))
            Bio3d(i,N(ng),iiIceNH4) = it(i,j,nstp,idice(3))
# elif defined BERING_10K
            Bio3d(i,N(ng),iiIcePhL) = IcePhL(i,j,nstp)
            Bio3d(i,N(ng),iiIceNO3) = IceNO3(i,j,nstp)
            Bio3d(i,N(ng),iiIceNH4) = IceNH4(i,j,nstp)
# endif

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

        ! Save a copy of the original biomass

        Bio_bak = Bio2d

        ! If any biomass is negtive, replace with 0 for all source/sink
        ! calculations

        Bio2d = max(0.0_r8, Bio2d)
        Bio3d = max(0.0_r8, Bio3d)

        extrabio = Bio2d - Bio_bak ! biomass added to make non-negative

#ifdef ICE_BIO
        ! Move tracers between surface water layer and ice skeletal layer
        ! if ice appeared or disappeared.  Note that the freezing rate
        ! diagnostics Frz_X_IX assume this process takes place over a
        ! full time step, rather than instantaneously, just to maintain
        ! some comparability to the other fluxes

        Frz_PhL_IPhL = 0.0_r8
        Frz_NO3_INO3 = 0.0_r8
        Frz_NH4_INH4 = 0.0_r8

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

            Frz_PhL_IPhL(i,N(ng)) = (Bio3d(i,N(ng),iiIcePhl) * aidz -  Bio2d(i,N(ng),iiIcePhL))/(dt(ng)*sec2day)
            Frz_NO3_INO3(i,N(ng)) = (Bio3d(i,N(ng),iiIceNO3) * aidz -  Bio2d(i,N(ng),iiIceNO3))/(dt(ng)*sec2day)
            Frz_NH4_INH4(i,N(ng)) = (Bio3d(i,N(ng),iiIceNH4) * aidz -  Bio2d(i,N(ng),iiIceNH4))/(dt(ng)*sec2day)

            Bio2d(i,N(ng),iiPhL)    = Bio3d(i,N(ng),iiPhl) * Hz(i,j,N(ng))
            Bio2d(i,N(ng),iiNO3)    = Bio3d(i,N(ng),iiNO3) * Hz(i,j,N(ng))
            Bio2d(i,N(ng),iiNH4)    = Bio3d(i,N(ng),iiNH4) * Hz(i,j,N(ng))

            Bio2d(i,N(ng),iiIcePhL) = Bio3d(i,N(ng),iiIcePhl) * aidz
            Bio2d(i,N(ng),iiIceNO3) = Bio3d(i,N(ng),iiIceNO3) * aidz
            Bio2d(i,N(ng),iiIceNH4) = Bio3d(i,N(ng),iiIceNH4) * aidz

          elseif (ice_status(i,j) .eq. -1.0) then

            ! If ice disappeared, biomass that was in the ice gets dumped
            ! into the water surface layer

            Frz_PhL_IPhL(i,N(ng)) = -Bio2d(i,N(ng),iiIcePhL)/(dt(ng)*sec2day)
            Frz_NO3_INO3(i,N(ng)) = -Bio2d(i,N(ng),iiIceNO3)/(dt(ng)*sec2day)
            Frz_NH4_INH4(i,N(ng)) = -Bio2d(i,N(ng),iiIceNH4)/(dt(ng)*sec2day)

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

          cff1 = 2.0_r8 * pi * ( yday-1.0_r8 ) / 365.0_r8
          cff2 = 0.006918_r8 - 0.399912_r8*cos(cff1)                    &
     &           + 0.070257_r8*sin(cff1) - 0.006758_r8*cos(2*cff1)      &
     &           + 0.000907_r8*sin(2*cff1) - 0.002697_r8*cos(3*cff1)    &
     &           + 0.00148_r8*sin(3*cff1) ! Solar declination from Oberhuber (1988) (COADS documentation)
          cff3 = lat(i,j) * pi /180.0_r8
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
          ! TODO: this option is definitely deprecated, and doesn't
          ! actually calculate PAR.  Should be removed soon.

          cff2 = 41.0_r8 - 35.0_r8                                      &
     &           * COS( ( 12.0_r8 + yday) * 2.0_r8 * pi / 365.0_r8 )
#else

          ! Calculate PAR at the surface
          ! Shortwave radiation input (in deg C m/s) is converted to
          ! W/m^2 (assuming standard density rho0=1025 kg/m^3 and heat
          ! capacity Cp=3985 J/kg/degC of seawater), then to photon flux
          ! (E/m^2/day) via
          !
          ! Eq = I * lambda/(h*c*A) * s2d
          !   where
          !
          !   Eq = photon flux (E/m^2/d)
          !   lambda = wavelength (m)
          !   I = irradiance (W/m^2)
          !   h = Plank's constant (Js)
          !   c = speed of light (m/s)
          !   A = Avogadro's number (mol^-1)
          !   s2d = 86400 s/d
          !
          ! The constant used here (copied from GOANPZ) implies an
          ! wavelength of ~547 nm... prob derives from a average across
          ! chl range (~400-700nm).

          PARs(i) = PARfrac(ng) * srflx(i,j) * rho0 * Cp * 0.394848_r8
#endif
        END DO

        !------------------------------------------
        ! Calculate light decay in the water column
        !------------------------------------------

        DO i=Istr,Iend
          DO k=1,N(ng)
            PAR(i,k) = Ifrac(i,j,k) * PARs(i)
          END DO
        END DO

! #ifdef NEWSHADE
!         ! Georgina Gibsons version after Morel 1988 (in Loukos 1997)
!
!         DO i=Istr,Iend
!
!           cff10=h(i,j)
!           k_extV= k_ext+2.00_r8*exp(-cff10*.05)
!
!           cff0=PARs(i)
!
!           DO k=N(ng),1,-1
!
!             dz=0.5_r8*(z_w(i,j,k)-z_w(i,j,k-1))
!             cff5=(Bio3d(i,k,iiPhS)/ccr)+ (Bio3d(i,k,iiPhL)/ccrPhL)
!             cff2 = (k_chl*(cff5)**(0.428_r8))
!             !cff2 = min(0.05_r8,max(0.0067_r8,(k_chl*(cff5)**(-0.428_r8))))
!             PAR(i,k) = cff0 * EXP(-(k_extV+cff2)*dz)
!             cff0=cff0 * EXP(-(k_extV+cff2)*dz*2.0_r8)
!
!           END DO
!         END DO
!
! #elif defined COKELET
!         ! Version from Ned Cokelet
!
!         DO i=Istr,Iend
!
!           cff10=h(i,j)
!           !  k_extV= k_ext+k_extZ*exp(-cff10*.05)
!           k_extV= k_ext
!           cff0=PARs(i)
!
!           DO k=N(ng),1,-1
!
!             dz=0.5_r8*(z_w(i,j,k)-z_w(i,j,k-1))
!             cff5=(Bio3d(i,k,iiPhS)/ccr)+ (Bio3d(i,k,iiPhL)/ccrPhL)
!             cff2 = (k_chlA*(cff5)**(k_chlB))
!
!             PAR(i,k) = cff0 * EXP(-(k_extV+cff2)*dz)
!             cff0=cff0 * EXP(-(k_extV+cff2)*dz*2.0_r8)
!
!           END DO
!         END DO
! #else
!         ! Version from Sarah Hinckley old C code
!         ! (KAK: probably wrong... do not use)
!
!         DO k=N(ng),1,-1
!           DO i=Istr,Iend
!             cff3 = z_r(i,j,k)+2.5_r8
!             IF ( cff3 .gt. -71.0_r8 ) THEN
!               cff1 = k_ext + k_chl *                                 &
!      &                  ( Bio3d(i,k,iiPhS) + Bio3d(i,k,iiPhL) ) / ccr
!             ELSE
!                 cff1 = 0.077_r8
!             END IF
!             PAR(i,k) = PARfrac(ng) * cff2 * exp( cff1 * cff3 )
!
!           END DO
!         END DO
! #endif


        !================================================================
        !  Begin time loop (if BioIter > 1, this divides the main time
        !  step into smaller steps for biological calculations)
        !================================================================

        ITER_LOOP: DO Iter=1,BioIter(ng)

          ! Initialize the rate of change, dB/dt, to 0 for all elements.
          ! Same for all intermediate flux arrays.  Note that these
          ! fluxes will hold the 2D equivalent of all the fluxes (i.e.
          ! per area, rather than per volume); this makes it easier to
          ! keep track of things that are moving between different-sized
          ! layers (e.g. ice to surface layer, or benthos to water
          ! column)

          DBio = 0.0_r8 ! Initializes entire array to 0

          Gpp_NO3_PhS    = 0.0_r8
          Gpp_NO3_PhL    = 0.0_r8
          Gpp_NH4_PhS    = 0.0_r8
          Gpp_NH4_PhL    = 0.0_r8
          Gra_PhS_MZL    = 0.0_r8
          Gra_PhL_MZL    = 0.0_r8
          Ege_MZL_Det    = 0.0_r8
          Gra_PhS_Cop    = 0.0_r8
          Gra_PhL_Cop    = 0.0_r8
          Gra_MZL_Cop    = 0.0_r8
          Gra_IPhL_Cop   = 0.0_r8
          Ege_Cop_DetF   = 0.0_r8
          Gra_PhS_NCaS   = 0.0_r8
          Gra_PhL_NCaS   = 0.0_r8
          Gra_MZL_NCaS   = 0.0_r8
          Gra_IPhL_NCaS  = 0.0_r8
          Ege_NCaS_DetF  = 0.0_r8
          Gra_PhS_NCaO   = 0.0_r8
          Gra_PhL_NCaO   = 0.0_r8
          Gra_MZL_NCaO   = 0.0_r8
          Gra_IPhL_NCaO  = 0.0_r8
          Ege_NCaO_DetF  = 0.0_r8
          Gra_PhS_EupS   = 0.0_r8
          Gra_PhL_EupS   = 0.0_r8
          Gra_MZL_EupS   = 0.0_r8
          Gra_Cop_EupS   = 0.0_r8
          Gra_IPhL_EupS  = 0.0_r8
          Gra_Det_EupS   = 0.0_r8
          Gra_DetF_EupS  = 0.0_r8
          Ege_EupS_DetF  = 0.0_r8
          Gra_PhS_EupO   = 0.0_r8
          Gra_PhL_EupO   = 0.0_r8
          Gra_MZL_EupO   = 0.0_r8
          Gra_Cop_EupO   = 0.0_r8
          Gra_IPhL_EupO  = 0.0_r8
          Gra_Det_EupO   = 0.0_r8
          Gra_DetF_EupO  = 0.0_r8
          Ege_EupO_DetF  = 0.0_r8
          Gra_Cop_Jel    = 0.0_r8
          Gra_EupS_Jel   = 0.0_r8
          Gra_EupO_Jel   = 0.0_r8
          Gra_NCaS_Jel   = 0.0_r8
          Gra_NCaO_Jel   = 0.0_r8
          Ege_Jel_DetF   = 0.0_r8
          Mor_PhS_Det    = 0.0_r8
          Mor_PhL_Det    = 0.0_r8
          Mor_MZL_Det    = 0.0_r8
          Mor_Cop_DetF   = 0.0_r8
          Mor_NCaS_DetF  = 0.0_r8
          Mor_EupS_DetF  = 0.0_r8
          Mor_NCaO_DetF  = 0.0_r8
          Mor_EupO_DetF  = 0.0_r8
          Mor_Jel_DetF   = 0.0_r8
          Res_PhS_NH4    = 0.0_r8
          Res_PhL_NH4    = 0.0_r8
          Res_MZL_NH4    = 0.0_r8
          Res_Cop_NH4    = 0.0_r8
          Res_NCaS_NH4   = 0.0_r8
          Res_NCaO_NH4   = 0.0_r8
          Res_EupS_NH4   = 0.0_r8
          Res_EupO_NH4   = 0.0_r8
          Res_Jel_NH4    = 0.0_r8
          Rem_Det_NH4    = 0.0_r8
          Rem_DetF_NH4   = 0.0_r8
          Nit_NH4_NO3    = 0.0_r8
          Gra_Det_Ben    = 0.0_r8
          Gra_DetF_Ben   = 0.0_r8
          Gra_PhS_Ben    = 0.0_r8
          Gra_PhL_Ben    = 0.0_r8
          Gra_DetBen_Ben = 0.0_r8
          Exc_Ben_NH4    = 0.0_r8
          Exc_Ben_DetBen = 0.0_r8
          Res_Ben_NH4    = 0.0_r8
          Mor_Ben_DetBen = 0.0_r8
          Rem_DetBen_NH4 = 0.0_r8
          Gpp_INO3_IPhL  = 0.0_r8
          Gpp_INH4_IPhL  = 0.0_r8
          Res_IPhL_INH4  = 0.0_r8
          Mor_IPhL_INH4  = 0.0_r8
          Nit_INH4_INO3  = 0.0_r8
          Twi_IPhL_PhL   = 0.0_r8
          Twi_INO3_NO3   = 0.0_r8
          Twi_INH4_NH4   = 0.0_r8
          Ver_PhS_DetBen = 0.0_r8
          Ver_PhS_Out    = 0.0_r8
          Ver_PhL_DetBen = 0.0_r8
          Ver_PhL_Out    = 0.0_r8
          Ver_Det_DetBen = 0.0_r8
          Ver_Det_Out    = 0.0_r8
          Ver_DetF_DetBen= 0.0_r8
          Ver_DetF_Out   = 0.0_r8
          Ver_NCaO_DetBen= 0.0_r8
          Ver_NCaS_DetF  = 0.0_r8
          Ver_NCaS_DetBen = 0.0_r8
          Ver_NCaS_DetBen = 0.0_r8


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

              IronLimS = min(1.0_r8, eps + Bio3d(i,k,iiFe)/(kfePhS + Bio3d(i,k,iiFe))*(kfePhS + FeCritPS)/FeCritPS)
              IronLimL = min(1.0_r8, eps + Bio3d(i,k,iiFe)/(kfePhL + Bio3d(i,k,iiFe))*(kfePhL + FeCritPL)/FeCritPL)
#endif

              ! Light limitation

#ifdef DENMAN
              LightLimS = TANH(alphaPhSv * PAR(i,k) / PmaxsS)
              LightLimL = TANH(alphaPhLv * PAR(i,k) / PmaxsL)
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

              if (DiS.le.0.0_r8) THEN
                Gpp_NO3_PhS(i,k) = 0; ! mg C m^-2 d^-1
                Gpp_NH4_PhS(i,k) = 0; ! mg C m^-2 d^-1
              endif
              if (DiL.le.0.0_r8) THEN
                Gpp_NO3_PhL(i,k) = 0; ! mg C m^-2 d^-1
                Gpp_NH4_PhL(i,k) = 0; ! mg C m^-2 d^-1
              endif

#ifdef STATIONARY

              ! Save limitation terms for output

              st(i,j,k,nstp,1) = LightLimS
              st(i,j,k,nstp,2) = LightLimL
              st(i,j,k,nstp,3) = NOLimS
              st(i,j,k,nstp,4) = NOLimL
              st(i,j,k,nstp,5) = NHLimS
              st(i,j,k,nstp,6) = NHLimL
              st(i,j,k,nstp,7) = IronLimS
              st(i,j,k,nstp,8) = IronLimL
#endif
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

              cff2 = eEup * Bio3d(i,k,iiEupS) / (fEup + cff1 + cff0)
              cff3 = Q10Eup ** ((Temp(i,k)-Q10EupT) / 10.0_r8)

#ifdef DEPTHLIMITER
              cff4 = 1.0_r8 - (0.5_r8 + 0.5_r8*tanh((h(i,j) - 200_r8)/.3_r8)) ! depth-limiter, stops growth if they move deeper than 200m
#else
              cff4 = 1.0_r8
#endif

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
     &                              (1.0_r8 - 0.3_r8)   * cff0) *       &
     &                              cff2 * cff3 * cff4

              ! Off-shelf euphausiids

              cff0 = fpDetEupO * Bio3d(i,k,iiDet)**2                    &
     &             + fpDetEupO * Bio3d(i,k,iiDetF)**2 ! detrital food

              cff2 = eEup * Bio3d(i,k,iiEupO) / (fEup + cff1 + cff0)

#ifdef DEPTHLIMITER
              cff4 = 1.0_r8 - (0.5_r8 + 0.5_r8*tanh((200_r8 - h(i,j))/.3_r8)) ! depth-limiter, stops growth if they move shallower than 200m
#else
              cff4 = 1.0_r8
#endif

              Gra_PhS_EupO(i,k)  = fpPhSEup * Bio3d(i,k,iiPhS)**2  * cff2 * cff3 * cff4
              Gra_PhL_EupO(i,k)  = fpPhLEup * Bio3d(i,k,iiPhL)**2  * cff2 * cff3 * cff4
              Gra_MZL_EupO(i,k)  = fpMZLEup * Bio3d(i,k,iiMZL)**2  * cff2 * cff3 * cff4
              Gra_Cop_EupO(i,k)  = fpCopEup * Bio3d(i,k,iiCop)**2  * cff2 * cff3 * cff4
              Gra_IPhL_EupO(i,k) = fpPhLEup * (IcePhlAvail)**2     * cff2 * cff3 * cff4
              Gra_Det_EupO(i,k)  = fpDetEup * Bio3d(i,k,iiDet)**2  * cff2 * cff3 * cff4
              Gra_DetF_EupO(i,k) = fpDetEup * Bio3d(i,k,iiDetF)**2 * cff2 * cff3 * cff4

              Ege_EupO_DetF(i,k) = ((1.0_r8 - gammaEup) * cff1 +        &
     &                              (1.0_r8 - 0.3_r8)   * cff0) *       &
     &                              cff2 * cff3 * cff4

              ! Jellyfish

              cff1 = fpCopJel * Bio3d(i,k,iiCop)**2 +                   &
     &               fpNCaJel * Bio3d(i,k,iiNCaS)**2 +                  &
     &               fpNCaJel * Bio3d(i,k,iiNCaO)**2 +                  &
     &               fpEupJel * Bio3d(i,k,iiEupS)**2 +                  &
     &               fpEupJel * Bio3d(i,k,iiEupO)**2

              cff2 = eJel * Bio3d(i,k,iiJel) / (fJel + cff1)
              cff3= Q10Jele ** ((Temp(i,k)-Q10JelTe) / 10.0_r8)

              Gra_Cop_Jel(i,k)  = fpCopJel * Bio3d(i,k,iiCop)**2  * cff2 * cff3
              Gra_NCaS_Jel(i,k) = fpNCaJel * Bio3d(i,k,iiNCaS)**2 * cff2 * cff3
              Gra_NCaO_Jel(i,k) = fpNCaJel * Bio3d(i,k,iiNCaO)**2 * cff2 * cff3
              Gra_EupS_Jel(i,k) = fpEupJel * Bio3d(i,k,iiEupS)**2 * cff2 * cff3
              Gra_EupO_Jel(i,k) = fpEupJel * Bio3d(i,k,iiEupO)**2 * cff2 * cff3

              ! Note: mentioned in docs that gammaJel can be >1 to allow
              ! for an outside food source.  However, GG's code doesn't
              ! specify how the flux to detritus might change in that
              ! case (as written currently, that extra would come out of
              ! the DetF biomass via a negative egestion flux) TODO: Do
              ! we want to allow gammaJel>1, and if so, how should we
              ! handle egestion?

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

          DO k=1,N(ng)
            DO i=Istr,Iend

              ! Phytoplankton (linear senescence)

              Mor_PhS_Det(i,k) = mPhS * Bio3d(i,k,iiPhS)
              Mor_PhL_Det(i,k) = mPhL * Bio3d(i,k,iiPhL)

              ! Microzooplankton (quadratic mortality, with option for
              ! linear)

#ifdef MZLM0LIN
              Mor_MZL_Det(i,k) = mMZL*Bio3d(i,k,iiMZL)          ! linear
#else
              Mor_MZL_Det(i,k) = mpredMZL*Bio3d(i,k,iiMZL)**2   ! quadratic
#endif

#ifdef fixedPRED
              ! TODO: original DBio(i,k,iXXX) = DBio(i,k,iXXX) - 0.5*Hz(i,j,k)/dtdays
              ! Implies coefficient units of mg C * day * m^-4???  Typo?
              ! Supposed to be constant rate, or maybe constant fraction
              ! of biomass?  Assuming the former for now.  This option
              ! seems deprecated and should probably be removed.
              Mor_Cop_DetF(i,k)  = 0.5
              Mor_NCaS_DetF(i,k) = 0.5
              Mor_EupS_DetF(i,k) = 1.0
              Mor_NCaO_DetF(i,k) = 0.5
              Mor_EupO_DetF(i,k) = 1.0
#else
              TFEup = Q10Eup ** ((Temp(i,k)-Q10EupT) / 10.0_r8)
# ifdef FEAST
              ! Mesozooplankton (quadratic predation closure).  FEAST
              ! predation only affects zooplankton within a specific
              ! region (mostly EBS).  This term (modified spatially by
              ! forcing input variable zoop_force) tries to smooth out
              ! some boundary issues that can come of that. Note that
              ! this term does *not* include predation by FEAST fish yet;
              ! that loss term will be added in feast_step.h.

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
          ! Nitrification and
          ! remineralization
          !------------------------------

          DO k=1,N(ng)
            DO i=Istr,Iend

              ! Detrital remineralization

              PON = Bio3d(i,k,iiDet)*xi  ! Particulate organic nitrogen in Det
              Rem_Det_NH4(i,k) = (Pv0 * exp(PvT*Temp(i,k)) * PON) ! mmol N m^-3 d^-1

              PON = Bio3d(i,k,iiDetF)*xi  ! Particulate organic nitrogen in DetF
              Rem_DetF_NH4(i,k) = (Pv0 * exp(PvT*Temp(i,k)) * PON) ! mmol N m^-3 d^-1

              ! Nitrification

              NitrifMax = Nitr0 * exp(-ktntr*(Temp(i,k) - ToptNtr)**2)     ! Arhonditsis 2005 temperature dependence

              ParW = PAR(i,k)/0.394848_r8 ! convert to W
              DLNitrif = (1 - MAX(0.0_r8, (ParW - tI0)/(KI + ParW - tI0))) ! Fennel light dependence
              DLNitrif = 1.0_r8  ! No light/depth dependence (overrides previous line)

              cff1 = Bio3d(i,k,iiNH4)/(KNH4Nit + Bio3d(i,k,iiNH4))         ! Arhonditsis saturation

              Nit_NH4_NO3(i,k) = NitrifMax * Bio3d(i,k,iiNH4) * DLNitrif * cff1 !  mmol N m^-3 d^-1

            END DO
          END DO

          ! Convert fluxes from volumetric to integrated over layer, and
          ! from N to C for consistency with other fluxes

          DO k=1,N(ng)
            DO i=Istr,Iend
              Rem_Det_NH4(i,k)  = Rem_Det_NH4(i,k)  * Hz(i,j,k)/xi
              Rem_DetF_NH4(i,k) = Rem_DetF_NH4(i,k) * Hz(i,j,k)/xi
              Nit_NH4_NO3(i,k)  = Nit_NH4_NO3(i,k)  * Hz(i,j,k)/xi
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
              mfromlayer(k) = cff1
              frac1(k) = cff1/Hz(i,j,k)

              ! Food available to benthos

              totD  = totD  + Bio2d(i,k,iiDet) *frac1(k)
              totDF = totDF + Bio2d(i,k,iiDetF)*frac1(k)
              totPS = totPS + Bio2d(i,k,iiPhS) *frac1(k)
              totPL = totPL + Bio2d(i,k,iiPhL) *frac1(k)

              cff2 = cff2 + Hz(i,j,k)
            END DO

            ! Fraction each layer contributing to each food source

            DO k = 1,N(ng)
              frac2(k,1) = mfromlayer(k)*Bio3d(i,k,iiDet) /totD
              frac2(k,2) = mfromlayer(k)*Bio3d(i,k,iiDetF)/totDF
              frac2(k,3) = mfromlayer(k)*Bio3d(i,k,iiPhS) /totPS
              frac2(k,4) = mfromlayer(k)*Bio3d(i,k,iiPhL) /totPL
            END DO
            if (totD .le. 0) then
              frac2(:,1) = 0
            endif
            if (totDF .le. 0) then
              frac2(:,2) = 0
            endif
            if (totPS .le. 0) then
              frac2(:,3) = 0
            endif
            if (totPL .le. 0) then
              frac2(:,4) = 0
            endif

            ! Potential food available from water column

            cff1 = (prefD *totD /((prefD *totD )+LupP))*prefD *totD
            cff2 = (prefD *totDF/((prefD *totDF)+LupP))*prefD *totDF
            cff3 = (prefPS*totPS/((prefPS*totPS)+LupP))*prefPS*totPS
            cff4 = (prefPL*totPL/((prefPL*totPL)+LupP))*prefPL*totPL

            cff6 = cff1+cff2+cff3+cff4 ! Total pelagic food

            ! Potential food available from  sea floor

            totBD = Bio2d(i,1,iiDetBen)
            cff5 = (prefD *totBD/((prefD *totBD)+LupD))*prefD *totBD

            ! Temperature mediation (for feeding and mortality)

            cff0 = q10r**((Temp(i,1)-T0benr)/10.0_r8)

            ! Total uptake of each food category
            ! TODO: Unit mismatch in the part... cff1 is mC/m^2 and
            ! (cff0*cff1*Bio2d(i,1,iiBen)*Rup/(cff6+KupP)) is mgC/m^2/d
            ! Is this supposed to be a zero-trap?  Should be cff1/dtdays
            ! (i.e. highest rate that would keep losses positive?)  Of
            ! course, that assumes no fluxes into the layer to possibly
            ! balance out a seemingly too-high loss rate.

!             cff7  = min(cff1,(cff0*cff1*Bio2d(i,1,iiBen)*Rup/(cff6+KupP))) ! D
!             cff8  = min(cff2,(cff0*cff2*Bio2d(i,1,iiBen)*Rup/(cff6+KupP))) ! DF
!             cff9  = min(cff3,(cff0*cff3*Bio2d(i,1,iiBen)*Rup/(cff6+KupP))) ! PS
!             cff10 = min(cff4,(cff0*cff4*Bio2d(i,1,iiBen)*Rup/(cff6+KupP))) ! PL
!             cff11 = min(cff5,(cff0*cff5*Bio2d(i,1,iiBen)*Rup/(cff5+KupD))) ! DetBen

            cff7  = cff0*cff1*Bio2d(i,1,iiBen)*Rup/(cff6+KupP) ! D
            cff8  = cff0*cff2*Bio2d(i,1,iiBen)*Rup/(cff6+KupP) ! DF
            cff9  = cff0*cff3*Bio2d(i,1,iiBen)*Rup/(cff6+KupP) ! PS
            cff10 = cff0*cff4*Bio2d(i,1,iiBen)*Rup/(cff6+KupP) ! PL
            cff11 = cff0*cff5*Bio2d(i,1,iiBen)*Rup/(cff5+KupD) ! DetBen

            ! Distribute pelagic feeding losses to appropriate water
            ! column layers

            DO k = 1,N(ng)

              Gra_Det_Ben(i,k)  = cff7  * frac2(k,1) ! mg C m^-2 d^-1
              Gra_DetF_Ben(i,k) = cff8  * frac2(k,1)
              Gra_PhS_Ben(i,k)  = cff9  * frac2(k,1)
              Gra_PhL_Ben(i,k)  = cff10 * frac2(k,1)

            END DO

            ! Benthic feeding takes place in bottom layer for bookkeeping
            ! purposes

            Gra_DetBen_Ben(i,1) = cff11 ! mg C m^-2 d^-1

            ! Assume all excretion occurs in the bottom layer too.  Half
            ! goes to NH4 and half to DetBen

            Exc_Ben_DetBen(i,1) = (eexD * (cff7 + cff8 + cff11) +       &
     &                             eex  * (cff9 + cff10)) * 0.5_r8
            Exc_Ben_NH4(i,1) = Exc_Ben_DetBen(i,1)

            ! Respiration (also takes place in bottom layer)

            cff3 = cff0 * Bio2d(i,1,iiBen) * Rres
            cff4 = ((1_r8 - eexD) * (cff7 + cff8 + cff11) +             &
     &              (1_r8 - eex)  * (cff9 + cff10)) * Qres

            Res_Ben_NH4(i,1) = cff3 + cff4 ! mg C m^-2 d^-1

            ! Mortality (linear senescence and quadratic predation closure)

            Mor_Ben_DetBen(i,1) = cff0*rmort  *Bio2d(i,1,iiBen) +         &
     &                            cff0*BenPred*Bio2d(i,1,iiBen)**2  ! mg C m^-2 d^-1

            ! Benthic remineralization: assumes only the top 25% is
            ! available to remineralize to NH4 (in bottom layer)

            PON = Bio3d(i,1,iiDetBen)*0.25*xi  ! Benthic Particulate organic nitrogen
            cff1 = Pv0*exp(PvT*Temp(i,1))*PON  ! Kawamiya 2000, mmol N m^-3

            Rem_DetBen_NH4(i,1) = cff1*Hz(i,j,1)/xi ! mg C m^-2 d^-1

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
              Par1  = PARs(i)/0.394848_r8  ! surface light, W m^-2

              aiceIfrac = (1-exp(-alphaIb*Par1))*exp(-betaI*Par1) ! light limitation

              cff1 = Bio3d(i,N(ng),iiIceNO3)/(ksnut1 + Bio3d(i,N(ng),iiIceNO3)) ! NO3 limitation
              cff2 = Bio3d(i,N(ng),iiIceNH4)/(ksnut2 + Bio3d(i,N(ng),iiIceNH4)) ! NH4 limitation
              aiceNfrac = cff1*exp(-inhib*Bio3d(i,N(ng),iiIceNH4)) + cff2       ! N limitation
              fNO3      = cff1*exp(-inhib*Bio3d(i,N(ng),iiIceNH4))/aiceNfrac    ! f-ratio
              if (fNO3 /= fNO3) then ! catch NaN if aiceNfrac=0
                fNO3 = 0
              end if

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

              gesi = max(0.0_r8, (1.1e-2                                &
     &                          + 3.012e-2*sb                           &
     &                          + 1.0342e-3*sb**2                       &
     &                          - 4.6033e-5*sb**3                       &
     &                          + 4.926e-7*sb**4                        &
     &                          - 1.659e-9*sb**5              ))

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

              Nit_INH4_INO3(i,N(ng)) = (annit*Bio3d(i,N(ng),iiIceNH4)/xi)*aidz ! mg C m^-2 d^-1

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

          DBio(:,:,iiNO3   ) = (Nit_NH4_NO3                             &
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
     &                       +  Rem_Det_NH4                             &
     &                       +  Rem_DetF_NH4                            &
     &                       +  Exc_Ben_NH4                             &
     &                       +  Res_Ben_NH4                             &
     &                       +  Rem_DetBen_NH4                          &
     &                       +  Twi_INH4_NH4                            &
     &                       -  Gpp_NH4_PhS                             &
     &                       -  Gpp_NH4_PhL                             &
     &                       -  Nit_NH4_NO3)*xi*dtdays ! NH4: mmol N m^-2

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
     &                        -  Rem_Det_NH4                            &
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
     &                        -  Rem_DetF_NH4                           &
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
     &                        +  Gra_DetBen_Ben                         &
     &                        -  Exc_Ben_NH4                            &
     &                        -  Exc_Ben_DetBen                         &
     &                        -  Res_Ben_NH4                            &
     &                        -  Mor_Ben_DetBen)*dtdays ! Ben: mg C m^-2

          DBio(:,:,iiDetBen)  = (Exc_Ben_DetBen                         &
     &                        +  Mor_Ben_DetBen                         &
     &                        -  Gra_DetBen_Ben                         &
     &                        -  Rem_DetBen_NH4)*dtdays ! DetBen: mg C m^-2

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

          DBio(:,:,iiIceNO3)  = (Nit_INH4_INO3                          &
     &                        -  Gpp_INO3_IPhL                          &
     &                        -  Twi_INO3_NO3)*xi*dtdays ! IceNO3: mmol N m^-2

          DBio(:,:,iiIceNH4)  = (Res_IPhL_INH4                          &
     &                        +  Mor_IPhL_INH4                          &
     &                        -  Gpp_INH4_IPhL                          &
     &                        -  Nit_INH4_INO3                          &
     &                        -  Twi_INH4_NH4)*xi*dtdays ! IceNH4: mmol N m^-2


          ! Net production rates (for diagnostics)

          prod_PhS           = Gpp_NO3_PhS                              &
     &                       + Gpp_NH4_PhS                              &
     &                       - Res_PhS_NH4     ! PhS: mg C m^-2 d^-1

          prod_PhL           = Gpp_NO3_PhL                              &
     &                       + Gpp_NH4_PhL                              &
     &                       - Res_PhL_NH4     ! PhL: mg C m^-2 d^-1

          prod_MZL           = Gra_PhS_MZL                              &
     &                       + Gra_PhL_MZL                              &
     &                       - Ege_MZL_Det                              &
     &                       - Res_MZL_NH4     ! MZL: mg C m^-2 d^-1

          prod_Cop           = Gra_PhS_Cop                              &
     &                       + Gra_PhL_Cop                              &
     &                       + Gra_MZL_Cop                              &
     &                       + Gra_IPhL_Cop                             &
     &                       - Ege_Cop_DetF                             &
     &                       - Res_Cop_NH4     ! Cop: mg C m^-2 d^-1

          prod_NCaS          = Gra_PhS_NCaS                             &
     &                       + Gra_PhL_NCaS                             &
     &                       + Gra_MZL_NCaS                             &
     &                       + Gra_IPhL_NCaS                            &
     &                       - Ege_NCaS_DetF                            &
     &                       - Res_NCaS_NH4    ! NCaS: mg C m^-2 d^-1

          prod_EupS          = Gra_PhS_EupS                             &
     &                       + Gra_PhL_EupS                             &
     &                       + Gra_MZL_EupS                             &
     &                       + Gra_Cop_EupS                             &
     &                       + Gra_IPhL_EupS                            &
     &                       + Gra_Det_EupS                             &
     &                       + Gra_DetF_EupS                            &
     &                       - Ege_EupS_DetF                            &
     &                       - Res_EupS_NH4    ! EupS: mg C m^-2 d^-1

          prod_NCaO          = Gra_PhS_NCaO                             &
     &                       + Gra_PhL_NCaO                             &
     &                       + Gra_MZL_NCaO                             &
     &                       + Gra_IPhL_NCaO                            &
     &                       - Ege_NCaO_DetF                            &
     &                       - Res_NCaO_NH4    ! NCaO: mg C m^-2 d^-1

          prod_EupO          = Gra_PhS_EupO                             &
     &                       + Gra_PhL_EupO                             &
     &                       + Gra_MZL_EupO                             &
     &                       + Gra_Cop_EupO                             &
     &                       + Gra_IPhL_EupO                            &
     &                       + Gra_Det_EupO                             &
     &                       + Gra_DetF_EupO                            &
     &                       - Ege_EupO_DetF                            &
     &                       - Res_EupO_NH4    ! EupO: mg C m^-2 d^-1

          prod_Jel            = Gra_Cop_Jel                             &
     &                        + Gra_EupS_Jel                            &
     &                        + Gra_EupO_Jel                            &
     &                        + Gra_NCaS_Jel                            &
     &                        + Gra_NCaO_Jel                            &
     &                        - Ege_Jel_DetF                            &
     &                        - Mor_Jel_DetF                            &
     &                        - Res_Jel_NH4    ! Jel: mg C m^-2 d^-1

          prod_Ben            = Gra_Det_Ben                             &
     &                        + Gra_DetF_Ben                            &
     &                        + Gra_PhS_Ben                             &
     &                        + Gra_PhL_Ben                             &
     &                        + Gra_DetBen_Ben                          &
     &                        - Exc_Ben_NH4                             &
     &                        - Exc_Ben_DetBen                          &
     &                        - Res_Ben_NH4                             &
     &                        - Mor_Ben_DetBen ! Ben: mg C m^-2 d^-1

          prod_IcePhL         = Gpp_INO3_IPhL                           &
     &                        + Gpp_INH4_IPhL                           &
     &                        - Res_IPhL_INH4  ! IcePhL: mg C m^-2 d^-1

          ! Add DBio terms to existing biomass

          Bio2d = Bio2d + DBio

          ! Infauna (Ben) group can receive flux from water column layers.
          ! Move these additions to the bottom layer now, consistent with
          ! the initial setup of the Bio2d and Bio3d arrays.

          Bio2d(:,1,iiBen) = sum(Bio2d(:,:,iiBen), DIM=2)
          Bio2d(:,2:N(ng),iiBen) = 0.0_r8

          ! TODO: Eliminate negatives?  Hopefully processes are
          ! formulated to prevent any, but very fast overturning might
          ! result in numerical issues.  Brute force zero traps will
          ! eliminate conservation of mass, so I'd prefer to look into
          ! increasing BioIter if this is a problem

          Bio2d = max(Bio2d, 0.0_r8)

          ! TODO: Save this "input flux" as diagnostic

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
          ! detritus, large copepod seasonal diapause, and euphausiid
          ! diel vertical migration.

          ! Fraction used: 20% of what hits the bottom becomes
          ! biologically unavailable, and 1% is lost to denitrification.
          ! These fractions apply to PhS, PhL, Det, and DetF.

          fracUsed = 0.79_r8

          ! Initialize temporary arrays to 0

          dBtmp = 0.0_r8
          flxtmp = 0.0_r8

          ! Small phytoplankton: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wPhS, Bio3d(i,:,iiPhS), dBtmp,         &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiPhS) = Bio3d(i,1:N(ng),iiPhS) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen) + flxtmp*fracUsed

            Ver_PhS_DetBen(i,1) = flxtmp*fracUsed/dtdays ! mg C m^-2 d^-1
            Ver_PhS_Out(i,1)    = flxtmp*(1_r8-fracUsed)/dtdays ! mg C m^-2 d^-1
          END DO


          ! Large phytoplankton: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wPhL, Bio3d(i,:,iiPhL), dBtmp,         &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiPhL) = Bio3d(i,1:N(ng),iiPhL) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen) + flxtmp*fracUsed

            Ver_PhL_DetBen(i,1) = (flxtmp*fracUsed)/dtdays ! mg C m^-2 d^-1
            Ver_PhL_Out(i,1)    = (flxtmp*(1_r8-fracUsed))/dtdays ! mg C m^-2 d^-1

          END DO

          ! Slow-sinking detritus: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wDet, Bio3d(i,:,iiDet), dBtmp,         &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiDet) = Bio3d(i,1:N(ng),iiDet) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen) + flxtmp*fracUsed

            Ver_Det_DetBen(i,1) = flxtmp*fracUsed/dtdays ! mg C m^-2 d^-1
            Ver_Det_Out(i,1)    = flxtmp*(1_r8-fracUsed)/dtdays ! mg C m^-2 d^-1

          END DO

          ! Fast-sinking detritus: sinks, and 79% of what sinks out of the
          ! bottom goes to benthic detritus

          DO i=Istr,Iend

            call BioVert(N(ng), -wDetF, Bio3d(i,:,iiDetF), dBtmp,       &
     &                   Hz(i,j,:), dtdays, z_w(i,j,:),                 &
     &                   z_w(i,j,N(ng))+10, flxtmp)
            Bio3d(i,1:N(ng),iiDetF) = Bio3d(i,1:N(ng),iiDetF) + dBtmp(1,1:N(ng))
            Bio2d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen) + flxtmp*fracUsed

            Ver_DetF_DetBen(i,1) = flxtmp*fracUsed/dtdays ! mg C m^-2 d^-1
            Ver_DetF_Out(i,1)    = flxtmp*(1_r8-fracUsed)/dtdays ! mg C m^-2 d^-1

          END DO

#ifdef DIAPAUSE
          ! On-shelf large copepods (NCaS i.e. CM): Move up and down
          ! based on dates set in input file.  Downward movement is
          ! stopped at 200 m or halfway through the bottom layer,
          ! whichever is shallower; upward movement is stopped halfway
          ! through the top layer.

          DO i=Istr,Iend

            if (downwardCM) then

# ifdef DEPTHLIMITER
              if (z_w(i,j,1) < -200.0_r8) then

                ! If water is deeper than 200m, no flux boundary imposed.
                ! Die and go to DetF when cross 200m (this assumes there
                ! was nothing deeper than 200m to begin with... if there
                ! was, that also gets killed off)

                call BioVert(N(ng), -wNCsink, Bio3d(i,:,iiNCaS), dBtmp, &
     &                       Hz(i,j,:), dtdays, z_w(i,j,:),             &
     &                       z_w(i,j,N(ng))+10, flxtmp)

                Bio3d(i,1:N(ng),iiNCaS) = Bio3d(i,1:N(ng),iiNCaS) + dBtmp(1,1:N(ng))
                do k = 1,N(ng)
                  if (z_w(i,j,k) < -200.0_r8) then
                    Ver_NCaS_DetF(i,k) = (Bio3d(i,k,iiNCaS)*Hz(i,j,k))/dtdays ! mg C m^-2 d^-1
                    Bio3d(i,k,iiDetF) = Bio3d(i,k,iiDetF) + Bio3d(i,k,iiNCaS)
                    Bio3d(i,k,iiNCaS) = 0;
                  endif
                enddo

                ! If water depth is just over the 200m limit,
                ! there might be some flux across the bottom boundary.
                ! Send this to DetBen.

                Bio2d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen) + flxtmp

                Ver_NCaS_DetBen(i,1) = flxtmp/dtdays ! TODO: should this be subject to the 21% loss too?
              else

                ! In shallow water, flux boundary at bottom.

                call BioVert(N(ng), -wNCsink, Bio3d(i,:,iiNCaS), dBtmp, &
     &                       Hz(i,j,:), dtdays, z_w(i,j,:),             &
     &                       max((z_w(i,j,0)+z_w(i,j,1))/2, -200.0_r8), &
     &                       flxtmp)
                Bio3d(i,1:N(ng),iiNCaS) = Bio3d(i,1:N(ng),iiNCaS) + dBtmp(1,1:N(ng))

              endif

# else

              call BioVert(N(ng), -wNCsink, Bio3d(i,:,iiNCaS), dBtmp,   &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     max((z_w(i,j,0)+z_w(i,j,1))/2, -200.0_r8), flxtmp)
              Bio3d(i,1:N(ng),iiNCaS) = Bio3d(i,1:N(ng),iiNCaS) + dBtmp(1,1:N(ng))
# endif

            else if (upwardCM) then

              call BioVert(N(ng), wNCrise, Bio3d(i,:,iiNCaS), dBtmp,    &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     (z_w(i,j,N(ng)-1)+z_w(i,j,N(ng)))/2, flxtmp)
              Bio3d(i,1:N(ng),iiNCaS) = Bio3d(i,1:N(ng),iiNCaS) + dBtmp(1,1:N(ng))

            end if

          END DO

          ! Off-shelf large copepods (NCaO i.e. NC): Move up and down
          ! based on dates set in input file.  Downward movement is
          ! stopped at 400 m ; upward movement is stopped halfway
          ! through the top layer. If the water depth is less than 400 m,
          ! copepods are assumed to die when they hit the bottom, and
          ! biomass is transferred to benthic detritus.

          DO i=Istr,Iend

            if (downwardNC) then

              call BioVert(N(ng), -wNCsink, Bio3d(i,:,iiNCaO), dBtmp,   &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     -400.0_r8, flxtmp)
              Bio3d(i,1:N(ng),iiNCaO) = Bio3d(i,1:N(ng),iiNCaO) + dBtmp(1,1:N(ng))
              Bio2d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen) + flxtmp

              Ver_NCaO_DetBen(i,1) = flxtmp/dtdays ! TODO: should this be subject to the 21% loss too?

            else if (upwardNC) then

              call BioVert(N(ng), wNCrise, Bio3d(i,:,iiNCaO), dBtmp,    &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     (z_w(i,j,N(ng)-1)+z_w(i,j,N(ng)))/2, flxtmp)
              Bio3d(i,1:N(ng),iiNCaO) = Bio3d(i,1:N(ng),iiNCaO) + dBtmp(1,1:N(ng))

            end if

          END DO
#endif

#ifdef EUPDIEL

          ! Euphausiid migration controlled by light level.  To achieve
          ! the effect of swimming towards a target depth, I run two
          ! passes, first with sinking and then with rising.

          DO i=Istr,Iend

            ! Still testing... if entire water column is bright or dark,
            ! no migration movement.  If not, target depth is midpoint of
            ! highest layer where PAR < 0.5 E/m^2/d

            if (      (ANY(PAR(i,:) < 0.5_r8)) .and.                    &
     &          (.not. ALL(PAR(i,:) < 0.5_r8))) then

              do k = 1,N(ng)
                if (PAR(i,k) < 0.5_r8) then
                  targetdepth = (z_w(i,j,k-1) + z_w(i,j,k))/2
                endif
              end do

              call BioVert(N(ng), -100.0_r8, Bio3d(i,:,iiEupS), dBtmp,  &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     targetdepth, flxtmp)

              Bio3d(i,1:N(ng),iiEupS) = Bio3d(i,1:N(ng),iiEupS) + dBtmp(1,1:N(ng))

              call BioVert(N(ng), 100.0_r8, Bio3d(i,:,iiEupS), dBtmp,   &
     &                     Hz(i,j,:), dtdays, z_w(i,j,:),               &
     &                     targetdepth, flxtmp)

              Bio3d(i,1:N(ng),iiEupS) = Bio3d(i,1:N(ng),iiEupS) + dBtmp(1,1:N(ng))
            endif
          END DO
#endif

          ! Sync 2d arrays to 3d again

          DO i=Istr,Iend
            ! Sync pelagic (3D modified in sinking portion of code)
            DO k = 1,N(ng)
              DO itrc = 1,NBT ! Pelagic
                Bio2d(i,k,itrc) = Bio3d(i,k,itrc)*Hz(i,j,k)
              END DO
            END DO
            ! Sync benthic (2D modified in sinking portion of code)
            Bio3d(i,1,iiDetBen) = Bio2d(i,1,iiDetBen)/Hz(i,j,1)
          END DO


#ifdef STATIONARY

          DO i = Istr,Iend
            DO k = 1,N(ng)

              ! Intermediate fluxes

              st(i,j,k,nstp,  9) = Gpp_NO3_PhS(i,k)
              st(i,j,k,nstp, 10) = Gpp_NO3_PhL(i,k)
              st(i,j,k,nstp, 11) = Gpp_NH4_PhS(i,k)
              st(i,j,k,nstp, 12) = Gpp_NH4_PhL(i,k)
              st(i,j,k,nstp, 13) = Gra_PhS_MZL(i,k)
              st(i,j,k,nstp, 14) = Gra_PhL_MZL(i,k)
              st(i,j,k,nstp, 15) = Ege_MZL_Det(i,k)
              st(i,j,k,nstp, 16) = Gra_PhS_Cop(i,k)
              st(i,j,k,nstp, 17) = Gra_PhL_Cop(i,k)
              st(i,j,k,nstp, 18) = Gra_MZL_Cop(i,k)
              st(i,j,k,nstp, 19) = Gra_IPhL_Cop(i,k)
              st(i,j,k,nstp, 20) = Ege_Cop_DetF(i,k)
              st(i,j,k,nstp, 21) = Gra_PhS_NCaS(i,k)
              st(i,j,k,nstp, 22) = Gra_PhL_NCaS(i,k)
              st(i,j,k,nstp, 23) = Gra_MZL_NCaS(i,k)
              st(i,j,k,nstp, 24) = Gra_IPhL_NCaS(i,k)
              st(i,j,k,nstp, 25) = Ege_NCaS_DetF(i,k)
              st(i,j,k,nstp, 26) = Gra_PhS_NCaO(i,k)
              st(i,j,k,nstp, 27) = Gra_PhL_NCaO(i,k)
              st(i,j,k,nstp, 28) = Gra_MZL_NCaO(i,k)
              st(i,j,k,nstp, 29) = Gra_IPhL_NCaO(i,k)
              st(i,j,k,nstp, 30) = Ege_NCaO_DetF(i,k)
              st(i,j,k,nstp, 31) = Gra_PhS_EupS(i,k)
              st(i,j,k,nstp, 32) = Gra_PhL_EupS(i,k)
              st(i,j,k,nstp, 33) = Gra_MZL_EupS(i,k)
              st(i,j,k,nstp, 34) = Gra_Cop_EupS(i,k)
              st(i,j,k,nstp, 35) = Gra_IPhL_EupS(i,k)
              st(i,j,k,nstp, 36) = Gra_Det_EupS(i,k)
              st(i,j,k,nstp, 37) = Gra_DetF_EupS(i,k)
              st(i,j,k,nstp, 38) = Ege_EupS_DetF(i,k)
              st(i,j,k,nstp, 39) = Gra_PhS_EupO(i,k)
              st(i,j,k,nstp, 40) = Gra_PhL_EupO(i,k)
              st(i,j,k,nstp, 41) = Gra_MZL_EupO(i,k)
              st(i,j,k,nstp, 42) = Gra_Cop_EupO(i,k)
              st(i,j,k,nstp, 43) = Gra_IPhL_EupO(i,k)
              st(i,j,k,nstp, 44) = Gra_Det_EupO(i,k)
              st(i,j,k,nstp, 45) = Gra_DetF_EupO(i,k)
              st(i,j,k,nstp, 46) = Ege_EupO_DetF(i,k)
              st(i,j,k,nstp, 47) = Gra_Cop_Jel(i,k)
              st(i,j,k,nstp, 48) = Gra_EupS_Jel(i,k)
              st(i,j,k,nstp, 49) = Gra_EupO_Jel(i,k)
              st(i,j,k,nstp, 50) = Gra_NCaS_Jel(i,k)
              st(i,j,k,nstp, 51) = Gra_NCaO_Jel(i,k)
              st(i,j,k,nstp, 52) = Ege_Jel_DetF(i,k)
              st(i,j,k,nstp, 53) = Mor_PhS_Det(i,k)
              st(i,j,k,nstp, 54) = Mor_PhL_Det(i,k)
              st(i,j,k,nstp, 55) = Mor_MZL_Det(i,k)
              st(i,j,k,nstp, 56) = Mor_Cop_DetF(i,k)
              st(i,j,k,nstp, 57) = Mor_NCaS_DetF(i,k)
              st(i,j,k,nstp, 58) = Mor_EupS_DetF(i,k)
              st(i,j,k,nstp, 59) = Mor_NCaO_DetF(i,k)
              st(i,j,k,nstp, 60) = Mor_EupO_DetF(i,k)
              st(i,j,k,nstp, 61) = Mor_Jel_DetF(i,k)
              st(i,j,k,nstp, 62) = Res_PhS_NH4(i,k)
              st(i,j,k,nstp, 63) = Res_PhL_NH4(i,k)
              st(i,j,k,nstp, 64) = Res_MZL_NH4(i,k)
              st(i,j,k,nstp, 65) = Res_Cop_NH4(i,k)
              st(i,j,k,nstp, 66) = Res_NCaS_NH4(i,k)
              st(i,j,k,nstp, 67) = Res_NCaO_NH4(i,k)
              st(i,j,k,nstp, 68) = Res_EupS_NH4(i,k)
              st(i,j,k,nstp, 69) = Res_EupO_NH4(i,k)
              st(i,j,k,nstp, 70) = Res_Jel_NH4(i,k)
              st(i,j,k,nstp, 71) = Rem_Det_NH4(i,k)
              st(i,j,k,nstp, 72) = Rem_DetF_NH4(i,k)
              st(i,j,k,nstp, 73) = Nit_NH4_NO3(i,k)
              st(i,j,k,nstp, 74) = Gra_Det_Ben(i,k)
              st(i,j,k,nstp, 75) = Gra_DetF_Ben(i,k)
              st(i,j,k,nstp, 76) = Gra_PhS_Ben(i,k)
              st(i,j,k,nstp, 77) = Gra_PhL_Ben(i,k)
              st(i,j,k,nstp, 78) = Gra_DetBen_Ben(i,k)
              st(i,j,k,nstp, 79) = Exc_Ben_NH4(i,k)
              st(i,j,k,nstp, 80) = Exc_Ben_DetBen(i,k)
              st(i,j,k,nstp, 81) = Res_Ben_NH4(i,k)
              st(i,j,k,nstp, 82) = Mor_Ben_DetBen(i,k)
              st(i,j,k,nstp, 83) = Rem_DetBen_NH4(i,k)
              st(i,j,k,nstp, 84) = Gpp_INO3_IPhL(i,k)
              st(i,j,k,nstp, 85) = Gpp_INH4_IPhL(i,k)
              st(i,j,k,nstp, 86) = Res_IPhL_INH4(i,k)
              st(i,j,k,nstp, 87) = Mor_IPhL_INH4(i,k)
              st(i,j,k,nstp, 88) = Nit_INH4_INO3(i,k)
              st(i,j,k,nstp, 89) = Twi_IPhL_PhL(i,k)
              st(i,j,k,nstp, 90) = Twi_INO3_NO3(i,k)
              st(i,j,k,nstp, 91) = Twi_INH4_NH4(i,k)
              st(i,j,k,nstp, 92) = Ver_PhS_DetBen(i,k)
              st(i,j,k,nstp, 93) = Ver_PhS_Out(i,k)
              st(i,j,k,nstp, 94) = Ver_PhL_DetBen(i,k)
              st(i,j,k,nstp, 95) = Ver_PhL_Out(i,k)
              st(i,j,k,nstp, 96) = Ver_Det_DetBen(i,k)
              st(i,j,k,nstp, 97) = Ver_Det_Out(i,k)
              st(i,j,k,nstp, 98) = Ver_DetF_DetBen(i,k)
              st(i,j,k,nstp, 99) = Ver_DetF_Out(i,k)
              st(i,j,k,nstp,100) = Ver_NCaO_DetBen(i,k)
              st(i,j,k,nstp,101) = Ver_NCaS_DetF(i,k)
              st(i,j,k,nstp,102) = Ver_NCaS_DetBen(i,k)
              st(i,j,k,nstp,103) = Frz_PhL_IPhL(i,k)
              st(i,j,k,nstp,104) = Frz_NO3_INO3(i,k)
              st(i,j,k,nstp,105) = Frz_NH4_INH4(i,k)

              DO itrc=9,105
                if (st(i,j,k,nstp,itrc) /= st(i,j,k,nstp,itrc)) then
                  write(*, '(A23,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1)') "NaN in flux array: st(", i, ",", j, ",", k, ",", nstp, ",", itrc, ")"
                  exit_flag = 1
                end if
              END DO

              ! Net production

              st(i,j,k,nstp,106) = prod_PhS(i,k)
              st(i,j,k,nstp,107) = prod_PhL(i,k)
              st(i,j,k,nstp,108) = prod_MZL(i,k)
              st(i,j,k,nstp,109) = prod_Cop(i,k)
              st(i,j,k,nstp,110) = prod_NCaS(i,k)
              st(i,j,k,nstp,111) = prod_EupS(i,k)
              st(i,j,k,nstp,112) = prod_NCaO(i,k)
              st(i,j,k,nstp,113) = prod_EupO(i,k)
              st(i,j,k,nstp,114) = prod_Jel(i,k)
              st(i,j,k,nstp,115) = prod_Ben(i,k)
              st(i,j,k,nstp,116) = prod_IcePhL(i,k)

            END DO
          END DO
#endif

        END DO ITER_LOOP

        !=============================================
        !  Update global tracer variables (m Tunits).
        !=============================================

        ! Calculate rate of change due to biogeochemical processes, based
        ! on difference between Bio2d now and Bio_bak from beginning of
        ! routine, and add this to the values in the predictor time step.
        ! (recall that everything in the t(nnew) step is in transport
        ! units, i.e. m*Tunits)

        ! TODO: Georgina's code has a max(t,0) applied to all tracers...
        ! maybe add?

        DO i=Istr,Iend
          DO k = 1,N(ng)

            t(i,j,k,nnew,iNO3 ) = t(i,j,k,nnew,iNO3 ) + (Bio2d(i,k,iiNO3 ) - Bio_bak(i,k,iiNO3 ))
            t(i,j,k,nnew,iNH4 ) = t(i,j,k,nnew,iNH4 ) + (Bio2d(i,k,iiNH4 ) - Bio_bak(i,k,iiNH4 ))
            t(i,j,k,nnew,iPhS ) = t(i,j,k,nnew,iPhS ) + (Bio2d(i,k,iiPhS ) - Bio_bak(i,k,iiPhS ))
            t(i,j,k,nnew,iPhL ) = t(i,j,k,nnew,iPhL ) + (Bio2d(i,k,iiPhL ) - Bio_bak(i,k,iiPhL ))
            t(i,j,k,nnew,iMZL ) = t(i,j,k,nnew,iMZL ) + (Bio2d(i,k,iiMZL ) - Bio_bak(i,k,iiMZL ))
            t(i,j,k,nnew,iCop ) = t(i,j,k,nnew,iCop ) + (Bio2d(i,k,iiCop ) - Bio_bak(i,k,iiCop ))
            t(i,j,k,nnew,iNCaS) = t(i,j,k,nnew,iNCaS) + (Bio2d(i,k,iiNCaS) - Bio_bak(i,k,iiNCaS))
            t(i,j,k,nnew,iEupS) = t(i,j,k,nnew,iEupS) + (Bio2d(i,k,iiEupS) - Bio_bak(i,k,iiEupS))
            t(i,j,k,nnew,iNCaO) = t(i,j,k,nnew,iNCaO) + (Bio2d(i,k,iiNCaO) - Bio_bak(i,k,iiNCaO))
            t(i,j,k,nnew,iEupO) = t(i,j,k,nnew,iEupO) + (Bio2d(i,k,iiEupO) - Bio_bak(i,k,iiEupO))
            t(i,j,k,nnew,iDet ) = t(i,j,k,nnew,iDet ) + (Bio2d(i,k,iiDet ) - Bio_bak(i,k,iiDet ))
            t(i,j,k,nnew,iDetF) = t(i,j,k,nnew,iDetF) + (Bio2d(i,k,iiDetF) - Bio_bak(i,k,iiDetF))
            t(i,j,k,nnew,iJel ) = t(i,j,k,nnew,iJel ) + (Bio2d(i,k,iiJel ) - Bio_bak(i,k,iiJel ))
            t(i,j,k,nnew,iFe  ) = t(i,j,k,nnew,iFe  ) + (Bio2d(i,k,iiFe  ) - Bio_bak(i,k,iiFe  ))
            t(i,j,k,nnew,iMZS ) = t(i,j,k,nnew,iMZS ) + 0.0_r8

            ! Check for negatives and NaNs (for debugging)

!             do itrc = iNO3,size(t,5)
!               if (t(i,j,k,nnew,itrc) < 0) then
!                 write(*, '(A19,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1)') "Negative tracer: t(", i, ",", j, ",", k, ",", nnew, ",", itrc, ")"
!               end if
!
!               if (t(i,j,k,nnew,itrc) /= t(i,j,k,nnew,itrc)) then
!                 write(*, '(A23,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1)') "NaN in tracer array: t(", i, ",", j, ",", k, ",", nnew, ",", itrc, ")"
!               end if
!             end do

#ifdef STATIONARY
            ! Add 2D tracer values to diagnostic array to allow
            ! calculation of advective-diffusive fluxes.  A bit silly to
            ! have these as extra outputs, but sticking the values here
            ! is a lot easier than creating and allocating entirely new
            ! arrays for it.
            ! (Note: step3d_t.F copies the nstp values to nnew, so when
            ! bestnpz.h is next called, the nnew values will hold the
            ! previous step... I think)

            st(i,j,k,nstp,117) = t(i,j,k,nnew,iNO3 )
            st(i,j,k,nstp,118) = t(i,j,k,nnew,iNH4 )
            st(i,j,k,nstp,119) = t(i,j,k,nnew,iPhS )
            st(i,j,k,nstp,120) = t(i,j,k,nnew,iPhL )
            st(i,j,k,nstp,121) = t(i,j,k,nnew,iMZL )
            st(i,j,k,nstp,122) = t(i,j,k,nnew,iCop )
            st(i,j,k,nstp,123) = t(i,j,k,nnew,iNCaS)
            st(i,j,k,nstp,124) = t(i,j,k,nnew,iEupS)
            st(i,j,k,nstp,125) = t(i,j,k,nnew,iNCaO)
            st(i,j,k,nstp,126) = t(i,j,k,nnew,iEupO)
            st(i,j,k,nstp,127) = t(i,j,k,nnew,iDet )
            st(i,j,k,nstp,128) = t(i,j,k,nnew,iDetF)
            st(i,j,k,nstp,129) = t(i,j,k,nnew,iJel )
            st(i,j,k,nstp,130) = t(i,j,k,nnew,iFe  )
#endif

#ifdef TS_MPDATA
            t(i,j,k,3,iNO3 ) = t(i,j,k,nnew,iNO3 ) * Hz_inv(i,k)
            t(i,j,k,3,iNH4 ) = t(i,j,k,nnew,iNH4 ) * Hz_inv(i,k)
            t(i,j,k,3,iPhS ) = t(i,j,k,nnew,iPhS ) * Hz_inv(i,k)
            t(i,j,k,3,iPhL ) = t(i,j,k,nnew,iPhL ) * Hz_inv(i,k)
            t(i,j,k,3,iMZL ) = t(i,j,k,nnew,iMZL ) * Hz_inv(i,k)
            t(i,j,k,3,iCop ) = t(i,j,k,nnew,iCop ) * Hz_inv(i,k)
            t(i,j,k,3,iNCaS) = t(i,j,k,nnew,iNCaS) * Hz_inv(i,k)
            t(i,j,k,3,iEupS) = t(i,j,k,nnew,iEupS) * Hz_inv(i,k)
            t(i,j,k,3,iNCaO) = t(i,j,k,nnew,iNCaO) * Hz_inv(i,k)
            t(i,j,k,3,iEupO) = t(i,j,k,nnew,iEupO) * Hz_inv(i,k)
            t(i,j,k,3,iDet ) = t(i,j,k,nnew,iDet ) * Hz_inv(i,k)
            t(i,j,k,3,iDetF) = t(i,j,k,nnew,iDetF) * Hz_inv(i,k)
            t(i,j,k,3,iJel ) = t(i,j,k,nnew,iJel ) * Hz_inv(i,k)
            t(i,j,k,3,iFe  ) = t(i,j,k,nnew,iFe  ) * Hz_inv(i,k)
            t(i,j,k,3,iMZS ) = 0
#endif
          END DO
        END DO

#ifdef BENTHIC

        DO i=Istr,Iend
          bt(i,j,1,nnew,iBen   ) = bt(i,j,1,nnew,iBen   ) + (Bio2d(i,1,iiBen   ) - Bio_bak(i,1,iiBen   ))
          bt(i,j,1,nnew,iDetBen) = bt(i,j,1,nnew,iDetBen) + (Bio2d(i,1,iiDetBen) - Bio_bak(i,1,iiDetBen))

# ifdef MASKING
          bt(i,j,1,nnew,iBen)    = bt(i,j,1,nnew,iBen   )*rmask(i,j)
          bt(i,j,1,nnew,iDetBen) = bt(i,j,1,nnew,iDetBen)*rmask(i,j)
# endif
        END DO
#endif

#ifdef ICE_BIO
        ! Note: Ice variables aren't subject to the same tranport
        ! equations as pelagic, so these use tracer units rather than
        ! tranport units in the nnew timestep.
# if defined CLIM_ICE_1D

        DO i=Istr,Iend
          it(i,j,nnew,iIceNO3) = it(i,j,nnew,iIceNO3) + (Bio2d(i,N(ng),iiIceNO3) - Bio_bak(i,N(ng),iiIceNO3))/aidz
          it(i,j,nnew,iIceNH4) = it(i,j,nnew,iIceNH4) + (Bio2d(i,N(ng),iiIceNH4) - Bio_bak(i,N(ng),iiIceNH4))/aidz
          it(i,j,nnew,iIcePhL) = it(i,j,nnew,iIcePhL) + (Bio2d(i,N(ng),iiIcePhL) - Bio_bak(i,N(ng),iiIcePhL))/aidz

          itL(i,j,nnew,iIceLog) = itL(i,j,nstp,iIceLog)

#  ifdef MASKING
          it( i,j,nnew,iIcePhL) = it( i,j,nnew,iIcePhL)*rmask(i,j)
          it( i,j,nnew,iIceNH4) = it( i,j,nnew,iIceNH4)*rmask(i,j)
          it( i,j,nnew,iIceNO3) = it( i,j,nnew,iIceNO3)*rmask(i,j)
          itL(i,j,nnew,iIceLog) = itL(i,j,nnew,iIceLog)*rmask(i,j)
#  endif
        END DO
# else
        DO i=Istr,Iend
          IceNO3(i,j,nnew) = IceNO3(i,j,nnew) + (Bio2d(i,N(ng),iiIceNO3) - Bio_bak(i,N(ng),iiIceNO3))/aidz
          IceNH4(i,j,nnew) = IceNH4(i,j,nnew) + (Bio2d(i,N(ng),iiIceNH4) - Bio_bak(i,N(ng),iiIceNH4))/aidz
          IcePhL(i,j,nnew) = IcePhL(i,j,nnew) + (Bio2d(i,N(ng),iiIcePhL) - Bio_bak(i,N(ng),iiIcePhL))/aidz

          IceLog(i,j,nnew) = IceLog(i,j,nstp) ! TODO: Current step value now in both positions... doublecheck that this is correct (real update in =happens in ice_limit.F))

#  ifdef MASKING
          IcePhL(i,j,nnew) = IcePhL(i,j,nnew)*rmask(i,j)
          IceNO3(i,j,nnew) = IceNO3(i,j,nnew)*rmask(i,j)
          IceNH4(i,j,nnew) = IceNH4(i,j,nnew)*rmask(i,j)
#  endif
        END DO

# endif
#endif

      END DO J_LOOP

      !=============================================
      !  FEAST
      !=============================================

#ifdef FEAST
      ! This block has something to do with updating the ghost points
      ! along the tile edges prior to running the FEAST advection code
      ! (which defines fish movement, independent from water advection).
      ! The history is a little murky...
      !
      ! TODO: maybe separate the fish movement code from the source/sinks
      ! code?  So the movement doesn't require this extra kludginess?
      ! (Copy/pasting the physics here seems like a recipe for disaster
      ! if we ever want to move to a more recent version of ROMS where
      ! these things may change slightly).

# if defined EW_PERIODIC || defined NS_PERIODIC

      ! Apply periodic boundary conditions.

      DO itrc=1,NBT
        ibio=idbio(itrc)

        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,ibio))
      END DO
# endif

# ifdef DISTRIBUTE
      ! Exchange boundary data.

      !ajh
      !added block on this for passives when feast is present
      !NOTE will need to exchange when passives/fish are changed elsewhere
#  ifdef FEAST_NOEXCHANGE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NAT,         &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,:))

      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NAT+NPT+1, NT(ng),                            &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,:))

#  else
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    t(:,:,:,nnew,:))
#  endif
!
#  ifdef STATIONARY
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NBTS,        &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    st(:,:,:,nnew,:))
#  endif
# endif

# if defined ICE_BIO
#  ifdef BERING_10K
#   ifndef MATLABCOMPILE

      CALL IcePhLbc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp, nnew,                                       &
     &                ui, vi, IcePhL)
      CALL IceNO3bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp, nnew,                                       &
     &                ui, vi, IceNO3)
      CALL IceNH4bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                nstp, nnew,                                       &
     &                ui, vi, IceNH4)

#    if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    IcePhL(:,:,nnew), IceNO3(:,:,nnew),           &
     &                    IceNH4(:,:,nnew))
#    endif
#   endif
#  endif
# endif

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









