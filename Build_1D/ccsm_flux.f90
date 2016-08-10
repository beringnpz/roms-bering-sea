      MODULE ccsm_flux_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!       CCSM's flux_mod.F90 by B. Kauffman                             !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC  :: ccsm_flux
contains
!
!***********************************************************************
      SUBROUTINE ccsm_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private storage
!  arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J)- and MAX(I,J)-directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 17)
      CALL ccsm_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
     &                     MIXING(ng) % alpha,                          &
     &                     MIXING(ng) % beta,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
     &                     FORCES(ng) % stflx,                          &
     &                     GRID(ng) % latr,                             &
     &                     FORCES(ng) % evap,                           &
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
      CALL wclock_off (ng, iNLM, 17)
      RETURN
      END SUBROUTINE ccsm_flux
!
!***********************************************************************
      SUBROUTINE ccsm_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
     &                           alpha, beta, rho, t, u, v,             &
     &                           Hair, Pair, Tair, Uwind, Vwind,        &
     &                           rain, lhflx, lrflx, shflx,             &
     &                           srflx, stflx,                          &
     &                           latr,                                  &
     &                           evap,                                  &
     &                           sustr, svstr)
!***********************************************************************
!     call srfflx_ao( nloc    ,   z   ,   u   ,   v   ,ptem   ,         &
!     &              shum     ,dens   ,uocn   ,vocn   ,                 &
!     &              tocn     ,mask   , shflx   , lhflx   ,lwup   ,     &
!     &              evap     ,taux   , tauy  ,tref   ,qref   ,         &
!     &              duu10n                                     )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: evap(LBi:,LBj:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k
      integer :: IterMax
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
! ratio of molecular weight of water to dry air
      real(r8), parameter :: epsilon = 0.622_r8
      real(r8) :: Bf, Ce, Ch, Hl, Hlw, Hscale, Hscale2, Hs, Hsr, IER
      real(r8) :: PairM,  RH, Taur, cff1, cff2
      real(r8) :: Wspeed, ZQoL, ZToL, cff
      real(r8) :: caw_curve
! ------------------------------------------------------------------
! !DESCRIPTION:
!     wrapper to atm/ocn flux calculation
!
! !REMARKS:
!     All data must be on the ocean domain (note: a domain includes a
!     particular decomposition).
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - first version
!
! ------------------------------------------------------------------
!
! !IROUTINE: srfflx_ao -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - brought in from cpl5.
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!
! !INTERFACE: ------------------------------------------------------------------
!SUBROUTINE srfflx_ao_tile(blk_ZX  ,Uwind  ,Vwind  ,Tair ,              &
!          &         Hair  ,rhoAir  ,u    ,v    ,                       &
!          &         t    ,rmask  ,shflx   ,lhflx   ,lwup  ,            &
!          &         evap  ,taux  ,tauy  ,tref  ,qref  ,                &
!          &         duu10n                                )
! !INPUT/OUTPUT PARAMETERS:
   !--- input arguments --------------------------------
!  integer(IN),intent(in) :: rmask(imax) ! ocn domain mask 0
!  real(R8)   ,intent(in) :: blk_ZX (imax) ! atm level height      (m)
!  real(R8)   ,intent(in) :: Uwind (imax) ! atm u wind            (m/s)
!  real(R8)   ,intent(in) :: Vwind (imax) ! atm v wind            (m/s)
!  real(R8)   ,intent(in) :: Tair(imax) ! atm potential T       (K)
!  real(R8)   ,intent(in) :: Hair (imax) ! atm specific humidity (kg/kg)
!  real(R8)   ,intent(in) :: rhoAir (imax) ! atm air density       (kg/m^3)
!  real(R8)   ,intent(in) :: u   (imax) ! ocn u-velocity        (m/s)
!  real(R8)   ,intent(in) :: v   (imax) ! ocn v-velocity        (m/s)
!  real(R8)   ,intent(in) :: t   (imax) ! ocn temperature       (K)
   !--- output arguments -------------------------------
!  real(R8),intent(out)  ::  shflx  (imax) ! heat flux: sensible    (W/m^2)
!  real(R8),intent(out)  ::  lhflx  (imax) ! heat flux: latent      (W/m^2)
!  real(R8),intent(out)  ::  lwup (imax) ! heat flux: lw upward   (W/m^2)
!  real(R8),intent(out)  ::  evap (imax) ! water flux: evap  ((kg/s)/m^2)
!  real(R8),intent(out)  ::  taux (imax) ! surface stress, zonal      (N)
!  real(R8),intent(out)  ::  tauy (imax) ! surface stress, maridional (N)
!  real(R8),intent(out)  ::  tref (imax) ! diag:  2m ref height T     (K)
!  real(R8),intent(out)  ::  qref (imax) ! diag:  2m ref humidity (kg/kg)
!  real(R8),intent(out)  :: duu10n(imax) ! diag: 10m wind speed squared (m/s)^2
!EOP
   !--- local constants --------------------------------
      real(R8),parameter :: umin  =  0.5    ! minimum wind speed       (m/s)
      real(R8),parameter :: zref  = 10.0    ! reference height           (m)
      real(R8),parameter :: ztref =  2.0    ! reference height for air T (m)
   !--- local variables --------------------------------
      real(r8)    :: lwup   ! upward longwave radiation
      real(r8)    :: vmag   ! surface wind magnitude   (m/s)
      real(r8)    :: thvbot ! virtual temperature      (K)
      real(r8)    :: ssq    ! sea surface humidity     (kg/kg)
!ajh
      real(r8)    :: TseaK  ! SST in Kelvins           (K)
      real(r8)    :: delt   ! potential T difference   (K)
      real(r8)    :: delq   ! humidity difference      (kg/kg)
      real(r8)    :: stable ! stability factor
      real(r8)    :: rdn    ! sqrt of neutral exchange coeff (momentum)
      real(r8)    :: rhn    ! sqrt of neutral exchange coeff (heat)
      real(r8)    :: ren    ! sqrt of neutral exchange coeff (water)
      real(r8)    :: rd     ! sqrt of exchange coefficient (momentum)
      real(r8)    :: re     ! sqrt of exchange coefficient (water)
      real(r8)    :: ustar  ! ustar
      real(r8)    :: qstar  ! qstar
      real(r8)    :: tstar  ! tstar
      real(r8)    :: hol    ! H (at blk_ZX) over L
      real(r8)    :: xsq    ! ?
      real(r8)    :: xqq    ! ?
      real(r8)    :: psimh  ! stability function at blk_ZX (momentum)
      real(r8)    :: psixh  ! stability function at blk_ZX (heat and water)
      real(r8)    :: psix2  ! stability function at ztref reference height
      real(r8)    :: alz    ! ln(blk_ZX/zref)
      real(r8)    :: al2    ! ln(zref/ztref)
      real(r8)    :: u10n   ! 10m neutral wind
      real(r8)    :: tau    ! stress at blk_ZX
      real(r8)    :: cpair     ! specific heat of moist air
      real(r8)    :: bn     ! exchange coef funct for interpolation
      real(r8)    :: bh     ! exchange coef funct for interpolation
      real(r8)    :: fac    ! vertical interpolation factor
      real(r8)    :: Hlv    ! Latent heat of vaporization at sea surface
      real(r8), dimension(IminS:ImaxS) :: rhoAir
      real(r8), dimension(IminS:ImaxS) :: TairK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SRad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy
   !--- local functions --------------------------------
      real(r8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
      real(r8)    :: cdn    ! function: neutral drag coeff at 10m
      real(r8)    :: psimhu ! function: unstable part of psimh
      real(r8)    :: psixhu ! function: unstable part of psimx
      real(r8)    :: tsk    ! water temperature (K)
      real(r8)    :: Tk     ! dummy arg ~ air temperature (K)
      real(r8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
      real(r8)    :: xd     ! dummy arg ~ ?
!  Stefan-Boltzmann constant ~ W/m^2/K^4
      real(r8),parameter :: shr_const_stebol  = 5.67e-8_r8
! Boltzmann's constant ~ J/K/molecule
      real(r8),parameter :: shr_const_boltz   = 1.38065e-23_r8
! Avogadro's number ~ molecules/kmole
      real(r8),parameter :: shr_const_avogad  = 6.02214e26_r8
! Universal gas constant ~ J/K/kmole
      real(r8),parameter :: shr_const_rgas =                            &
     &         shr_const_avogad*shr_const_boltz
! molecular weight dry air ~ kg/kmole
      real(r8),parameter :: shr_const_mwdair  = 28.966_r8
! molecular weight water vapor
      real(r8),parameter :: shr_const_mwwv = 18.016_r8
! Water vapor gas constant ~ J/K/kg
      real(r8),parameter :: shr_const_rwv =                             &
     &         shr_const_rgas/shr_const_mwwv
! Dry air gas constant     ~ J/K/kg
      real(r8),parameter :: shr_const_rdair   =                         &
     &         shr_const_rgas/shr_const_mwdair
      real(r8),parameter :: shr_const_zvir    =                         &
     &         (shr_const_rwv/shr_const_rdair)-1.0_r8
! Von Karman constant
      real(r8),parameter :: shr_const_karman  = 0.4_r8
! specific heat of dry air   ~ J/kg/K
      real(r8),parameter :: shr_const_cpdair  = 1.00464e3_r8
! specific heat of water vap ~ J/kg/K
      real(r8),parameter :: shr_const_cpwv    = 1.810e3_R8
! CPWV/CPDAIR - 1.0
      real(r8),parameter :: shr_const_cpvir   =                         &
     &         (shr_const_cpwv/shr_const_cpdair)-1.0_r8
! latent heat of evaporation ~ J/kg
      real(r8),parameter :: shr_const_latvap  = 2.501e6_r8
      qsat(Tk)   = 640380.0 / exp(5107.4/Tk)
      cdn(Umps)  = 0.0027 / Umps + 0.000142 + 0.0000764 * Umps
      psimhu(xd) = log((1.0+xd*(2.0+xd))*(1.0+xd*xd)/8.0) - 2.0*atan(xd) + 1.571
      psixhu(xd) = 2.0 * log((1.0 + xd*xd)/2.0)
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!-------------------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES:
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + shflx + lhflx
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!
! ASSUMPTIONS:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-------------------------------------------------------------------------------
      al2 = log(zref/ztref)
      Hscale=rho0*Cp
      Hscale2=1.0_r8/(rho0*Cp)
! Note that this loop needs to be cleaned of all global arrays for
! OpenMP.
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          !--- compute some needed quantities ---
          vmag   = max(umin, sqrt( (Uwind(i,j)-u(i,j,N(ng),nrhs))**2 +  &
     &                             (Vwind(i,j)-v(i,j,N(ng),nrhs))**2) )
          TairK(i) = Tair(i,j) + 273.16_r8
!ajh
          TseaK = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          thvbot = TairK(i) * (1.0 + shr_const_zvir * Hair(i,j)) ! virtual temp (K)
          rhoAir(i) = Pair(i,j)*100.0_r8/(blk_Rgas*TairK(i)*            &
     &                    (1.0_r8+0.61_r8*Hair(i,j))) !  Moist air density (kg/m3).
!ajh make Kate's fix
!          ssq    = 0.98_r8 * qsat(TairK(i)) / rhoAir(i)    ! sea surf hum (kg/kg)
          ssq    = 0.98_r8 * qsat(TseaK) / rhoAir(i)    ! sea surf hum (kg/kg)
          delt   = Tair(i,j) - t(i,j,N(ng),nrhs,itemp)                  ! pot temp diff (K)
          delq   = Hair(i,j) - ssq                     ! spec hum dif (kg/kg)
          alz    = log(blk_ZW(ng)/zref)
          cpair     = shr_const_cpdair*(1.0 + shr_const_cpvir*ssq)
          !------------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !------------------------------------------------------------
          !--- neutral coefficients, z/L = 0.0 ---
          stable = 0.5_r8 + sign(0.5_r8 , delt)
          rdn    = sqrt(cdn(vmag))
          rhn    = (1.0_r8-stable) * 0.0327_r8 + stable * 0.018_r8
          ren    = 0.0346_r8
          !--- ustar, tstar, qstar ---
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq
          !--- compute stability & evaluate all stability functions ---
          hol  = shr_const_karman*g*blk_ZT(ng)*                         &
     &          (tstar/thvbot+qstar/(1.0_r8/shr_const_zvir+Hair(i,j)))/ &
     &                     ustar**2
          hol  = sign( min(abs(hol),10.0_r8), hol )
          stable = 0.5_r8 + sign(0.5_r8 , hol)
          xsq    = max(sqrt(abs(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
          xqq    = sqrt(xsq)
          psimh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psimhu(xqq)
          psixh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psixhu(xqq)
          !--- shift wind speed using old coefficient ---
          rd   = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          u10n = vmag * rd / rdn
          !--- update transfer coeffs at 10m and neutral stability ---
          rdn = sqrt(cdn(u10n))
          ren = 0.0346_r8
          rhn = (1.0_r8-stable)*0.0327_r8 + stable * 0.018_r8
          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          rh = rhn / (1.0_r8 + rhn/shr_const_karman*(alz-psixh))
          re = ren / (1.0_r8 + ren/shr_const_karman*(alz-psixh))
          !--- update ustar, tstar, qstar using updated, shifted coeffs --
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq
          !------------------------------------------------------------
          ! iterate to converge on Z/L, ustar, tstar and qstar
          !------------------------------------------------------------
          !--- compute stability & evaluate all stability functions ---
          hol  = shr_const_karman*g*blk_ZQ(ng)*                         &
     &            (tstar/thvbot+qstar/(1.0/shr_const_zvir+Hair(i,j)))/  &
     &                         ustar**2
          hol  = sign( min(abs(hol),10.0_r8), hol )
          stable = 0.5_r8 + sign(0.5_r8 , hol)
          xsq    = max(sqrt(abs(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
          xqq    = sqrt(xsq)
          psimh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psimhu(xqq)
          psixh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psixhu(xqq)
          !--- shift wind speed using old coeffs ---
          rd   = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          u10n = vmag * rd/rdn
          !--- update transfer coeffs at 10m and neutral stability ---
          rdn = sqrt(cdn(u10n))
          ren = 0.0346_r8
          rhn = (1.0_r8 - stable)*0.0327_r8 + stable * 0.018_r8
          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          rh = rhn / (1.0_r8 + rhn/shr_const_karman*(alz-psixh))
          re = ren / (1.0_r8 + ren/shr_const_karman*(alz-psixh))
          !--- update ustar, tstar, qstar using updated, shifted coeffs ---
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq
          !------------------------------------------------------------
          ! compute the fluxes
          !------------------------------------------------------------
          tau = rhoAir(i) * ustar * ustar
	  SRad(i,j) = srflx(i,j) * Hscale
!g	   print*,'i=',i,'j=',j,'srflx(i,j)=',srflx(i,j)
          !--- momentum flux ---
          Taux(i,j) = tau * (Uwind(i,j)-u(i,j,N(ng),nrhs)) / vmag
          Tauy(i,j) = tau * (Vwind(i,j)-v(i,j,N(ng),nrhs)) / vmag
          !--- heat flux ---
          SHeat (i,j) = cpair * tau * tstar / ustar
          LHeat (i,j) = shr_const_latvap * tau * qstar / ustar
        END DO
      END DO
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!ajh apr 23 2011 try increasing latent and sensible heat flux
          shflx (i,j) = 1.25*SHeat(i,j)
!          shflx (i,j) = SHeat(i,j)
          lhflx (i,j) = 1.25*LHeat(i,j)
!          lhflx (i,j) = LHeat(i,j)
          tsk = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          lwup = -shr_const_stebol * tsk**4
          caw_curve = 1._r8-(0.069_r8 - 0.011_r8*cos(2*deg2rad*56.877_r8))
          srflx(i,j) = srflx(i,j)*caw_curve
          lrflx(i,j) = lwup*Hscale2 + lrflx(i,j)
          lhflx(i,j) = lhflx(i,j)*Hscale2
          shflx(i,j) = shflx(i,j)*Hscale2
          stflx(i,j,itemp) = (srflx(i,j)+lrflx(i,j)+                    &
     &                      lhflx(i,j)+shflx(i,j))
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
          Hlv = (2.501_r8-0.00237_r8*t(i,j,N(ng),nrhs,itemp))*1.0e+6_r8
          evap(i,j) = -LHeat(i,j)/Hlv
          stflx(i,j,isalt) = cff*(evap(i,j)-rain(i,j))
! Limit shortwave radiation to be used in computing penetrative 
! radiation by presence of ice and biology
!1D case when ice is read from climatology
!--- water flux ---
!         evap(i,j) = lhflx(i,j)/shr_const_latvap
          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
!         hol = hol*ztref/blk_ZQ(ng)
!         xsq = max( 1.0, sqrt(abs(1.0-16.0*hol)) )
!         xqq = sqrt(xsq)
!         psix2   = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
!         fac     = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
!         tref(i,j) = TairK - delt*fac
!         tref(i,j) = tref(i,j) - 0.01*ztref   ! pot temp to temp correction
!         fac     = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
!         qref(i,j) =  Hair(i,j) - delq*fac
!         duu10n(i,j) = u10n*u10n ! 10m wind speed squared
        END DO
      END DO
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      Hscale = 1./rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=Taux(i,j)*Hscale
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=Tauy(i,j)*Hscale
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        lrflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        lhflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        shflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        stflx(:,:,itemp))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        evap)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        stflx(:,:,isalt))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        sustr)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        svstr)
      RETURN
      END SUBROUTINE ccsm_flux_tile
      END module ccsm_flux_mod
