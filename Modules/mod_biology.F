#include "cppdefs.h"
      MODULE mod_biology
#ifdef BIOLOGY
!
!svn $Id: mod_biology.F 966 2009-03-28 00:29:06Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================

# ifdef BIO_FENNEL
!                                                                      !
!  Parameters for Fennel et al (2006) model:                           !
!                                                                      !
!   AttSW    Light attenuation due to sea water [1/m].                 !
!   AttChl   Light attenuation by Chlorophyll [1/(mg_Chl m2)].         !
!   BioIter  Maximum number of iterations to achieve convergence       !
!              of the nonlinear solution.                              !
!   Chl2C_m  Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].        !
!   ChlMin   Chlorophill minimum threshold value [mg_Chl/m3].          !
!   CoagR    Coagulation rate: agregation rate of SDeN + Phyt ==> LDeN !
!              [1/day].                                                !
!   D_p5NH4  Half-saturation radiation for nitrification inhibition    !
!              [Watts/m2].                                             !
!   I_thNH4  Radiation threshold for nitrification inhibition          !
!              [Watts/m2].                                             !
!   K_NH4    Inverse half-saturation for Phytoplankton NH4 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_NO3    Inverse half-saturation for Phytoplankton NO3 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_Phy    Zooplankton half-saturation, squared constant for         !
!              ingestion [mmol_N/m3]^2.                                !
!   LDeRR    Large Detrital re-mineralization rate [1/day].            !
!   NitriR   Nitrification rate: oxidation of NH4 to NO3 [1/day].      !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   PhyCN    Phytoplankton Carbon:Nitrogen ratio [mol_C/mol_N].        !
!   PhyIP    Phytoplankton NH4 inhibition parameter [1/(mmol_N)].      !
!   PhyIS    Phytoplankton, initial slope of the P-I curve             !
!              [mg_C/(mg_Chl W m-2 day)].                              !
!   ZooMin   Phytoplankton minimum threshold value [mmol_N/m3].        !
!   PhyMR    Phytoplankton mortality rate [1/day] to small detritus.   !
!   SDeAR    Small detritus aggregation rate into Large detritus       !
!              [1/day].                                                !
!   SDeBR    Small Detrital breakdown to NH4 rate [1/day].             !
!   SDeRR    Large Detrital re-mineralization rate [1/day].            !
!   Vp0      Eppley temperature-limited and light-limited growth       !
!              tuning parameter [nondimensional].                      !
!   wLDet    Vertical sinking velocities for Large Detritus            !
!              fraction [m/day].                                       !
!   wPhy     Vertical sinking velocity for Phytoplankton               !
!              fraction [m/day].                                       !
!   wSDet    Vertical sinking velocities for Small Detritus            !
!              fraction [m/day].                                       !
!   ZooAE_N  Zooplankton nitrogen assimilation efficiency fraction     !
!              [nondimensional].                                       !
!   ZooBM    Zooplankton basal metabolism [1/day].                     !
!   ZooCN    Zooplankton Carbon:Nitrogen ratio [mol_C/mol_N].          !
!   ZooER    Zooplankton specific excretion rate [1/day].              !
!   ZooGR    Zooplankton maximum growth rate [1/day].                  !
!   ZooMin   Zooplankton minimum threshold value [mmol_N/m3].          !
!   ZooMR    Zooplankton mortality to Detritus [1/day].                !
!   pCO2air  CO2 partial pressure in the air [ppmv].                   !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: AttChl        ! 1/(mg_Chl m2)
        real(r8), dimension(Ngrids) :: Chl2C_m       ! mg_Chl/mg_C
        real(r8), dimension(Ngrids) :: ChlMin        ! mg_Chl/m3
        real(r8), dimension(Ngrids) :: CoagR         ! 1/day
        real(r8), dimension(Ngrids) :: D_p5NH4       ! Watts/m2
        real(r8), dimension(Ngrids) :: I_thNH4       ! Watts/m2
        real(r8), dimension(Ngrids) :: K_NH4         ! m3/mmol_N
        real(r8), dimension(Ngrids) :: K_NO3         ! m3/mmol_N
        real(r8), dimension(Ngrids) :: K_Phy         ! (mmol_N/m3)^2
        real(r8), dimension(Ngrids) :: LDeRRN        ! 1/day
        real(r8), dimension(Ngrids) :: LDeRRC        ! 1/day
        real(r8), dimension(Ngrids) :: NitriR        ! 1/day
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: PhyCN         ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: PhyIP         ! 1/mmol_N
        real(r8), dimension(Ngrids) :: PhyIS         ! 1/(Watts m-2 day)
        real(r8), dimension(Ngrids) :: PhyMin        ! mmol_N/m3
        real(r8), dimension(Ngrids) :: PhyMR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeAR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeBR         ! 1/day
        real(r8), dimension(Ngrids) :: SDeRRN        ! 1/day
        real(r8), dimension(Ngrids) :: SDeRRC        ! 1/day
        real(r8), dimension(Ngrids) :: Vp0           ! nondimensional
        real(r8), dimension(Ngrids) :: wLDet         ! m/day 
        real(r8), dimension(Ngrids) :: wPhy          ! m/day
        real(r8), dimension(Ngrids) :: wSDet         ! m/day
        real(r8), dimension(Ngrids) :: ZooAE_N       ! nondimensional
        real(r8), dimension(Ngrids) :: ZooBM         ! 1/day
        real(r8), dimension(Ngrids) :: ZooCN         ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: ZooER         ! 1/day
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMin        ! mmol_N/m3
        real(r8), dimension(Ngrids) :: ZooMR         ! 1/day
        real(r8), dimension(Ngrids) :: pCO2air       ! ppmv

# endif
# ifdef BIO_UMaine
!                                                                      !
!  Parameters for UMaine CoSiNE model:                     !
!                                                                      !
!   reg1     Microzooplankton excretion rate to ammonium [1/day].      !
!   reg2     Mesozooplankton excretion rate to ammonium [1/day].       !
!   gmaxs1   Maximum specific growth rate of small phytoplankton       !
!              [1/day]                                                 !
!   gmaxs2   Maximum specific growth rate of diatom [1/day]            !
!   beta1    Microzooplankton maximum grazing rate [1/day]             !
!   beta2    Mesozooplankton maximum grazing rate [1/day]              !
!   akz1     Half saturation constant for microzooplankton grazing     !
!              [mmol_N/m3]                                             !
!   akz2     Half saturation constant for mesozooplankton grazing      !
!              [mmol_N/m3]                                             !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   amaxs2   Initial slope of P-I curve of small phytoplankton         !
!              [1/(Watts/m2)/day]                                      !
!   parsats1 PAR saturation onset parameter of small phytoplankton     !
!              [Watts/m2]                                              !
!   parsats2 PAR saturation onset parameter of diatom [Watts/m2]       !
!              [Watts/m2]                                              !
!   pis1     Ammonium inhibition parameter for small phytoplankton     !
!              [mmol_N/m3]                                             !
!   pis2     Ammonium inhibition parameter for diatom [mmol_N/m3]      !
!   akno3s1  Half saturation concentration for nitrate uptake by       !
!              small phytoplankton [mmol_N/m3].                        !
!   akno3s2  Half saturation concentration for nitrate uptake by       !
!              diatom [mmol_N/m3].                                     !
!   aknh4s1  Half saturation concentration for ammonium uptake by      !
!              small phytoplankton [mmol_N/m3].                        !
!   aknh4s2  Half saturation concentration for ammonium uptake by      !
!              diatom [mmol_N/m3].                                     !
!   akpo4s1  Half saturation concentration for phosphate uptake by     !
!              small phytoplankton [mmol_P/m3].                        !
!   akpo4s2  Half saturation concentration for phosphate uptake by     !
!              diatom [mmol_P/m3].                                     !
!   akco2s1  Half saturation concentration for co2 uptake by           !
!              small phytoplankton [mmol_C/m3].                        !
!   akco2s2  Half saturation concentration for co2 uptake by           !
!              diatom [mmol_C/m3].                                     !
!   aksio4s2 Half saturation constant for silicate uptake by           !
!              diatom [mmol_Si/m3].                                     !
!   ak1      Light attenuation coefficient of water [1/m]              !
!   ak2      Specific light attenuation coefficient for                !
!              phytoplankton [1/m/(mmol_N/m3)].                        !
!   bgamma0   Mesozooplankton specific mortality rate [1/day].         !
!   bgamma1   Grazing efficiency of microzooplankton [nondimensional]. !
!   bgamma2   Grazing efficiency of mesozooplankton [nondimensional].  !
!   bgamma3   Death rate of small phytoplankton [1/day].               !
!   bgamma4   Death rate of large phytoplankton [1/day].               !
!   bgamma5   Decay rate of detritus [1/day].                          !
!   bgamma6                                                            !
!   bgamma7   Nitrafication rate [1/day].                              !
!   wsd      Sinking velocity of detritus [m/day].                     !
!   wsdsi    Sinking velocity of detritus silicate [m/day].            !
!   wsp      Sinking velocity of large phytoplankton [m/day].          !
!   pco2a    Air pCO2 [ppmv].                                          !
!   si2n     Silicate to nitrogen ratio [mol_Si/mol_N].                !
!   p2n      Phosphorus to nitrogen ratio [mol_P/mol_N].               !
!   o2no     Oxygen to nitrate ratio [mol_O2/mol_NO3].                 !
!   o2nh     Oxygen to ammonium ratio [mol_O2/mol_NH4].                !
!   c2n      Carbon to nitrogen ratio [mol_C/mol_N].                   !
!   ro5      Grazing preference for diatom [nondimensional].           !
!   ro6      Grazing preference for mesozooplankton [nondimensional]   !
!   ro7      Grazing preference for detritus [nondimensional].         !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: reg1          ! 1/day
        real(r8), dimension(Ngrids) :: reg2          ! 1/day
        real(r8), dimension(Ngrids) :: gmaxs1        ! 1/day
        real(r8), dimension(Ngrids) :: gmaxs2        ! 1/day
        real(r8), dimension(Ngrids) :: beta1         ! 1/day
        real(r8), dimension(Ngrids) :: beta2         ! 1/day
        real(r8), dimension(Ngrids) :: akz1          ! mmol_N/m3
        real(r8), dimension(Ngrids) :: akz2          ! mmol_N/m3
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: amaxs2        ! 1/(Watts/m2)/day
        real(r8), dimension(Ngrids) :: pis1          ! m3/mmol_N
        real(r8), dimension(Ngrids) :: pis2          ! m3/mmol_N
        real(r8), dimension(Ngrids) :: akno3s1       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: aknh4s1       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: akpo4s1       ! mmol_P/m3
        real(r8), dimension(Ngrids) :: akco2s1       ! mmol_C/m3
        real(r8), dimension(Ngrids) :: akno3s2       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: aknh4s2       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: aksio4s2      ! mmol_Si/m3
        real(r8), dimension(Ngrids) :: akpo4s2       ! mmol_P/m3
        real(r8), dimension(Ngrids) :: akco2s2       ! mmol_C/m3
        real(r8), dimension(Ngrids) :: ak1           ! 1/m
        real(r8), dimension(Ngrids) :: ak2           ! 1/m/(mmol_N/m3)
        real(r8), dimension(Ngrids) :: parsats1      ! Watts/m2
        real(r8), dimension(Ngrids) :: parsats2      ! Watts/m2
        real(r8), dimension(Ngrids) :: bgamma0       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma1       ! [nondimensional]
        real(r8), dimension(Ngrids) :: bgamma2       ! [nondimensional]
        real(r8), dimension(Ngrids) :: bgamma3       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma4       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma5       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma6       ! 
        real(r8), dimension(Ngrids) :: bgamma7       ! 1/day
        real(r8), dimension(Ngrids) :: wsd           ! m/day
        real(r8), dimension(Ngrids) :: wsdsi         ! m/day
        real(r8), dimension(Ngrids) :: wsp           ! m/day
        real(r8), dimension(Ngrids) :: si2n          ! mol_Si/mol_N
        real(r8), dimension(Ngrids) :: pco2a         ! ppmv
        real(r8), dimension(Ngrids) :: p2n           ! mol_P/mol_N
        real(r8), dimension(Ngrids) :: o2no          ! mol_O2/mol_NO3
        real(r8), dimension(Ngrids) :: o2nh          ! mol_O2/mol_NH4
        real(r8), dimension(Ngrids) :: c2n           ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: ro5           ! nondimensional
        real(r8), dimension(Ngrids) :: ro6           ! nondimensional
        real(r8), dimension(Ngrids) :: ro7           ! nondimensional

# endif
# ifdef BIO_LIMADONEY
!                                                                      !
!  Parameters for Lima and Doney (2004) model:                         !
!                                                                      !
!   AttSW     Light attenuation due to seawater [1/m].                 !
!   AttChl    Light attenuation by chlorophyll [1/(mg_Chl m2)].        !
!   PARfrac   Photosynthetically active fraction of shortwave          !
!               radiation.                                             !
!   alphaPI   Initial slope of P_I curve (in carbon).                  !
!   thetaN0   Maximum Chlorophyll to phytoplankton Nitrogen ratio.     !
!   Tref      Reference temperature.                                   !
!   e_a       Temperature coefficient.                                 !
!   p_cref    Maximum phytoplankton C-spec. growth rate at Tref.       !
!   z_umax1   Maximum grazing rate on small phytoplankton at Tref.     !
!   z_umax2   Maximum grazing rate on large phytoplankton at Tref.     !
!   p1_excr   Small phytoplankton excretion rate.                      !
!   p2_excr   Large phytoplankton excretion rate.                      !
!   p1_mort2  Small phytoplankton quadratic mortality rate.            !
!   p2_mort2  Large phytoplankton quadratic mortality rate.            !
!   z_mort    Zooplankton mortality rate.                              !
!   z_mort2   zooplankton quadratic mortality rate.                    !
!   d1_rem    Small detrital remineralization rate.                    !
!   d2_rem    Large detrital remineralization rate.                    !
!   wDet2     Large detritial sinking rate.                            !
!   p1_kno3   Small phytoplankton NO3 uptake coefficient.              !
!   p2_kno3   Large phytoplankton NO3 uptake coefficient.              !
!   p1_knh4   Small phytoplankton NH4 uptake coefficient.              !
!   p2_knh4   Large phytoplankton NH4 uptake coefficient.              !
!   p2_si     Large phytoplankton Si uptake coefficient.               !
!   p_nitr    Nitrification rate.                                      !
!   falloc0   Minimum recycling allocation factor.                     !
!   si_alloc  Silica recycling allocation factor.                      !
!   z_grz     Grazing coefficient.                                     !
!   zegest0   Zooplankton egestion coefficient.                        !
!   r_ref     Respiration rate at temperature Tref.                    !
!   q_min     Minimum phytoplankton N:C ratio (mol/mol).               !
!   q_max     Maximum phytoplankton N:C ratio (mol/mol).               !
!   q_simin   Minimum large phytoplankton Si:C ratio (mol/mol).        !
!   q_simax   Maximum large phytoplankton Si:C ratio (mol/mol).        !
!   vncref    Value of VNC at Tref.                                    !
!   biolambda cost of biosynthesis (mmolC/mmolN) GD98.                 !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: AttChl        ! 1(mg_Chl m2)
        real(r8), dimension(Ngrids) :: PARfrac
        real(r8), dimension(Ngrids) :: alphaPI
        real(r8), dimension(Ngrids) :: thetaN0
        real(r8), dimension(Ngrids) :: tref
        real(r8), dimension(Ngrids) :: e_a
        real(r8), dimension(Ngrids) :: p_cref
        real(r8), dimension(Ngrids) :: z_umax1
        real(r8), dimension(Ngrids) :: z_umax2
        real(r8), dimension(Ngrids) :: p1_excr
        real(r8), dimension(Ngrids) :: p2_excr
        real(r8), dimension(Ngrids) :: p1_mort2
        real(r8), dimension(Ngrids) :: p2_mort2
        real(r8), dimension(Ngrids) :: z_mort
        real(r8), dimension(Ngrids) :: z_mort2
        real(r8), dimension(Ngrids) :: d1_rem
        real(r8), dimension(Ngrids) :: d2_rem
        real(r8), dimension(Ngrids) :: wDet2         ! m/day    
        real(r8), dimension(Ngrids) :: p1_kno3
        real(r8), dimension(Ngrids) :: p2_kno3
        real(r8), dimension(Ngrids) :: p1_knh4
        real(r8), dimension(Ngrids) :: p2_knh4
        real(r8), dimension(Ngrids) :: p2_si
        real(r8), dimension(Ngrids) :: p_nitr
        real(r8), dimension(Ngrids) :: falloc0
        real(r8), dimension(Ngrids) :: si_alloc
        real(r8), dimension(Ngrids) :: z_grz
        real(r8), dimension(Ngrids) :: zegest0
        real(r8), dimension(Ngrids) :: r_ref
        real(r8), dimension(Ngrids) :: q_min
        real(r8), dimension(Ngrids) :: q_max
        real(r8), dimension(Ngrids) :: q_simin
        real(r8), dimension(Ngrids) :: q_simax
        real(r8), dimension(Ngrids) :: vncref
        real(r8), dimension(Ngrids) :: biolambda

# endif
# ifdef NEMURO
!
!  Parameters for Nemuro ecosystem model:                              ! 
!                                                                      !
!  AlphaPS     Small Phytoplankton photochemical reaction coefficient: !
!                initial slope (low light) of the P-I curve,           !
!                [1/(W/m2) 1/day].                                     !
!  AlphaPL     Large Phytoplankton photochemical reaction coefficient: !
!                initial slope (low light) of the P-I curve,           !
!                [1/(W/m2) 1/day].                                     !
!  AlphaZL     Large Zooplankton assimilation efficiency,              !
!                [nondimemsional].                                     !
!  AlphaZP     Predator Zooplankton assimilation efficiency,           !
!                [nondimemsional].                                     !
!  AlphaZS     Small Zooplankton assimilation efficiency,              !
!                [nondimemsional].                                     !
!  AttPL       Light attenuation due to Large Phytoplankton, self-     !
!                shading coefficient, [m2/millimole_N].                !
!  AttPS       Light attenuation due to Small Phytoplankton, self-     !
!                shading coefficient, [m2/millimole_N].                !
!  AttSW       Light attenuation due to sea water, [1/m].              !
!  BetaPL      Large Phytoplankton photoinhibition coefficient,        !
!                [1/(W/m2) 1/day].                                     !
!  BetaPS      Small Phytoplankton photoinhibition coefficient,        !
!                [1/(W/m2) 1/day].                                     !
!  BetaZL      Large Zooplankton growth efficiency [nondimensional].   !
!  BetaZP      Predator Zooplankton growth efficiency [nondimensional].!
!  BetaZS      Small Zooplankton growth efficiency [nondimensional].   !
!  BioIter     Maximum number of iterations to achieve convergence of  !
!                the nonlinear solution.                               !
!  GammaL      Large Phytoplankton ratio of extracellular excretion to !
!                photosynthesis [nondimensional].                      !
!  GammaS      Small Phytoplankton ratio of extracellular excretion to !
!                photosynthesis [nondimensional].                      !
!  GRmaxLpl    Large Zooplankton maximum grazing rate on Large         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxLps    Large Zooplankton maximum grazing rate on Small         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxLzs    Small Zooplankton maximum grazing rate on Small         !
!                Zooplankton at 0 Celsius, [1/day].                    !
!  GRmaxPpl    Predator Zooplankton maximum grazing rate on Large      !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxPzl    Predator Zooplankton maximum grazing rate on Large      !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxPzs    Predator Zooplankton maximum grazing rate on Small      !
!                Zooplankton at 0 Celsius, [1/day].                    !
!  GRmaxSps    Small Zooplankton maximum grazing rate on Small         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxSpl    Small Zooplankton maximum grazing rate on Large         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  KD2N        Temperature coefficient for DON to NH4 decomposition,   !
!                [1/Celsius].                                          !
!  KGppL       Large Phytoplankton temperature coefficient for         !
!                photosynthetic rate, [1/Celsius].                     !
!  KGppS       Small Phytoplankton temperature coefficient for         !
!                photosynthetic rate, [1/Celsius].                     !
!  KGraL       Large Zooplankton temperature coefficient for grazing,  !
!                [1/Celsius].                                          !
!  KGraP       Predator Zooplankton temperature coefficient for        !
!                grazing,[1/Celsius].                                  !
!  KGraS       Small Zooplankton temperature coefficient for grazing,  !
!                [1/Celsius].                                          !
!  KMorPL      Large Phytoplankton temperature coefficient for         !
!                mortality, [1/Celsius].                               !
!  KMorPS      Small Phytoplankton temperature coefficient for         !
!                mortality, [1/Celsius].                               !
!  KMorZL      Large Zooplankton temperature coefficient for           !
!                mortality, [1/Celsius].                               !
!  KMorZP      Predator Zooplankton temperature coefficient for        !
!                mortality, [1/Celsius].                               !
!  KMorZS      Small Zooplankton temperature coefficient for           !
!                mortality, [1/Celsius].                               !
!  KNit        Temperature coefficient for nitrification (NH4 to NO3)  !
!                decomposition, [1/Celsius].                           !
!  KNH4L       Large Phytoplankton half satuation constant for NH4,    !
!                [millimole_N/m3].                                     !
!  KNH4S       Small Phytoplankton half satuation constant for NH4,    !
!                [millimole_N/m3].                                     !
!  KNO3L       Large Phytoplankton half satuation constant for NO3,    !
!                [millimole_N/m3].                                     !
!  KNO3S       Small Phytoplankton half satuation constant for NO3,    !
!                [millimole_N/m3].                                     !
!  KO2S        Temperature coefficient for Opal to SiOH4 decomposition,!
!                [1/Celsius].                                          !
!  KP2D        Temperature coefficient for PON to DON decomposition,   !
!                [1/Celsius].                                          !
!  KP2N        Temperature coefficient for PON to NH4 decomposition,   !
!                [1/Celsius].                                          !
!  KPL2ZS      Small Zooplankton half-saturation coefficient for       !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPL2ZL      Large Zooplankton half-saturation coefficient for       !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPL2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPS2ZL      Larg Zooplankton half-saturation coefficient for        !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KPS2ZS      Small Zooplankton half-saturation coefficient for       !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KResPL      Large Phytoplankton temperature coefficient for         !
!                respiration, [1/Celsius].                             !
!  KResPS      Small Phytoplankton temperature coefficient for         !
!                respiration, [1/Celsius].                             !
!  KSiL        Large Phytoplankton half satuation constant for SiOH4,  !
!                [millimole_Si/m3].                                    !
!  KZL2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Large Zooplankton [millimole_N/m3]^2.    !
!  KZS2ZL      Large Zooplankton half-saturation coefficient for       !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KZS2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Small Zooplankton [millimole_N/m3]^2.    !
!  LamL        Large Zooplankton Ivlev constant, [m3/millimole_N].     !
!  LamP        Predator Zooplankton Ivlev constant, [m3/millimole_N].  !
!  LamS        Small Zooplankton Ivlev constant, [m3/millimole_N].     !
!  MorPL0      Large Phytoplankton mortality rate at 0 Celsius,        !
!                [m3/millimole_N 1/day].                               !
!  MorPS0      Small Phytoplankton mortality rate at 0 Celsius,        !
!                [m3/millimole_N 1/day].                               !
!  MorZL0      Large Zooplankton mortality rate at 0 Celsius,          !
!                [m3/millimole_N 1/day].                               !
!  MorZP0      Predator Zooplankton mortality rate at 0 Celsius,       !
!                [m3/millimole_N 1/day].                               !
!  MorZS0      Small Zooplankton mortality rate at 0 Celsius,          !
!                [m3/millimole_N 1/day].                               !
!  Nit0        Nitrification (NH4 to NO3) rate at 0 Celsius, [1/day].  !
!  PARfrac     Fraction of shortwave radiation that is available for   !
!                photosyntesis [nondimensional].                       !
!  PL2ZSstar   Small Zooplankton threshold value for grazing on        !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PL2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PL2ZPstar   Predator Zooplankton threshold value for grazing on     !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PS2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  PS2ZSstar   Small Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  PusaiL      Large Phytoplankton Ammonium inhibition coefficient,    !
!                [m3/millimole_N].                                     !
!  PusaiPL     Predator Zooplankton grazing on Large Phytoplankton     !
!                inhibition coefficient, [m3/millimole_N].             !
!  PusaiS      Small Phytoplankton Ammonium inhibition coefficient,    !
!                [m3/millimole_N].                                     !
!  PusaiZS     Predator Zooplankton grazing on Small Zooplankton       !
!                inhibition coefficient, [m3/millimole_N].             !
!  ResPL0      Large Phytoplankton respiration rate at 0 Celsius,      !
!                [1/day].                                              !
!  ResPS0      Small Phytoplankton respiration rate at 0 Celsius,      !
!                [1/day].                                              !
!  RSiN        Si:N ratio [millimole_Si/millimole_N].                  !
!  setVOpal    Opal Settling (sinking) velocity [m/day].               !
!  setVPON     PON Settling (sinking) velocity [m/day].                !
!  VD2N0       DON to NH4 decomposition rate at 0 Celsius, [1/day].    !
!  VmaxL       Maximum Large Phytoplankton photosynthetic rate [1/day] !
!                in the absence of photoinhibition under optimal light.!
!  VmaxS       Maximum Small Phytoplankton photosynthetic rate [1/day] !
!                in the absence of photoinhibition under optimal light.!
!  VO2S0       Opal to Silicate decomposition rate at 0 Celsius,       !
!                [1/day].                                              !
!  VP2D0       PON to DON decomposition rate at 0 Celsius, [1/day].    !
!  VP2N0       PON to NH4 decomposition rate at 0 Celsius, [1/day].    !
!  ZL2ZPstar   Small Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  ZS2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Small Zooplankton, [millimole_N/m3].                  !
!  ZS2ZPstar   Predator Zooplankton threshold value for grazing on     !
!                Small Zooplankton, [millimole_N/m3].                  !
!                                                                      !
! Parameters for iron limitation                                       !
!  T_Fe      Iron uptake timescale, [day].                             !
!  A_Fe      Empirical Fe:C power, [-].                                !
!  B_Fe      Empirical Fe:C coefficient, [1/M-C].                      !
!  SK_FeC    Small phytoplankton Fe:C at F=0.5, [muM-Fe/M-C].          !
!  LK_FeC    Large phytoplankton Fe:C at F=0.5, [muM-Fe/M-C].          !
!  FeRR      Fe remineralization rate, [1/day].                        !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: AlphaPL       ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: AlphaPS       ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: AlphaZL       ! nondimensional
        real(r8), dimension(Ngrids) :: AlphaZP       ! nondimensional
        real(r8), dimension(Ngrids) :: AlphaZS       ! nondimensional
        real(r8), dimension(Ngrids) :: AttPL         ! m2/mmole_N
        real(r8), dimension(Ngrids) :: AttPS         ! m2/mmole_N
        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: BetaPL        ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: BetaPS        ! 1/(W/m2) 1/day
        real(r8), dimension(Ngrids) :: BetaZS        ! nondimensional
        real(r8), dimension(Ngrids) :: BetaZL        ! nondimensional
        real(r8), dimension(Ngrids) :: BetaZP        ! nondimensional
        real(r8), dimension(Ngrids) :: GammaL        ! nondimensional
        real(r8), dimension(Ngrids) :: GammaS        ! nondimensional
        real(r8), dimension(Ngrids) :: GRmaxLpl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxLps      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxLzs      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxPpl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxPzl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxPzs      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxSpl      ! 1/day
        real(r8), dimension(Ngrids) :: GRmaxSps      ! 1/day
        real(r8), dimension(Ngrids) :: KD2N          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGppL         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGppS         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGraL         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGraP         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KGraS         ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorPL        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorPS        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorZL        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorZP        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KMorZS        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KNH4L         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KNH4S         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KNit          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KNO3L         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KNO3S         ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KO2S          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KP2D          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KP2N          ! 1/Celsius
        real(r8), dimension(Ngrids) :: KPL2ZL        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPL2ZS        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPS2ZL        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPS2ZS        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KPL2ZP        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KResPL        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KResPS        ! 1/Celsius
        real(r8), dimension(Ngrids) :: KSiL          ! mmole_Si/m3
        real(r8), dimension(Ngrids) :: KZL2ZP        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KZS2ZL        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: KZS2ZP        ! mmole_N/m3
        real(r8), dimension(Ngrids) :: LamL          ! m3/mmole_N
        real(r8), dimension(Ngrids) :: LamP          ! m3/mmole_N
        real(r8), dimension(Ngrids) :: LamS          ! m3/mmole_N
        real(r8), dimension(Ngrids) :: MorPL0        ! m3/mmole_N/day
        real(r8), dimension(Ngrids) :: MorPS0        ! m3/mmole_N/day
        real(r8), dimension(Ngrids) :: MorZL0        ! m3/mmole_N 1/day
        real(r8), dimension(Ngrids) :: MorZP0        ! m3/mmole_N 1/day
        real(r8), dimension(Ngrids) :: MorZS0        ! m3/mmole_N 1/day
        real(r8), dimension(Ngrids) :: Nit0          ! 1/day
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: PusaiL        ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PusaiPL       ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PusaiS        ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PusaiZS       ! m3/mmole_N
        real(r8), dimension(Ngrids) :: PL2ZLstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PL2ZPstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PL2ZSstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PS2ZLstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: PS2ZSstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: ResPL0        ! 1/day
        real(r8), dimension(Ngrids) :: ResPS0        ! 1/day
        real(r8), dimension(Ngrids) :: RSiN          ! mmole_Si/mmole_N
        real(r8), dimension(Ngrids) :: setVOpal      ! m/day
        real(r8), dimension(Ngrids) :: setVPON       ! m/day
        real(r8), dimension(Ngrids) :: VD2N0         ! 1/day
        real(r8), dimension(Ngrids) :: VmaxL         ! 1/day
        real(r8), dimension(Ngrids) :: VmaxS         ! 1/day
        real(r8), dimension(Ngrids) :: VO2S0         ! 1/day
        real(r8), dimension(Ngrids) :: VP2D0         ! 1/day
        real(r8), dimension(Ngrids) :: VP2N0         ! 1/day
        real(r8), dimension(Ngrids) :: ZL2ZPstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: ZS2ZLstar     ! mmole_N/m3
        real(r8), dimension(Ngrids) :: ZS2ZPstar     ! mmole_N/m3
#  ifdef IRON_LIMIT
        real(r8), dimension(Ngrids) :: T_Fe          ! day
        real(r8), dimension(Ngrids) :: A_Fe          ! nondimensional
        real(r8), dimension(Ngrids) :: B_Fe          ! 1/M-C
        real(r8), dimension(Ngrids) :: SK_FeC        ! muM-Fe/M-C
        real(r8), dimension(Ngrids) :: LK_FeC        ! muM-Fe/M-C
        real(r8), dimension(Ngrids) :: FeRR          ! 1/day
#  endif
# endif

# ifdef NPZD_FRANKS
!
!  Parameters for Franks et al. (1986) type model:                     !
!                                                                      !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  DetRR     Detritus remineraliztion rate, [1/day].                   !
!  K_ext     Light extinction coefficient, [1/m].                      !
!  K_NO3     Inverse half-saturation for phytoplankton nitrate uptake  !
!              [1/(millimole_N m-3)].                                  !
!  K_Phy     Phytoplankton saturation coefficient, [millimole_N m-3].  !
!  PhyMR     Phytoplankton senescence/mortality rate, [1/day].         !
!  Vm_NO3    Nitrate uptake rate, [1/day].                             !
!  wDet      Detrital sinking rate, [m/day].                           !
!  ZooGR     Zooplankton maximum growth rate, [1/day].                 !
!  ZooMR     Zooplankton mortality rate, [1/day].                      !
!  ZooMD     Zooplankton death bits rate, [1/day].                     !
!  ZooGA     Zooplankton grazing inefficiency, [nondimensional].       !
!  ZooEC     Zooplankton excreted fraction, [nondimensional].          !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

#  ifdef ANA_BIOLOGY
        real(r8), allocatable :: BioIni(:,:)
#  endif
        real(r8), dimension(Ngrids) :: DetRR         ! 1/day
        real(r8), dimension(Ngrids) :: K_ext         ! 1/m
        real(r8), dimension(Ngrids) :: K_NO3         ! 1/(mmol/m3)
        real(r8), dimension(Ngrids) :: K_Phy         ! mmol/m3
        real(r8), dimension(Ngrids) :: PhyMR         ! 1/day
        real(r8), dimension(Ngrids) :: Vm_NO3        ! 1/day
        real(r8), dimension(Ngrids) :: wDet          ! m/day
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wDet
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wDet
#  endif
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMD         ! 1/day
        real(r8), dimension(Ngrids) :: ZooGA         ! nondimensional
        real(r8), dimension(Ngrids) :: ZooEC         ! nondimensional
# endif
# ifdef NPZD_POWELL
!
!  Parameters for Franks et al. (1986) type model:                     !
!                                                                      !
!  AttPhy    Light attenuation due to phytoplankton (self-shading      !
!              coefficient), [m2/millimole_N].                         !
!  AttSW     Light attenuation due to sea water, [1/m].                !
!  BioIter   Maximum number of iterations to achieve convergence of    !
!              the nonlinear solution.                                 !
!  BioIni    Initial concentration for analytical initial (uniform)    !
!              conditions.                                             !
!  DetRR     Detritus remineraliztion rate, [1/day].                   !
!  K_NO3     Half-saturation for phytoplankton nitrate uptake          !
!              [millimole_N m-3].                                      !
!  Ivlev     Ivlev constant for zooplankton grazing parameterization,  !
!              [nondimensional].                                       !
!  PARfrac   Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!  PhyIS     Phytoplankton, initial slope of the P-I curve [m2/W].     !
!  PhyMRD    Phytoplankton mortality rate to the Detritus pool,        !
!              [1/day].                                                !
!  PhyMRN    Phytoplankton mortality rate to the Nitrogen pool,        !
!              [1/day].                                                !
!  Vm_NO3    Nitrate uptake rate, [1/day].                             !
!  wDet      Detrital sinking rate, [m/day].                           !
!  wPhy      Phytoplankton sinking rate, [m/day].                      !
!  ZooEED    Zooplankton excretion efficiency to Detritus pool,        !
!              {nondimensional].                                       !
!  ZooEEN    Zooplankton excretion efficiency to Nitrogen pool,        !
!              {nondimensional].                                       !
!  ZooGR     Zooplankton grazing rate, [1/day].                        !
!  ZooMRD    Zooplankton mortality rate to Detritus pool, [1/day].     !
!  ZooMRN    Zooplankton mortality rate to Nitrogen pool, [1/day].     !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

#  ifdef ANA_BIOLOGY
        real(r8), allocatable :: BioIni(:,:)
#  endif
        real(r8), dimension(Ngrids) :: AttPhy        ! m2/mmole
        real(r8), dimension(Ngrids) :: AttSW         ! 1/m
        real(r8), dimension(Ngrids) :: DetRR         ! 1/day
        real(r8), dimension(Ngrids) :: K_NO3         ! mmol/m3
        real(r8), dimension(Ngrids) :: Ivlev         ! nondimensional
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_PARfrac
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_PARfrac
#  endif
        real(r8), dimension(Ngrids) :: PhyIS         ! m2/W
        real(r8), dimension(Ngrids) :: PhyMRD        ! 1/day
        real(r8), dimension(Ngrids) :: PhyMRN        ! 1/day
        real(r8), dimension(Ngrids) :: Vm_NO3        ! 1/day
        real(r8), dimension(Ngrids) :: wDet          ! m/day
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wDet
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wDet
#  endif
        real(r8), dimension(Ngrids) :: wPhy          ! m/day
#  ifdef TANGENT
        real(r8), dimension(Ngrids) :: tl_wPhy
#  endif
#  ifdef ADJOINT
        real(r8), dimension(Ngrids) :: ad_wPhy
#  endif
        real(r8), dimension(Ngrids) :: ZooEED        ! nondimensional
        real(r8), dimension(Ngrids) :: ZooEEN        ! nondimensional
        real(r8), dimension(Ngrids) :: ZooGR         ! 1/day
        real(r8), dimension(Ngrids) :: ZooMRD        ! 1/day
        real(r8), dimension(Ngrids) :: ZooMRN        ! 1/day
# endif
# ifdef ECOSIM
!
!  EcoSim Model Phytoplaknton Parameters:                              !
!                                                                      !
!                                                                      !
!  HsNO3          Half-saturation for phytoplankton NO3 uptake         !
!                   (micromole_NO3/liter).                             !
!  HsNH4          Half-saturation for phytoplankton NH4 uptake         !
!                   (micromole_NH4/liter).                             !
!  HsSiO          Half-saturation for phytoplankton SiO uptake         !
!                   (micromole_SiO/liter).                             !
!  HsPO4          Half-saturation for phytoplankton PO4 uptake         !
!                   (micromole_PO4/liter).                             !
!  HsFe           Half-saturation for phytoplankton Fe uptake          !
!                  (micromole_Fe/liter).                               !
!  GtALG_max      Maximum 24 hour growth rate (1/day).                 !
!  PhyTbase       Phytoplankton temperature base for exponential       !
!                   response to temperature (Celsius).                 !
!  PhyTfac        Phytoplankton exponential temperature factor         !
!                   (1/Celsius).                                       !
!  BET_           Nitrate uptake inhibition for NH4 (l/micromole).     !
!  maxC2nALG      Maximum phytoplankton C:N ratio                      !
!                   (micromole_C/micromole_N).                         !
!  minC2nALG      Balanced phytoplankton C:N ratio                     !
!                   (micromole_C/micromole_N).                         !
!  C2nALGminABS   Absolute minimum phytoplankton C:N ratio             !
!                   (micromole_C/micromole_N).                         !
!  maxC2SiALG     Maximum phytoplankton C:Si ratio                     !
!                   (micromole_C/micromole_Si).                        !
!  minC2SiALG     Balanced phytoplankton C:Si ratio                    !
!                   (micromole_C/micromole_Si).                        !
!  C2SiALGminABS  Absolute minimum phytoplankton C:Si ratio            !
!                  (micromole_C/micromole_Si).                         !
!  maxC2pALG      Maximum phytoplankton C:P ratio                      !
!                   (micromole_C/micromole_P).                         !
!  minC2pALG      Balanced phytoplankton C:P ratio                     !
!                   (micromole_C/micromole_P).                         !
!  C2pALGminABS   Absolute minimum phytoplankton C:P ratio             !
!                   (micromole_C/micromole_P).                         !
!  maxC2FeALG     Maximum phytoplankton C:Fe ratio                     !
!                   (micromole_C/micromole_Fe).                        !
!  minC2FeALG     Balanced phytoplankton C:Fe ratio                    !
!                   (micromole_C/micromole_Fe).                        !
!  C2FeALGminABS  Absolute minimum phytoplankton C:Fe ratio            !
!                   (micromole_C/micromole_Fe).                        !
!  qu_yld         Maximum quantum yield                                !
!                   (micromole_C/micromole_quanta).                    !
!  E0_comp        Compensation light level (micromole_quanta).         !
!  E0_inhib       Light level for onset of photoinhibition             !
!                   (micromole_quanta).                                !
!  inhib_fac      Exponential decay factor for light limited growth    !
!                   (1/micromole_quanta).                              !
!  C2Chl_max      Maximum lighted limited (nutrient replete) C:Chl     !
!                   ratio (microgram_C/microgram_Chl).                 !
!  mxC2Cl         Rate of change in the lighted limited C:Chl ratio    !
!                   (microgram_C/microgram_Chl/micromole_quanta).      !
!  b_C2Cl         Minimum lighted limited (nutrient replete) C:Chl     !
!                   ratio (microgram_C/microgram_Chl).                 !
!  mxC2Cn         Rate of change in the nutrient limited C:Chl ratio   !
!                   [(ug_C/ug_Chl)/(umole_C/umole_N)].                 !
!  b_C2Cn         Minimum nutrient limited C:Chl ratio                 !
!                   (microgram_C/microgram_Chl).                       !
!  mxPacEff       Rate of change in package effect                     !
!                   [1/(microgram_C/microgram_Chl)].                   !
!  b_PacEff       Maximum package effect                               !
!                   [1/(microgram_C/microgram_Chl)].                   !
!  mxChlB         Rate of change in the Chl_b:Chl_a ratio              !
!                   [(ug_Chl_b/ug_Chl_a)/(ug_C/ug_Chl_a)].             !
!  b_ChlB         Maximum Chl_b:Chl_a ratio                            !
!                   (microgram_Chl_b/microgram_Chl_a).                 !
!  mxChlC         Rate of change in the Chl_c:Chl_a ratio              !             
!                   [(ug_Chl_c/ug_Chl_a)/(ug_C/ug_Chl_a)].             !
!  b_ChlC         Maximum Chl_c:Chl_a ratio                            !
!                   (microgram_Chl_c/microgram_Chl_a).                 !
!  mxPSC          Rate of change in the PSC:Chl_a ratio                !
!                   [(ug_PSC/ug_Chl_a)/(ug_C/ug_Chl_a)].               !
!  b_PSC          Maximum PSC:Chl_a ratio                              !
!                  (microgram_Chl_c/microgram_Chl_a).                  !
!  mxPPC          Rate of change in the PPC:Chl_a ratio                !
!                   [(ug_PPC/ug_Chl_a)/(ug_C/ug_Chl_a)].               !
!  b_PPC          Maximum PPC:Chl_a ratio                              !
!                  (microgram_Chl_c/microgram_Chl_a).                  !
!  mxLPUb         Rate of change in the LPUb:Chl_a ratio               !
!                   [(ug_LPUb/ug_Chl_a)/(ug_C/ug_Chl_a)].              !
!  b_LPUb         Maximum LPUb:Chl_a ratio                             !
!                   (microgram_HPUb/microgram_Chl_a).                  !
!  mxHPUb         Rate of change in the HPUb:Chl_a ratio               !
!                   [(ug_HPUb/ug_Chl_a)/(ug_C/ug_Chl_a)].              !
!  b_HPUb         Maximum HPUb:Chl_a ratio                             !
!                   (microgram_HPUb/microgram_Chl_a).                  !
!  FecDOC         Proportion of grazing stress which is apportioned    !
!                   to DOM (nondimensional).                           !
!  FecPEL         Proportion of grazing stress which is apportioned    !
!                   to fecal pellets (nondimesional).                  !
!  FecCYC         Proportion of grazing stress which is apportioned    !
!                   to direct remineralization (nondimensional).       !
!  ExALG          Proportion of daily production that is lost to       !
!                   excretion (nondimensional).                        !
!  WS             Phytoplankton sinking speed (meters/day).            !
!  HsGRZ          Phytoplankton grazing parameter (nondimensional).    !
!  MinRefuge      Refuge Phytoplankton population (micromole_C/liter). !
!  RefugeDep      Maximum Refuge Phytoplankton depth (meters).         !
!  Norm_Vol       Normalized Volume factor (nondimensional).           !
!  Norm_Surf      Normalized surface area factor (nondimensional).     !
!  HsDOP          Half Saturation Constant for DOP uptake              !
!                   (micromole_DOP/liter).                             !
!  C2pALKPHOS     C:P ratio where DOP uptake begins                    !
!                   (micromole_C/micromole_P).                         !
!  HsDON          Half Saturation Constant for DON uptake              !
!                   (micromole_DON/liter).                             !
!  C2nNupDON      C:N ratio where DON uptake begins                    !
!                   (micromole_C/micromole_N).                         !
!                                                                      !
! Bacteria Parameters:                                                 !
!                                                                      !
!  HsDOC_ba       Half saturation constant for bacteria DOC uptake     !
!                   (micromole_DOC/liter).                             !
!  GtBAC_max      Maximum 24 hour bacterial growth rate (1/day).       !
!  BacTbase       Phytoplankton temperature base for exponential       !
!                   response to temperature, (Celsius).                !
!  BacTfac        Phytoplankton exponential temperature factor         !
!                   (1/Celsius).                                       !
!  C2nBAC         Carbon to Nitrogen ratio of Bacteria                 !
!                   (micromole_C/micromole_N).                         !
!  C2pBAC         Carbon to Phosphorus ratio of Bacteria               !
!                   (micromole_C/micromole_P).                         !
!  C2FeBAC        Carbon to Iron ratio of Bacteria                     !
!                   (micromole_C/micromole_Fe)                         !
!  BacDOC         Proportion of bacterial grazing stress which is      !
!                   apportioned to DOM (nondimensional).               !
!  BacPEL         Proportion of bacterial grazing stress which is      !
!                   apportioned to fecal pellets (nondimensional).     !
!  BacCYC         Proportion of bacterial grazing stress which is      !
!                   apportioned to direct remineralization             !
!                   (nondimensional).                                  !
!  ExBAC_c        Bacterial recalcitrant carbon excretion as a         !
!                   proportion of uptake (nondimensional)              !
!  ExBacC2N       Bacterial recalcitrant excretion carbon to nitrogen  !
!                   ratio (micromole_C/micromole_N).                   !
!  Bac_Ceff       Bacterial gross growth carbon efficiency             !
!                   (nondimensional).                                  !
!  RtNIT          Maximum bacterial nitrification rate (1/day).        !
!  HsNIT          Half saturation constant for bacterial nitrification !
!                   (micromole NH4/liter)                              !
!                                                                      !
! Dissolved Organic Matter Parameters:                                 !
!                                                                      !
!  cDOCfrac_c     Colored fraction of DOC production from              !
!                   phytoplankton and bacterial losses                 !
!                   (nondimensional).                                  !
!  RtUVR_DIC      UV degradation of DOC into DIC at 410 nanometers     !
!                   (micromole/meter/liter/hour).                      !
!  RtUVR_DIC      UV degradation of DOC into colorless labile DOC at   !
!                   410 nanometers (micromole/meter/liter/hour).       !
!                                                                      !
! Fecal and detritus Parameters:                                       !
!                                                                      !
!  WF             Fecal sinking flux (meters/day).                     !
!  RegTbase       Fecal regeneration temperature base for exponential  !
!                   response to temperature (Celsius).                 !
!  RegTfac        Fecal regeneration exponential temperature factor    !
!                   (1/Celsius).                                       !
!  RegCR          Fecal carbon regeneration rate (1/day).              !
!  RegNR          Fecal nitrogen regeneration rate (1/day).            !
!  RegSR          Fecal silica regeneration rate (1/day).              !
!  RegPR          Fecal phosphorus regeneration rate (1/day).          !
!  RegFR          Fecal iron regeneration rate (1/day).                !
!                                                                      !
!======================================================================!
!
      USE mod_param
      USE mod_eclight
      USE mod_scalars

      implicit none
!
!-----------------------------------------------------------------------
!  Standard input parameters.
!-----------------------------------------------------------------------
!
!  Number of biological iterations.
!
      integer, dimension(Ngrids) :: BioIter
!
!  Control flags.
!
      logical, dimension(Ngrids) :: RtUVR_flag
      logical, dimension(Ngrids) :: NFIX_flag
      logical, dimension(Ngrids) :: Regen_flag
!
!  Phytoplankton parameters.
!
      real(r8), dimension(Nphy,Ngrids) :: HsNO3
      real(r8), dimension(Nphy,Ngrids) :: HsNH4
      real(r8), dimension(Nphy,Ngrids) :: HsSiO
      real(r8), dimension(Nphy,Ngrids) :: HsPO4
      real(r8), dimension(Nphy,Ngrids) :: HsFe
      real(r8), dimension(Nphy,Ngrids) :: GtALG_max
      real(r8), dimension(Nphy,Ngrids) :: PhyTbase
      real(r8), dimension(Nphy,Ngrids) :: PhyTfac
      real(r8), dimension(Nphy,Ngrids) :: BET_
      real(r8), dimension(Nphy,Ngrids) :: maxC2nALG
      real(r8), dimension(Nphy,Ngrids) :: minC2nALG
      real(r8), dimension(Nphy,Ngrids) :: C2nALGminABS
      real(r8), dimension(Nphy,Ngrids) :: maxC2SiALG
      real(r8), dimension(Nphy,Ngrids) :: minC2SiALG
      real(r8), dimension(Nphy,Ngrids) :: C2SiALGminABS
      real(r8), dimension(Nphy,Ngrids) :: maxC2pALG
      real(r8), dimension(Nphy,Ngrids) :: minC2pALG
      real(r8), dimension(Nphy,Ngrids) :: C2pALGminABS
      real(r8), dimension(Nphy,Ngrids) :: maxC2FeALG
      real(r8), dimension(Nphy,Ngrids) :: minC2FeALG
      real(r8), dimension(Nphy,Ngrids) :: C2FeALGminABS
      real(r8), dimension(Nphy,Ngrids) :: qu_yld
      real(r8), dimension(Nphy,Ngrids) :: E0_comp
      real(r8), dimension(Nphy,Ngrids) :: E0_inhib
      real(r8), dimension(Nphy,Ngrids) :: inhib_fac
      real(r8), dimension(Nphy,Ngrids) :: C2CHL_max
      real(r8), dimension(Nphy,Ngrids) :: mxC2Cl
      real(r8), dimension(Nphy,Ngrids) :: b_C2Cl
      real(r8), dimension(Nphy,Ngrids) :: mxC2Cn
      real(r8), dimension(Nphy,Ngrids) :: b_C2Cn
      real(r8), dimension(Nphy,Ngrids) :: mxPacEff
      real(r8), dimension(Nphy,Ngrids) :: b_PacEff
      real(r8), dimension(Nphy,Ngrids) :: mxChlB
      real(r8), dimension(Nphy,Ngrids) :: b_ChlB
      real(r8), dimension(Nphy,Ngrids) :: mxChlC
      real(r8), dimension(Nphy,Ngrids) :: b_ChlC
      real(r8), dimension(Nphy,Ngrids) :: mxPSC
      real(r8), dimension(Nphy,Ngrids) :: b_PSC
      real(r8), dimension(Nphy,Ngrids) :: mxPPC
      real(r8), dimension(Nphy,Ngrids) :: b_PPC
      real(r8), dimension(Nphy,Ngrids) :: mxLPUb
      real(r8), dimension(Nphy,Ngrids) :: b_LPUb
      real(r8), dimension(Nphy,Ngrids) :: mxHPUb
      real(r8), dimension(Nphy,Ngrids) :: b_HPUb
      real(r8), dimension(Nphy,Ngrids) :: FecDOC
      real(r8), dimension(Nphy,Nfec,Ngrids) :: FecPEL
      real(r8), dimension(Nphy,Ngrids) :: FecCYC
      real(r8), dimension(Nphy,Ngrids) :: ExALG
      real(r8), dimension(Nphy,Ngrids) :: WS
      real(r8), dimension(Nphy,Ngrids) :: HsGRZ
      real(r8), dimension(Nphy,Ngrids) :: MinRefuge
      real(r8), dimension(Nphy,Ngrids) :: RefugeDep
      real(r8), dimension(Nphy,Ngrids) :: Norm_Vol
      real(r8), dimension(Nphy,Ngrids) :: Norm_Surf
      real(r8), dimension(Nphy,Ngrids) :: HsDOP
      real(r8), dimension(Nphy,Ngrids) :: C2pALKPHOS
      real(r8), dimension(Nphy,Ngrids) :: HsDON
      real(r8), dimension(Nphy,Ngrids) :: C2nNupDON
!
!  Bacteria parameters.
!
      real(r8), dimension(Nbac,Ngrids) :: HsDOC_ba
      real(r8), dimension(Nbac,Ngrids) :: GtBAC_max
      real(r8), dimension(Nbac,Ngrids) :: BacTbase
      real(r8), dimension(Nbac,Ngrids) :: BacTfac
      real(r8), dimension(Ngrids) :: C2nBAC
      real(r8), dimension(Ngrids) :: C2pBAC
      real(r8), dimension(Ngrids) :: C2FeBAC
      real(r8), dimension(Ngrids) :: BacDOC
      real(r8), dimension(Ngrids) :: BacPEL
      real(r8), dimension(Ngrids) :: BacCYC
      real(r8), dimension(Ngrids) :: ExBAC_c
      real(r8), dimension(Ngrids) :: ExBacC2N
      real(r8), dimension(Ngrids) :: Bac_Ceff
      real(r8), dimension(Ngrids) :: RtNIT
      real(r8), dimension(Ngrids) :: HsNIT  
!
!  DOM parameters.
!
      real(r8), dimension(Ndom,Ngrids) :: cDOCfrac_c
      real(r8), dimension(Ngrids) :: RtUVR_DIC
      real(r8), dimension(Ngrids) :: RtUVR_DOC
!
!  Fecal parameters.
!
      real(r8), dimension(Nfec,Ngrids) :: WF
      real(r8), dimension(Nfec,Ngrids) :: RegTbase
      real(r8), dimension(Nfec,Ngrids) :: RegTfac
      real(r8), dimension(Nfec,Ngrids) :: RegCR
      real(r8), dimension(Nfec,Ngrids) :: RegNR
      real(r8), dimension(Nfec,Ngrids) :: RegSR
      real(r8), dimension(Nfec,Ngrids) :: RegPR
      real(r8), dimension(Nfec,Ngrids) :: RegFR
!
!-----------------------------------------------------------------------
!  Internal parameters.
!-----------------------------------------------------------------------
!    
!  Spectral band width used in light calculations.

      real(r8), parameter :: DLAM  = 5.0_r8
!
!  Flags used for testing purposes.
!
      real(r8), parameter :: SMALL  = 1.0e-6_r8
      real(r8), parameter :: VSMALL = 1.0e-14_r8
      real(r8), parameter :: LARGE  = 1.0e+10_r8
      real(r8), parameter :: VLARGE = 1.0e+50_r8
!
!  Array indexes for frequently used constituents.
!
      integer, parameter :: ilab=1    ! labile index for DOC.
      integer, parameter :: irct=2    ! relict index for DOC.
      integer, parameter :: ichl=1    ! pigment index for chlorophyll-a
      integer, parameter :: isfc=1    ! slow fecal pellet index
      integer, parameter :: iffc=2    ! fast fecal pellet index
!
!  Phytoplankton calculated paramters.
!
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2nALG   ! inverse C2nALG
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2SiALG  ! inverse C2SiALG
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2pALG   ! inverse C2pALG
      real(r8), dimension(Nphy,Ngrids) :: ImaxC2FeALG  ! inverse C2FeALG
!
!  Bacteria calculated parameters.
!
      real(r8), dimension(Ngrids) :: N2cBAC
      real(r8), dimension(Ngrids) :: P2cBAC
      real(r8), dimension(Ngrids) :: Fe2cBAC
      real(r8), dimension(Nbac,Ngrids) :: HsNH4_ba
      real(r8), dimension(Nbac,Ngrids) :: HsPO4_ba
      real(r8), dimension(Nbac,Ngrids) :: HsFe_ba
      real(r8), dimension(Ngrids) :: R_ExBAC_c
      real(r8), dimension(Ngrids) :: ExBAC_n
      real(r8), dimension(Ngrids) :: Frac_ExBAC_n
      real(r8), dimension(Ngrids) :: I_Bac_Ceff
!
!  Absorption parameters.
!
      real(r8), dimension(NBands) :: wavedp   ! a and b factor
      real(r8), dimension(Ndom) :: aDOC410    ! CDM absorption at 410
      real(r8), dimension(Ndom) :: aDOC300    ! CDM absorption at 300

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine initializes several parameters in module "mod_biology" !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: ibac, iband, ifec, iphy, ng
!
!-----------------------------------------------------------------------
!  Derived parameters.
!-----------------------------------------------------------------------
!
!  Convert rates from day-1 to second-1.
!
      DO ng=1,Ngrids
        DO iphy=1,Nphy
          GtALG_max(iphy,ng)=GtALG_max(iphy,ng)*sec2day
          ExALG(iphy,ng)=ExALG(iphy,ng)*sec2day
          HsGRZ(iphy,ng)=HsGRZ(iphy,ng)*sec2day
          WS(iphy,ng)=WS(iphy,ng)*sec2day
        END DO
        DO ibac=1,Nbac
          GtBAC_max(ibac,ng)=GtBAC_max(ibac,ng)*sec2day
        END DO
        DO ifec=1,Nfec
          WF(ifec,ng)=WF(ifec,ng)*sec2day
        END DO
        RtNIT(ng)=RtNIT(ng)*sec2day
      END DO
!
!  Calculated reciprocal phytoplankton parameters.
!
      DO ng=1,Ngrids
        DO iphy=1,Nphy
          IF (maxC2nALG(iphy,ng).gt.SMALL) THEN
            ImaxC2nALG(iphy,ng)=1.0_r8/maxC2nALG(iphy,ng)
          ELSE
            ImaxC2nALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2SiALG(iphy,ng).gt.SMALL) THEN
            ImaxC2SiALG(iphy,ng)=1.0_r8/maxC2SiALG(iphy,ng)
          ELSE
            ImaxC2SiALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2pALG(iphy,ng).gt.SMALL) THEN
            ImaxC2pALG(iphy,ng)=1.0_r8/maxC2pALG(iphy,ng)
          ELSE
            ImaxC2pALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2FeALG(iphy,ng).gt.SMALL) THEN
            ImaxC2FeALG(iphy,ng)=1.0_r8/maxC2FeALG(iphy,ng)
          ELSE
            ImaxC2FeALG(iphy,ng)=0.0_r8
          END IF
        END DO
      END DO
!
!  Calculated bacterial parameters.
!
      DO ng=1,Ngrids
        DO ibac=1,Nbac
          HsNH4_ba(ibac,ng)=HsDOC_ba(ibac,ng)/C2nBAC(ng)
          HsPO4_ba(ibac,ng)=HsDOC_ba(ibac,ng)/C2pBAC(ng)
          HsFe_ba (ibac,ng)=HsDOC_ba(ibac,ng)/C2FeBAC(ng)
        END DO
      END DO
!
!  Inverse parameters for computational efficiency.
!
      DO ng=1,Ngrids
        N2cBAC(ng)=1.0_r8/C2nBAC(ng)
        P2cBAC(ng)=1.0_r8/C2pBAC(ng)
        Fe2cBAC(ng)=1.0_r8/C2FeBAC(ng)
        I_Bac_Ceff(ng)=1.0_r8/Bac_Ceff(ng) 
      END DO
!
!  Reciprocal of non baterial recalcitran carbon excretion.
!
      DO ng=1,Ngrids
        R_ExBAC_c(ng)=1.0_r8/(1.0_r8-ExBAC_c(ng))
      END DO
!
!  Bacterial recalcitrant nitrogen excretion as a function of uptake.
!
      DO ng=1,Ngrids
        ExBAC_n(ng)=ExBAC_c(ng)*C2nBAC(ng)/ExBacC2N(ng)
        Frac_ExBAC_n(ng)=1.0_r8-ExBAC_n(ng)
      END DO
!
!  Scale UV degradation parameters.
!
      DO ng=1,Ngrids
        RtUVR_DIC(ng)=RtUVR_DIC(ng)/3600.0_r8
        RtUVR_DOC(ng)=RtUVR_DOC(ng)/3600.0_r8 
      END DO
!
!  If applicable, zero-out fecal regeneration parameters.
!
      DO ng=1,Ngrids
        IF (Regen_flag(ng)) THEN
          DO ifec=1,Nfec
            RegCR(ifec,ng)=RegCR(ifec,ng)*sec2day
            RegNR(ifec,ng)=RegNR(ifec,ng)*sec2day
            RegPR(ifec,ng)=RegPR(ifec,ng)*sec2day
            RegFR(ifec,ng)=RegFR(ifec,ng)*sec2day
            RegSR(ifec,ng)=RegSR(ifec,ng)*sec2day
          END DO
        ELSE
          DO ifec=1,Nfec
            RegCR(ifec,ng)=0.0_r8
            RegNR(ifec,ng)=0.0_r8
            RegPR(ifec,ng)=0.0_r8
            RegFR(ifec,ng)=0.0_r8
            RegSR(ifec,ng)=0.0_r8
          END DO
        END IF
      END DO
!
!  Spectral dependency for scattering and backscattering.
!
      DO iband=1,NBands
        wavedp(iband)=(550.0_r8/(397.0_r8+REAL(iband,r8)*DLAM))
      END DO
!
!  Calculated IOP parameter values.
!
      aDOC410(ilab)=aDOC(ilab,1)*EXP(0.014_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC410(irct)=aDOC(irct,1)*EXP(0.025_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC300(ilab)=EXP(0.0145_r8*(410.0_r8-300.0_r8))
      aDOC300(irct)=EXP(0.0145_r8*(410.0_r8-300.0_r8))

      RETURN
      END SUBROUTINE initialize_biology
# endif
#ifdef BEST_NPZ
!=======================================================================
!   Parameters for BEST NPZ Biological Model
!      Georgina Gibson 
!=======================================================================

      USE mod_param

        implicit none

        integer, dimension(Ngrids) :: BioIter
        real(r8) :: VertMixIncr
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
!  Bio- conversions
        real(r8) :: xi, ccr, ccrPhL
!  extinction coefficients
        real(r8) :: k_ext, k_chl
!  PhS growth params
        real(r8) :: DiS, DpS,aPS
        real(r8) :: alphaPhS
        real(r8) :: psiPhS
        real(r8) :: k1PhS, k2PhS
!  PhL growth params
        real(r8) :: DiL,aPL
        real(r8) :: DpL
        real(r8) :: alphaPhL
        real(r8) :: psiPhL
        real(r8) :: k1PhL
        real(r8) :: k2PhL
!  MZS preference params
        real(r8) :: fpPhSMZS, fpPhLMZS
!  MZS growth and feeding params
        real(r8) :: eMZS
        real(r8) :: Q10MZS
        real(r8) :: Q10MZST
        real(r8) :: fMZS
        real(r8) :: kMZS
        real(r8) :: gammaMZS
!  MZL preferences params
        real(r8) :: fpPhSMZL, fpPhLMZL, fpMZSMZL
!  MZL growth and feeding params
        real(r8) :: eMZL
        real(r8) :: Q10MZL
        real(r8) :: Q10MZLT
        real(r8) :: fMZL
        real(r8) :: kMZL
        real(r8) :: gammaMZL
!  Cop preference params
        real(r8) :: fpPhSCop, fpPhLCop, fpMZSCop, fpMZLCop
!  Cop growth and feeding params
        real(r8) :: eCop
        real(r8) :: Q10Cop
        real(r8) :: Q10CopT
        real(r8) :: fCop
        real(r8) :: gammaCop
        real(r8) :: kCop
!  NCa preference params
        real(r8) :: fpPhSNCa, fpPhLNCa, fpMZSNCa, fpMZLNCa
!  NCa growth and feeding params
        real(r8) :: eNCa
        real(r8) :: Q10NCa
        real(r8) :: Q10NCaT
        real(r8) :: fNCa
        real(r8) :: gammaNCa
        real(r8) :: kNCa
!  Eup preference params
        real(r8) :: fpPhSEup, fpPhLEup, fpMZSEup, fpMZLEup, fpCopEup
!  Eup growth and feeding params
        real(r8) :: eEup
        real(r8) :: Q10Eup
        real(r8) :: Q10EupT
        real(r8) :: fEup
        real(r8) :: gammaEup
        real(r8) :: kEup
# ifdef JELLY
!  Jellyfish Parameters
        real(r8) :: eJel, gammaJel,fJel
        real(r8) :: respJel,mpredJel
        real(r8) :: fpCopJel, fpNCaJel, fpEupJel
        real(r8) :: Q10Jelr, Q10JelTr,Q10Jele, Q10JelTe
        real(r8) :: bmJ,ktbmJ,TrefJ

# endif

!  Phytoplankton senescence
        real(r8) :: mPhS, maxmPhS, NcritPhS,minmPhS
        real(r8) :: mPhL, maxmPhL, NcritPhL, minmPhL
!  Zoopkankton mortality
        real(r8) :: mMZS, mMZL, mCop, mNCa, mEup
!  predation closure
        real(r8) :: mpredCop, mpredNCa, mpredEup
        real(r8) :: mpredMZS, mpredMZL
!  sinking
        real(r8) :: wPhS, wPhL, wDet, wDetF
!  Terms to define the Iron climatology field
        real(r8) :: Feinlo, Feinhi, Feinh, Feofflo, Feoffhi, Feoffh

!  Terms to define respiration
        real(r8) :: respPhS, respPhL, respMZS, respMZL
        real(r8) :: respCop, respNCa, respEup
        real(r8) :: TmaxPhS,TminPhS, Topt_PhS, KtBm_PhS
        real(r8) :: TmaxPhL, TminPhL, Topt_PhL, KtBm_PhL
        real(r8) :: TmaxMZS, KtBm_MZS, TmaxMZL, KtBm_MZL
        real(r8) :: ktbmC,TrefC
        real(r8) :: ktbmN,TrefN
        real(r8) :: ktbmE,TrefE
!  Detrital Remineralization and Nitrification
        real(r8) :: regen, dgrad
        real(r8) :: Pv0, PvT
        real(r8) :: KnT, Nitr0,ToptNtr,ktntr,KNH4Nit
        real(r8) :: tI0,KI

!  Iron limitation terms
        real(r8) :: kfePhS, kfePhL, FeC
!  Diapause
        real(r8) :: NCmaxz
        real(r8) :: wNCrise,wNCsink
        real(r8) :: RiseStart, RiseEnd, SinkStart, SinkEnd
# ifdef BENTHIC
        real(r8) :: bmB,ktbmB,TrefB
        real(r8) :: iremin
        real(r8) :: q10,q10r
        real(r8) :: Rup,KupP,LupD, LupP,KupD
        real(r8) :: Qres,Rres,rmort,eex,eexD,BenPred
        real(r8) :: prefD,prefPS,prefPL,T0ben,T0benr
# endif
# ifdef ICE_BIO
        real(r8) :: alphaIb, betaI,  inhib
        real(r8) :: ksnut1,ksnut2,mu0, R0i
        real(r8) :: rg0,rg,annit
        real(r8) :: aidz
# endif

# ifdef FEAST

  !  FEAST parameters
 
!!MOD_PARAM   ! Local set constant parameters (e.g. #defines in C)
!!MOD_PARAM     real(r8),    parameter :: EPSILON      = 1e-20
!!MOD_PARAM     !integer, parameter :: ROMS_LAYERS    = 60  ! MUST MATCH UBk
!!MOD_PARAM     integer, parameter :: LAYERS         =  3  !// Number of vert layers for fish
!!MOD_PARAM     !integer, parameter :: FEAST_PASSIVES = 21  ! FEAST passive tracers
!!MOD_PARAM
!!MOD_PARAM!integer, parameter :: NUM_ROMS_DEPTHS       = 60
!!MOD_PARAM!integer, parameter :: NUM_PLANKTON          =  4 ! COP, NCA, EUP, BEN
!!MOD_PARAMinteger, parameter :: NUM_PLANKTON          =  6 ! COP, NCAS, NCA0, EUPS, EUPO, BEN
!!MOD_PARAMinteger, parameter :: MAX_LCLASS            = 60 !// Maximum length classes
!!MOD_PARAMinteger, parameter :: NUM_AGES              = 11 ! nac
!!MOD_PARAMinteger, parameter :: NUM_AGED_LENGTHS      = 14 ! nlc
!!MOD_PARAMinteger, parameter :: NUM_NOAGE_LENGTHS     = 20 !nlct
!!MOD_PARAMinteger, parameter :: NUM_AGED_SPECIES      = 3 ! POL, COD, ATF
!!MOD_PARAMinteger, parameter :: NUM_LENGTHED_SPECIES  = 7 ! HER, CAP, EUL, SAN, MYC, SAL1, SAL2 
!!MOD_PARAMinteger, parameter :: NUM_SIMPLE_SPECIES    = 5 ! SHR, SQU, EPI, CRA, OTH 
!!MOD_PARAMinteger, parameter :: NUM_GEARS             = 8 ! COD_CVPOT, COD_CVLGL, COD_CVTWL, COD_CPPOT, COD_CPLGL, COD_CPTWL, POL_CVTWL, POL_CPTWL
!!MOD_PARAMinteger, parameter :: NUM_TVAR              = 20
!!MOD_PARAM
!!MOD_PARAMinteger, parameter :: nfvaral = 3          !Numbers, CF, Cal
!!MOD_PARAMinteger, parameter :: nfvarl  = 3          !Numbers, CF, Cal
!!MOD_PARAMinteger, parameter :: nfvar   = 2          !Biomass, Cal
!!MOD_PARAMinteger, parameter :: nrates  = 5          !Catch,happiness,dN,dCF,dCAL
!!MOD_PARAM
!!MOD_PARAMinteger, parameter :: TOT_LENGTHED = NUM_AGED_SPECIES + NUM_LENGTHED_SPECIES
!!MOD_PARAMinteger, parameter :: TOT_FEAST    = NUM_AGED_SPECIES + NUM_LENGTHED_SPECIES + NUM_SIMPLE_SPECIES
!!MOD_PARAMinteger, parameter :: TOT_LINKS = NUM_AGED_SPECIES * NUM_AGES * NUM_AGED_LENGTHS + NUM_NOAGE_LENGTHS * NUM_LENGTHED_SPECIES + 1
                                     
  !     TYPE species_params   !struct species_params{
          !  Identifiers
   	         !character(len=20) Name 
             !integer, dimension(NUM_SPECIES)  :: Guild 
          !  Parameters         
        integer ALL_LINKS
        integer, dimension(TOT_LENGTHED)   :: fsh_Lsize     !// width of length bins
        integer, dimension(NUM_AGED_SPECIES,NUM_AGES) :: fsh_age_offset
        real(r8), dimension(TOT_LENGTHED)  :: fsh_a_T, fsh_b_T, fsh_c_T
        real(r8), dimension(TOT_LENGTHED)  :: fsh_A_L, fsh_B_L 
        !real(r8), dimension(TOT_LENGTHED)  :: fsh_H, fsh_theta 
        real(r8), dimension(TOT_FEAST,TOT_FEAST+NUM_PLANKTON) :: fsh_q_G
				integer, dimension(TOT_LENGTHED)   :: fsh_base_prey      
				real(r8), dimension(TOT_LENGTHED,TOT_FEAST+NUM_PLANKTON)   :: fsh_alpha_G, fsh_beta_G  
        real(r8), dimension(TOT_LENGTHED)  :: fsh_A_S, fsh_B_S 
        !real(r8), dimension(TOT_LENGTHED)  :: fsh_a_C, fsh_b_C !fsh_CQ
        !real(r8), dimension(TOT_LENGTHED)  :: fsh_A_Q, fsh_B_Q 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_ED_m, fsh_ED_b 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_A_R, fsh_B_R 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_F_A, fsh_U_A, fsh_SDA 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_R_TM, fsh_R_T0, fsh_R_Q 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_C_TM, fsh_C_T0, fsh_C_Q 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_omega, fsh_mu, fsh_zeta, fsh_smallomega 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_g_W 
        real(r8), dimension(TOT_LENGTHED)  :: fsh_lambda, fsh_phi, fsh_D, fsh_kappa, fsh_chi 
      ! recruitment added 9/6/2009 KYA                
        real(r8), dimension(TOT_LENGTHED)  :: fsh_z_muL, fsh_z_sdL
        real(r8), dimension(TOT_LENGTHED)   :: fsh_z_sday, fsh_z_eday, fsh_sp_sday, fsh_sp_eday
      ! recruitment of 12/23/2101
         real(r8), dimension(TOT_LENGTHED)  :: fsh_mat_a, fsh_mat_b, fsh_fec_a, fsh_fec_b
         real(r8), dimension(TOT_LENGTHED)  :: fsh_fec_max, fsh_fem_prop, fsh_rec_prop
      ! variable activity added 10/28/11
			    real(r8), dimension(TOT_LENGTHED) :: fsh_A_enc, fsh_B_enc
			    real(r8), dimension(TOT_LENGTHED) :: fsh_Bv_min, fsh_B_Lzero, fsh_B_Lone, fsh_B_Lpow
					!real(r8), dimension(TOT_LENGTHED) :: fsh_minVel, fsh_A_v, fsh_B_v, fsh_minACT 
          
					                   !fsh_minVel(sp) = (fsh_A_v(sp)/fsh_minACT(sp)) **(-1.0/fsh_B_v(sp))	            
       ! feeding for simple and zoops
       real(r8), dimension(NUM_SIMPLE_SPECIES) :: fsh_simple_len, fsh_simple_wt, fsh_simple_JperG 
       real(r8), dimension(NUM_PLANKTON) :: fsh_zoop_len, fsh_zoop_wt, fsh_zoop_JperG 
       
      ! movement
        real(r8), dimension(TOT_FEAST)     :: fsh_max_speed, fsh_happy_01, fsh_happy_99  !bodL/day
      ! derived quantities
        real(r8), dimension(TOT_LENGTHED) :: fsh_RTM_RT0, fsh_CTM_CT0 
        real(r8), dimension(TOT_LENGTHED) :: fsh_CTMoverCTM_CT0, fsh_RTMoverRTM_RT0 
        real(r8), dimension(TOT_LENGTHED) :: fsh_X_R, fsh_X_C, fsh_P_j 
        real(r8) :: fpredCop, fpredNcaS, fpredNcaO, fpredEupS, fpredEupO
      ! parameter vectors      
        !real, dimension(TOT_LENGTHED,MAX_LCLASS) :: fsh_LL, fsh_AL_LL_BL, fsh_AS_LL_BS
        !real, dimension(TOT_LENGTHED,MAX_LCLASS) :: fsh_L_deriv, fsh_ED, fsh_zeros
     !  Gear vectors
     	 real, dimension(NUM_GEARS,TOT_FEAST) :: fsh_catch_sel, fsh_catch_01, fsh_catch_99
     	 !real, dimension(NUM_GEARS,TOT_FEAST,MAX_LCLASS) :: fsh_catch_CPUF
     	 !real(r8) fsh_qq(NUM_GEARS,NUM_AGED_SPECIES,NUM_AGES,NUM_AGED_LENGTHS)
     	 
  integer rlook_fal(NUM_AGED_SPECIES,NUM_AGES,NUM_AGED_LENGTHS)
  integer rlook_fl(NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS)
  integer, dimension(TOT_LINKS) :: sp_all,age_all,lc_all,all_order,rev_order
  integer, dimension(TOT_LINKS) :: prev_len,next_len,type_all
  real(r8), dimension(TOT_LINKS) :: all_AL_LL_BL, all_AS_LL_BS
	real(r8), dimension(TOT_LINKS) :: all_Aenc_LL_Benc, all_Ar_v_LL 
	real(r8), dimension(TOT_LINKS) :: all_L_deriv,all_zeros,L_all,CALbase_all, L_width
  real(r8), dimension(TOT_LINKS) :: swim_all, swim_selm, swim_selk
  real(r8), dimension(TOT_LINKS) :: all_mat_rec, sum_preference
  
  real(r8) :: all_QQ(NUM_GEARS,TOT_LINKS), all_preference(TOT_LINKS,TOT_LINKS)
  real(r8) :: zoop_preference(TOT_LINKS,NUM_PLANKTON)
  real(r8) :: simple_preference(TOT_LINKS,NUM_SIMPLE_SPECIES)
  integer :: tspal(nfvaral,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
  integer :: tspl(nfvarl,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS )
  integer :: tsp(nfvar,NUM_SIMPLE_SPECIES)
  integer :: kspal(nfvaral,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
  integer :: kspl(nfvarl,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS )
  integer :: ksp(nfvar,NUM_SIMPLE_SPECIES)	
	integer :: teggs(TOT_LENGTHED), tspawn(TOT_LENGTHED)
	integer :: keggs(TOT_LENGTHED), kspawn(TOT_LENGTHED)
	      	 
  !    END TYPE   
  !     real(r8), dimension(NUM_PLANKTON) ::  fsh_Zoop_JperG 

  real(r8) :: feast_test
  real(r8) :: force_ration
  real(r8) :: force_vonB
  integer :: dye_update
  integer :: dye_loaded = 0 
  integer :: feast_mixed
  integer :: feast_coupled
  integer :: feast_mort
  integer :: feast_fishing
  integer :: feast_growth
  integer :: feast_movement
  integer :: feast_recruitment       ! flag to do or not do recruitment
  integer :: feast_promotion         ! flag to do promotion or not
  integer :: feast_time = 0          ! Set the time of day at which FEAST is updated
  integer :: feast_sp_view  = 1      ! species to output detailed statistics
  integer :: feast_age_view = 1      ! age to output detailed statistics
  !real(r8), dimension(TOT_FEAST) :: age_zeros


# endif
!#  if defined BIOFLUX
!        real(r8), dimension(NAT+NBT,NAT+NBT) :: bflx
!#  endif
# endif
       




# ifdef BIO_GOANPZ
!=======================================================================
!   Parameters for Gulf of Alaska Biological Model
!      Sarah Hinckley and Elizabeth Dobbins
!=======================================================================
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!=======================================================================
        USE mod_param

        implicit none

        integer, dimension(Ngrids) :: BioIter
        real(r8) :: VertMixIncr
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
!  Bio- conversions
        real(r8) :: xi, ccr, ccrPhL
!  extinction coefficients
        real(r8) :: k_ext, k_chl
!  PhS growth params
        real(r8) :: DiS, DpS
        real(r8) :: alphaPhS
        real(r8) :: psiPhS
        real(r8) :: k1PhS, k2PhS
!  PhL growth params
        real(r8) :: DiL
        real(r8) :: DpL
        real(r8) :: alphaPhL
        real(r8) :: psiPhL
        real(r8) :: k1PhL
        real(r8) :: k2PhL
!  MZS preference params
        real(r8) :: fpPhSMZS, fpPhLMZS
!  MZS growth and feeding params
        real(r8) :: eMZS
        real(r8) :: Q10MZS
        real(r8) :: Q10MZST
        real(r8) :: fMZS
        real(r8) :: kMZS
        real(r8) :: gammaMZS
!  MZL preferences params
        real(r8) :: fpPhSMZL, fpPhLMZL, fpMZSMZL
!  MZL growth and feeding params
        real(r8) :: eMZL
        real(r8) :: Q10MZL
        real(r8) :: Q10MZLT
        real(r8) :: fMZL
        real(r8) :: kMZL
        real(r8) :: gammaMZL
!  Cop preference params
        real(r8) :: fpPhSCop, fpPhLCop, fpMZSCop, fpMZLCop
!  Cop growth and feeding params
        real(r8) :: eCop
        real(r8) :: Q10Cop
        real(r8) :: Q10CopT
        real(r8) :: fCop
        real(r8) :: gammaCop
        real(r8) :: kCop
!  NCa preference params
        real(r8) :: fpPhSNCa, fpPhLNCa, fpMZSNCa, fpMZLNCa
!  NCa growth and feeding params
        real(r8) :: eNCa
        real(r8) :: Q10NCa
        real(r8) :: Q10NCaT
        real(r8) :: fNCa
        real(r8) :: gammaNCa
        real(r8) :: kNCa
!  Eup preference params
        real(r8) :: fpPhSEup, fpPhLEup, fpMZSEup, fpMZLEup, fpCopEup
!  Eup growth and feeding params
        real(r8) :: eEup
        real(r8) :: Q10Eup
        real(r8) :: Q10EupT
        real(r8) :: fEup
        real(r8) :: gammaEup
        real(r8) :: kEup
!  Phytoplankton senescence
        real(r8) :: minmPhS, maxmPhS, NcritPhS
        real(r8) :: minmPhL, maxmPhL, NcritPhL
!  Zoopkankton mortality
        real(r8) :: mMZS, mMZL, mCop, mNCa, mEup
!  predation closure
        real(r8) :: mpredCop, mpredNCa, mpredEup
	real(r8) :: mpredMZS, mpredMZL
!  sinking and regeneration terms
        real(r8) :: regen, dgrad
        real(r8) :: wPhS, wPhL, wDet, terms
!  Terms to define the Iron climatology field
        real(r8) :: Feinlo, Feinhi, Feinh, Feofflo, Feoffhi, Feoffh
!  Terms to define respiration
        real(r8) :: respPhS, respPhL, respMZS, respMZL
	real(r8) :: respCop, respNCa, respEup
!  Iron limitation terms
        real(r8) :: kfePhS, kfePhL, FeC
!  Diapause
        real(r8) :: NCmaxz
        real(r8) :: wNCrise,wNCsink
        real(r8) :: RiseStart, RiseEnd, SinkStart, SinkEnd
#  if defined BIOFLUX
        real(r8), dimension(NAT+NBT,NAT+NBT) :: bflx = 0.0_r8
#  endif
# endif
        integer :: dummy_bio
        
#ifdef FEAST
      CONTAINS
!!-----------------------------------------------------------------------------
          subroutine QsortC(A, index)
            real, intent(in out), dimension(:) :: A
            integer, intent(in out), dimension(:) :: index
            integer :: L,i,j,tind
            real :: temp
            
            L=size(A)
            do i = 1,L-1
              do j = i+1,L
                 if (A(i) < A(j)) then
                  temp = A(i)
                  A(i) = A(j)
                  A(j) = temp
                  tind = index(i)
                  index(i) = index(j)
                  index(j) = tind          
                 end if
              end do
            end do
          end subroutine QsortC

!!------------------------------------------------------------------------------
SUBROUTINE initialize_feast
      
  USE mod_param
          
  implicit none
  integer :: isp, sp, gr, ac, lc, py, ng, lastT, lastK
  integer :: fvaral,fvarl,fvar,spal,spl,spalindex,splindex,spindex
  integer :: spalindexmax,splindexmax,spindexmax          
  integer :: NL
  real(r8) :: Z,Y, mu, sd, lp, selk, selm, sels, Smax, Ratio
  real(r8) :: rec,mat
  integer:: i,L, NN, SS, spn, spy, YY
  !integer :: pcount(TOT_FEAST,200)
  real(r8), dimension(TOT_LINKS) :: prob 
  real(r8), dimension(TOT_FEAST) :: probtot
  real, dimension(TOT_LINKS) :: tmp_order
! initialize paramters for fish

  real(r8),    parameter :: PIPI = 3.14159265358979
    
  ! print *, 'Initializing FEAST variables'
  open (unit = 5, file = "init.txt", status='REPLACE')

! derived and stored bioenergetics quantities    
  DO isp=1,TOT_LENGTHED
     !fsh_minVel(isp) = (fsh_A_v(isp)/fsh_minACT(isp)) **(-1.0/fsh_B_v(isp))
    
     fsh_CTM_CT0(isp) = fsh_C_TM(isp) - fsh_C_T0(isp)    
     fsh_CTMoverCTM_CT0(isp) = fsh_C_TM(isp) / fsh_CTM_CT0(isp)
     Z = log(fsh_C_Q(isp)) * fsh_CTM_CT0(isp)
     Y = log(fsh_C_Q(isp)) * (fsh_CTM_CT0(isp) + 2.0)  			 
     fsh_X_C(isp) = ( Z*Z * (1.0 + ( (1.0 + 40.0/Y) ** 0.5)) ** 2.0)/400.0 
		 		  
     fsh_RTM_RT0(isp) = fsh_R_TM(isp) - fsh_R_T0(isp)    
     fsh_RTMoverRTM_RT0(isp) = fsh_R_TM(isp) / fsh_RTM_RT0(isp)   
     Z = log(fsh_R_Q(isp)) * fsh_RTM_RT0(isp)
     Y = log(fsh_R_Q(isp)) * (fsh_RTM_RT0(isp) + 2.0)  			 
     fsh_X_R(isp) = ( Z*Z * (1.0 + ( (1.0 + 40.0/Y) ** 0.5)) ** 2.0)/400.0 

     fsh_P_j(isp) = (fsh_F_A(isp) + (fsh_U_A(isp) + fsh_SDA(isp))*(1.0-fsh_F_A(isp))) 
  END DO
          
          NN = 0
          all_QQ=0.0
          
          			rlook_fal= 0
          			rlook_fl=0
                sp_all= 0
          			age_all= 0
          			lc_all= 0
          			all_order= 0
                rev_order=0 
          			!labs_all= 0
          			
          do sp = 1,NUM_AGED_SPECIES
            do ac = 1,NUM_AGES
              do lc = 1,NUM_AGED_LENGTHS
                NN = NN+1
          			rlook_fal(sp,ac,lc) = NN
                sp_all(NN)   = sp
          			age_all(NN)  = ac
          			lc_all(NN)   = lc
          			all_order(NN) = NN 
          			!labs_all(NN) = lc + fsh_age_offset(sp,ac)
                type_all(NN)  = 1
          		end do
           end do
          end do
          !write(5,"(I5 I5)") NN, ALL_LINKS
          do sp = 1,NUM_LENGTHED_SPECIES
             do lc = 1,NUM_NOAGE_LENGTHS
                NN=NN+1
          			rlook_fl(sp,lc)=NN
                sp_all(NN)   = sp + NUM_AGED_SPECIES
                age_all(NN)  = 0
                lc_all(NN)   = lc
                !labs_all(NN) = lc
                all_order(NN) = NN
                type_all(NN)  = 2 
             end do
          end do
          
          
          ALL_LINKS = NN
          
            !do SS = 1,ALL_LINKS
            !   spn = all_order(SS)
            !   write (5,*),SS,spn,L_all(spn),sp_all(spn),age_all(spn),lc_all(spn),labs_all(spn)
            !end do
          !write(5,"(I5 I5)") NN, ALL_LINKS
          probtot = 0.0
          prob    = 0.0
          
          prev_len = 0
          next_len = 0
          
          do NN = 1,ALL_LINKS
             sp = sp_all(NN)
             lc = lc_all(NN)
             ac = age_all(NN)
             if (ac>0) then
                if (lc>1)                  prev_len(NN) = rlook_fal(sp,ac,lc-1) ! else 0
                if (lc<NUM_AGED_LENGTHS)   next_len(NN) = rlook_fal(sp,ac,lc+1) ! else 0
                !if (ac>1)                  prev_age(NN) = rlook_fal(sp,ac-1,lc) ! else 0
                !if (ac<NUM_AGES)           next_age(NN) = rlook_fal(sp,ac+1,lc) ! else 0
             else
                if (lc>1)                  prev_len(NN) = rlook_fl(sp-NUM_AGED_SPECIES,lc-1)  ! else 0
                if (lc<NUM_NOAGE_LENGTHS)  next_len(NN) = rlook_fl(sp-NUM_AGED_SPECIES,lc+1)	 ! else 0
          	 end if 
						    ! Age 0s and noaged species are half the length set for adults
						 
						 if (ac>0) then
						   if (ac<2) then
								L_width(NN) = fsh_Lsize(sp_all(NN))/2.0
								L_all(NN)   = lc * L_width(NN) - L_width(NN)/2.0
						   else            
								L_width(NN) = fsh_Lsize(sp_all(NN))
							  L_all(NN)   = (lc+fsh_age_offset(sp,ac))*L_width(NN) - L_width(NN)/2.0
               end if
						 else
								L_width(NN) = fsh_Lsize(sp_all(NN))
							  L_all(NN)   = lc*L_width(NN) - L_width(NN)/2.0
						 end if	
						 					    
						 !print *,sp,ac,lc,L_all(NN) 
                CALbase_all(NN)  = fsh_ED_b(sp) + fsh_ED_m(sp)*L_all(NN)
                all_AL_LL_BL(NN) = fsh_A_L(sp) * (L_all(NN) ** fsh_B_L(sp))    
                all_AS_LL_BS(NN) = fsh_A_S(sp) * (L_all(NN) ** fsh_B_S(sp))
          			all_Aenc_LL_Benc(NN) = fsh_A_enc(sp) * (L_all(NN) ** fsh_B_enc(sp))
          			all_Ar_v_LL(NN) = fsh_Bv_min(sp)+fsh_B_Lzero(sp) / &
								                & (1.0 +(fsh_B_Lzero(sp)/fsh_B_Lone(sp)-1.0) * &
																& (L_all(NN) ** fsh_B_Lpow(sp)))
								
								!write(*,'(A10,I4,4G14.6)'),"stuff",sp,L_all(NN),fsh_A_enc(sp),fsh_B_enc(sp),all_Aenc_LL_Benc(NN) 
								!write(*,'(A10,I4,2G14.6)'),"Lpow",sp,L_all(NN),all_Ar_v_LL(NN)
								!write(*,'(A10,G14.6)'),"Bv_min",fsh_Bv_min(sp)
								!write(*,'(A10,G14.6)'),"Lzero",fsh_B_Lzero(sp)
								!write(*,'(A10,G14.6)'),"Lone",fsh_B_Lone(sp)
								!write(*,'(A10,G14.6)'),"Lpow",fsh_B_Lpow(sp)
																								
								all_L_deriv(NN)  = fsh_A_L(sp) * fsh_B_L(sp) * L_all(NN) ** (fsh_B_L(sp) - 1.0)
          			!all_ED(NN)       = fsh_ED_b(sp) + fsh_ED_m(sp) * all_AS_LL_BS(NN)
          			
          		  swim_all(NN)  = L_all(NN) * fsh_max_speed(sp) / 100.0  ! Units are m/s
								swim_selk(NN) = -log((1.0-0.99)/0.99)/((fsh_happy_99(sp)-fsh_happy_01(sp))/2.0)  !=-LN((1-0.99)/0.99)/((L_99-L_01)/2) 
          			swim_selm(NN) = fsh_happy_01(sp)+(fsh_happy_99(sp)-fsh_happy_01(sp))/2.0 !=L_01+(L_99-L_01)/2

                do gr = 1, NUM_GEARS
                   if (fsh_catch_sel(gr,sp)>0.0) then
             		      selk  = -log((1.0-0.99)/0.99)/((fsh_catch_99(gr,sp)-fsh_catch_01(gr,sp))/2.0)  !=-LN((1-0.99)/0.99)/((L_99-L_01)/2) 
                      selm  = fsh_catch_01(gr,sp)+(fsh_catch_99(gr,sp)-fsh_catch_01(gr,sp))/2.0 !=L_01+(L_99-L_01)/2
          		        sels  = selm-(log((1.0/0.001)-1.0))/selk
          		        !write (5,*) selk,selm,sels,L_all(NN),1.0/(1.0 + exp(-selk*(L_all(NN)-selm)))
                      if (L_all(NN)>=sels) then
                       		all_QQ(gr,NN) = fsh_catch_sel(gr,sp) * 1.0/(1.0 + exp(-selk*(L_all(NN)-selm)))
                      end if
                      !write (5,*) selk,selm,sels,L_all(NN),1.0/(1.0 + exp(-selk*(L_all(NN)-selm))),fsh_catch_sel(gr,sp),all_QQ(gr,NN)
                    end if
                end do
                ! recruitment
                  mat =             1.0 / (1.0 + exp( -(fsh_mat_a(sp) + fsh_mat_b(sp) * L_all(NN))))
									if (mat < 0.05) then
									   mat = 0.0
									end if
                  rec = fsh_fec_max(sp) / (1.0 + exp( -(fsh_fec_a(sp) + fsh_fec_b(sp) * L_all(NN))))  
									all_mat_rec(NN) = mat * rec  * fsh_fem_prop(sp)  ! * fsh_rec_prop(sp)
                  
                  if (ac.eq.1 .or. ac.eq.0) then
                    mu = fsh_z_muL(sp)  
                    sd = fsh_z_sdL(sp)
									  lp = log(mu/L_all(NN))
									  prob(NN) = (1.0/sqrt(2.0*PIPI*sd*sd))*exp(-(lp*lp)/(2.0*sd*sd))
                    probtot(sp) = probtot(sp) + prob(NN)
                  end if
          end do
          
      ! 10/27/11 Made preference a normalized sum
          sum_preference = 0.0
            !#prefbase <-  ifelse(pynum==BASE_PREY,
            !#          ifelse(logratio <= py_Smax,  py_q,  py_q*((logratio/py_Smax)^(py_alpha-1))*exp(-(logratio-py_Smax)/py_beta) ),
            !#          ifelse(logratio <= 0.0,       0.0,  py_q*((logratio/py_Smax)^(py_alpha-1))*exp(-(logratio-py_Smax)/py_beta) )
					  !#  )
          !pcount = 0
					do spn = 1,ALL_LINKS
             sp=sp_all(spn)
             do spy = 1,ALL_LINKS
                py=sp_all(spy)
          		! 11/3/11 Made cutoff 80% of length, not 100%, maybe replace 0.80
          		! with species-specific parameter at some point
								if (L_all(spn)*0.80 <= L_all(spy)) THEN 
          			   all_preference(spn,spy)=0.0
          			else
          		   Smax  = (fsh_alpha_G(sp,py)-1.0) * fsh_beta_G(sp,py)
                 Ratio = log (L_all(spn)/L_all(spy))
								 ! uncomment this if-else for saturation
								 !if (Ratio .ge. Smax) then
								 !	 all_preference(spn,spy) = fsh_q_G(sp,py)
								 !else				                 
                   all_preference(spn,spy) = fsh_q_G(sp,py) *                  & 
          					              &((Ratio/Smax) ** (fsh_alpha_G(sp,py)-1.0))  &
          						            & * exp(-(Ratio-Smax)/fsh_beta_G(sp,py))	 	            
          			 !end if
								end if
						    sum_preference(spn) = sum_preference(spn) + all_preference(spn,spy)
						    !if (pcount(sp,int(L_all(spn))) .eq. 0) then
                !write(*,'(I4,F10.2,I4,F10.2,G14.6)'),sp,real(L_all(spn)),py,real(L_all(spy)),all_preference(spn,spy)
								!end if
             end do
             !pcount(sp,int(L_all(spn)))=1
          end do
          
          !pcount = 0
          do spn = 1,ALL_LINKS
             sp=sp_all(spn)
             do spy = 1,NUM_SIMPLE_SPECIES
              py=spy + NUM_AGED_SPECIES + NUM_LENGTHED_SPECIES
          		if (fsh_base_prey(sp) == py) then
          		  Smax  = (fsh_alpha_G(sp,py)-1.0) * fsh_beta_G(sp,py)
          		  Ratio = log (L_all(spn)/fsh_simple_len(spy))
          		  if (Ratio .le. Smax ) then
							     simple_preference(spn,spy) = fsh_q_G(sp,py)
							  else
                   simple_preference(spn,spy) = fsh_q_G(sp,py) *               & 
          					              &((Ratio/Smax) ** (fsh_alpha_G(sp,py)-1.0))  &
          						            & * exp(-(Ratio-Smax)/fsh_beta_G(sp,py))							  
							  end if
							else					  
          			if (L_all(spn) <= fsh_simple_len(spy)) THEN 
          			   simple_preference(spn,spy)=0.0
          			else
                 Smax  = (fsh_alpha_G(sp,py)-1.0) * fsh_beta_G(sp,py)
								 Ratio = log (L_all(spn)/fsh_simple_len(spy))
								 !if (Ratio .ge. Smax) then
								 !  simple_preference(spn,spy) = fsh_q_G(sp,py)
								 !else						                 
                   simple_preference(spn,spy) = fsh_q_G(sp,py) *               & 
          					              &((Ratio/Smax) ** (fsh_alpha_G(sp,py)-1.0))  &
          						            & * exp(-(Ratio-Smax)/fsh_beta_G(sp,py))
          			 !end if
								end if    
          		end if
          		sum_preference(spn) = sum_preference(spn) + simple_preference(spn,spy)
						    !if (pcount(sp,int(L_all(spn))) .eq. 0) then
              !write(*,'(I4,F10.2,I4,F10.2,G14.6)'),sp,real(L_all(spn)),py,real(fsh_simple_len(spy)),simple_preference(spn,spy)

								!end if
             end do
             !pcount(sp,int(L_all(spn)))=1
          end do
          
          !pcount = 0     
          do spn = 1,ALL_LINKS
             sp=sp_all(spn)
             !write(0,*),sp,fsh_base_prey(sp)
             do spy = 1,NUM_PLANKTON
               py=spy + NUM_AGED_SPECIES + NUM_LENGTHED_SPECIES + NUM_SIMPLE_SPECIES
               if (fsh_base_prey(sp) == py) then 
          		    Smax  = (fsh_alpha_G(sp,py)-1.0) * fsh_beta_G(sp,py)
          		    Ratio = log (L_all(spn)/fsh_zoop_len(spy))
          		    if (Ratio .le. Smax ) then
                     zoop_preference(spn,spy) = fsh_q_G(sp,py)
                  else
									 zoop_preference(spn,spy) = fsh_q_G(sp,py) *                 & 
          					              &((Ratio/Smax) ** (fsh_alpha_G(sp,py)-1.0))  &
          						            & * exp(-(Ratio-Smax)/fsh_beta_G(sp,py))
									end if
									!write(0,*),sp,py,Ratio,Smax,zoop_preference(spn,spy)                 
							 else
          			if (L_all(spn) <= fsh_zoop_len(spy)) THEN 
          			   zoop_preference(spn,spy)=0.0
          			else
          		    Smax  = (fsh_alpha_G(sp,py)-1.0) * fsh_beta_G(sp,py)
                  Ratio = log (L_all(spn)/fsh_zoop_len(spy))						                 
                  !if (Ratio .ge. Smax) then
                  ! zoop_preference(spn,spy) = fsh_q_G(sp,py)
									!else   
									 zoop_preference(spn,spy) = fsh_q_G(sp,py) *                 & 
          					              &((Ratio/Smax) ** (fsh_alpha_G(sp,py)-1.0))  &
          						            & * exp(-(Ratio-Smax)/fsh_beta_G(sp,py))
									!end if	 
          			end if
          		 end if	
							sum_preference(spn)  = sum_preference(spn) + zoop_preference(spn,spy) 
						    !if (pcount(sp,int(L_all(spn))) .eq. 0) then
          		  !  write(*,'(I4,F10.2,I4,F10.2,G14.6)'),sp,real(L_all(spn)),py,real(fsh_zoop_len(spy)),zoop_preference(spn,spy)
								!end if 
             end do
            !pcount(sp,int(L_all(spn)))=1
          end do 
 
         ! Regardless of above, set BASE_PREY preference to 1.0

				   
 
         ! Use following if preference is normalized
          !do spn = 1,ALL_LINKS
          !  if (sum_preference(spn) > 0.0) then
          !      do spy = 1,ALL_LINKS
          !         all_preference(spn,spy) = all_preference(spn,spy)/sum_preference(spn)
          !      end do
          !      do spy = 1,NUM_SIMPLE_SPECIES
          !         simple_preference(spn,spy) = simple_preference(spn,spy)/sum_preference(spn)
          !      end do
          !      do spy = 1,NUM_PLANKTON
          !         zoop_preference(spn,spy) = zoop_preference(spn,spy)/sum_preference(spn)
          !      end do
          !  end if
          !end do
                
          do NN = 1,ALL_LINKS
             sp = sp_all(NN)
             all_zeros(NN) = prob(NN)/probtot(sp) 
             if (all_zeros(NN)<1e-20) then
                all_zeros(NN) = 0.0
             end if
          end do
          
          
          !write(5,"(I5 I5)") NN, ALL_LINKS
          
          
          !probtot = 0.0
          !mu = fsh_z_muL(isp)  
          !sd = fsh_z_sdL(isp)  
          !
          ! DO L=1,MAX_LCLASS
          !    lp = fsh_LL(isp,L)
          !    prob(L) = (1.0/sqrt(2.0*PIPI*sd*sd)) * exp ( -(log(mu/lp)*log(mu/lp)) / (2.0*sd*sd) )
          !    probtot  = probtot+prob(L)
          !    !print *,'INIT,',isp,lp,mu,sd,PIPI,prob(L),probtot,(1.0/sqrt(2.0*PIPI*sd*sd))*exp(-(log(mu/lp)*log(mu/lp))/(2.0*sd*sd))
          ! END DO
          !             
          ! DO L=1,MAX_LCLASS  
          !    prob(L) = prob(L)/ probtot  
          !    fsh_zeros(isp,L) = prob(L) * fsh_z_N(isp)
          !    if (fsh_zeros(isp,L)<1e-20) then
          !       fsh_zeros(isp,L) = 0.0
          !    end if
          !    ! print *,'INIT,',isp,L,mu,sd,prob(L),fsh_z_N,fsh_zeros(L)
          ! end DO
          
          
          
          !    do gr = 1, NUM_GEARS
          !		 selk  = -log((1.0-0.99)/0.99)/((fsh_catch_99(gr,sp)-fsh_catch_01(gr,sp))/2.0)  !=-LN((1-0.99)/0.99)/((L_99-L_01)/2) 
          !     selm  = fsh_catch_01(gr,sp)+(fsh_catch_99(gr,sp)-fsh_catch_01(gr,sp))/2.0 !=L_01+(L_99-L_01)/2
          !		 sels  = selm-(log((1.0/0.001)-1.0))/selk
          !     if (fsh_catch_sel(gr,sp)>0.0) then
          !          do ac = 1,NUM_AGES
          !            do lc = 1,NUM_AGED_LENGTHS
          !                !NS   = rlook(sp,ac,lc)
          !								lp = (lc + fsh_age_offset(sp,ac)) * fsh_Lsize(sp) + fsh_Lsize(sp)/2
          !								if (lp>=sels) then
          !							     fsh_qq(gr,sp,ac,lc) = fsh_catch_sel(gr,sp) * 1.0/(1.0 + exp(-selk*(lp-selm)))
          !                else
          !								   fsh_qq(gr,sp,ac,lc) = 0.0
          !								end if
          !            end do
          !          end do
          !      end if
          !    end do
            tmp_order = L_all
            call QsortC(tmp_order,all_order)
            do NN = 1,ALL_LINKS
               rev_order(all_order(NN))=NN
            end do
            
				 write(5,'(3A4,A6,A4,A4,A14)')," pdS"," pdA"," PdL", &
				      & "  Prey","pyA","PyL","Pref"              						
						!do SS = 1,ALL_LINKS
						!
            !   spn = all_order(SS)
            !   DO py = 1,NUM_PLANKTON
						!      write(5,'(I4,I4,F4.0,A4,I2.2,I4,F4.0,G14.6)'), &
						!			sp_all(spn), age_all(spn), real(L_all(spn)), &
						!	    &	'ZZ_',py, 0, 0.0, &
						!	    & real(zoop_preference(spn,py))
						!	 end do
            !  DO py = 1,NUM_SIMPLE_SPECIES
						!      write(5,'(I4,I4,F4.0,A4,I2.2,I4,F4.0,G14.6)'), &
						!			sp_all(spn), age_all(spn), real(L_all(spn)), &
						!	    &	'SS_',py, 0, 0.0, &
						!	    & real(simple_preference(spn,py))
						!	 end do							 
						!	 do YY = 1,ALL_LINKS
						!	    spy = all_order(YY)							    
						!     write(5,'(I4,I4,F4.0,A4,I2.2,I4,F4.0,G14.6)'), &
						!			sp_all(spn), age_all(spn), real(L_all(spn)), &
						!	    &	'FF_',sp_all(spy), age_all(spy), L_all(spy), &
						!	    & real(all_preference(spn,spy))
						!	 end do			
						!	 !write (5,*),sp_all(spn),age_all(spn),L_all(spn),all_AL_LL_BL(spn),all_mat_rec(spn)
            !end do
            
            close(5)

!subroutine fish_indices(tspal,tspl,tsp,kspal,kspl,ksp,NL)

!USE mod_param
!USE mod_feast

!integer, intent(inout) :: tspal(nfvaral,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
!integer, intent(inout) :: tspl(nfvarl,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS )
!integer, intent(inout) :: tsp(nfvar,NUM_SIMPLE_SPECIES)
!integer, intent(inout) :: kspal(nfvaral,NUM_AGED_SPECIES,NUM_AGED_LENGTHS,NUM_AGES)
!integer, intent(inout) :: kspl(nfvarl,NUM_LENGTHED_SPECIES,NUM_NOAGE_LENGTHS )
!integer, intent(inout) :: ksp(nfvar,NUM_SIMPLE_SPECIES)
ng = 1
NL = N(ng) !ROMS_LAYERS !integer :: NL
!print *,NAT,NPT,NL

spalindexmax = nfvaral*NUM_AGED_SPECIES*NUM_AGED_LENGTHS*NUM_AGES
splindexmax = nfvarl*NUM_LENGTHED_SPECIES*NUM_NOAGE_LENGTHS 
spindexmax = nfvar*NUM_SIMPLE_SPECIES

open (unit = 5, file = "dyes.txt", status='REPLACE')
write (5,'(A6,6A4)'),"typ","dye","lev","var","sp","lc","ac"

do fvaral=1,nfvaral
do spal=1,NUM_AGED_SPECIES
do lc=1,NUM_AGED_LENGTHS
do ac=1,NUM_AGES
spalindex = (fvaral-1)*NUM_AGED_SPECIES*NUM_AGED_LENGTHS*NUM_AGES + (spal-1)*NUM_AGED_LENGTHS*NUM_AGES + (lc-1)*NUM_AGES +ac
tspal(fvaral,spal,lc,ac) = int((spalindex-1)/NL) + 1 + NAT
kspal(fvaral,spal,lc,ac) = mod((spalindex-1),NL) + 1


write (5,'(A6,6I4)'),"age",tspal(fvaral,spal,lc,ac)-NAT,kspal(fvaral,spal,lc,ac),fvaral,spal,lc,ac
!print*,"fvaral,spal,lc,ac,spalindex,&
!           &tspal(fvaral,spal,lc,ac),&
!           &kspal(fvaral,spal,lc,ac)",&
!           &fvaral,spal,lc,ac,spalindex,&
!           &tspal(fvaral,spal,lc,ac),&
!           &kspal(fvaral,spal,lc,ac)

end do
end do
end do
end do

do fvarl=1,nfvarl
do spl=1,NUM_LENGTHED_SPECIES
do lc=1,NUM_NOAGE_LENGTHS 
splindex = spalindexmax + (fvarl-1)*NUM_LENGTHED_SPECIES*NUM_NOAGE_LENGTHS  + (spl-1)*NUM_NOAGE_LENGTHS  +lc
tspl(fvarl,spl,lc) = int((splindex-1)/NL) + 1 + NAT
kspl(fvarl,spl,lc) = mod((splindex-1),NL) + 1

write (5,'(A6,6I4)'),"len",tspl(fvarl,spl,lc)-NAT,kspl(fvarl,spl,lc),fvarl,spl,lc,0
!print*,"fvarl,spl,lc,splindex,&
!           &tspl(fvarl,spl,lc),&
!           &kspl(fvarl,spl,lc)",&
!           &fvarl,spl,lc,splindex,&
!           &tspl(fvarl,spl,lc),&
!           &kspl(fvarl,spl,lc)


end do
end do
end do

do fvar=1,nfvar
do sp=1,NUM_SIMPLE_SPECIES
spindex = spalindexmax + splindexmax + (fvar-1)*NUM_SIMPLE_SPECIES + sp
tsp(fvar,sp) = int((spindex-1)/NL) + 1 + NAT
ksp(fvar,sp) = mod((spindex-1),NL) + 1

write (5,'(A6,6I4)'),"sim",tsp(fvar,sp)-NAT,ksp(fvar,sp),fvar,sp,0,0

lastT = tsp(fvar,sp)
lastK = ksp(fvar,sp)
end do
end do

do sp=1,TOT_LENGTHED
   lastK = lastK + 1
   if (lastK > NL) then
       lastK = 1
       lastT = lastT + 1
   end if
   teggs(sp) = lastT
   keggs(sp) = lastK
   write (5,'(A6,6I4)'),"egg",teggs(sp)-NAT,keggs(sp),0,sp,0,0
   
end do

do sp=1,TOT_LENGTHED
   lastK = lastK + 1
   if (lastK > NL) then
       lastK = 1
       lastT = lastT + 1
   end if
   tspawn(sp) = lastT
   kspawn(sp) = lastK
   write (5,'(A6,6I4)'),"spn",tspawn(sp)-NAT,kspawn(sp),0,sp,0,0
end do

 close(5)
!print *,"last tspawn",tspawn(TOT_LENGTHED)
!print *,"last kspawn",kspawn(TOT_LENGTHED)

!return
!END SUBROUTINE FISH_INDICES


      
      end SUBROUTINE initialize_feast
#endif        
        
#endif

      END MODULE mod_biology
