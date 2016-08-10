      MODULE mod_param
!
!svn $Id: mod_param.F 1076 2009-09-25 23:18:43Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Grid parameters:                                                    !
!                                                                      !
!  Im         Number of global grid points in the XI-direction         !
!               for each nested grid.                                  !
!  Jm         Number of global grid points in the ETA-direction        !
!               for each nested grid.                                  !
!  Lm         Number of interior grid points in the XI-direction       !
!               for each nested grid.                                  !
!  Mm         Number of internal grid points in the ETA-direction.     !
!               for each nested grid.                                  !
!  N          Number of vertical levels for each nested grid.          !
!  Ngrids     Number of nested and/or connected grids to solve.        !
!  NtileI     Number of XI-direction tiles or domain partitions for    !
!               each nested grid. Values used to compute tile ranges.  !
!  NtileJ     Number of ETA-direction tiles or domain partitions for   !
!               each nested grid. Values used to compute tile ranges.  !
!  NtileX     Number of XI-direction tiles or domain partitions for    !
!               each nested grid. Values used in parallel loops.       !
!  NtileE     Number of ETA-direction tiles or domain partitions for   !
!               each nested grid. Values used in parallel loops.       !
!  HaloSizeI  Maximum halo size, in grid points, in XI-direction.      !
!  HaloSizeJ  Maximum halo size, in grid points, in ETA-direction.     !
!  TileSide   Maximun tile side length in XI- or ETA-directions.       !
!  TileSize   Maximum tile size.                                       !
!                                                                      !
!  Configuration parameters:                                           !
!                                                                      !
!  Nbico      Number of balanced SSH elliptic equation iterations.     !
!  Nfloats    Number of floats tracjectories.                          !
!  Nstation   Number of output stations.                               !
!  MTC        Maximum number of tidal components.                      !
!  MTY        Maximum number of tidal nodal correction years           !
!  NSA        Number of state array for error covariance.              !
!  NSV        Number of model state variables.                         !
!                                                                      !
!  Tracer parameters:                                                  !
!                                                                      !
!  NAT        Number of active tracer type variables (usually,         !
!               NAT=2 for potential temperature and salinity).         !
!  NBT        Number of biological tracer type variables.              !
!  NST        Number of sediment tracer type variables (NCS+NNS).      !
!  NPT        Number of extra passive tracer type variables to         !
!               advect and diffuse only (dyes, etc).                   !
!  NT         Total number of tracer type variables.                   !
!  MT         Maximum number of tracer type variables.                 !
!                                                                      !
!  Nbed       Number of sediment bed layers.                           !
!  NCS        Number of cohesive (mud) sediment tracers.               !
!  NNS        Number of non-cohesive (sand) sediment tracers.          !
!                                                                      !
!  NBands     Number of spectral irradiance bands in EcoSim.           !
!  Nbac       Number of bacteria constituents in EcoSim.               !
!  Ndom       Number of dissolved matter constituents in EcoSim.       !
!  Nfec       Number of fecal matter constituents in EcoSim.           !
!  Nphy       Number of phytoplankton constituents in EcoSim.          !
!  Npig       Number of pigment constituents in EcoSim.                !
!  PHY        EcoSim indices of phytoplankton species considered.      !
!  PIG        EcoSim phytoplankton-pigment matrix.                     !
!                                                                      !
!  Diagnostic fields parameters:                                       !
!                                                                      !
!  NDbio2d    Number of diagnostic 2D biology fields.                  !
!  NDbio3d    Number of diagnostic 3D biology fields.                  !
!  NDT        Number of diagnostic tracer fields.                      !
!  NDM2d      Number of diagnostic 2D momentum fields.                 !
!  NDM3d      Number of diagnostic 3D momentum fields.                 !
!  NDrhs      Number of diagnostic 3D right-hand-side fields.          !
!
!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none
!
!-----------------------------------------------------------------------
!  Number of nested and/or connected grids to solve.
!-----------------------------------------------------------------------
!
!  Because of ROMS design, the Ngrids parameter is assigned during
!  C-preprocessing before compilation in the "makefile" or "build"
!  script. This is the only way that can be done.
!
        integer, parameter :: Ngrids = 1
!
!-----------------------------------------------------------------------
!  Lower and upper bounds indices per domain partition for all grids.
!-----------------------------------------------------------------------
!
!  All the 1D array indices are of size -1:NtileI(ng)*NtileJ(ng)-1. The
!  -1 index include the values for the full (no partitions) grid.
!
!  Notice that the starting (Imin, Jmin) and ending (Imax, Jmax) indices
!  for I/O processing are 3D arrays. The first dimension (1:4) is for
!  1=PSI, 2=RHO, 3=u, 4=v points; the second dimension (0:1) is number
!  of ghost points (0: no ghost points, 1: Nghost points), and the
!  the third dimension is for 0:NtileI(ng)*NtileJ(ng)-1.
!
        TYPE T_BOUNDS
          integer, pointer :: tile(:)  ! tile partition
          integer, pointer :: LBi(:)   ! lower bound I-dimension
          integer, pointer :: UBi(:)   ! upper bound I-dimension
          integer, pointer :: LBj(:)   ! lower bound J-dimension
          integer, pointer :: UBj(:)   ! upper bound J-dimension
          integer :: LBij              ! lower bound MIN(I,J)-dimension
          integer :: UBij              ! upper bound MAX(I,J)-dimension
          integer :: edge(4,4)         ! boundary edges I- or J-indices
          integer, pointer :: Istr(:)  ! starting tile I-direction
          integer, pointer :: Iend(:)  ! ending   tile I-direction
          integer, pointer :: Jstr(:)  ! starting tile J-direction
          integer, pointer :: Jend(:)  ! ending   tile J-direction
          integer, pointer :: IstrR(:) ! starting tile I-direction (RHO)
          integer, pointer :: IstrT(:) ! starting nest I-direction (RHO)
          integer, pointer :: IstrU(:) ! starting tile I-direction (U)
          integer, pointer :: IendR(:) ! ending   tile I-direction (RHO)
          integer, pointer :: IendT(:) ! ending   nest I-direction (RHO)
          integer, pointer :: JstrR(:) ! starting tile J-direction (RHO)
          integer, pointer :: JstrT(:) ! starting nest J-direction (RHO)
          integer, pointer :: JstrV(:) ! starting tile J-direction (V)
          integer, pointer :: JendR(:) ! ending   tile J-direction (RHO)
          integer, pointer :: JendT(:) ! ending   nest J-direction (RHO)
          integer, pointer :: Imin(:,:,:)  ! starting ghost I-direction
          integer, pointer :: Imax(:,:,:)  ! ending   ghost I-direction
          integer, pointer :: Jmin(:,:,:)  ! starting ghost J-direction
          integer, pointer :: Jmax(:,:,:)  ! ending   ghost J-direction
        END TYPE T_BOUNDS
        TYPE (T_BOUNDS), allocatable :: BOUNDS(:)
!
!-----------------------------------------------------------------------
!  Lower and upper bounds in NetCDF files.
!-----------------------------------------------------------------------
!
        TYPE T_IOBOUNDS
          integer :: ILB_psi       ! I-direction lower bound (PSI)
          integer :: IUB_psi       ! I-direction upper bound (PSI)
          integer :: JLB_psi       ! J-direction lower bound (PSI)
          integer :: JUB_psi       ! J-direction upper bound (PSI)
          integer :: ILB_rho       ! I-direction lower bound (RHO)
          integer :: IUB_rho       ! I-direction upper bound (RHO)
          integer :: JLB_rho       ! J-direction lower bound (RHO)
          integer :: JUB_rho       ! J-direction upper bound (RHO)
          integer :: ILB_u         ! I-direction lower bound (U)
          integer :: IUB_u         ! I-direction upper bound (U)
          integer :: JLB_u         ! J-direction lower bound (U)
          integer :: JUB_u         ! J-direction upper bound (U)
          integer :: ILB_v         ! I-direction lower bound (V)
          integer :: IUB_v         ! I-direction upper bound (V)
          integer :: JLB_v         ! J-direction lower bound (V)
          integer :: JUB_v         ! J-direction upper bound (V)
          integer :: IorJ          ! number of MAX(I,J)-direction points
          integer :: xi_psi        ! number of I-direction points (PSI)
          integer :: xi_rho        ! number of I-direction points (RHO)
          integer :: xi_u          ! number of I-direction points (U)
          integer :: xi_v          ! number of I-direction points (V)
          integer :: eta_psi       ! number of J-direction points (PSI)
          integer :: eta_rho       ! number of J-direction points (RHO)
          integer :: eta_u         ! number of I-direction points (U)
          integer :: eta_v         ! number of I-direction points (V)
        END TYPE T_IOBOUNDS
        TYPE (T_IOBOUNDS) :: IOBOUNDS(Ngrids)
!
!-----------------------------------------------------------------------
!  Tiles minimum and maximum fractional grid coordinates.
!-----------------------------------------------------------------------
!
        TYPE T_DOMAIN
          real(r8), pointer :: Xmin_psi(:)
          real(r8), pointer :: Xmax_psi(:)
          real(r8), pointer :: Ymin_psi(:)
          real(r8), pointer :: Ymax_psi(:)
          real(r8), pointer :: Xmin_rho(:)
          real(r8), pointer :: Xmax_rho(:)
          real(r8), pointer :: Ymin_rho(:)
          real(r8), pointer :: Ymax_rho(:)
          real(r8), pointer :: Xmin_u(:)
          real(r8), pointer :: Xmax_u(:)
          real(r8), pointer :: Ymin_u(:)
          real(r8), pointer :: Ymax_u(:)
          real(r8), pointer :: Xmin_v(:)
          real(r8), pointer :: Xmax_v(:)
          real(r8), pointer :: Ymin_v(:)
          real(r8), pointer :: Ymax_v(:)
        END TYPE T_DOMAIN
        TYPE (T_DOMAIN), allocatable :: DOMAIN(:)
!
!-----------------------------------------------------------------------
!  Model grid(s) parameters.
!-----------------------------------------------------------------------
!
!  Number of interior RHO-points in the XI- and ETA-directions. The
!  size of models state variables (C-grid) at input and output are:
!
!    RH0-type variables:  [0:Lm+1, 0:Mm+1]        ----v(i,j+1)----
!    PSI-type variables:  [1:Lm+1, 1:Mm+1]        |              |
!      U-type variables:  [1:Lm+1, 0:Mm+1]     u(i,j)  r(i,j)  u(i+1,j)
!      V-type variables:  [0:Lm+1, 1:Mm+1]        |              |
!                                                 -----v(i,j)-----
        integer, dimension(Ngrids) :: Lm
        integer, dimension(Ngrids) :: Mm
!
!  Global horizontal size of model arrays including padding.  All the
!  model state arrays are of same size to facilitate parallelization.
!
        integer, dimension(Ngrids) :: Im
        integer, dimension(Ngrids) :: Jm
!
!  Number of vertical levels. The vertical ranges of model state
!  variables are:
!                                                 -----W(i,j,k)-----
!    RHO-, U-, V-type variables: [1:N]            |                |
!              W-type variables: [0:N]            |    r(i,j,k)    |
!                                                 |                |
!                                                 ----W(i,j,k-1)----
        integer, dimension(Ngrids) :: N
       integer, parameter, dimension(Ngrids) :: NBL  =  1
!-----------------------------------------------------------------------
!  Tracers parameters.
!-----------------------------------------------------------------------
!
!  Total number of tracer type variables, NT(:) = NAT + NBT + NPT + NST.
!  The MT corresponds to the maximum number of tracers between all
!  nested grids.
!
        integer, dimension(Ngrids) :: NT
        integer :: MT
!
!  Number of active tracers. Usually, NAT=2 for potential temperature
!  and salinity.
!
        integer :: NAT = 0
!
!  Total number of inert passive tracers to advect and diffuse only
!  (like dyes, etc). This parameter is independent of the number of
!  biological and/or sediment tracers.
!
        integer :: NPT = 0
!
!-----------------------------------------------------------------------
!  Sediment tracers parameters.
!-----------------------------------------------------------------------
!
!  Number of sediment bed layes.
!
        integer :: Nbed = 0
!
!  Total number of sediment tracers, NST = NCS + NNS.
!
        integer :: NST = 0
!
!  Number of cohesive (mud) sediments.
!
        integer :: NCS
!
!  Number of non-cohesive (sand) sediments.
!
        integer :: NNS
!
!-----------------------------------------------------------------------
!  Biological tracers parameters.
!-----------------------------------------------------------------------
!
!  Number of tracers for the Bering Sea ecosystem model.
!
        integer, parameter :: NBT = 14
        integer,parameter :: NBEN=2
        integer :: MBT
        integer, dimension(Ngrids) :: NBeT
!
!-----------------------------------------------------------------------
!  Maximum number of tidal constituents to process.
!-----------------------------------------------------------------------
!
!  Number of tidal constituents to process.
!
        integer :: MTC
!
!-----------------------------------------------------------------------
!  Model state parameters.
!-----------------------------------------------------------------------
!
!  Number of model state variables.
!
        integer, dimension(Ngrids) :: NSV
!
!  Set nonlinear, tangent linear, and adjoint models identifiers.
!
        integer :: iNLM = 1
        integer :: iTLM = 2
        integer :: iRPM = 3
        integer :: iADM = 4
!
!-----------------------------------------------------------------------
!  Domain partition parameters.
!-----------------------------------------------------------------------
!
!  Number of tiles or domain partitions in the XI- and ETA-directions.
!  These values are used to compute tile ranges [Istr:Iend, Jstr:Jend].
!
        integer, dimension(Ngrids) :: NtileI
        integer, dimension(Ngrids) :: NtileJ
!
!  Number of tiles or domain partitions in the XI- and ETA-directions.
!  These values are used to parallel loops to differentiate between
!  shared-memory and distributed-memory.  Notice that in distributed
!  memory both values are set to one.
!
        integer, dimension(Ngrids) :: NtileX
        integer, dimension(Ngrids) :: NtileE
!
!  Maximum number of points in the halo region in the XI- and
!  ETA-directions.
!
        integer, dimension(Ngrids) :: HaloSizeI
        integer, dimension(Ngrids) :: HaloSizeJ
!
!  Maximum tile side length in XI- or ETA-directions.
!
        integer, dimension(Ngrids) :: TileSide
!
!  Maximum number of points in a tile partition.
!
        integer, dimension(Ngrids) :: TileSize
!
!  Set number of ghost-points in the halo region.  It is only used
!  in distributed-memory applications.
!
        integer :: NghostPoints = 2
      CONTAINS
        SUBROUTINE initialize_param
!
!=======================================================================
!                                                                      !
!  This routine initializes several parameters in module "mod_param"   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
        integer :: I_padd, J_padd, ng
!
!-----------------------------------------------------------------------
!  Derived dimension parameters.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          I_padd=(Lm(ng)+2)/2-(Lm(ng)+1)/2
          J_padd=(Mm(ng)+2)/2-(Mm(ng)+1)/2
          Im(ng)=Lm(ng)+I_padd
          Jm(ng)=Mm(ng)+J_padd
          NT(ng)=NAT+NBT+NST+NPT
          NSV(ng)=NT(ng)+5
        END DO
!
!  Set maximum number of tracer between all nested grids.
!
        MT=MAX(2,MAXVAL(NT))
!---------------------------------------------
!Adding stationary and production tracers to model
!---------------------------------------------
      DO ng=1,Ngrids
           NBeT(ng) =NBEN
      END DO
        RETURN
        END SUBROUTINE initialize_param
      END MODULE mod_param
