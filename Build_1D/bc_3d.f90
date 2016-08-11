      MODULE bc_3d_mod
!
!svn $Id: bc_3d.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient or periodic boundary conditions for   !
!  generic 3D fields.                                                  !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r3d_tile        Boundary conditions for field at RHO-points    !
!    bc_u3d_tile        Boundary conditions for field at U-points      !
!    bc_v3d_tile        Boundary conditions for field at V-points      !
!    bc_w3d_tile        Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
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
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
      RETURN
      END SUBROUTINE bc_r3d_tile
!
!***********************************************************************
      SUBROUTINE bc_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_u3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
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
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
      RETURN
      END SUBROUTINE bc_u3d_tile
!
!***********************************************************************
      SUBROUTINE bc_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_v3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
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
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
      RETURN
      END SUBROUTINE bc_v3d_tile
!
!***********************************************************************
      SUBROUTINE bc_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
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
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
      RETURN
      END SUBROUTINE bc_w3d_tile
      END MODULE bc_3d_mod
