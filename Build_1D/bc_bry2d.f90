      MODULE bc_bry2d_mod
!
!svn $Id: bc_bry2d.F 933 2009-02-24 19:25:01Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient conditions for generic 2D boundary    !
!  fields.                                                             !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r2d_bry_tile    Boundary conditions for field at RHO-points    !
!    bc_u2d_bry_tile    Boundary conditions for field at U-points      !
!    bc_v2d_bry_tile    Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
! 
!***********************************************************************
      SUBROUTINE bc_r2d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij,                           &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij
      real(r8), intent(inout) :: A(LBij:)
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
!  Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.iwest) THEN
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
      IF (boundary.eq.ieast) THEN
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((Jstr.eq.1).and.(Istr.eq.1)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((Jstr.eq.1).and.(Iend.eq.Lm(ng))) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF
      IF (boundary.eq.inorth) THEN
        IF ((Jend.eq.Mm(ng)).and.(Istr.eq.1)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((Jend.eq.Mm(ng)).and.(Iend.eq.Lm(ng))) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF
      RETURN
      END SUBROUTINE bc_r2d_bry_tile
! 
!***********************************************************************
      SUBROUTINE bc_u2d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij,                           &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij
      real(r8), intent(inout) :: A(LBij:)
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
!  Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.iwest) THEN
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
      IF (boundary.eq.ieast) THEN
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((Jstr.eq.1).and.(Istr.eq.1)) THEN
          A(IstrU-1)=A(IstrU)
        END IF
        IF ((Jstr.eq.1).and.(Iend.eq.Lm(ng))) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF
      IF (boundary.eq.inorth) THEN
        IF ((Jend.eq.Mm(ng)).and.(Istr.eq.1)) THEN
          A(IstrU-1)=A(IstrU)
        END IF
        IF ((Jend.eq.Mm(ng)).and.(Iend.eq.Lm(ng))) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF
      RETURN
      END SUBROUTINE bc_u2d_bry_tile
! 
!***********************************************************************
      SUBROUTINE bc_v2d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij,                           &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij
      real(r8), intent(inout) :: A(LBij:)
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
!  Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.iwest) THEN
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(JstrV-1)=A(JstrV)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
      IF (boundary.eq.ieast) THEN
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(JstrV-1)=A(JstrV)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((Jstr.eq.1).and.(Istr.eq.1)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((Jstr.eq.1).and.(Iend.eq.Lm(ng))) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF
      IF (boundary.eq.inorth) THEN
        IF ((Jend.eq.Mm(ng)).and.(Istr.eq.1)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((Jend.eq.Mm(ng)).and.(Iend.eq.Lm(ng))) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF
      RETURN
      END SUBROUTINE bc_v2d_bry_tile
      END MODULE bc_bry2d_mod
