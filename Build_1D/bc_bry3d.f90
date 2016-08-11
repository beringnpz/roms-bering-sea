      MODULE bc_bry3d_mod
!
!svn $Id: bc_bry3d.F 933 2009-02-24 19:25:01Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient conditions for generic 3D boundary    !
!  fields.                                                             !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r3d_bry_tile    Boundary conditions for field at RHO-points    !
!    bc_u3d_bry_tile    Boundary conditions for field at U-points      !
!    bc_v3d_bry_tile    Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
! 
!***********************************************************************
      SUBROUTINE bc_r3d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij, LBk, UBk,                 &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk
      real(r8), intent(inout) :: A(LBij:,LBk:)
!
!  Local variable declarations.
!
      integer :: k
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
          DO k=LBk,UBk
            A(Jstr-1,k)=A(Jstr,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Jend+1,k)=A(Jend,k)
          END DO
        END IF
      END IF
      IF (boundary.eq.ieast) THEN
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Jstr-1,k)=A(Jstr,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Jend+1,k)=A(Jend,k)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((Jstr.eq.1).and.(Istr.eq.1)) THEN
          DO k=LBk,UBk
            A(Istr-1,k)=A(Istr,k)
          END DO
        END IF
        IF ((Jstr.eq.1).and.(Iend.eq.Lm(ng))) THEN
          DO k=LBk,UBk
            A(Iend+1,k)=A(Iend,k)
          END DO
        END IF
      END IF
      IF (boundary.eq.inorth) THEN
        IF ((Jend.eq.Mm(ng)).and.(Istr.eq.1)) THEN
          DO k=LBk,UBk
            A(Istr-1,k)=A(Istr,k)
          END DO
        END IF
        IF ((Jend.eq.Mm(ng)).and.(Iend.eq.Lm(ng))) THEN
          DO k=LBk,UBk
            A(Iend+1,k)=A(Iend,k)
          END DO
        END IF
      END IF
      RETURN
      END SUBROUTINE bc_r3d_bry_tile
! 
!***********************************************************************
      SUBROUTINE bc_u3d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij, LBk, UBk,                 &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk
      real(r8), intent(inout) :: A(LBij:,LBk:)
!
!  Local variable declarations.
!
      integer :: k
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
          DO k=LBk,UBk
            A(Jstr-1,k)=A(Jstr,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Jend+1,k)=A(Jend,k)
          END DO
        END IF
      END IF
      IF (boundary.eq.ieast) THEN
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Jstr-1,k)=A(Jstr,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Jend+1,k)=A(Jend,k)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((Jstr.eq.1).and.(Istr.eq.1)) THEN
          DO k=LBk,UBk
            A(IstrU-1,k)=A(IstrU,k)
          END DO
        END IF
        IF ((Jstr.eq.1).and.(Iend.eq.Lm(ng))) THEN
          DO k=LBk,UBk
            A(Iend+1,k)=A(Iend,k)
          END DO
        END IF
      END IF
      IF (boundary.eq.inorth) THEN
        IF ((Jend.eq.Mm(ng)).and.(Istr.eq.1)) THEN
          DO k=LBk,UBk
            A(IstrU-1,k)=A(IstrU,k)
          END DO
        END IF
        IF ((Jend.eq.Mm(ng)).and.(Iend.eq.Lm(ng))) THEN
          DO k=LBk,UBk
            A(Iend+1,k)=A(Iend,k)
          END DO
        END IF
      END IF
      RETURN
      END SUBROUTINE bc_u3d_bry_tile
! 
!***********************************************************************
      SUBROUTINE bc_v3d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij, LBk, UBk,                 &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk
      real(r8), intent(inout) :: A(LBij:,LBk:)
!
!  Local variable declarations.
!
      integer :: k
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
          DO k=LBk,UBk
            A(JstrV-1,k)=A(JstrV,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Jend+1,k)=A(Jend,k)
          END DO
        END IF
      END IF
      IF (boundary.eq.ieast) THEN
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(JstrV-1,k)=A(JstrV,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Jend+1,k)=A(Jend,k)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((Jstr.eq.1).and.(Istr.eq.1)) THEN
          DO k=LBk,UBk
            A(Istr-1,k)=A(Istr,k)
          END DO
        END IF
        IF ((Jstr.eq.1).and.(Iend.eq.Lm(ng))) THEN
          DO k=LBk,UBk
            A(Iend+1,k)=A(Iend,k)
          END DO
        END IF
      END IF
      IF (boundary.eq.inorth) THEN
        IF ((Jend.eq.Mm(ng)).and.(Istr.eq.1)) THEN
          DO k=LBk,UBk
            A(Istr-1,k)=A(Istr,k)
          END DO
        END IF
        IF ((Jend.eq.Mm(ng)).and.(Iend.eq.Lm(ng))) THEN
          DO k=LBk,UBk
            A(Iend+1,k)=A(Iend,k)
          END DO
        END IF
      END IF
      RETURN
      END SUBROUTINE bc_v3d_bry_tile
      END MODULE bc_bry3d_mod
