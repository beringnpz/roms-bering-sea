      MODULE exchange_3d_mod
!
!svn $Id: exchange_3d.F 895 2009-01-12 21:06:20Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This package contains periodic boundary conditions and parallel     !
!  exchage (distributed-memory only) routines for 3D variables.        !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p3d_tile    periodic conditions/exchange at PSI-points   !
!    exchange_r3d_tile    periodic conditions/exchange at RHO-points   !
!    exchange_u3d_tile    periodic conditions/exchange at U-points     !
!    exchange_v3d_tile    periodic conditions/exchange at V-points     !
!    exchange_w3d_tile    periodic conditions/exchange at W-points     !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
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
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_p3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
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
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_r3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
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
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_u3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
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
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_v3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
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
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(Lm(ng)+1,j,k)=A(1,j,k)
              A(Lm(ng)+2,j,k)=A(2,j,k)
            END DO
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              A(-2,j,k)=A(Lm(ng)-2,j,k)
              A(-1,j,k)=A(Lm(ng)-1,j,k)
              A( 0,j,k)=A(Lm(ng)  ,j,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,Mm(ng)+1,k)=A(i,1,k)
              A(i,Mm(ng)+2,k)=A(i,2,k)
            END DO
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              A(i,-2,k)=A(i,Mm(ng)-2,k)
              A(i,-1,k)=A(i,Mm(ng)-1,k)
              A(i, 0,k)=A(i,Mm(ng)  ,k)
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
            A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
            A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
            A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          DO k=LBk,UBk
            A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
            A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
            A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
            A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
            A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
            A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
          END DO
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
            A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
            A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
            A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
            A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
            A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
          END DO
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          DO k=LBk,UBk
            A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
            A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
            A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
            A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
            A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
            A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
            A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
            A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
            A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_w3d_tile
      END MODULE exchange_3d_mod
