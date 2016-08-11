      MODULE exchange_2d_mod
!
!svn $Id: exchange_2d.F 895 2009-01-12 21:06:20Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This package contains periodic boundary conditions routines for 2D  !
!  variables.                                                          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p2d_tile    periodic conditions/exchange at PSI-points   !
!    exchange_r2d_tile    periodic conditions/exchange at RHO-points   !
!    exchange_u2d_tile    periodic conditions/exchange at U-points     !
!    exchange_v2d_tile    periodic conditions/exchange at V-points     !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
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
          DO j=Jstr,Jend
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=Jstr,Jend
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO i=Istr,Iend
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO i=Istr,Iend
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      RETURN
      END SUBROUTINE exchange_p2d_tile
!
!***********************************************************************
      SUBROUTINE exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
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
          DO j=Jstr,Jend
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=Jstr,Jend
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO i=Istr,Iend
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO i=Istr,Iend
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      RETURN
      END SUBROUTINE exchange_r2d_tile
!
!***********************************************************************
      SUBROUTINE exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
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
          DO j=Jstr,Jend
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=Jstr,Jend
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO i=Istr,Iend
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO i=Istr,Iend
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      RETURN
      END SUBROUTINE exchange_u2d_tile
!
!***********************************************************************
      SUBROUTINE exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
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
          DO j=Jstr,Jend
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=Jstr,Jend
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Jstr.eq.1) THEN
          DO i=Istr,Iend
            A(i,Mm(ng)+1)=A(i,1)
            A(i,Mm(ng)+2)=A(i,2)
          END DO
        END IF
        IF (Jend.eq.Mm(ng)) THEN
          DO i=Istr,Iend
            A(i,-2)=A(i,Mm(ng)-2)
            A(i,-1)=A(i,Mm(ng)-1)
            A(i, 0)=A(i,Mm(ng)  )
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((Istr.eq.1).and.(Jstr.eq.1)) THEN
          A(Lm(ng)+1,Mm(ng)+1)=A(1,1)
          A(Lm(ng)+1,Mm(ng)+2)=A(1,2)
          A(Lm(ng)+2,Mm(ng)+1)=A(2,1)
          A(Lm(ng)+2,Mm(ng)+2)=A(2,2)
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jstr.eq.1)) THEN
          A(-2,Mm(ng)+1)=A(Lm(ng)-2,1)
          A(-1,Mm(ng)+1)=A(Lm(ng)-1,1)
          A( 0,Mm(ng)+1)=A(Lm(ng)  ,1)
          A(-2,Mm(ng)+2)=A(Lm(ng)-2,2)
          A(-1,Mm(ng)+2)=A(Lm(ng)-1,2)
          A( 0,Mm(ng)+2)=A(Lm(ng)  ,2)
        END IF
        IF ((Istr.eq.1).and.(Jend.eq.Mm(ng))) THEN
          A(Lm(ng)+1,-2)=A(1,Mm(ng)-2)
          A(Lm(ng)+1,-1)=A(1,Mm(ng)-1)
          A(Lm(ng)+1, 0)=A(1,Mm(ng)  )
          A(Lm(ng)+2,-2)=A(2,Mm(ng)-2)
          A(Lm(ng)+2,-1)=A(2,Mm(ng)-1)
          A(Lm(ng)+2, 0)=A(2,Mm(ng)  )
        END IF
        IF ((Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          A(-2,-2)=A(Lm(ng)-2,Mm(ng)-2)
          A(-2,-1)=A(Lm(ng)-2,Mm(ng)-1)
          A(-2, 0)=A(Lm(ng)-2,Mm(ng)  )
          A(-1,-2)=A(Lm(ng)-1,Mm(ng)-2)
          A(-1,-1)=A(Lm(ng)-1,Mm(ng)-1)
          A(-1, 0)=A(Lm(ng)-1,Mm(ng)  )
          A( 0,-2)=A(Lm(ng)  ,Mm(ng)-2)
          A( 0,-1)=A(Lm(ng)  ,Mm(ng)-1)
          A( 0, 0)=A(Lm(ng)  ,Mm(ng)  )
        END IF
      RETURN
      END SUBROUTINE exchange_v2d_tile
      END MODULE exchange_2d_mod
