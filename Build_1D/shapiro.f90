      MODULE shapiro_mod
!
!svn $Id: shapiro.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group         Kate Hedstrom   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package contains shapiro filter routines for order 2 and       !
!  reduced order at the boundary and mask edges.                       !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    shapirp2d_tile       Shapiro filter for 2D fields.                !
!    shapirp3d_tile       Shapiro filter for 3D fields.                !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE shapiro2d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk2
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
!  Shapiro filter requested 2D field.
!-----------------------------------------------------------------------
!
!  This subroutine will apply a Shapiro filter of order 2 (defined
!  as twice the order in Shapiro (1970), with N even) to an array, A.
!  The order of the filter is reduced at the boundaries and at the
!  mask edges, if any.
!
!  Initialize filter in the Y-direction.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
          Awrk1(i,j)=0.25_r8*                                           &
     &             (A(i,j-1)+A(i,j+1)-2.0_r8*A(i,j))
        END DO
      END DO
!
!  Add the changes to the field.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
          Awrk2(i,j)=A(i,j)+Awrk1(i,j)
        END DO
      END DO
!
!  Initialize filter in the X-direction.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Awrk1(i,j)=0.25_r8*                                           &
     &               (Awrk2(i-1,j)+Awrk2(i+1,j)-2.0_r8*Awrk2(i,j))
        END DO
      END DO
!
!  Add changes to field.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          A(i,j)=Awrk2(i,j)+Awrk1(i,j)
        END DO
      END DO
      RETURN
      END SUBROUTINE shapiro2d_tile
!
!***********************************************************************
      SUBROUTINE shapiro3d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk2
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
!  Shapiro filter requested 3D field.
!-----------------------------------------------------------------------
!
!  This subroutine will apply a Shapiro filter of order 2 (defined
!  as twice the order in Shapiro (1970), with N even) to an array, A.
!  The order of the filter is reduced at the boundaries and at the
!  mask edges, if any.
!
!  Initialize filter in the Y-direction.
!
      DO k=LBk,UBk
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (A(i,j-1,k)+A(i,j+1,k)-2.0_r8*A(i,j,k))
          END DO
        END DO
!
!  Add the changes to the field.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            Awrk2(i,j)=A(i,j,k)+Awrk1(i,j)
          END DO
        END DO
!
!  Initialize filter in the X-direction.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (Awrk2(i-1,j)+Awrk2(i+1,j)-2.0_r8*Awrk2(i,j))
          END DO
        END DO
!
!  Add changes to field.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            A(i,j,k)=Awrk2(i,j)+Awrk1(i,j)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE shapiro3d_tile
      END MODULE shapiro_mod
