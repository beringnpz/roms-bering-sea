      MODULE state_dotprod_mod
!
!svn $Id: state_dotprod.F 999 2009-06-09 23:48:31Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the dot product between two model states:     !
!                                                                      !
!      DotProd(0:NstateVars) = < s1, s2 >                              !
!                                                                      !
!  where                                                               !
!                                                                      !
!      DotProd(0)           All state variable dot product             !
!      DotProd(isUvel)      3D U-momentum contribution                 !
!      DotProd(isVvel)      3D V-momentum contribution                 !
!      DotProd(isTvar(:))   Tracer-type variables contribution         !
!      DotProd(isFsur)      Free-surface contribution                  !
!                                                                      !
!=======================================================================
!
      implicit none
      PUBLIC  :: state_dotprod
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_dotprod (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          NstateVars, DotProd,                    &
     &                          s1_t, s2_t,                             &
     &                          s1_u, s2_u,                             &
     &                          s1_v, s2_v,                             &
     &                          s1_zeta, s2_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: NstateVars
!
      real(r8), intent(in) :: s1_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s1_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)
!
      real(r8), intent(out), dimension(0:NstateVars) :: DotProd
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      integer :: ir, it
      real(r8) :: cff
      real(r8), dimension(0:NstateVars) :: my_DotProd
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
!  Compute dot product between S1 and S2 model state trajectories.
!-----------------------------------------------------------------------
!
      DO i=0,NstateVars
        my_DotProd(i)=0.0_r8
      END DO
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=s1_zeta(i,j)*s2_zeta(i,j)
          my_DotProd(0)=my_DotProd(0)+cff
          my_DotProd(isFsur)=my_DotProd(isFsur)+cff
        END DO
      END DO
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=s1_u(i,j,k)*s2_u(i,j,k)
            my_DotProd(0)=my_DotProd(0)+cff
            my_DotProd(isUvel)=my_DotProd(isUvel)+cff
          END DO
        END DO
      END DO
!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff=s1_v(i,j,k)*s2_v(i,j,k)
            my_DotProd(0)=my_DotProd(0)+cff
            my_DotProd(isVvel)=my_DotProd(isVvel)+cff
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              cff=s1_t(i,j,k,it)*s2_t(i,j,k,it)
              my_DotProd(0)=my_DotProd(0)+cff
              my_DotProd(isTvar(it))=my_DotProd(isTvar(it))+cff
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform parallel global reduction operations.
!-----------------------------------------------------------------------
!
      IF ((Istr.eq.1).and.(Jstr.eq.1).and.                              &
     &    (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVars
          DotProd(i)=0.0_r8
        END DO
      END IF
      DO i=0,NstateVars
        DotProd(i)=DotProd(i)+my_DotProd(i)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
      END IF
      RETURN
      END SUBROUTINE state_dotprod
      END MODULE state_dotprod_mod
