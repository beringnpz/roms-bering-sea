      MODULE state_addition_mod
!
!svn $Id: state_addition.F 999 2009-06-09 23:48:31Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the following model state addition:           !
!                                                                      !
!      s1_var(...,Lout) = fac1 * s1_var(...,Lin1) +                    !
!                         fac2 * s2_var(...,Lin2)                      !
!                                                                      !
!  where fac1 and fac2 are scalars.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
      PUBLIC  :: state_addition
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_addition (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           Lin1, Lin2, Lout,                      &
     &                           fac1, fac2,                            &
     &                           s1_t, s2_t,                            &
     &                           s1_u, s2_u,                            &
     &                           s1_v, s2_v,                            &
     &                           s1_zeta, s2_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lin1, Lin2, Lout
!
      real(r8), intent(in) :: fac1, fac2
!
      real(r8), intent(in) :: s2_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: s2_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: s1_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s1_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      integer :: ib, ir, it
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
!  Compute the following operation between S1 and S2 model state
!  trajectories:
!                 S1(Lout) = fac1 * S1(Lin1) + fac2 * S2(Lin2)
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          s1_zeta(i,j,Lout)=fac1*s1_zeta(i,j,Lin1)+                     &
     &                      fac2*s2_zeta(i,j,Lin2)
        END DO
      END DO
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s1_u(i,j,k,Lout)=fac1*s1_u(i,j,k,Lin1)+                     &
     &                       fac2*s2_u(i,j,k,Lin2)
          END DO
        END DO
      END DO
!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s1_v(i,j,k,Lout)=fac1*s1_v(i,j,k,Lin1)+                     &
     &                       fac2*s2_v(i,j,k,Lin2)
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
              s1_t(i,j,k,Lout,it)=fac1*s1_t(i,j,k,Lin1,it)+             &
     &                            fac2*s2_t(i,j,k,Lin2,it)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE state_addition
      END MODULE state_addition_mod
