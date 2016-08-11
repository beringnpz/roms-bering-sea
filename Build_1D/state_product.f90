      MODULE state_product_mod
!
!svn $Id: state_product.F 1020 2009-07-10 23:10:30Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the product of each element of two model      !
!  states:                                                             !
!                                                                      !
!    s3(...) = s1(...) * s2(...)                                       !
!                                                                      !
!=======================================================================
!
      implicit none
      PUBLIC  :: state_product
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_product (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          s1_t, s2_t, s3_t,                       &
     &                          s1_u, s2_u, s3_u,                       &
     &                          s1_v, s2_v, s3_v,                       &
     &                          s1_zeta, s2_zeta, s3_zeta)
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
!
      real(r8), intent(in) :: s1_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s1_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: s3_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s3_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: s3_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)
      real(r8), intent(inout) :: s3_zeta(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      integer :: ir, it
      real(r8) :: cff
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
!  Compute product between S1 and S2 model state trajectories and save 
!  in S3.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff=s1_zeta(i,j)*s2_zeta(i,j)
          s3_zeta(i,j)=cff
        END DO
      END DO
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff=s1_u(i,j,k)*s2_u(i,j,k)
            s3_u(i,j,k)=cff
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
            s3_v(i,j,k)=cff
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
              s3_t(i,j,k,it)=cff
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE state_product
      END MODULE state_product_mod
