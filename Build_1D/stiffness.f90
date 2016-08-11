      MODULE stiffness_mod
!
!svn $Id: stiffness.F 1038 2009-08-11 22:29:40Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine surveys the 3D grid in order to determine maximum      !
!  grid stiffness ratio:                                               !
!                                                                      !
!             z(i,j,k)-z(i-1,j,k)+z(i,j,k-1)-z(i-1,j,k-1)              !
!      r_x = ---------------------------------------------             !
!             z(i,j,k)+z(i-1,j,k)-z(i,j,k-1)-z(i-1,j,k-1)              !
!                                                                      !
!  This is done for diagnostic purposes and it does not affect the     !
!  computations.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: stiffness
      CONTAINS
!
!***********************************************************************
      SUBROUTINE stiffness (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private storage
!  arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J)- and MAX(I,J)-directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL stiffness_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % omn,                              &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_w,                              &
     &                     OCEAN(ng)% zeta)
      RETURN
      END SUBROUTINE stiffness
!
!***********************************************************************
      SUBROUTINE stiffness_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           h, omn,                                &
     &                           Hz, z_w,                               &
     &                           zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      real(r8) :: cff, ratio
      real(r8) :: my_rx0 = 0.0_r8
      real(r8) :: my_rx1 = 0.0_r8
      real(r8) :: my_volume0 = 0.0_r8
      real(r8) :: my_volume1 = 1.0E+20_r8
      real(r8) :: my_volume2 = 0.0_r8
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
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i-1,j,0))/            &
     &                            (z_w(i,j,0)+z_w(i-1,j,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i-1,j,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i-1,j,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i-1,j,k-1))))
            END DO
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
            my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i,j-1,0))/            &
     &                            (z_w(i,j,0)+z_w(i,j-1,0))))
            DO k=1,N(ng)
              my_rx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i,j-1,k  )+       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))/      &
     &                              (z_w(i,j,k  )+z_w(i,j-1,k  )-       &
     &                               z_w(i,j,k-1)-z_w(i,j-1,k-1))))
            END DO
        END DO
      END DO
!
!-------------------------------------------------------------------------
!  Compute initial basin volume and grid cell minimum and maximum volumes.
!-------------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
              cff=omn(i,j)*Hz(i,j,k)
              my_volume0=my_volume0+cff
              my_volume1=MIN(my_volume1,cff)
              my_volume2=MAX(my_volume2,cff)
          END DO
        END DO
      END DO
!
!-------------------------------------------------------------------------
!  Compute global values.
!-------------------------------------------------------------------------
!
      IF ((Istr.eq.1).and.(Jstr.eq.1).and.                              &
     &    (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
      IF (tile_count.eq.0) THEN
        rx0=my_rx0
        rx1=my_rx1
        TotVolume=my_volume0
        MinVolume=my_volume1
        MaxVolume=my_volume2
      ELSE
        rx0=MAX(rx0,my_rx0)
        rx1=MAX(rx1,my_rx1)
        TotVolume=my_volume0
        MinVolume=MIN(MinVolume,my_volume1)
        MaxVolume=MAX(MaxVolume,my_volume2)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
        IF (Master) THEN
          WRITE (stdout,10) rx0, rx1
  10      FORMAT (/,' Maximum grid stiffness ratios:  rx0 = ',1pe14.6,  &
     &              ' (Beckmann and Haidvogel)',/,t34,'rx1 = ',1pe14.6, &
     &              ' (Haney)',/)
          IF (MinVolume.ne.0.0_r8) THEN
            ratio=MaxVolume/MinVolume
          ELSE
            ratio=0.0_r8
          END IF
          WRITE (stdout,20) TotVolume, MinVolume, MaxVolume, ratio
  20      FORMAT (/,' Initial basin volumes: TotVolume = ',1p,e17.10,0p,&
     &            ' m3',/,t25,'MinVolume = ',1p,e17.10,0p,' m3',        &
     &            /,t25,'MaxVolume = ',1p,e17.10,0p,' m3',              &
     &            /,t25,'  Max/Min = ',1p,e17.10,0p,/)
        END IF
      END IF
      RETURN
      END SUBROUTINE stiffness_tile
      END MODULE stiffness_mod
