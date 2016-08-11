      MODULE set_zeta_mod
!
!svn $Id: set_zeta.F 1027 2009-07-23 15:50:57Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface to its fast-time averaged value.     !
!                                                                      !
!=======================================================================
!
!
      implicit none
      PRIVATE
      PUBLIC  :: set_zeta
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_zeta (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
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
      CALL wclock_on (ng, iNLM, 12)
      CALL set_zeta_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    COUPLING(ng) % Zt_avg1,                       &
     &                    OCEAN(ng) % zeta)
      CALL wclock_off (ng, iNLM, 12)
      RETURN
      END SUBROUTINE set_zeta
!
!***********************************************************************
      SUBROUTINE set_zeta_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Zt_avg1,                                &
     &                          zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      real(r8) :: cff
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
!
!-----------------------------------------------------------------------
!  Prepare to time-step 2D equations:  set initial free-surface
!  to its fast-time averaged values (which corresponds to the time
!  step "n").
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,1)=Zt_avg1(i,j)
          zeta(i,j,2)=Zt_avg1(i,j)
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,1))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,2))
      RETURN
      END SUBROUTINE set_zeta_tile
      END MODULE set_zeta_mod
