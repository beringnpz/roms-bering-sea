#include "cppdefs.h"
      MODULE ice_limit_mod
#if defined ICE_ADVECT || defined ICE_THERMO
!
!============================================== W. Paul Budgell =======!
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!============================================== Hernan G. Arango ======!
!                                                                      !
!  Limit ice tracer variables, for instance to not allow negative      !
!  ice concentration or thickness.                                     !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC ice_limit

      CONTAINS

      SUBROUTINE ice_limit (ng, tile)

      USE mod_param
      USE mod_ice
      USE mod_stepping

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!


      CALL ice_limit_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     liold(ng), linew(ng),                        &
     &                     ICE(ng) % ui,                                &
     &                     ICE(ng) % vi,                                &
     &                     ICE(ng) % ai,                                &
     &                     ICE(ng) % hi,                                &
     &                     ICE(ng) % hsn,                               &
     &                     ICE(ng) % ti,                                &
     &                     ICE(ng) % enthalpi,                          &
     &                     ICE(ng) % sfwat                              &
# if defined ICE_BIO
# ifdef BERING_10K
     &                     ,ICE(ng) % IcePhL                            &
     &                     ,ICE(ng) % IceNO3                            &
     &                     ,ICE(ng) % IceNH4                            &
     &                     ,ICE(ng) % IceLog                            &
# endif
# endif
     &                     )
!
      RETURN
      END SUBROUTINE ice_limit
!
!***********************************************************************
      SUBROUTINE ice_limit_tile (ng, tile,                              &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, linew,                             &
     &                        ui, vi, ai, hi, hsn, ti, enthalpi, sfwat  &
# if defined ICE_BIO
#  ifdef BERING_10K
     &                        ,IcePhL,IceNO3,IceNH4,IceLog              &
#  endif
# endif

     &                        )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE aibc_mod, ONLY : aibc_tile
      USE hibc_mod, ONLY : hibc_tile
# if defined ICE_BIO
#  ifdef BERING_10K

      USE IcePhLbc_mod, ONLY : IcePhLbc_tile
      USE IceNO3bc_mod, ONLY : IceNO3bc_tile
      USE IceNH4bc_mod, ONLY : IceNH4bc_tile
      USE mod_biology, ONLY: aidz
#  endif
# endif
      USE hsnbc_mod, ONLY : hsnbc_tile

      USE tibc_mod, ONLY : tibc_tile
      USE sfwatbc_mod, ONLY : sfwatbc_tile
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(inout) :: liold, linew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
#  if defined ICE_BIO
#   ifdef BERING_10K
      real(r8), intent(inout) :: IcePhL(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNO3(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNH4(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceLog(LBi:,LBj:,:)
#   endif
#  endif
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
      real(r8), intent(inout) :: sfwat(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
#  if defined ICE_BIO
#   ifdef BERING_10K
      real(r8), intent(inout) :: IcePhL(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNO3(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNH4(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceLog(LBi:UBi,LBj:UBj,2)
#   endif
#  endif
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sfwat(LBi:UBi,LBj:UBj,2)
# endif

! Local variable definitions
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j

# include "set_bounds.h"
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ai(i,j,linew) = MIN(ai(i,j,linew),max_a(ng))
          ai(i,j,linew) = MAX(ai(i,j,linew),0.0_r8)
          hi(i,j,linew) = MAX(hi(i,j,linew),0.0_r8)
          hsn(i,j,linew) = MAX(hsn(i,j,linew),0.0_r8)

          sfwat(i,j,linew) = MAX(sfwat(i,j,linew),0.0_r8)
          ti(i,j,linew) = MAX(ti(i,j,linew),-70.0_r8)
          IF (hi(i,j,linew) .le. min_h(ng)) ai(i,j,linew) = min_a(ng)
          IF (ai(i,j,linew) .le. min_a(ng)) THEN
            ai(i,j,linew) = min_a(ng)
            hi(i,j,linew) = min_h(ng)
            sfwat(i,j,linew) = 0.0_r8
            hsn(i,j,linew) = 0.0_r8

          END IF
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
# if defined ICE_BIO
#  ifdef BERING_10K

          IcePhL(i,j,linew) = MAX(IcePhL(i,j,linew),0.0_r8)
          IceNO3(i,j,linew) = MAX(IceNO3(i,j,linew),0.0_r8)
          IceNH4(i,j,linew) = MAX(IceNH4(i,j,linew),0.0_r8)

          if (hi(i,j,linew).le.aidz .or. ai(i,j,linew).le.0.5) THEN
            IceLog(i,j,linew) = -1.0_r8
          else
            IceLog(i,j,linew) = 1.0_r8 !min(IceLog(i,j,linew),1.0_r8)
          endif

#  endif
# endif
        END DO
      END DO
!
      CALL aibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, ai)
      CALL hibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, hi)
      CALL hsnbc_tile (ng, tile,                                        &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, hsn)
# if defined ICE_BIO
#  ifdef BERING_10K

      CALL IcePhLbc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, IcePhL)
      CALL IceNO3bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, IceNO3)
      CALL IceNH4bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, IceNH4)

#  endif
# endif

      CALL tibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew, min_h(ng),                          &
     &                ui, vi, hi, ti, enthalpi)

      CALL sfwatbc_tile (ng, tile,                                      &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, sfwat)

# if defined EW_PERIODIC || defined NS_PERIODIC

      CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ai(:,:,linew))

      CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hi(:,:,linew))

      CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hsn(:,:,linew))
#  if defined ICE_BIO
#   ifdef BERING_10K


     CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                    &
     &                          IcePhL(:,:,linew))
     CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                    &
     &                          IceNO3(:,:,linew))
     CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                    &
     &                          IceNH4(:,:,linew))

     CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                    &
     &                          IceLog(:,:,linew))

#   endif
#  endif

      CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ti(:,:,linew))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sfwat(:,:,linew))
# endif

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ai(:,:,linew), hi(:,:,linew),                 &
     &                    hsn(:,:,linew),ti(:,:,linew))

      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sfwat(:,:,linew))

#  if defined ICE_BIO
#   ifdef BERING_10K

      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    IcePhL(:,:,linew), IceNO3(:,:,linew),         &
     &                    IceNH4(:,:,linew),IceLog(:,:,linew))

#   endif
#  endif
# endif
      RETURN
      END SUBROUTINE ice_limit_tile
#endif
      END MODULE ice_limit_mod
