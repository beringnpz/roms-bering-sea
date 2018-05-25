      SUBROUTINE ana_tclima (ng, tile, model)
!
!! svn $Id: ana_tclima.h 223 2007-04-27 04:10:48Z arango $
!!======================================================================
!! Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets analytical tracer climatology fields.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tclima_tile (ng, model, tile,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      CLIMA(ng) % tclm)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(33)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tclima
!
!***********************************************************************
      SUBROUTINE ana_tclima_tile (ng, model, tile,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            tclm)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
#if defined BEST_NPZ && defined IRON_LIMIT
      USE mod_biology
#endif
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tclm(LBi:,LBj:,:,:)
#else
      real(r8), intent(out) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: i, itrc, j, k
#if defined BEST_NPZ && defined IRON_LIMIT
      real(r8) val3,val1,val2
#endif
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#if defined BEST_NPZ && defined IRON_LIMIT

! Iron - linear from surface value to value at 100m and increase onshore
      DO i=IstrR,IendR
        DO j=JstrR,JendR
          val3 = MAX(0._r8,MIN(1._r8,(1._r8*GRID(ng)%h(i,j)-Feinh)/(Feoffh-Feinh)))
          val1 = Feinlo + val3*(Feofflo-Feinlo)
          val2 = Feinhi + val3*(Feoffhi-Feinhi)
          val3 = (val2-val1) / 300._r8
          DO k=1,N(ng)
            if(GRID(ng)%h(i,j).gt.Feoffh.and.GRID(ng)%z_r(i,j,k).ge.-50_r8)THEN
              tclm(i,j,k,iFe) =  val1  
            else 
              tclm(i,j,k,iFe) = MIN(val2, val1 - GRID(ng)%z_r(i,j,k)*val3 )
            endif
                 
            !ajh nitrate
            tclm(i,j,k,iNO3) = MIN(30._r8, 18._r8 - GRID(ng)%z_r(i,j,k)*(30._r8-18._r8)/200._r8)
            !ajh nh4
            tclm(i,j,k,iNH4) = 1.d-16
          END DO
        END DO
      END DO


#else
      print*,'ana_tclima_user.h: No values provided for tclm.'
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      DO itrc=1,NAT
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tclm(:,:,:,itrc))
      END DO
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NAT,         &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tclm)
#endif

      RETURN
      END SUBROUTINE ana_tclima_tile
