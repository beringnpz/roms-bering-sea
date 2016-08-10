      MODULE state_scale_mod
!
!svn $Id: state_scale.F 999 2009-06-09 23:48:31Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine scales a state variable by a constant:                 !
!                                                                      !
!     s_var(...,Lout) = fac * s_var(...,Linp)                          !
!                                                                      !
!=======================================================================
!
      implicit none
      PUBLIC  :: state_scale
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_scale (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        Linp, Lout, fac,                          &
     &                        s_t, s_u, s_v,                            &
     &                        s_zeta)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Linp, Lout
!
      real(r8), intent(in) :: fac
!
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
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
!  Scale model state variable by a constant.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          s_zeta(i,j,Lout)=fac*s_zeta(i,j,Linp)
        END DO
      END DO
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s_u(i,j,k,Lout)=fac*s_u(i,j,k,Linp)
          END DO
        END DO
      END DO
!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s_v(i,j,k,Lout)=fac*s_v(i,j,k,Linp)
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
              s_t(i,j,k,Lout,it)=fac*s_t(i,j,k,Linp,it)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE state_scale
      END MODULE state_scale_mod
