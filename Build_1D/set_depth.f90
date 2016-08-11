      MODULE set_depth_mod
!
!svn $Id: set_depth.F 1012 2009-07-07 20:52:45Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!  Currently, two vertical coordinate transformations are available    !
!  with various possible vertical stretching, C(s), functions, (see    !
!  routine "set_scoord.F" for details).                                !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: set_depth, set_depth_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
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
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp(ng), nnew(ng),                          &
     &                     GRID(ng) % h,                                &
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
      RETURN
      END SUBROUTINE set_depth
!
!***********************************************************************
      SUBROUTINE set_depth_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp, nnew,                            &
     &                           h,                                     &
     &                           Zt_avg1,                               &
     &                           Hz, z_r, z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: C2_r, C2_w, hh2, vert_n1, vert_a, vert_h0, vert_s0
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
!  Original formulation: Compute vertical depths (meters, negative) at
!                        RHO- and W-points, and vertical grid
!  thicknesses. Various stretching functions are possible.
!
!         z_w(x,y,s,t) = Zo_w + zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
!
!                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!-----------------------------------------------------------------------
!
      IF (Vtransform(ng).eq.3) THEN
        vert_s0 = 90./120.
        vert_n1 = 2./3. !3./5.
        vert_h0 = 60./(vert_s0)**(vert_n1)
        vert_a = 1./(vert_s0-1)
      END IF
      IF (Vtransform(ng).eq.1) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N(ng)
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrR,IendR
              hwater=h(i,j)
              hinv=1.0_r8/hwater
              z_w0=cff_w+cff1_w*hwater
              z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
              z_r0=cff_r+cff1_r*hwater
              z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
              Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
            END DO
          END DO
        END DO
      ELSE IF (Vtransform(ng).eq.3) THEN
        DO j=JstrR,JendR
          DO k=1,N(ng)
            DO i=IstrR,IendR
              hwater=h(i,j)
              hinv=1.0_r8/hwater
              hh2=(min(vert_h0,hwater))*hinv
              IF (SCALARS(ng)%sc_w(k).gt.-vert_s0) THEN 
                C2_w=-hh2*(-SCALARS(ng)%sc_w(k))**(vert_n1) 
                C2_r=-hh2*(-SCALARS(ng)%sc_r(k))**(vert_n1)
                cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-C2_w)
                cff1_w=C2_w
                cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-C2_r)
                cff1_r=C2_r
                z_w0=cff_w+cff1_w*hwater
                z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
                z_r0=cff_r+cff1_r*hwater
                z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
              ELSE
                C2_w=-hh2*(-SCALARS(ng)%sc_w(k))**(vert_n1)             &
     &               -(1-hh2)*(vert_a*(SCALARS(ng)%sc_w(k)+vert_s0))**2
                C2_r=-hh2*(-SCALARS(ng)%sc_r(k))**(vert_n1)             &
     &               -(1-hh2)*(vert_a*(SCALARS(ng)%sc_r(k)+vert_s0))**2
                cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-C2_w)
                cff1_w=C2_w
                cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-C2_r)
                cff1_r=C2_r
                z_w0=cff_w+cff1_w*hwater
                z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
                z_r0=cff_r+cff1_r*hwater
                z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
              END IF
              Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  New formulation: Compute vertical depths (meters, negative) at
!                   RHO- and W-points, and vertical grid thicknesses.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!-----------------------------------------------------------------------
!
      ELSE IF (Vtransform(ng).eq.2) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N(ng)
            cff_r=hc(ng)*SCALARS(ng)%sc_r(k)
            cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrR,IendR
              hwater=h(i,j)
              hinv=1.0_r8/(hc(ng)+hwater)
              cff2_r=(cff_r+cff1_r*hwater)*hinv
              cff2_w=(cff_w+cff1_w*hwater)*hinv
              z_w(i,j,k)=Zt_avg1(i,j)+(Zt_avg1(i,j)+hwater)*cff2_w
              z_r(i,j,k)=Zt_avg1(i,j)+(Zt_avg1(i,j)+hwater)*cff2_r
              Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        h)
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        z_w)
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        z_r)
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Hz)
      RETURN
      END SUBROUTINE set_depth_tile
      END MODULE set_depth_mod
