      MODULE step3d_uv_mod
!
!svn $Id: step3d_uv.F 1060 2009-09-12 00:25:38Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This subroutine time-steps the nonlinear  horizontal  momentum      !
!  equations. The vertical viscosity terms are time-stepped using      !
!  an implicit algorithm.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: step3d_uv
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_uv (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_mixing
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
      CALL wclock_on (ng, iNLM, 34)
      CALL step3d_uv_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng), nstp(ng), nnew(ng),                &
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     MIXING(ng) % Akv,                            &
     &                     COUPLING(ng) % DU_avg1,                      &
     &                     COUPLING(ng) % DV_avg1,                      &
     &                     COUPLING(ng) % DU_avg2,                      &
     &                     COUPLING(ng) % DV_avg2,                      &
     &                     OCEAN(ng) % ru,                              &
     &                     OCEAN(ng) % rv,                              &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % ubar,                            &
     &                     OCEAN(ng) % vbar,                            &
     &                     GRID(ng) % Huon,                             &
     &                     GRID(ng) % Hvom)
      CALL wclock_off (ng, iNLM, 34)
      RETURN
      END SUBROUTINE step3d_uv
!
!***********************************************************************
      SUBROUTINE step3d_uv_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs, nstp, nnew,                      &
     &                           om_v, on_u, pm, pn,                    &
     &                           Hz, z_r, z_w,                          &
     &                           Akv,                                   &
     &                           DU_avg1, DV_avg1,                      &
     &                           DU_avg2, DV_avg2,                      &
     &                           ru, rv,                                &
     &                           u, v,                                  &
     &                           ubar, vbar,                            &
     &                           Huon, Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: rv(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: Hvom(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, idiag, is, j, k
      real(r8) :: cff, cff1, cff2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: AK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: oHz
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
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          AK(i,0)=0.5_r8*(Akv(i-1,j,0)+                                 &
     &                    Akv(i  ,j,0))
          DO k=1,N(ng)
            AK(i,k)=0.5_r8*(Akv(i-1,j,k)+                               &
     &                      Akv(i  ,j,k))
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
          END DO
        END DO
!
!  Time step right-hand-side terms.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.25_r8*dt(ng)
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
        ELSE
          cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
        END IF
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+                                &
     &                    DC(i,0)*ru(i,j,k,nrhs)
          END DO
        END DO
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
!
        cff=-lambda*dt(ng)/0.5_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                   &
     &                   z_r(i,j,k  )-z_r(i-1,j,k  ))
            FC(i,k)=cff*cff1*AK(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
        END DO
!
!  Solve the tridiagonal system.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
          CF(i,1)=cff*FC(i,1)
          DC(i,1)=cff*DC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FC(i,k)
            DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
          END DO
        END DO
!
!  Compute new solution by back substitution.
!
        DO i=IstrU,Iend
          DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/        &
     &                (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
          u(i,j,N(ng),nnew)=DC(i,N(ng))
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            u(i,j,k,nnew)=DC(i,k)
          END DO
        END DO
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          cff1=1.0_r8/(CF(i,0)*on_u(i,j))
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j))*cff1      ! recursive
        END DO
!
!  Couple and update new solution.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew)-DC(i,0)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            AK(i,0)=0.5_r8*(Akv(i,j-1,0)+                               &
     &                      Akv(i,j  ,0))
            DO k=1,N(ng)
              AK(i,k)=0.5_r8*(Akv(i,j-1,k)+                             &
     &                        Akv(i,j  ,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
            END DO
          END DO
!
!  Time step right-hand-side terms.
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.25_r8*dt(ng)
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
          ELSE
            cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
          END IF
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs)
            END DO
          END DO
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
          cff=-lambda*dt(ng)/0.5_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                 &
     &                     z_r(i,j,k  )-z_r(i,j-1,k  ))
              FC(i,k)=cff*cff1*AK(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Solve the tridiagonal system.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
     &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
            v(i,j,N(ng),nnew)=DC(i,N(ng))
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              v(i,j,k,nnew)=DC(i,k)
            END DO
          END DO
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
!
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            cff1=1.0_r8/(CF(i,0)*om_v(i,j))
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j))*cff1    ! recursive
          END DO
!
!  Couple and update new solution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)-DC(i,0)
            END DO
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
! Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp, nnew,                                      &
     &                 u)
      CALL v3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp, nnew,                                      &
     &                 v)
!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!
!
!  Couple velocity component in the XI-direction.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
          FC(i,0)=0.0_r8
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.5_r8*on_u(i,j)
            DC(i,k)=cff*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+                                            &
     &              DC(i,k)*u(i,j,k,nnew)
          END DO
        END DO
        DO i=Istr,Iend
          cff1=DC(i,0)                                  ! intermediate
          DC(i,0)=1.0_r8/DC(i,0)                        ! recursive
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))        ! recursive
          ubar(i,j,1)=DC(i,0)*DU_avg1(i,j)
          ubar(i,j,2)=ubar(i,j,1)
        END DO
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice.
!
!
!  Compute correct mass flux, Hz*u/n.
!
        DO k=N(ng),1,-1
          DO i=Istr,Iend
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
            FC(i,0)=FC(i,0)+Huon(i,j,k)
          END DO
        END DO
        DO i=Istr,Iend
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))        ! recursive
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=Istr,Iend
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
            FC(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=0.5_r8*om_v(i,j)
              DC(i,k)=cff*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+                                          &
     &                DC(i,k)*v(i,j,k,nnew)
            END DO
          END DO
          DO i=Istr,Iend
            cff1=DC(i,0)                                 ! Intermediate
            DC(i,0)=1.0_r8/DC(i,0)                       ! recursive
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j))       ! recursive
            vbar(i,j,1)=DC(i,0)*DV_avg1(i,j)
            vbar(i,j,2)=vbar(i,j,1)
          END DO
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice.
!  
!
!  Compute correct mass flux, Hz*v/m.
!
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))      ! recursive
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nnew))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nnew))
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Huon)
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Hvom)
!
      DO k=1,2
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,k))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,k))
      END DO
      RETURN
      END SUBROUTINE step3d_uv_tile
      END MODULE step3d_uv_mod
