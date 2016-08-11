      MODULE rhs3d_mod
!
!svn $Id: rhs3d.F 895 2009-01-12 21:06:20Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine evaluates right-hand-side terms for 3D momentum     !
!  and tracers equations.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: rhs3d
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rhs3d (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE pre_step3d_mod, ONLY : pre_step3d
      USE prsgrd_mod, ONLY : prsgrd
      USE t3dmix_mod, ONLY : t3dmix2
      USE uv3dmix_mod, ONLY : uv3dmix2
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
!-----------------------------------------------------------------------
!  Initialize computations for new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
!
      CALL pre_step3d (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
!-----------------------------------------------------------------------
!
      CALL prsgrd (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
      CALL t3dmix2 (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for the 3D momentum equations.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 21)
      CALL rhs3d_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nrhs(ng),                                        &
     &                 GRID(ng) % Hz,                                   &
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
     &                 GRID(ng) % fomn,                                 &
     &                 GRID(ng) % om_u,                                 &
     &                 GRID(ng) % om_v,                                 &
     &                 GRID(ng) % on_u,                                 &
     &                 GRID(ng) % on_v,                                 &
     &                 GRID(ng) % pm,                                   &
     &                 GRID(ng) % pn,                                   &
     &                 FORCES(ng) % bustr,                              &
     &                 FORCES(ng) % bvstr,                              &
     &                 FORCES(ng) % sustr,                              &
     &                 FORCES(ng) % svstr,                              &
     &                 OCEAN(ng) % u,                                   &
     &                 OCEAN(ng) % v,                                   &
     &                 OCEAN(ng) % W,                                   &
     &                 COUPLING(ng) % rufrc,                            &
     &                 COUPLING(ng) % rvfrc,                            &
     &                 OCEAN(ng) % ru,                                  &
     &                 OCEAN(ng) % rv)
      CALL wclock_off (ng, iNLM, 21)
!
!-----------------------------------------------------------------------
!  Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
!
      CALL uv3dmix2 (ng, tile)
      RETURN
      END SUBROUTINE rhs3d
!
!***********************************************************************
      SUBROUTINE rhs3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs,                                      &
     &                       Hz, Huon, Hvom,                            &
     &                       fomn,                                      &
     &                       om_u, om_v, on_u, on_v, pm, pn,            &
     &                       bustr, bvstr,                              &
     &                       sustr, svstr,                              &
     &                       u, v, W,                                   &
     &                       rufrc, rvfrc,                              &
     &                       ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
      real(r8), intent(out) :: rufrc(LBi:,LBj:)
      real(r8), intent(out) :: rvfrc(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: Gadv = -0.25_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: fac, fac1, fac2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
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
      K_LOOP : DO k=1,N(ng)
!
!-----------------------------------------------------------------------
!  Add in Coriolis and curvilinear transformation terms, if any.
!-----------------------------------------------------------------------
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Hz(i,j,k)*(                                      &
     &          fomn(i,j)                                               &
     &          )
            UFx(i,j)=cff*(v(i,j  ,k,nrhs)+                              &
     &                    v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i  ,j,k,nrhs)+                              &
     &                    u(i+1,j,k,nrhs))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
!  of tensor of momentum flux due to horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend+1
            uxx(i,j)=u(i-1,j,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
     &               u(i+1,j,k,nrhs)
            Huxx(i,j)=Huon(i-1,j,k)-2.0_r8*Huon(i,j,k)+Huon(i+1,j,k)
          END DO
        END DO
!
!  Third-order, upstream bias u-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff1=u(i  ,j,k,nrhs)+                                       &
     &           u(i+1,j,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
            UFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Huon(i  ,j,k)+                                    &
     &                Huon(i+1,j,k)+                                    &
     &                Gadv*0.5_r8*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
        DO j=Jstr-1,Jend+1
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
     &               u(i,j+1,k,nrhs)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
           Hvxx(i,j)=Hvom(i-1,j,k)-2.0_r8*Hvom(i,j,k)+Hvom(i+1,j,k)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            cff1=u(i,j  ,k,nrhs)+                                       &
     &           u(i,j-1,k,nrhs)
            cff2=Hvom(i,j,k)+Hvom(i-1,j,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=uee(i,j-1)
            ELSE
              cff=uee(i,j)
            END IF
            UFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Hvxx(i  ,j)+                    &
     &                                  Hvxx(i-1,j)))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr-1,Iend+1
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
     &               v(i+1,j,k,nrhs)
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
           Huee(i,j)=Huon(i,j-1,k)-2.0_r8*Huon(i,j,k)+Huon(i,j+1,k)
          END DO
        END DO
!
!  Third-order, upstream bias v-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            cff1=v(i  ,j,k,nrhs)+                                       &
     &           v(i-1,j,k,nrhs)
            cff2=Huon(i,j,k)+Huon(i,j-1,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=vxx(i-1,j)
            ELSE
              cff=vxx(i,j)
            END IF
            VFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Huee(i,j  )+                    &
     &                                  Huee(i,j-1)))
          END DO
        END DO
        DO j=JstrV-1,Jend+1
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
     &               v(i,j+1,k,nrhs)
            Hvee(i,j)=Hvom(i,j-1,k)-2.0_r8*Hvom(i,j,k)+Hvom(i,j+1,k)
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1=v(i,j  ,k,nrhs)+                                       &
     &           v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF
            VFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Hvom(i,j  ,k)+                                    &
     &                Hvom(i,j+1,k)+                                    &
     &                Gadv*0.5_r8*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
!
!  Add in horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=UFx(i,j)-UFx(i-1,j)+                                    &
     &          UFe(i,j+1)-UFe(i,j)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=VFx(i+1,j)-VFx(i,j)+                                    &
     &          VFe(i,j)-VFe(i,j-1)
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
          END DO
        END DO
      END DO K_LOOP
!
      J_LOOP : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Add in vertical advection.
!-----------------------------------------------------------------------
!
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of u-momentum.
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=cff1*(Hz(i  ,j,k)+                                  &
     &                    Hz(i-1,j,k))-                                 &
     &              cff2*(Hz(i+1,j,k)+                                  &
     &                    Hz(i-2,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.0_r8*(u(i,j,k+1,nrhs)-                       &
     &                           u(i,j,k  ,nrhs))-                      &
     &                   DC(i,k)*CF(i,k-1))
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          END DO
        END DO
!
! Compute spline-interpolated, vertical advective u-momentum flux.
!
        cff3=1.0_r8/3.0_r8
        cff4=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))*                                &
     &              (u(i,j,k,nrhs)+                                     &
     &               DC(i,k)*(cff3*CF(i,k  )+                           &
     &                        cff4*CF(i,k-1)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=FC(i,k)-FC(i,k-1)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
          END DO
        END DO
        IF (j.ge.JstrV) THEN
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of v-momentum.
!
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=(cff1*(Hz(i,j  ,k)+                               &
     &                       Hz(i,j-1,k))-                              &
     &                 cff2*(Hz(i,j+1,k)+                               &
     &                       Hz(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.0_r8*(v(i,j,k+1,nrhs)-                     &
     &                             v(i,j,k  ,nrhs))-                    &
     &                     DC(i,k)*CF(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            END DO
          END DO
!
! Compute spline-interpolated, vertical advective v-momentum flux.
!
          cff3=1.0_r8/3.0_r8
          cff4=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=(cff1*(W(i,j  ,k)+                                &
     &                       W(i,j-1,k))-                               &
     &                 cff2*(W(i,j+1,k)+                                &
     &                       W(i,j-2,k)))*                              &
     &                (v(i,j,k,nrhs)+                                   &
     &                 DC(i,k)*(cff3*CF(i,k  )+                         &
     &                          cff4*CF(i,k-1)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=FC(i,k)-FC(i,k-1)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        DO i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1,nrhs)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k,nrhs)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=om_u(i,j)*on_u(i,j)
          cff1= sustr(i,j)*cff
          cff2=-bustr(i,j)*cff
          rufrc(i,j)=rufrc(i,j)+cff1+cff2
        END DO
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            rvfrc(i,j)=rv(i,j,1,nrhs)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k,nrhs)
            END DO
          END DO
          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
            cff1= svstr(i,j)*cff
            cff2=-bvstr(i,j)*cff
            rvfrc(i,j)=rvfrc(i,j)+cff1+cff2
          END DO
        END IF
      END DO J_LOOP
      RETURN
      END SUBROUTINE rhs3d_tile
      END MODULE rhs3d_mod
