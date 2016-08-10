      MODULE ini_fields_mod
!
!svn $Id: ini_fields.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for 2D fields. It also   !
!  couples 3D and 2D momentum equations:  it initializes 2D momentum   !
!  (ubar,vbar) to the vertical integral of initial 3D momentum (u,v).  !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC :: ini_fields
      PUBLIC :: ini_zeta
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_coupling
      USE mod_ocean
      USE mod_stepping
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
      CALL ini_fields_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      kstp(ng), krhs(ng), knew(ng),               &
     &                      nstp(ng), nnew(ng),                         &
     &                      GRID(ng) % h,                               &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      COUPLING(ng) % Zt_avg1,                     &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE ini_fields
!
!***********************************************************************
      SUBROUTINE ini_fields_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            kstp, krhs, knew,                     &
     &                            nstp, nnew,                           &
     &                            h,                                    &
     &                            Hz,                                   &
     &                            z_r, z_w,                             &
     &                            Zt_avg1,                              &
     &                            t, u, v,                              &
     &                            ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE set_depth_mod, ONLY : set_depth_tile
      USE t3dbc_mod, ONLY : t3dbc_tile
      USE u3dbc_mod, ONLY : u3dbc_tile
      USE v3dbc_mod, ONLY : v3dbc_tile
      USE u2dbc_mod, ONLY : u2dbc_tile
      USE v2dbc_mod, ONLY : v2dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
      integer, intent(in) :: nstp, nnew
!
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k, kbed
      real(r8) :: cff1, cff2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
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
!  Initialize other 3D momentum time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff1=u(i,j,k,nstp)
              u(i,j,k,nstp)=cff1
              u(i,j,k,nnew)=cff1
            END DO
          END DO
          IF (j.ge.JstrV) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff2=v(i,j,k,nstp)
                v(i,j,k,nstp)=cff2
                v(i,j,k,nnew)=cff2
              END DO
            END DO
          END IF
        END DO
!
!  Apply boundary conditions.
!
        CALL u3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nstp,                                    &
     &                   u)
        CALL v3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nstp,                                    &
     &                   v)
        CALL u3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   u)
        CALL v3dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng),                     &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   v)
      END IF
!
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nstp))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nstp))
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        u(:,:,:,nnew))
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        v(:,:,:,nnew))
!
!-----------------------------------------------------------------------
!  Initialize other tracers time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff1=t(i,j,k,nstp,itrc)
                t(i,j,k,nstp,itrc)=cff1
                t(i,j,k,nnew,itrc)=cff1
              END DO
            END DO
          END DO
!
!  Apply boundary conditions.
!
          CALL t3dbc_tile (ng, tile, itrc,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nstp,                                  &
     &                     t)
          CALL t3dbc_tile (ng, tile, itrc,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     t)
        END DO
      END IF
!
      DO itrc=1,NT(ng)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nstp,itrc))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,itrc))
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial depths and thicknesses.
!-----------------------------------------------------------------------
!
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     h,                                           &
     &                     Zt_avg1,                                     &
     &                     Hz, z_r, z_w)
!
!-----------------------------------------------------------------------
!  Compute vertically integrated momentum (ubar,vbar) from initial 3D
!  momentum (u,v).
!-----------------------------------------------------------------------
!
!  Compute adjoint 2D velocity component in the XI-direction.  Here
!  DC(i,1:N) are the thicknesses of U-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrU,Iend
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
            END DO
          END DO
          DO i=IstrU,Iend
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
            ubar(i,j,kstp)=cff2
            ubar(i,j,knew)=cff2
          END DO
!
!  Compute adjoint 2D velocity component in the ETA-direction.  Here
!  DC(i,1:N) are the thicknesses of V-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
          IF (j.ge.Jstr) THEN
            DO i=Istr,Iend
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
                CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
              END DO
            END DO
            DO i=Istr,Iend
              cff1=1.0_r8/DC(i,0)
              cff2=CF(i,0)*cff1
              vbar(i,j,kstp)=cff2
              vbar(i,j,knew)=cff2
            END DO
          END IF   
        END DO
!
!  Apply boundary conditions.
!
        CALL u2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, kstp,                              &
     &                   ubar, vbar, zeta)
        CALL u2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, knew,                              &
     &                   ubar, vbar, zeta)
        CALL v2dbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs, kstp, knew,                              &
     &                   ubar, vbar, zeta)
      END IF
!
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,kstp))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,kstp))
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubar(:,:,knew))
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbar(:,:,knew))
      RETURN
      END SUBROUTINE ini_fields_tile
!
!***********************************************************************
      SUBROUTINE ini_zeta (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_coupling
      USE mod_ocean
      USE mod_stepping
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
      CALL ini_zeta_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    kstp(ng), krhs(ng), knew(ng),                 &
     &                    COUPLING(ng) % Zt_avg1,                       &
     &                    OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE ini_zeta
!
!***********************************************************************
      SUBROUTINE ini_zeta_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          kstp, krhs, knew,                       &
     &                          Zt_avg1,                                &
     &                          zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE zetabc_mod, ONLY : zetabc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
!
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff1, cff2
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
!  Initialize other free-surface time-levels.
!-----------------------------------------------------------------------
!
      IF (.not.PerfectRST(ng)) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff1=zeta(i,j,kstp)
            zeta(i,j,kstp)=cff1
            zeta(i,j,knew)=cff1
          END DO
        END DO
!
!  Apply boundary conditions.
!
        CALL zetabc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, kstp,                             &
     &                    zeta)
        CALL zetabc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    zeta)
      END IF
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,kstp))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,knew))
      IF (PerfectRST(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,krhs))
      END IF
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Zt_avg1(i,j)=zeta(i,j,kstp)
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Zt_avg1)
      RETURN
      END SUBROUTINE ini_zeta_tile
      END MODULE ini_fields_mod
