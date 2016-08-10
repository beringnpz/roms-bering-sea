      MODULE step3d_t_mod
!
!svn $Id: step3d_t.F 991 2009-05-28 23:37:09Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations. Notice that advective     !
!  and diffusive terms are time-stepped differently.                   !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: step3d_t
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
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
      CALL wclock_on (ng, iNLM, 35)
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
     &                    CLIMA(ng) % Tnudgcof,                         &
     &                    CLIMA(ng) % tclm,                             &
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % t                                 &
     &                    ,OCEAN(ng) % bt                               &
     &                               )
      CALL wclock_off (ng, iNLM, 35)
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
     &                          Tnudgcof, tclm,                         &
     &                          Akt,                                    &
     &                          W,                                      &
     &                          t                                       &
     &                          ,bt                                     &
     &                           )
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Tnudgcof(LBi:,LBj:,:)
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
     real(r8), intent(inout) :: bt(LBi:UBi,LBj:UBj,NBL(ng),3,NBeT(ng))
!
!  Local variable declarations.
!
      integer :: i, ibt, is, itrc, j, k, ltrc
      integer :: idiag
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
      real(r8) :: my_maxbio(15)
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
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng) 
        K_LOOP : DO k=1,N(ng)
!
!  Fourth-order, centered differences horizontal advective fluxes.
!  
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+2
              FX(i,j)=t(i  ,j,k,3,itrc)-                                &
     &                t(i-1,j,k,3,itrc)
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc)-                               &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
            END DO
          END DO
!
          DO j=Jstr-1,Jend+2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,3,itrc)-                                &
     &                t(i,j-1,k,3,itrc)
            END DO
          END DO
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc)-                               &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
            END DO
          END DO
!
!  Time-step horizontal advection term.
!
!       print*,'before horz advect','t=',t(1,1,19,2,1)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j)+                              &
     &                  FE(i,j+1)-FE(i,j))
!g	   IF ( itrc.gt.2 ) THEN
!g               
!g		cff1=0.0_r8
!g		
!g           END IF
	   t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
!Constraining nutrients in inlets that are only 1 grid cell wide 	   
!Cook Inlet	   
!	 if(i.ge.171.and.i.le.179.and.j.ge.18.and.j.le.27)THEN
!	  
!	  if(itrc.eq.iFe)THEN
!	   t(i,j,k,nnew,itrc)=min(2.0_r8,t(i,j,k,nnew,itrc))
!	   else if(itrc.eq.iNH4)THEN
!	   t(i,j,k,nnew,itrc)=min(5.0_r8,t(i,j,k,nnew,itrc))
!	    else if(itrc.eq.iNO3)THEN
!	   t(i,j,k,nnew,itrc)=min(15.0_r8,t(i,j,k,nnew,itrc))
!	  endif
!	 endif 
!Hagemeister Strait   
!	  if(i.ge.121.and.i.le.124.and.j.ge.56.and.j.le.59)THEN
!	  if(itrc.eq.iFe)THEN
!	   t(i,j,k,nnew,itrc)=min(2.0_r8,t(i,j,k,nnew,itrc))
!	   else if(itrc.eq.iNH4)THEN
!	   t(i,j,k,nnew,itrc)=min(5.0_r8,t(i,j,k,nnew,itrc))
!	    else if(itrc.eq.iNO3)THEN
!	   t(i,j,k,nnew,itrc)=min(15.0_r8,t(i,j,k,nnew,itrc))
!	  endif
!	 endif 
!            endif
            END DO
          END DO
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
!      print*,'before vert advect','t=',t(1,1,19,2,1)
      DO j=Jstr,Jend
        DO itrc=1,NT(ng) 
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,3,itrc)+                         &
     &                       t(i,j,k+1,3,itrc))-                        &
     &                 cff3*(t(i,j,k-1,3,itrc)+                         &
     &                       t(i,j,k+2,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,3,itrc)+                              &
     &               cff2*t(i,j,2,3,itrc)-                              &
     &               cff3*t(i,j,3,3,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step vertical advection term.
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
!g	   IF ( itrc.gt.2 ) THEN
!g                
!g		cff1=0.0_r8
!g		
!g           END IF
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng) 
          ltrc=MIN(NAT,itrc)
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
          cff=-dt(ng)*lambda
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              DC(i,k)=t(i,j,k,nnew,itrc)
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
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
     &                   (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
!g      IF (itrc.eq.iNCaS.or.itrc.eq.iNCaO.or.itrc.eq.iEupS.or.         &
!g     &        itrc.eq.iEupO.or.itrc.eq.iJel) THEN
!g      t(i,j,N(ng),nnew,itrc)=t(i,j,N(ng),nnew,itrc)      
!g      else    
      t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
!g             END IF
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!g        IF (itrc.eq.iNCaS.or.itrc.eq.iNCaO.or.itrc.eq.iEupS.or.       &
!g     &   itrc.eq.iEupO.or.itrc.eq.iJel) THEN
!g              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
!g           else
              t(i,j,k,nnew,itrc)=DC(i,k)
!g              END IF
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng) 
!
!  Set lateral boundary conditions.
!
        CALL t3dbc_tile (ng, tile, itrc,                                &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)
!
!  Nudge towards tracer climatology.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
!              
            END DO
          END DO
        END DO
!
!  Apply periodic boundary conditions.
!
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,itrc))
      END DO
        do i=IstrR,IendR
        do j=JstrR,JendR
            do k=1,NBL(ng)
                  bt(i,j,k,nstp,1) = bt(i,j,k,nnew,1)
                  bt(i,j,k,nstp,2) = bt(i,j,k,nnew,2)
!          if(i.eq.144.and.j.eq.144) THEN
!	  print*,'Step3d'
!	  print*,'nstp=',nstp,'nnew=',nnew
!          print*,'bt(i,j,k,nstp,1)=', bt(i,j,k,nstp,1)
!         
!          print*,'bt(i,j,k,nnew,1)=', bt(i,j,k,nnew,1)
!          print*,''
!          end if
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE step3d_t_tile
      END MODULE step3d_t_mod
