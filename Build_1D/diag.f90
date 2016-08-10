      MODULE diag_mod
!
!svn $Id: diag.F 916 2009-02-05 00:38:29Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes various diagnostic fields.                    !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: diag
      CONTAINS
!
!***********************************************************************
      SUBROUTINE diag (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
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
      CALL wclock_on (ng, iNLM, 7)
      CALL diag_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                nstp(ng), krhs(ng),                               &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % omn,                                   &
     &                GRID(ng) % Hz,                                    &
     &                GRID(ng) % z_r,                                   &
     &                GRID(ng) % z_w,                                   &
     &                OCEAN(ng) % rho,                                  &
     &                OCEAN(ng) % u,                                    &
     &                OCEAN(ng) % v,                                    &
     &                OCEAN(ng) % t,                                    &
     &                OCEAN(ng) % ubar,                                 &
     &                OCEAN(ng) % vbar,                                 &
     &                OCEAN(ng) % zeta)
      CALL wclock_off (ng, iNLM, 7)
      RETURN
      END SUBROUTINE diag
!
!***********************************************************************
      SUBROUTINE diag_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, krhs,                                 &
     &                      h, omn,                                     &
     &                      Hz, z_r, z_w,                               &
     &                      rho, u, v, t,                               &
     &                      ubar, vbar, zeta)
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
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, krhs
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k, trd, ibt, itrc
      integer :: my_threadnum      
      real(r8) :: cff, my_avgke, my_avgpe, my_volume
      real(r8) :: my_maxspeed, u2v2
      real(r8) :: my_maxrho
      real(r8) :: my_maxbio(NBT)
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ke2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: pe2d
      character (len=8) :: kechar, pechar
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
!  Compute and report out volume averaged kinetic, potential
!  total energy, and volume.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng)-1,ninfo(ng)).eq.0) THEN
        my_maxspeed=0.0_r8
        my_maxrho=-1.0E+20_r8
        DO ibt=1,NBT
          my_maxbio(ibt)=-1.0E+20_r8
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            ke2d(i,j)=0.0_r8
            pe2d(i,j)=0.5_r8*g*z_w(i,j,N(ng))*z_w(i,j,N(ng))
          END DO
          cff=g/rho0
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              u2v2=u(i  ,j,k,nstp)*u(i  ,j,k,nstp)+                     &
     &             u(i+1,j,k,nstp)*u(i+1,j,k,nstp)+                     &
     &             v(i,j  ,k,nstp)*v(i,j  ,k,nstp)+                     &
     &             v(i,j+1,k,nstp)*v(i,j+1,k,nstp)
              ke2d(i,j)=ke2d(i,j)+                                      &
     &                  Hz(i,j,k)*0.25_r8*u2v2
              pe2d(i,j)=pe2d(i,j)+                                      &
     &                  cff*Hz(i,j,k)*(rho(i,j,k)+1000.0_r8)*           &
     &                  (z_r(i,j,k)-z_w(i,j,0))
              my_maxspeed=MAX(my_maxspeed,SQRT(0.5_r8*u2v2))
              my_maxrho=MAX(my_maxrho,rho(i,j,k))
            END DO
          END DO
          DO k=N(ng),1,-1
            DO ibt=1,NBT
              itrc=idbio(ibt)
              DO i=Istr,Iend
                my_maxbio(ibt)=MAX(my_maxbio(ibt),t(i,j,k,nstp,itrc))
              END DO
            END DO
          END DO
        END DO
!
!  Integrate horizontally within one tile. In order to reduce the
!  round-off errors, the summation is performed in two stages. First,
!  the index j is collapsed and then the accumulation is carried out
!  along index i. In this order, the partial sums consist on much
!  fewer number of terms than in a straightforward two-dimensional
!  summation. Thus, adding numbers which are orders of magnitude
!  apart is avoided.
!
        DO i=Istr,Iend
          pe2d(i,Jend+1)=0.0_r8
          pe2d(i,Jstr-1)=0.0_r8
          ke2d(i,Jstr-1)=0.0_r8
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+                              &
     &                     omn(i,j)*(z_w(i,j,N(ng))-z_w(i,j,0))
            pe2d(i,Jstr-1)=pe2d(i,Jstr-1)+omn(i,j)*pe2d(i,j)
            ke2d(i,Jstr-1)=ke2d(i,Jstr-1)+omn(i,j)*ke2d(i,j)
          END DO
        END DO
        my_volume=0.0_r8
        my_avgpe=0.0_r8
        my_avgke=0.0_r8
        DO i=Istr,Iend
          my_volume=my_volume+pe2d(i,Jend+1)
          my_avgpe =my_avgpe +pe2d(i,Jstr-1)
          my_avgke =my_avgke +ke2d(i,Jstr-1)
        END DO
!
!  Perform global summation: whoever gets first to the critical region
!  resets global sums before global summation starts; after the global
!  summation is completed, thread, which is the last one to enter the
!  critical region, finalizes the computation of diagnostics and prints
!  them out.
!
        IF ((Istr.eq.1).and.(Jstr.eq.1).and.                            &
     &      (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
        IF (tile_count.eq.0) THEN
          volume=0.0_r8
          avgke=0.0_r8
          avgpe=0.0_r8
          maxspeed(ng)=0.0_r8
          maxrho(ng)=0.0_r8
          DO ibt=1,NBT
            maxbio(ibt,ng)=0.0_r8
          END DO
        END IF
        volume=volume+my_volume
        avgke=avgke+my_avgke
        avgpe=avgpe+my_avgpe
        maxspeed(ng)=MAX(maxspeed(ng),my_maxspeed)
        maxrho(ng)=MAX(maxrho(ng),my_maxrho)
        DO ibt=1,NBT
          maxbio(ibt,ng)=MAX(maxbio(ibt,ng),my_maxbio(ibt))
        END DO
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
          trd=my_threadnum()
          avgke=avgke/volume
          avgpe=avgpe/volume
          avgkp=avgke+avgpe
          IF (first_time(ng).eq.0) THEN
            first_time(ng)=1
            IF (Master) WRITE (stdout,10) 'STEP', 'Day HH:MM:SS',       &
     &                                    'KINETIC_ENRG', 'POTEN_ENRG', &
     &                                    'TOTAL_ENRG', 'NET_VOLUME'
 10         FORMAT (/,3x,a,3x,a,2x,a,3x,a,4x,a,4x,a,/)
          END IF
          IF (Master) WRITE (stdout,20) iic(ng)-1, time_code(ng),       &
     &                                  avgke, avgpe, avgkp, volume
 20       FORMAT (i7,1x,a,4(1pe14.6))
!          IF (Master) print *, 'BIO in diag ', maxbio(:,ng)
          IF (Master) CALL my_flush (stdout)
!
!  If blowing-up, set exit_flag to stop computations.
!
          WRITE (kechar,'(1pe8.1)') avgke
          WRITE (pechar,'(1pe8.1)') avgpe
          DO i=1,8
            IF ((kechar(i:i).eq.'N').or.(pechar(i:i).eq.'N').or.        &
     &          (kechar(i:i).eq.'n').or.(pechar(i:i).eq.'n').or.        &
     &          (kechar(i:i).eq.'*').or.(pechar(i:i).eq.'*')) THEN
              exit_flag=1
            END IF
          END DO
!
!  Stop computations if exceeding maximum speed allowed.  This will be
!  useful during debugging to avoid NaNs in output NetCDF files.
!
          IF (maxspeed(ng).gt.max_speed) THEN
            exit_flag=1
	    IF (Master) print *, 'DIAG speed trouble ', maxspeed(ng)
         END IF
!
!  Stop computation if exceeding maximum density anomaly allowed.
!  Recall that  density is computed from potential temperature and
!  salinity.  This is a good way to screen for very bad values which
!  indicates that the model is blowing-up.
!
          IF (maxrho(ng).gt.max_rho) THEN
            exit_flag=1
	    IF (Master) print *, 'DIAG rho trouble ', maxrho(ng)
          END IF
!          DO ibt=1,NBT
!            IF (maxbio(ibt,ng).gt.max_bio(ibt)) THEN
!              IF (Master) print *, 'BIO trouble in DIAG ', maxbio(:,ng)
!              exit_flag=1
!            END IF
!          END DO
        END IF
      END IF
      RETURN
      END SUBROUTINE diag_tile
      END MODULE diag_mod
