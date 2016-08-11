      SUBROUTINE wrt_avg2 (ng)
!
!svn $Id: wrt_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged fields into averages     !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average2
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
        gfactor=1
!
!  Set time record index.
!
      tAVG2indx(ng)=tAVG2indx(ng)+1
      NrecAVG2(ng)=NrecAVG2(ng)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, AVG2name(ng),                     &
     &                      TRIM(Vname(idtime,ng)), AVG2time(ng:),      &
     &                      (/tAVG2indx(ng)/), (/1/),                   &
     &                      ncid = ncAVG2id(ng),                        &
     &                      varid = avg2Vid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Hout2(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idFsur,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgzeta)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUbar,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgu2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVbar,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgv2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUvel,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgu3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVvel,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgv3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout2(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Tid(itrc,ng),  &
     &                       tAVG2indx(ng), gtype,                      &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       AVERAGE2(ng) % avgt(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                        tAVG2indx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout2(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idDano,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgrho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out depth of surface boundary layer.
!
      IF (Hout2(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idHsbl,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avghsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net heat flux.
!
      IF (Hout2(idTsur(itemp),ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng),                      &
     &                     avg2Vid(idTsur(itemp),ng),                   &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgstf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Hout2(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng),                      &
     &                     avg2Vid(idTsur(isalt),ng),                   &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgswf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out latent heat flux.
!
      IF (Hout2(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idLhea,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avglhf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout2(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idShea,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgshf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout2(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idLrad,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avglrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-wind.
!
      IF (Hout2(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUair,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avguwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), tAVG2indx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-wind.
!
      IF (Hout2(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVair,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgvwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Hout2(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idevap,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgevap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout2(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idrain,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgrain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out shortwave radiation flux.
!
      IF (Hout2(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idSrad,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgsrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tAVG2indx(ng)
          END IF          
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-momentun stress.
!
      IF (Hout2(idUsms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUsms,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgsus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tAVG2indx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentun stress.
!
      IF (Hout2(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVsms,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgsvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentun stress.
!
      IF (Hout2(idUbms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUbms,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgbus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tAVG2indx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentun stress.
!
      IF (Hout2(idVbms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVbms,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE2(ng) % avgbvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, AVG2name(ng), ncAVG2id(ng))
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) tAVG2indx(ng)
!
  10  FORMAT (/,' WRT_AVG2 - error while writing variable: ',a,/,11x,   &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_AVG2   - wrote averaged fields into time ',       &
     &        'record =',t72,i7.7)
      RETURN
      END SUBROUTINE wrt_avg2
