      SUBROUTINE wrt_avg (ng)
!
!svn $Id: wrt_avg.F 1023 2009-07-20 21:45:24Z kate $
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
      USE mod_average
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
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
      SourceFile='wrt_avg.F'
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
      tAVGindx(ng)=tAVGindx(ng)+1
      NrecAVG(ng)=NrecAVG(ng)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, AVGname(ng),                      &
     &                      TRIM(Vname(idtime,ng)), AVGtime(ng:),       &
     &                      (/tAVGindx(ng)/), (/1/),                    &
     &                      ncid = ncAVGid(ng),                         &
     &                      varid = avgVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idFsur,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgzeta)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUbar,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgu2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVbar,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgv2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idUvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     AVERAGE(ng) % avgu3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idVvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     AVERAGE(ng) % avgv3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idOvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     AVERAGE(ng) % avgw3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idWvel,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     AVERAGE(ng) % avgwvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tAVGindx(ng)
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
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgTid(itrc,ng),    &
     &                       tAVGindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       AVERAGE(ng) % avgt(:,:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                        tAVGindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!-----------------------------
!  Write out benthic variables.
!-----------------------------
      DO itrc=1,NBeT(ng)
        IF (Hout(idBvar(itrc),ng)) THEN
          scale=1.0_r8
	    gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgBid(itrc,ng),      &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgbt(:,:,1,itrc))
! switch if have more than one benthic level	  
!          gtype=gfactor*r3dvar
!          status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgBid(itrc,ng),   &
!     &                       tAVGindx(ng), gtype,                      &
!     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,      &
!     &                       GRID(ng) % rmask(LBi,LBj),                &
!     &                       AVERAGE(ng) % avgbt(LBi,LBj,1,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBvar(itrc))),            &
     &                        tAVGindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!-----------------------------
!  Write out ice bio variables.
!-----------------------------
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idDano,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     AVERAGE(ng) % avgrho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out depth of surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idHsbl,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avghsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idVvis,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   AVERAGE(ng) % avgAKv)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idTdif,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   AVERAGE(ng) % avgAKt)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, ncAVGid(ng), avgVid(idSdif,ng),      &
     &                   tAVGindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   AVERAGE(ng) % avgAKs)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)), tAVGindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface net heat flux.
!
      IF (Hout(idTsur(itemp),ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng),                       &
     &                     avgVid(idTsur(itemp),ng),                    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgstf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Hout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng),                       &
     &                     avgVid(idTsur(isalt),ng),                    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgswf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idLhea,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avglhf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idShea,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgshf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idLrad,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avglrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-wind.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUair,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avguwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), tAVGindx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-wind.
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVair,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgvwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idevap,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgevap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idrain,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgrain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idSrad,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgsrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tAVGindx(ng)
          END IF          
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-momentun stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUsms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgsus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tAVGindx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentun stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVsms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgsvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tAVGindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentun stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idUbms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgbus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tAVGindx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentun stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVGid(ng), avgVid(idVbms,ng),    &
     &                     tAVGindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     AVERAGE(ng) % avgbvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tAVGindx(ng)
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
      CALL netcdf_sync (ng, iNLM, AVGname(ng), ncAVGid(ng))
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) tAVGindx(ng)
!
  10  FORMAT (/,' WRT_AVG - error while writing variable: ',a,/,11x,    &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_AVG   - wrote averaged fields into time ',        &
     &        'record =',t72,i7.7)
      RETURN
      END SUBROUTINE wrt_avg
