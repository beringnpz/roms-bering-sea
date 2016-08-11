      SUBROUTINE wrt_his (ng)
!
!svn $Id: wrt_his.F 1023 2009-07-20 21:45:24Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested model fields at requested levels      !
!  into history NetCDF file.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod, ONLY : scale_omega
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
      integer :: gfactor, gtype, status
      integer :: i, itrc, j, k, tile
      real(r8) :: scale
      real(r8), allocatable :: wrk(:,:,:)
!
      SourceFile='wrt_his.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      tHISindx(ng)=tHISindx(ng)+1
      NrecHIS(ng)=NrecHIS(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, HISname(ng),                      &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/tHISindx(ng)/), (/1/),                    &
     &                      ncid = ncHISid(ng),                         &
     &                      varid = hisVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idFsur,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     OCEAN(ng) % zeta(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbar,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     OCEAN(ng) % ubar(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbar,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     OCEAN(ng) % vbar(:,:,kstp(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idUvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     OCEAN(ng) % u(:,:,:,nrhs(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     OCEAN(ng) % v(:,:,:,nrhs(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tHISindx(ng)
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
        IF (.not.allocated(wrk)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj,0:N(ng)))
          wrk(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_r8
        gtype=gfactor*w3dvar
        DO tile=0,NtileX(ng)*NtileE(ng)-1
          CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),     &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % W,                              &
     &                      wrk)
        END DO
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idOvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     wrk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (wrk)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idWvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     OCEAN(ng) % wvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), tHISindx(ng)
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
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisTid(itrc,ng),    &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!--------------------------------------
!  Write out 2D production array
!--------------------------------------
!--------------------------------------
!  Write out 3D production array
!--------------------------------------
!------------------------------
!Write out benthic variables
!-----------------------------
              DO itrc=1,NBeT(ng)
              IF (Hout(idBvar(itrc),ng)) THEN
                 scale=1.0_r8
                 gtype=gfactor*r2dvar
             status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisBid(itrc,ng), &
     &                          tHISindx(ng), gtype,                    &
     &                          LBi, UBi, LBj, UBj, scale,              &
!Will need to switch this if have more than one depth level for benthos 
!      status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisBid(itrc,ng),       &
!     &                          tHISindx(ng), gtype,                   &
!     &                          LBi, UBi, LBj, UBj, 1, NBL(ng), scale, &
     &                       OCEAN(ng) % bt(:,:,1,nrhs(ng),itrc))
!             print*,OCEAN(ng) % bt(LBi,LBj,1,nrhs(ng),itrc)
!             print *, 'LBi=',LBi,'LBj=',LBj,'NOUT=',nrhs(ng),'itrc=',itrc
           IF (status.ne.nf90_noerr) THEN
           IF (Master) THEN
            WRITE (stdout,10)TRIM(Vname(1,idBvar(itrc))),               &
     &                       tHISindx(ng)
           END IF
              exit_flag=3
              ioerror=status
            RETURN
          END IF
        END IF
      END DO
!------------------------------
!Write out ice bio variables
!-----------------------------
!--------------------------------------
!  Write out stationary tracer variable
!--------------------------------------
!----------------------------
!
!  Write out density anomaly.
!----------------------------
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idDano,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     OCEAN(ng) % rho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idHsbl,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     MIXING(ng) % hsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVvis,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     MIXING(ng) % Akv)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idTdif,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     MIXING(ng) % Akt(:,:,:,itemp))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idSdif,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     MIXING(ng) % Akt(:,:,:,isalt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface active traces fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idTsur(itrc),ng),                   &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idLhea,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % lhflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tHISindx(ng)
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
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idShea,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % shflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tHISindx(ng)
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
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idLrad,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % lrflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-wind.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % Uwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-wind.
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % Vwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), tHISindx(ng)
          END IF 
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idevap,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % evap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tHISindx(ng)
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
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idrain,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % rain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tHISindx(ng)
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
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idSrad,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % srflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUsms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVsms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tHISindx(ng)
          END IF 
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % bustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     FORCES(ng) % bvstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HISname(ng), ncHISid(ng))
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) kstp(ng), nrhs(ng), tHISindx(ng)
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
      RETURN
      END SUBROUTINE wrt_his
