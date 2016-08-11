      SUBROUTINE def_avg (ng, ldef)
!
!svn $Id: def_avg.F 990 2009-05-28 00:54:01Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates averages NetCDF file, it defines its           !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
      integer, parameter :: Natt = 25
      integer :: i, itrc, j, model, nvd3, nvd4
      integer :: recdim, status
      integer :: DimIDs(31), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)
      integer :: def_dim
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
      real(r8) :: Aval(6)
      character (len=13) :: Prefix
      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
      SourceFile='def_avg.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=AVGname(ng)
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
      model=iNLM
!
!=======================================================================
!  Create a new averages file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), ncAVGid(ng))
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, model, ncAVGid(ng), ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
!
!  Define dimension vector for staggered w-momemtum type variables.
!
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        tAVGindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Set long name prefix string.
!
        Prefix='time-averaged'
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, ncAVGid(ng), ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,1x,a)') 'averaged', TRIM(Vname(2,idtime))
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idtime,ng),       &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idFsur,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUbar,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVbar,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUvel,ng),     &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVvel,ng),     &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOvel))
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idOvel,ng),     &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define "true" vertical momentum component.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvel))
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idWvel,ng),     &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, ncAVGid(ng), avgTid(itrc,ng),     &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!------------------
!benthic bio variables
!-------------------
        DO itrc=1,NBeT(ng)
         IF (Hout(idBvar(itrc),ng)) THEN
             Vinfo( 1)=Vname(1,idBvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idBvar(itrc)))
         Vinfo( 1)=Vname(1,idBvar(itrc))
         Vinfo( 2)=Vname(2,idBvar(itrc))
         Vinfo(3)=Vname(3,idBvar(itrc))
         Vinfo(14)=Vname(4,idBvar(itrc))
          Vinfo(22)='coordinates'
         Vinfo(16)=Vname(1,idtime)
	 Aval(5)=REAL(r2dvar,r8)
         status=def_var(ng, model,ncavgid(ng),avgBid(itrc,ng),NF_FOUT,  &
     &                     nvd3,t2dgrd,Aval,Vinfo,ncname)
              END IF
           END DO
!------------------
!ice bio variables
!-------------------
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idDano,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHsbl))
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idHsbl,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvis))
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idVvis,ng),       &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTdif))
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idTdif,ng),       &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idSdif,ng),       &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface net heat flux.
!
        IF (Hout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTsur(itemp),ng), NF_FOUT,             &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Hout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTsur(isalt),ng), NF_FOUT,             &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLhea))
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idLhea,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idShea))
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idShea,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLrad))
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idLrad,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTair))
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idTair,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface u-wind.
!
        IF (Hout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUair))
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUair,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface v-wind.
!
        IF (Hout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVair))
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVair,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idevap))
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idevap,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idrain,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSrad))
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idSrad,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface u-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUsms,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVsms,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUbms,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVbms,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, ncAVGid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, ncAVGid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVGname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open averages file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, ncAVGid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  average variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            avgVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            avgVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            avgVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            avgVid(idVbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            avgVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            avgVid(idVvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            avgVid(idOvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            avgVid(idWvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            avgVid(idDano,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            avgVid(idHsbl,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            avgVid(idVvis,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            avgVid(idTdif,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            avgVid(idSdif,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            avgVid(idTsur(itemp),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            avgVid(idTsur(isalt),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            avgVid(idLhea,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            avgVid(idShea,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            avgVid(idLrad,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            avgVid(idTair,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            avgVid(idUair,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            avgVid(idVair,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            avgVid(idevap,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            avgVid(idrain,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            avgVid(idSrad,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            avgVid(idUsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            avgVid(idVsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            avgVid(idUbms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            avgVid(idVbms,ng)=var_id(i)
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             avgTid(itrc,ng)=var_id(i)
            END IF
          END DO 
        END DO
!
!  Check if averages variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(itemp)).and.Hout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Hout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUair)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVair)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nAVG(ng)) THEN
          IF (ndefAVG(ng).gt.0) THEN
            tAVGindx(ng)=((ntstart(ng)-1)-                              &
     &                    ndefAVG(ng)*((ntstart(ng)-1)/ndefAVG(ng)))/   &
     &                   nAVG(ng)
          ELSE
            tAVGindx(ng)=(ntstart(ng)-1)/nAVG(ng)
          END IF
        ELSE
          tAVGindx(ng)=rec_size
        END IF
      END IF QUERY
!
!  Set initial average time. Notice that the value is offset by half
!  nAVG*dt so there is not a special case when computing its value
!  in "set_avg".
!
      IF (ntsAVG(ng).eq.1) THEN
        AVGtime(ng)=time(ng)-0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      ELSE
        AVGtime(ng)=time(ng)+REAL(ntsAVG(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (6x,'DEF_AVG   - creating average file: ',a)
  20  FORMAT (6x,'DEF_AVG   - inquiring average file: ',a)
  30  FORMAT (/,' DEF_AVG - unable to create averages NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG - unable to open averages NetCDF file: ',a)
  60  FORMAT (/,' DEF_AVG - unable to find variable: ',a,2x,            &
     &        ' in averages NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_avg
