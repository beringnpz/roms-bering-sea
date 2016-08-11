      SUBROUTINE def_rst (ng)
!
!svn $Id: def_rst.F 1023 2009-07-20 21:45:24Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates restart NetCDF file, it defines its            !
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
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)
      integer, parameter :: Natt = 25
      integer :: i, j, nvd3, nvd4, nvd5
      integer :: recdim, status, varid
      integer :: DimIDs(31)
      integer :: r2dgrd(4), ru2dgrd(4), rv2dgrd(4)
      integer :: sr2dgrd(3), su2dgrd(3), sv2dgrd(3)
      integer :: sr3dgrd(4), su3dgrd(4), sv3dgrd(4)
      integer :: t2dgrd(4), u2dgrd(4), v2dgrd(4)
      integer :: Vsize(4)
      integer :: def_dim
      integer :: itrc
      integer :: k3dgrd(5), t3dgrd(5)
      integer :: r3dgrd(4), ru3dgrd(5), rv3dgrd(5) 
      integer :: u3dgrd(5), v3dgrd(5), w3dgrd(4)
      real(r8) :: Aval(6)
      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
      SourceFile='def_rst.F'
!
!=======================================================================
!  Create a new restart NetCDF file.
!=======================================================================
!
!  Activate creation of restart NetCDF file.  Create a new restart
!  file if during a restart run, the restart filename "RSTname"
!  is different than the initial filename "INIname".
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=RSTname(ng)
      Ldefine=.FALSE.
      IF (((nrrec(ng).eq.0).and.(iic(ng).eq.ntstart(ng))).or.           &
     &    ((nrrec(ng).ne.0).and.                                        &
     &     (TRIM(ncname).ne.TRIM(INIname(ng))))) THEN
        Ldefine=.TRUE.
      END IF
!
      DEFINE : IF (Ldefine) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), ncRSTid(ng))
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'xi_rho',         &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'xi_u',           &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'xi_v',           &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'xi_psi',         &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'eta_rho',        &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'eta_u',          &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'eta_v',          &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'eta_psi',        &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'N',              &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 's_rho',          &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 's_w',            &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'tracer',         &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname, 'boundary',       &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, ncRSTid(ng), ncname,                   &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
        nvd5=5
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        sr2dgrd(1)=DimIDs( 1)
        sr2dgrd(2)=DimIDs( 5)
        sr2dgrd(3)=DimIDs(12)
        t2dgrd(3)=DimIDs(12)
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        r3dgrd(1)=DimIDs( 1)
        r3dgrd(2)=DimIDs( 5)
        r3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
        r3dgrd(4)=DimIDs(12)
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
!  Define dimension vector for sediment, radiation stress variables.
!
        su2dgrd(1)=DimIDs( 2)
        su2dgrd(2)=DimIDs( 6)
        su2dgrd(3)=DimIDs(12)
        sv2dgrd(1)=DimIDs( 3)
        sv2dgrd(2)=DimIDs( 7)
        sv2dgrd(3)=DimIDs(12)
        sr3dgrd(1)=DimIDs( 1)
        sr3dgrd(2)=DimIDs( 5)
        sr3dgrd(3)=DimIDs(16)
        sr3dgrd(4)=DimIDs(12)
        su3dgrd(1)=DimIDs( 2)
        su3dgrd(2)=DimIDs( 6)
        su3dgrd(3)=DimIDs( 9)
        su3dgrd(4)=DimIDs(12)
        sv3dgrd(1)=DimIDs( 3)
        sv3dgrd(2)=DimIDs( 7)
        sv3dgrd(3)=DimIDs( 9)
        sv3dgrd(4)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        tRSTindx(ng)=0
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
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, ncRSTid(ng), ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
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
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idtime,ng),        &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idFsur,ng),        &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)=Vname(3,idUbar)
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idUbar,ng),        &
     &                 NF_FRST, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)=Vname(3,idVbar)
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idVbar,ng),        &
     &                 NF_FRST, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D momentum component in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUvel)
        Vinfo( 2)=Vname(2,idUvel)
        Vinfo( 3)=Vname(3,idUvel)
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idUvel,ng),        &
     &                 NF_FRST, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D momentum component in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVvel)
        Vinfo( 2)=Vname(2,idVvel)
        Vinfo( 3)=Vname(3,idVvel)
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idVvel,ng),        &
     &                 NF_FRST, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, ncRSTid(ng), rstTid(itrc,ng),        &
     &                   NF_FRST, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!------------------------------------------
!  Define benthic tracer type variables.
!------------------------------------------   
         DO itrc=1,NBeT(ng)
          rstBid(itrc,ng) = var_id(i)
          IF (Hout(idBvar(itrc),ng)) THEN
          Vinfo( 1)=Vname(1,idBvar(itrc))
          Vinfo( 2)=Vname(2,idBvar(itrc))
          Vinfo(3)=Vname(3,idBvar(itrc))
          Vinfo(14)=Vname(4,idBvar(itrc))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng,iNLM,ncRSTid(ng),rstBid(itrc,ng),NF_FOUT,   &
     &                     nvd3,t2dgrd,Aval,Vinfo,ncname)
!Will need to switch this if have more than one depth level for benthos 
!         status=def_var(ncHISid(ng),hisBid(itrc,ng),NF_FOUT,           &
!    &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
               END IF
            END DO
!------------------------------------------
!  Define ice bio tracer type variables.
!------------------------------------------  
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idDano,ng),        &
     &                 NF_FRST, nvd4, r3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define depth of surface boundary layer.
!
        Vinfo( 1)=Vname(1,idHsbl)
        Vinfo( 2)=Vname(2,idHsbl)
        Vinfo( 3)=Vname(3,idHsbl)
        Vinfo(14)=Vname(4,idHsbl)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idHsbl,ng),        &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        Vinfo( 2)=Vname(2,idVvis)
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idVvis,ng),        &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        Vinfo( 2)=Vname(2,idTdif)
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idTdif,ng),        &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        Vinfo( 2)=Vname(2,idSdif)
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
        status=def_var(ng, iNLM, ncRSTid(ng), rstVid(idSdif,ng),        &
     &                 NF_FRST, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, ncRSTid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, ncRSTid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing restart file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.Ldefine) THEN
        ncname=RSTname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open restart file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, ncRSTid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,30) TRIM(ncname)
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
!  restart variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            rstVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            rstVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            rstVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            rstVid(idVbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            rstVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            rstVid(idVvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            rstVid(idDano,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            rstVid(idHsbl,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            rstVid(idVvis,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            rstVid(idTdif,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            rstVid(idSdif,ng)=var_id(i)
          END IF
          DO itrc=1,NT(ng)
           IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              rstTid(itrc,ng)=var_id(i)
	       END IF
          END DO
          DO itrc=1,NBeT(ng)
           IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBvar(itrc)))) THEN
              got_var(idBvar(itrc))=.TRUE.
              rstBid(itrc,ng)=var_id(i)
	       END IF
          END DO
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to current value.
!
        IF (LcycleRST(ng)) THEN
          tRSTindx(ng)=0
        ELSE
          tRSTindx(ng)=rec_size
        END IF
      END IF QUERY
!
  10  FORMAT (/,' DEF_RST - unable to create restart NetCDF file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_RST - unable to open restart NetCDF file: ',a)
  40  FORMAT (/,' DEF_RST - unable to find variable: ',a,2x,            &
     &        ' in restart NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_rst
