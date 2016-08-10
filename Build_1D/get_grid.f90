      SUBROUTINE get_grid (ng, model)
!
!svn $Id: get_grid.F 1056 2009-09-08 18:51:43Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine reads grid information from GRID NetCDF file.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE exchange_2d_mod
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile, LBi, UBi, LBj, UBj
      integer :: gtype, i, j, it, status, vindex
      integer :: Vsize(4)
      real(r8), parameter :: Fscl = 1.0_r8
      real(r8) :: Fmax, Fmin
      character (len=1 ) :: char1
      character (len=80) :: ncname
!
      SourceFile='get_grid.F'
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=GRDname(ng)
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'xl',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'el',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'spherical',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'h',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'f',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'f', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pm',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pn',vindex)) THEN
        IF (Master) WRITE (stdout,10) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Open grid NetCDF file for reading.
!
      IF (ncGRDid(ng).eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, ncGRDid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,20) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL netcdf_get_svar (ng, model, ncname, 'spherical',           &
     &                        char1,                                    &
     &                         ncid = ncGRDid(ng))
        IF (exit_flag.eq.NoError) THEN 
          IF ((char1.eq.'t').or.(char1.eq.'T')) THEN
            spherical=.TRUE.
          END IF
        ELSE
          WRITE (stdout,30) 'spherical', TRIM(ncname)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid variables.
!-----------------------------------------------------------------------
!
!  Set 2D arrays bounds.
!
      tile=-1
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Scan the variable list and read in needed variables.
!
      DO it=1,n_var
        SELECT CASE (TRIM(ADJUSTL(var_name(it))))
!
!  Read in basin X-length.
!
          CASE ('xl')
            CALL netcdf_get_fvar (ng, model, ncname, 'xl',              &
     &                            xl(ng),                               &
     &                            ncid = ncGRDid(ng))
            IF (exit_flag.ne.NoError) EXIT
!
!  Read in basin Y-length.
!
          CASE ('el')
            CALL netcdf_get_fvar (ng, model, ncname, 'el',              &
     &                            el(ng),                               &
     &                            ncid = ncGRDid(ng))
            IF (exit_flag.ne.NoError) EXIT
!
!  Read in bathymetry.
!
          CASE ('h')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, hmin(ng), hmax(ng),                 &
     &                        GRID(ng) % h)
            IF (status.eq.nf90_noerr) THEN
              DO j=LBj,UBj
                DO i=LBi,UBi
                  IF (GRID(ng)%h(i,j) .eq. 0.0) GRID(ng) % h(i,j) = -1.
                END DO
              END DO
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'h', TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % h)
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % f)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % f)
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pm)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % pm)
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pn)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % pn)
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % xp)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % yp)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % xr)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % yr)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % xu)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % yu)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % xv)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % yv)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % lonp)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % latp)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LonMin(ng), LonMax(ng),           &
     &                          GRID(ng) % lonr)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LatMin(ng), LatMax(ng),           &
     &                          GRID(ng) % latr)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % lonu)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % latu)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % lonv)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncname, ncGRDid(ng),         &
     &                          var_name(it), var_id(it),               &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % latv)
              IF (status.ne.nf90_noerr) THEN
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, ncGRDid(ng),           &
     &                        var_name(it), var_id(it),                 &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % angler)
            IF (status.ne.nf90_noerr) THEN
              exit_flag=2
              ioerror=status
              EXIT
            END IF
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              GRID(ng) % angler)
        END SELECT
      END DO
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,30) TRIM(var_name(it)), TRIM(ncname)
        RETURN
      END IF
!
! Close GRID NetCDF file.
!
      CALL netcdf_close (ng, model, ncGRDid(ng), ncname)
      IF (exit_flag.ne.NoError) RETURN
!
  10  FORMAT (/,' GET_GRID - unable to find grid variable: ',a,         &
     &        /,12x,'in grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID - unable to open grid NetCDF file: ',a)
  30  FORMAT (/,' GET_GRID - error while reading variable: ',a,         &
     &        /,12x,'in grid NetCDF file: ',a)
  40  FORMAT (/,' GET_GRID - Reading adjoint sensitivity scope arrays', &
     &        ' from file:',/12x,a)
      RETURN
      END SUBROUTINE get_grid
