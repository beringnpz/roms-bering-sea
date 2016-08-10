      SUBROUTINE close_io
!
!svn $Id: close_io.F 1020 2009-07-10 23:10:30Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This subroutine flushes and closes all IO files.                     !
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
      implicit none
!
!  Local variable declarations.
!
      logical :: First
      integer :: MyError, i, ng, status
!
      SourceFile='close_io.F'
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (ncRSTid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncRSTid(ng))
        END IF
        IF (ncHISid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncHISid(ng))
        END IF
        IF (ncAVGid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncAVGid(ng))
        END IF
        IF (ncAVG2id(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncAVG2id(ng))
        END IF
!
!  Report number of time records written.
!
        IF (Master) THEN
          WRITE (stdout,10) ng
          IF (NrecHIS(ng).gt.0) THEN
            WRITE (stdout,20) 'HISTORY', NrecHIS(ng)
          END IF
          IF (NrecRST(ng).gt.0) THEN
            IF (LcycleRST(ng)) THEN
              IF (NrecRST(ng).gt.1) THEN
                NrecRST(ng)=2
              ELSE
                NrecRST(ng)=1
              END IF
            END IF
            WRITE (stdout,20) 'RESTART', NrecRST(ng)
          END IF
          IF (NrecAVG(ng).gt.0) THEN
            WRITE (stdout,20) 'AVERAGE', NrecAVG(ng)
          END IF
          IF (NrecAVG2(ng).gt.0) THEN
            WRITE (stdout,20) 'AVERAGE', NrecAVG2(ng)
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report analytical header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,47
          IF (LEN_TRIM(ANANAME(i)).gt.0) THEN
            IF (First) THEN
              First=.FALSE.
              WRITE (stdout,30) ' Analytical header files used:'
            END IF
            WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(ANANAME(i)))
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Report biology model header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,4
          IF (LEN_TRIM(BIONAME(i)).gt.0) THEN
            IF (First) THEN
              First=.FALSE.
              WRITE (stdout,30) ' Biology model header files used:'
            END IF
            WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(BIONAME(i)))
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  If applicable, report internal exit errors.
!-----------------------------------------------------------------------
!
      IF (Master.and.(exit_flag.ne.NoError)) THEN
        WRITE (stdout,40) Rerror(exit_flag), exit_flag
      END IF
      IF (exit_flag.eq.NoError) THEN
        CALL get_date (date_str)
        IF (Master) WRITE (stdout,50) TRIM(date_str)
      ELSE IF ((exit_flag.eq.1).or.(blowup.ne.0)) THEN
        IF (Master) WRITE (stdout,60)
      ELSE IF (exit_flag.eq.2) THEN
        IF (Master) WRITE (stdout,70) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.3) THEN
        IF (Master) WRITE (stdout,80) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.4) THEN
        IF (Master) WRITE (stdout,90)
      ELSE IF (exit_flag.eq.5) THEN
        IF (Master) WRITE (stdout,100)
      ELSE IF (exit_flag.eq.6) THEN
        IF (Master) WRITE (stdout,110)
      ELSE IF (exit_flag.eq.7) THEN
        IF (Master) WRITE (stdout,120)
      ELSE IF (exit_flag.eq.8) THEN
        IF (Master) WRITE (stdout,130)
      END IF
!
 10   FORMAT (/,' ROMS/TOMS - Output NetCDF summary for Grid ',         &
     &        i2.2,':')
 20   FORMAT (13x,'number of time records written in ',                 &
     &        a,' file = ',i8.8)
 30   FORMAT (/,a,/)
 40   FORMAT (/,a,i3,/)
 50   FORMAT (/,' ROMS/TOMS: DONE... ',a)
 60   FORMAT (/,' MAIN: Abnormal termination: BLOWUP.')
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.',/,        &
     &          ' REASON: ',a)
 80   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.',/,       &
     &          ' REASON: ',a)
 90   FORMAT (/,' ERROR: I/O related problem.')
100   FORMAT (/,' ERROR: Illegal model configuration.')
110   FORMAT (/,' ERROR: Illegal domain partition.')
120   FORMAT (/,' ERROR: Illegal input parameter.')
130   FORMAT (/,' ERROR: Fatal algorithm result.')
      RETURN
      END SUBROUTINE close_io
