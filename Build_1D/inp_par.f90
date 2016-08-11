      SUBROUTINE inp_par (model)
!     
!svn $Id: inp_par.F 966 2009-03-28 00:29:06Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in input model parameters from standard input.   !
!  It also writes out these parameters to standard output.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE ran_state, ONLY: ran_seed
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      logical :: Lwrite
      integer :: Itile, Jtile, Nghost, Ntiles, tile
      integer :: Imin, Imax, Jmin, Jmax
      integer :: inp, out, itrc, ng, npts, sequence
      real(r8) :: cff
      real(r8), parameter :: spv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
!
!  Set input units.
!
      Lwrite=Master
      inp=stdinp
      out=stdout
!
!  Get current data.
!
      CALL get_date (date_str)
!
!-----------------------------------------------------------------------
!  Read in physical model input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) WRITE (out,10) version, TRIM(date_str)
 10   FORMAT (/,' Model Input Parameters:  ROMS/TOMS version ',a,/,     &
     &       26x,a,/,1x,77('-'))
      CALL read_PhyPar (model, inp, out, Lwrite)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set lower and upper bounds indices per domain partition for all
!  nested grids.
!-----------------------------------------------------------------------
!
!  Allocate structure.
!
      IF (.not.allocated(BOUNDS)) THEN
        allocate ( BOUNDS(Ngrids) )
        DO ng=1,Ngrids
          Ntiles=NtileI(ng)*NtileJ(ng)-1
          allocate ( BOUNDS(ng) % tile (-1:Ntiles) )
          allocate ( BOUNDS(ng) % LBi  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % UBi  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % LBj  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % UBj  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iend (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Istr (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jstr (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jend (-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrV(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Imin (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Imax (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Jmin (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Jmax (4,0:1,0:Ntiles) )
        END DO
      END IF
!
!  Set boundary edge I- or J-indices for each variable type.
!
      DO ng=1,Ngrids
        BOUNDS(ng) % edge(iwest ,p2dvar) = 1
        BOUNDS(ng) % edge(iwest ,r2dvar) = 0
        BOUNDS(ng) % edge(iwest ,u2dvar) = 1
        BOUNDS(ng) % edge(iwest ,v2dvar) = 0
        BOUNDS(ng) % edge(ieast ,p2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(ieast ,r2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(ieast ,u2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(ieast ,v2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(isouth,p2dvar) = 1
        BOUNDS(ng) % edge(isouth,u2dvar) = 0
        BOUNDS(ng) % edge(isouth,r2dvar) = 0
        BOUNDS(ng) % edge(isouth,v2dvar) = 1
        BOUNDS(ng) % edge(inorth,p2dvar) = Mm(ng)+1
        BOUNDS(ng) % edge(inorth,r2dvar) = Mm(ng)+1
        BOUNDS(ng) % edge(inorth,u2dvar) = Mm(ng)+1
        BOUNDS(ng) % edge(inorth,v2dvar) = Mm(ng)+1
      END DO
!
!  Set tile computational indices and arrays allocation bounds.
!
      Nghost=2
      DO ng=1,Ngrids
        BOUNDS(ng) % LBij = 0
        BOUNDS(ng) % UBij = MAX(Lm(ng)+1,Mm(ng)+1)
        DO tile=-1,NtileI(ng)*NtileJ(ng)-1
          BOUNDS(ng) % tile(tile) = tile
          CALL get_tile (ng, tile, Itile, Jtile,                        &
     &                   BOUNDS(ng) % Istr(tile),                       &
     &                   BOUNDS(ng) % Iend(tile),                       &
     &                   BOUNDS(ng) % Jstr(tile),                       &
     &                   BOUNDS(ng) % Jend(tile),                       &
     &                   BOUNDS(ng) % IstrR(tile),                      &
     &                   BOUNDS(ng) % IstrT(tile),                      &
     &                   BOUNDS(ng) % IstrU(tile),                      &
     &                   BOUNDS(ng) % IendR(tile),                      &
     &                   BOUNDS(ng) % IendT(tile),                      &
     &                   BOUNDS(ng) % JstrR(tile),                      &
     &                   BOUNDS(ng) % JstrT(tile),                      &
     &                   BOUNDS(ng) % JstrV(tile),                      &
     &                   BOUNDS(ng) % JendR(tile),                      &
     &                   BOUNDS(ng) % JendT(tile))
          CALL get_bounds (ng, tile, 0, Nghost, Itile, Jtile,           &
     &                     BOUNDS(ng) % LBi(tile),                      &
     &                     BOUNDS(ng) % UBi(tile),                      &
     &                     BOUNDS(ng) % LBj(tile),                      &
     &                     BOUNDS(ng) % UBj(tile))
        END DO
      END DO
!
!  Set I/O processing minimum (Imin, Jmax) and maximum (Imax, Jmax)
!  indices for non-overlapping (Nghost=0) and overlapping (Nghost>0)
!  tiles for each C-grid type variable.
!
      Nghost=2
      DO ng=1,Ngrids
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_bounds (ng, tile, p2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(1,0,tile),                 &
     &                     BOUNDS(ng) % Imax(1,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(1,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(1,0,tile))
          CALL get_bounds (ng, tile, p2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(1,1,tile),                 &
     &                     BOUNDS(ng) % Imax(1,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(1,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(1,1,tile))
          CALL get_bounds (ng, tile, r2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(2,0,tile),                 &
     &                     BOUNDS(ng) % Imax(2,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(2,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(2,0,tile))
          CALL get_bounds (ng, tile, r2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(2,1,tile),                 &
     &                     BOUNDS(ng) % Imax(2,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(2,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(2,1,tile))
          CALL get_bounds (ng, tile, u2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(3,0,tile),                 &
     &                     BOUNDS(ng) % Imax(3,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(3,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(3,0,tile))
          CALL get_bounds (ng, tile, u2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(3,1,tile),                 &
     &                     BOUNDS(ng) % Imax(3,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(3,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(3,1,tile))
          CALL get_bounds (ng, tile, v2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(4,0,tile),                 &
     &                     BOUNDS(ng) % Imax(4,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(4,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(4,0,tile))
          CALL get_bounds (ng, tile, v2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(4,1,tile),                 &
     &                     BOUNDS(ng) % Imax(4,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(4,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(4,1,tile))
        END DO
      END DO
!
!  Set NetCDF IO bounds.
!
      DO ng=1,Ngrids
        CALL get_iobounds (ng)
      END DO
!
!-----------------------------------------------------------------------
!  Check tile partition starting and ending (I,J) indices for illegal
!  domain decomposition parameters NtileI and NtileJ in standard input
!  file.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        DO ng=1,Ngrids
          WRITE (stdout,50) ng, Lm(ng), Mm(ng), N(ng),                  &
     &                      NtileI(ng), NtileJ(ng)
          Print *,'Number of tracers: ', NT(ng)
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            npts=(BOUNDS(ng)%Iend(tile)-                                &
     &            BOUNDS(ng)%Istr(tile)+1)*                             &
     &           (BOUNDS(ng)%Jend(tile)-                                &
     &            BOUNDS(ng)%Jstr(tile)+1)*N(ng)
            WRITE (stdout,70) tile,                                     &
     &                        BOUNDS(ng)%Istr(tile),                    &
     &                        BOUNDS(ng)%Iend(tile),                    &
     &                        BOUNDS(ng)%Jstr(tile),                    &
     &                        BOUNDS(ng)%Jend(tile),                    &
     &                        npts
            IF ((BOUNDS(ng)%Iend(tile)-                                 &
     &           BOUNDS(ng)%Istr(tile)+1).lt.2) THEN
              WRITE (stdout,80) ng, 'NtileI = ', NtileI(ng),            &
     &                              'Lm = ', Lm(ng),                    &
     &                              'Istr = ', BOUNDS(ng)%Istr(tile),   &
     &                              '  Iend = ', BOUNDS(ng)%Iend(tile), &
     &                              'NtileI'
              exit_flag=6
              RETURN
            END IF
            IF ((BOUNDS(ng)%Jend(tile)-                                 &
     &           BOUNDS(ng)%Jstr(tile)+1).lt.2) THEN
              WRITE (stdout,80) ng, 'NtileJ = ', NtileJ(ng),            &
     &                              'Mm = ', Mm(ng),                    &
     &                              'Jstr = ', BOUNDS(ng)%Jstr(tile),   &
     &                              '  Jend = ', BOUNDS(ng)%Jend(tile), &
     &                              'NtileJ'
              exit_flag=6
              RETURN
            END IF
          END DO
        END DO
 50     FORMAT (/,' Tile partition information for Grid ',i2.2,':',2x,  &
     &          i4.4,'x',i4.4,'x',i4.4,2x,'tiling: ',i3.3,'x',i3.3,/,/, &
     &          5x,'tile',5x,'Istr',5x,'Iend',5x,'Jstr',5x,'Jend',      &
     &          5x,'Npts',/)
 70     FORMAT (5(5x,i4),2x,i7)
 80     FORMAT (/,' INP_PAR - domain decomposition error in input ',    &
     &                        'script file for grid: ',i2,/,            &
     &          /,11x,'The domain partition parameter, ',a,i3,          &
     &          /,11x,'is incompatible with grid size, ',a,i4,          &
     &          /,11x,'because it yields too small tile, ',a,i3,a,i3,   &
     &          /,11x,'Decrease partition parameter: ',a)
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
!  Report tile minimum and maximum fractional grid coordinates.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,90)
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            WRITE (stdout,100) tile,                                    &
     &                         DOMAIN(ng)%Xmin_rho(tile),               &
     &                         DOMAIN(ng)%Xmax_rho(tile),               &
     &                         DOMAIN(ng)%Ymin_rho(tile),               &
     &                         DOMAIN(ng)%Ymax_rho(tile), 'RHO-points'
          END DO
          WRITE (stdout,'(1x)')
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            WRITE (stdout,100) tile,                                    &
     &                         DOMAIN(ng)%Xmin_u(tile),                 &
     &                         DOMAIN(ng)%Xmax_u(tile),                 &
     &                         DOMAIN(ng)%Ymin_u(tile),                 &
     &                         DOMAIN(ng)%Ymax_u(tile), '  U-points'
          END DO
          WRITE (stdout,'(1x)')
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            WRITE (stdout,100) tile,                                    &
     &                         DOMAIN(ng)%Xmin_v(tile),                 &
     &                         DOMAIN(ng)%Xmax_v(tile),                 &
     &                         DOMAIN(ng)%Ymin_v(tile),                 &
     &                         DOMAIN(ng)%Ymax_v(tile), '  V-points'
          END DO
 90       FORMAT (/,' Tile minimum and maximum fractional grid ',       &
     &            'coordinates:',/,                                     &
     &            '   (interior points only)',/,/,                      &
     &            5x,'tile',5x,'Xmin',5x,'Xmax',5x,'Ymin',5x,'Ymax',    &
     &            5x,'grid',/)
 100      FORMAT (5x,i4,4f9.2,2x,a)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Read in biological model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (15, FILE=TRIM(bparnam), FORM='formatted', STATUS='old')
      CALL read_BioPar (model, 15, out, Lwrite)
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options and definitions.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        CALL checkdefs
        CALL my_flush (out)
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Compute various constants.
!-----------------------------------------------------------------------
!
      gorho0=g/rho0
      DO ng=1,Ngrids
        dtfast(ng)=dt(ng)/REAL(ndtfast(ng),r8)
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
        visc4(ng)=SQRT(ABS(visc4(ng)))
        tkenu4(ng)=SQRT(ABS(tkenu4(ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
        IF (Znudg(ng).gt.0.0_r8) THEN
          Znudg(ng)=1.0_r8/(Znudg(ng)*86400.0_r8)
        ELSE
          Znudg(ng)=0.0_r8
        END IF
        IF (M2nudg(ng).gt.0.0_r8) THEN
          M2nudg(ng)=1.0_r8/(M2nudg(ng)*86400.0_r8)
        ELSE
          M2nudg(ng)=0.0_r8
        END IF
        IF (M3nudg(ng).gt.0.0_r8) THEN
          M3nudg(ng)=1.0_r8/(M3nudg(ng)*86400.0_r8)
        ELSE
          M3nudg(ng)=0.0_r8
        END IF
!
!  Convert momentum stresses and tracer flux scales to kinematic
!  Values. Recall, that all the model fluxes are kinematic.
!
        cff=1.0_r8/rho0
        Fscale(idUsms,ng)=cff*Fscale(idUsms,ng)
        Fscale(idVsms,ng)=cff*Fscale(idVsms,ng)
        Fscale(idUbms,ng)=cff*Fscale(idUbms,ng)
        Fscale(idVbms,ng)=cff*Fscale(idVbms,ng)
        Fscale(idUbrs,ng)=cff*Fscale(idUbrs,ng)
        Fscale(idVbrs,ng)=cff*Fscale(idVbrs,ng)
        Fscale(idUbws,ng)=cff*Fscale(idUbws,ng)
        Fscale(idVbws,ng)=cff*Fscale(idVbws,ng)
        Fscale(idUbcs,ng)=cff*Fscale(idUbcs,ng)
        Fscale(idVbcs,ng)=cff*Fscale(idVbcs,ng)
        cff=1.0_r8/(rho0*Cp)
        Fscale(idTsur(itemp),ng)=cff*Fscale(idTsur(itemp),ng)
        Fscale(idTbot(itemp),ng)=cff*Fscale(idTbot(itemp),ng)
        Fscale(idSrad,ng)=cff*Fscale(idSrad,ng)
        Fscale(idLdwn,ng)=cff*Fscale(idLdwn,ng)
        Fscale(idLrad,ng)=cff*Fscale(idLrad,ng)
        Fscale(idLhea,ng)=cff*Fscale(idLhea,ng)
        Fscale(idShea,ng)=cff*Fscale(idShea,ng)
        Fscale(iddQdT,ng)=cff*Fscale(iddQdT,ng)
      END DO
!
!-----------------------------------------------------------------------
!  Initialize random number sequence so we can get identical results
!  everytime that we run the same solution.
!-----------------------------------------------------------------------
!
      sequence=759
      CALL ran_seed (sequence)
      RETURN
      END SUBROUTINE inp_par
      SUBROUTINE read_PhyPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in physical model input parameters.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: inhere
      integer :: Lstr, Npts, Nval, i, itrc, k, ng, status
      integer :: decode_line, load_i, load_l, load_r
      logical, allocatable :: Ltracer(:,:)
      real(r8), allocatable :: Rtracer(:,:)
      real(r8), allocatable :: tracer(:,:)
      real(r8), dimension(300) :: Rval
      character (len=1 ), parameter :: blank = ' '
      character (len=19) :: ref_att
      character (len=40) :: KeyWord
      character (len=160) :: fname, line
      character (len=160), dimension(300) :: Cval
!
!-----------------------------------------------------------------------
!  Read in physical model parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'TITLE') THEN
            IF (Nval.eq.1) THEN
              title=TRIM(ADJUSTL(Cval(Nval)))
            ELSE
              WRITE(title,'(a,1x,a)') TRIM(ADJUSTL(title)),             &
     &                                TRIM(ADJUSTL(Cval(Nval)))
            END IF           
          ELSE IF (TRIM(KeyWord).eq.'MyAppCPP') THEN
            DO i=1,LEN(MyAppCPP)
              MyAppCPP(i:i)=blank
            END DO
            MyAppCPP=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'VARNAME') THEN
            DO i=1,LEN(varname)
              varname(i:i)=blank
            END DO
            varname=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'Lm') THEN
            Npts=load_i(Nval, Rval, Ngrids, Lm)
            DO ng=1,Ngrids
              IF (Lm(ng).le.0) THEN
                IF (Master) WRITE (out,300) 'Lm', ng,                   &
     &                                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Mm') THEN
            Npts=load_i(Nval, Rval, Ngrids, Mm)
            DO ng=1,Ngrids
              IF (Mm(ng).le.0) THEN
                IF (Master) WRITE (out,300) 'Mm', ng,                   &
     &                                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'N') THEN
            Npts=load_i(Nval, Rval, Ngrids, N)
            DO ng=1,Ngrids
              IF (N(ng).lt.0) THEN
                IF (Master) WRITE (out,300) 'N', ng,                    &
     &                                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'NAT') THEN
            Npts=load_i(Nval, Rval, 1, NAT)
            IF ((NAT.lt.1).or.(NAT.gt.2)) THEN
              IF (Master) WRITE (out,290) 'NAT = ', NAT,                &
     &                          'make sure that NAT is either 1 or 2.'
              exit_flag=5
              RETURN
            END IF
            IF (NAT.ne.2) THEN
              IF (Master) WRITE (out,290) 'NAT = ', NAT,                &
     &                          'make sure that NAT is equal to 2.'
              exit_flag=5
              RETURN
            END IF
          ELSE IF (TRIM(KeyWord).eq.'NtileI') THEN
            Npts=load_i(Nval, Rval, Ngrids, NtileI)
            NtileX(1:Ngrids)=NtileI(1:Ngrids)
          ELSE IF (TRIM(KeyWord).eq.'NtileJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, NtileJ)
            NtileE(1:Ngrids)=NtileJ(1:Ngrids)
            CALL initialize_param
            CALL initialize_scalars
            CALL initialize_ncparam
            IF (.not.allocated(Ltracer)) THEN
              allocate (Ltracer(NAT+NPT,Ngrids))
            END IF
            IF (.not.allocated(Rtracer)) THEN
              allocate (Rtracer(NAT+NPT,Ngrids))
            END IF
            IF (.not.allocated(tracer)) THEN
              allocate (tracer(MT,Ngrids))
            END IF
          ELSE IF (TRIM(KeyWord).eq.'NTIMES') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntimes)
          ELSE IF (TRIM(KeyWord).eq.'DT') THEN
            Npts=load_r(Nval, Rval, Ngrids, dt)
          ELSE IF (TRIM(KeyWord).eq.'NDTFAST') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndtfast)
          ELSE IF (TRIM(KeyWord).eq.'ERstr') THEN
            Npts=load_i(Nval, Rval, 1, ERstr)
          ELSE IF (TRIM(KeyWord).eq.'ERend') THEN
            Npts=load_i(Nval, Rval, 1, ERend)
          ELSE IF (TRIM(KeyWord).eq.'Nouter') THEN
            Npts=load_i(Nval, Rval, 1, Nouter)
          ELSE IF (TRIM(KeyWord).eq.'Ninner') THEN
            Npts=load_i(Nval, Rval, 1, Ninner)
          ELSE IF (TRIM(KeyWord).eq.'Nintervals') THEN
            Npts=load_i(Nval, Rval, 1, Nintervals)
          ELSE IF (TRIM(KeyWord).eq.'NRREC') THEN
            Npts=load_i(Nval, Rval, Ngrids, nrrec)
            DO ng=1,Ngrids
              IF (nrrec(ng).lt.0) THEN
                LastRec(ng)=.TRUE.
              ELSE
                LastRec(ng)=.FALSE.
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'LcycleRST') THEN
            Npts=load_l(Nval, Cval, Ngrids, LcycleRST)
          ELSE IF (TRIM(KeyWord).eq.'NRST') THEN
            Npts=load_i(Nval, Rval, Ngrids, nRST)
          ELSE IF (TRIM(KeyWord).eq.'NSTA') THEN
            Npts=load_i(Nval, Rval, Ngrids, nSTA)
          ELSE IF (TRIM(KeyWord).eq.'NFLT') THEN
            Npts=load_i(Nval, Rval, Ngrids, nFLT)
          ELSE IF (TRIM(KeyWord).eq.'NINFO') THEN	 
            Npts=load_i(Nval, Rval, Ngrids, ninfo)
          ELSE IF (TRIM(KeyWord).eq.'LDEFOUT') THEN
            Npts=load_l(Nval, Cval, Ngrids, ldefout)
          ELSE IF (TRIM(KeyWord).eq.'NHIS') THEN
            Npts=load_i(Nval, Rval, Ngrids, nHIS)
          ELSE IF (TRIM(KeyWord).eq.'NDEFHIS') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefHIS)
          ELSE IF (TRIM(KeyWord).eq.'NTSAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntsAVG)
          ELSE IF (TRIM(KeyWord).eq.'NAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, nAVG)
          ELSE IF (TRIM(KeyWord).eq.'NDEFAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefAVG)
          ELSE IF (TRIM(KeyWord).eq.'NTSAVG2') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntsAVG2)
          ELSE IF (TRIM(KeyWord).eq.'NAVG2') THEN
            Npts=load_i(Nval, Rval, Ngrids, nAVG2)
          ELSE IF (TRIM(KeyWord).eq.'NDEFAVG2') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefAVG2)
          ELSE IF (TRIM(KeyWord).eq.'NTSDIA') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntsDIA)
          ELSE IF (TRIM(KeyWord).eq.'NDIA') THEN
            Npts=load_i(Nval, Rval, Ngrids, nDIA)
          ELSE IF (TRIM(KeyWord).eq.'NDEFDIA') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefDIA)
          ELSE IF (TRIM(KeyWord).eq.'LcycleTLM') THEN
            Npts=load_l(Nval, Cval, Ngrids, LcycleTLM)
          ELSE IF (TRIM(KeyWord).eq.'NTLM') THEN
            Npts=load_i(Nval, Rval, Ngrids, nTLM)
          ELSE IF (TRIM(KeyWord).eq.'NDEFTLM') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefTLM)
          ELSE IF (TRIM(KeyWord).eq.'LcycleADJ') THEN
            Npts=load_l(Nval, Cval, Ngrids, LcycleADJ)
          ELSE IF (TRIM(KeyWord).eq.'NADJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, nADJ)
          ELSE IF (TRIM(KeyWord).eq.'NDEFADJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefADJ)
          ELSE IF (TRIM(KeyWord).eq.'NOBC') THEN
            Npts=load_i(Nval, Rval, Ngrids, nOBC)
          ELSE IF (TRIM(KeyWord).eq.'NSFF') THEN
            Npts=load_i(Nval, Rval, Ngrids, nSFF)
          ELSE IF (TRIM(KeyWord).eq.'LrstGST') THEN
            Npts=load_l(Nval, Cval, 1, LrstGST)
          ELSE IF (TRIM(KeyWord).eq.'MaxIterGST') THEN
            Npts=load_i(Nval, Rval, 1, MaxIterGST)
          ELSE IF (TRIM(KeyWord).eq.'NGST') THEN
            Npts=load_i(Nval, Rval, 1, nGST)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                tnu2(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                tnu4(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VISC2') THEN
            Npts=load_r(Nval, Rval, Ngrids, visc2)
          ELSE IF (TRIM(KeyWord).eq.'VISC4') THEN
            Npts=load_r(Nval, Rval, Ngrids, visc4)
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                Akt_bak(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKV_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
          ELSE IF (TRIM(KeyWord).eq.'AKK_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akk_bak)
          ELSE IF (TRIM(KeyWord).eq.'AKP_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akp_bak)
          ELSE IF (TRIM(KeyWord).eq.'TKENU2') THEN
            Npts=load_r(Nval, Rval, Ngrids, tkenu2)
          ELSE IF (TRIM(KeyWord).eq.'TKENU4') THEN
            Npts=load_r(Nval, Rval, Ngrids, tkenu4)
          ELSE IF (TRIM(KeyWord).eq.'GLS_P') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_p)
          ELSE IF (TRIM(KeyWord).eq.'GLS_M') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_m)
          ELSE IF (TRIM(KeyWord).eq.'GLS_N') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_n)
          ELSE IF (TRIM(KeyWord).eq.'GLS_Kmin') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_Kmin)
          ELSE IF (TRIM(KeyWord).eq.'GLS_Pmin') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_Pmin)
          ELSE IF (TRIM(KeyWord).eq.'GLS_CMU0') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_cmu0)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C1') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c1)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C2') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c2)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C3M') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c3m)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C3P') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c3p)
          ELSE IF (TRIM(KeyWord).eq.'GLS_SIGK') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_sigk)
          ELSE IF (TRIM(KeyWord).eq.'GLS_SIGP') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_sigp)
          ELSE IF (TRIM(KeyWord).eq.'CHARNOK_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, charnok_alpha)
          ELSE IF (TRIM(KeyWord).eq.'ZOS_HSIG_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, zos_hsig_alpha)
          ELSE IF (TRIM(KeyWord).eq.'SZ_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, sz_alpha)
          ELSE IF (TRIM(KeyWord).eq.'CRGBAN_CW') THEN
            Npts=load_r(Nval, Rval, Ngrids, crgban_cw)
          ELSE IF (TRIM(KeyWord).eq.'RDRG') THEN
            Npts=load_r(Nval, Rval, Ngrids, rdrg)
          ELSE IF (TRIM(KeyWord).eq.'RDRG2') THEN
            Npts=load_r(Nval, Rval, Ngrids, rdrg2)
          ELSE IF (TRIM(KeyWord).eq.'Zob') THEN
            Npts=load_r(Nval, Rval, Ngrids, Zob)
          ELSE IF (TRIM(KeyWord).eq.'Zos') THEN
            Npts=load_r(Nval, Rval, Ngrids, Zos)
          ELSE IF (TRIM(KeyWord).eq.'BLK_ZQ') THEN
            Npts=load_r(Nval, Rval, Ngrids, blk_ZQ)
          ELSE IF (TRIM(KeyWord).eq.'BLK_ZT') THEN
            Npts=load_r(Nval, Rval, Ngrids, blk_ZT)
          ELSE IF (TRIM(KeyWord).eq.'BLK_ZW') THEN
            Npts=load_r(Nval, Rval, Ngrids, blk_ZW)
          ELSE IF (TRIM(KeyWord).eq.'DCRIT') THEN
            Npts=load_r(Nval, Rval, Ngrids, Dcrit)
          ELSE IF (TRIM(KeyWord).eq.'WTYPE') THEN
            Npts=load_i(Nval, Rval, Ngrids, lmd_Jwt)
          ELSE IF (TRIM(KeyWord).eq.'LEVSFRC') THEN
            Npts=load_i(Nval, Rval, Ngrids, levsfrc)
          ELSE IF (TRIM(KeyWord).eq.'LEVBFRC') THEN
            Npts=load_i(Nval, Rval, Ngrids, levbfrc)
          ELSE IF (TRIM(KeyWord).eq.'Vtransform') THEN
            Npts=load_i(Nval, Rval, Ngrids, Vtransform)
            DO ng=1,Ngrids
              IF ((Vtransform(ng).lt.0).or.                             &
     &            (Vtransform(ng).gt.3)) THEN
                IF (Master) WRITE (out,260) 'Vtransform = ',            &
     &                                      Vtransform(ng),             &
     &                                      'Must be either 1, 2 or 3'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Vstretching') THEN
            Npts=load_i(Nval, Rval, Ngrids, Vstretching)
            DO ng=1,Ngrids
              IF ((Vstretching(ng).lt.0).or.                            &
     &            (Vstretching(ng).gt.3)) THEN
                IF (Master) WRITE (out,260) 'Vstretching = ',           &
     &                                      Vstretching(ng),            &
     &                                      'Must between 1 and 3'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'THETA_S') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_s)
          ELSE IF (TRIM(KeyWord).eq.'THETA_B') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_b)
          ELSE IF (TRIM(KeyWord).eq.'TCLINE') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcline)
            DO ng=1,Ngrids
              hc(ng)=Tcline(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'RHO0') THEN
            Npts=load_r(Nval, Rval, 1, rho0)
          ELSE IF (TRIM(KeyWord).eq.'BVF_BAK') THEN
            Npts=load_r(Nval, Rval, 1, bvf_bak)
          ELSE IF (TRIM(KeyWord).eq.'DSTART') THEN
            Npts=load_r(Nval, Rval, 1, dstart)
          ELSE IF (TRIM(KeyWord).eq.'TIDE_START') THEN
            Npts=load_r(Nval, Rval, 1, tide_start)
          ELSE IF (TRIM(KeyWord).eq.'TIME_REF') THEN
            Npts=load_r(Nval, Rval, 1, time_ref)
            r_text=ref_att(time_ref,r_date)
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                Tnudg(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZNUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, Znudg)
          ELSE IF (TRIM(KeyWord).eq.'M2NUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, M2nudg)
          ELSE IF (TRIM(KeyWord).eq.'M3NUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, M3nudg)
          ELSE IF (TRIM(KeyWord).eq.'OBCFAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, obcfac)
          ELSE IF (TRIM(KeyWord).eq.'R0') THEN
            Npts=load_r(Nval, Rval, Ngrids, R0)
            DO ng=1,Ngrids
              IF (R0(ng).lt.100.0_r8) R0(ng)=R0(ng)+1000.0_r8
            END DO
          ELSE IF (TRIM(KeyWord).eq.'T0') THEN
            Npts=load_r(Nval, Rval, Ngrids, T0)
          ELSE IF (TRIM(KeyWord).eq.'S0') THEN
            Npts=load_r(Nval, Rval, Ngrids, S0)
          ELSE IF (TRIM(KeyWord).eq.'TCOEF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcoef)
            DO ng=1,Ngrids
              Tcoef(ng)=ABS(Tcoef(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SCOEF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Scoef)
            DO ng=1,Ngrids
              Scoef(ng)=ABS(Scoef(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'GAMMA2') THEN
            Npts=load_r(Nval, Rval, Ngrids, gamma2)
          ELSE IF (TRIM(KeyWord).eq.'Hout(idFsur)') THEN
            IF (idFsur.eq.0) THEN
              IF (Master) WRITE (out,280) 'idFsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idFsur,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbar)') THEN
            IF (idUbar.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbar,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbar)') THEN
            IF (idVbar.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbar,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUvel)') THEN
            IF (idUvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvel)') THEN
            IF (idVvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWvel)') THEN
            IF (idWvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idOvel)') THEN
            IF (idOvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idOvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idOvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            IF (MAXVAL(idTvar).eq.0) THEN
              IF (Master) WRITE (out,280) 'idTvar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                i=idTvar(itrc)
                Hout(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUair)') THEN
            IF (idUair.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUair'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUair,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVair)') THEN
            IF (idVair.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVair'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVair,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUsms)') THEN
            IF (idUsms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUsms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUsms,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVsms)') THEN
            IF (idVsms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVsms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVsms,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbms)') THEN
            IF (idUbms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbms,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbms)') THEN
            IF (idVbms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbms,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbrs)') THEN
            IF (idUbrs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbrs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbrs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbrs)') THEN
            IF (idVbrs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbrs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbrs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbws)') THEN
            IF (idUbws.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbws'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbws,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbws)') THEN
            IF (idVbws.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbws'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbws,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbcs)') THEN
            IF (idUbcs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbcs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbcs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbcs)') THEN
            IF (idVbcs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbcs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbcs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbot)') THEN
            IF (idUbot.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbot'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbot,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbot)') THEN
            IF (idVbot.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbot'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbot,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbur)') THEN
            IF (idUbur.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbur,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbvr)') THEN
            IF (idVbvr.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbvr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbvr,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2xx)') THEN
            IF (idW2xx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW2xx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xx,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2xy)') THEN
            IF (idW2xy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW2xy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xy,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2yy)') THEN
            IF (idW2yy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW2yy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2yy,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3xx)') THEN
            IF (idW3xx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3xx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xx,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3xy)') THEN
            IF (idW3xy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3xy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xy,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3yy)') THEN
            IF (idW3yy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3yy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3yy,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3zx)') THEN
            IF (idW3zx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3zx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zx,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3zy)') THEN
            IF (idW3zy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3zy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zy,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU2rs)') THEN
            IF (idU2rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU2rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU2rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV2rs)') THEN
            IF (idV2rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV2rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV2rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU2Sd)') THEN
            IF (idU2Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU2Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU2Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV2Sd)') THEN
            IF (idV2Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV2Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV2Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU3rs)') THEN
            IF (idU3rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU3rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU3rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV3rs)') THEN
            IF (idV3rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV3rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV3rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU3Sd)') THEN
            IF (idU3Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU3Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU3Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV3Sd)') THEN
            IF (idV3Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV3Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV3Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWamp)') THEN
            IF (idWamp.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWamp'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWamp,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWlen)') THEN
            IF (idWlen.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWlen'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWlen,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdir)') THEN
            IF (idWdir.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWdir'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdir,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            IF (idTsur(itemp).eq.0) THEN
              IF (Master) WRITE (out,280) 'idTsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                i=idTsur(itrc)
                Hout(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLhea)') THEN
            IF (idLhea.eq.0) THEN
              IF (Master) WRITE (out,280) 'idLhea'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLhea,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idShea)') THEN
            IF (idShea.eq.0) THEN
              IF (Master) WRITE (out,280) 'idShea'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idShea,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLrad)') THEN
            IF (idLrad.eq.0) THEN
              IF (Master) WRITE (out,280) 'idLrad'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLrad,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSrad)') THEN
            IF (idSrad.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSrad'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idevap)') THEN
            IF (idevap.eq.0) THEN
              IF (Master) WRITE (out,280) 'idevap'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idevap,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idrain)') THEN
            IF (idrain.eq.0) THEN
              IF (Master) WRITE (out,280) 'idrain'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idrain,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idDano)') THEN
            IF (idDano.eq.0) THEN
              IF (Master) WRITE (out,280) 'idDano'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idDano,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvis)') THEN
            IF (idVvis.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVvis'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvis,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTdif)') THEN
            IF (idTdif.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTdif'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTdif,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSdif)') THEN
            IF (idSdif.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSdif'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSdif,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHsbl)') THEN
            IF (idHsbl.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHsbl'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHsbl,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHbbl)') THEN
            IF (idHbbl.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHbbl'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHbbl,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtke)') THEN
            IF (idMtke.eq.0) THEN
              IF (Master) WRITE (out,280) 'idMtke'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtke,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtls)') THEN
            IF (idMtls.eq.0) THEN
              IF (Master) WRITE (out,280) 'idMtls'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtls,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idFsur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idFsur,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idUbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idUbar,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idVbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idVbar,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idUvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idVvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idTvar)') THEN
            Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                i=idTvar(itrc)
                Hout2(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idUair)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idUair,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idVair)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idVair,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idUsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idUsms,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idVsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idVsms,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idWamp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idWamp,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idWlen)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idWlen,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idWdir)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idWdir,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idTsur)') THEN
            Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                i=idTsur(itrc)
                Hout2(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idLhea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idLhea,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idShea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idShea,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idLrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idLrad,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idSrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idSrad,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idevap)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idevap,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idrain)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idrain,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idDano)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idDano,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idHsbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idHsbl,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout2(idHbbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Hout2(idHbbl,:))
          ELSE IF (TRIM(KeyWord).eq.'NUSER') THEN
            Npts=load_i(Nval, Rval, 1, Nuser)
          ELSE IF (TRIM(KeyWord).eq.'USER') THEN
            Npts=load_r(Nval, Rval, MAX(1,Nuser), user)
          ELSE IF (TRIM(KeyWord).eq.'NC_SHUFFLE') THEN
            Npts=load_i(Nval, Rval, 1, shuffle)
          ELSE IF (TRIM(KeyWord).eq.'NC_DEFLATE') THEN
            Npts=load_i(Nval, Rval, 1, deflate)
          ELSE IF (TRIM(KeyWord).eq.'NC_DLEVEL') THEN
            Npts=load_i(Nval, Rval, 1, deflate_level)
          ELSE IF (TRIM(KeyWord).eq.'GSTNAME') THEN
            DO i=1,LEN(GSTname(Nval))
              GSTname(Nval)(i:i)=blank
            END DO
            GSTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'RSTNAME') THEN
            DO i=1,LEN(RSTname(Nval))
              RSTname(Nval)(i:i)=blank
            END DO
            RSTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'HISNAME') THEN
            DO i=1,LEN(HISname(Nval))
              HISname(Nval)(i:i)=blank
              HISbase(Nval)(i:i)=blank
            END DO
            HISname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            HISbase(Nval)=TRIM(ADJUSTL(HISname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'TLMNAME') THEN
            DO i=1,LEN(TLMname(Nval))
              TLMname(Nval)(i:i)=blank
              TLMbase(Nval)(i:i)=blank
            END DO
            TLMname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            TLMbase(Nval)=TRIM(ADJUSTL(TLMname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'TLFNAME') THEN
            DO i=1,LEN(TLMname(Nval))
              TLFname(Nval)(i:i)=blank
            END DO
            TLFname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ADJNAME') THEN
            DO i=1,LEN(ADJname(Nval))
              ADJname(Nval)(i:i)=blank
              ADJbase(Nval)(i:i)=blank
            END DO
            ADJname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            ADJbase(Nval)=TRIM(ADJUSTL(ADJname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'AVGNAME') THEN
            DO i=1,LEN(AVGname(Nval))
              AVGname(Nval)(i:i)=blank
              AVGbase(Nval)(i:i)=blank
            END DO
            AVGname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            AVGbase(Nval)=TRIM(ADJUSTL(AVGname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'AVG2NAME') THEN
            DO i=1,LEN(AVG2name(Nval))
              AVG2name(Nval)(i:i)=blank
              AVG2base(Nval)(i:i)=blank
            END DO
            AVG2name(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            AVG2base(Nval)=TRIM(ADJUSTL(AVG2name(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'DIANAME') THEN
            DO i=1,LEN(DIAname(Nval))
              DIAname(Nval)(i:i)=blank
              DIAbase(Nval)(i:i)=blank
            END DO
            DIAname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            DIAbase(Nval)=TRIM(ADJUSTL(DIAname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'STANAME') THEN
            DO i=1,LEN(STAname(Nval))
              STAname(Nval)(i:i)=blank
            END DO
            STAname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FLTNAME') THEN
            DO i=1,LEN(FLTname(Nval))
              FLTname(Nval)(i:i)=blank
            END DO
            FLTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'GRDNAME') THEN
            DO i=1,LEN(GRDname(Nval))
              GRDname(Nval)(i:i)=blank
            END DO
            GRDname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ININAME') THEN
            DO i=1,LEN(INIname(Nval))
              INIname(Nval)(i:i)=blank
            END DO
            INIname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'IRPNAME') THEN
            DO i=1,LEN(ITLname(Nval))
              IRPname(Nval)(i:i)=blank
            END DO
            IRPname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ITLNAME') THEN
            DO i=1,LEN(ITLname(Nval))
              ITLname(Nval)(i:i)=blank
            END DO
            ITLname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'IADNAME') THEN
            DO i=1,LEN(IADname(Nval))
              IADname(Nval)(i:i)=blank
            END DO
            IADname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NFFILES') THEN
            Npts=load_i(Nval, Rval, Ngrids, nFfiles)
            DO ng=1,Ngrids
              IF (nFfiles(ng).le.0) THEN
                IF (Master) WRITE (out,260) 'NFFILES',                  &
     &                            'Must be equal or greater than one.'
                exit_flag=4
                RETURN
              END IF
            END DO                
            Npts=MAXVAL(nFfiles)           
            allocate ( FRCids (Npts,Ngrids) )
            allocate ( FRCname(Npts,Ngrids) )
            FRCids(1:Npts,1:Ngrids)=-1
            DO ng=1,Ngrids
              DO k=1,Npts
                DO i=1,LEN(FRCname(k,ng))
                  FRCname(k,ng)(i:i)=blank
                END DO
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'FRCNAME') THEN
            DO ng=1,Ngrids
              DO i=1,nFfiles(ng)
                IF (Nval.eq.(i+(Ngrids-1)*nFfiles(ng))) THEN
                  FRCname(i,ng)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'CLMNAME') THEN
            DO i=1,LEN(CLMname(Nval))
              CLMname(Nval)(i:i)=blank
            END DO
            CLMname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'BRYNAME') THEN
            DO i=1,LEN(BRYname(Nval))
              BRYname(Nval)(i:i)=blank
            END DO
            BRYname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FWDNAME') THEN
            DO i=1,LEN(FWDname(Nval))
              FWDname(Nval)(i:i)=blank
              FWDbase(Nval)(i:i)=blank
            END DO
            FWDname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            FWDbase(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ADSNAME') THEN
            DO i=1,LEN(ADSname(Nval))
              ADSname(Nval)(i:i)=blank
            END DO
            ADSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'APARNAM') THEN
            DO i=1,LEN(aparnam)
              aparnam(i:i)=blank
            END DO
            aparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SPOSNAM') THEN
            DO i=1,LEN(sposnam)
              sposnam(i:i)=blank
            END DO
            sposnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FPOSNAM') THEN
            DO i=1,LEN(fposnam)
              fposnam(i:i)=blank
            END DO
            fposnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'BPARNAM') THEN
            DO i=1,LEN(bparnam)
              bparnam(i:i)=blank
            END DO
            bparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'IPARNAM') THEN
            DO i=1,LEN(iparnam)
              iparnam(i:i)=blank
            END DO
            iparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SPARNAM') THEN
            DO i=1,LEN(sparnam)
              sparnam(i:i)=blank
            END DO
            sparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'USRNAME') THEN
            DO i=1,LEN(USRname)
              USRname(i:i)=blank
            END DO
            USRname=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!
!  Set switch to create NetCDF file.
!
      DO ng=1,Ngrids
        DO i=1,NV
          IF (Hout(i,ng)) LdefHIS(ng)=.TRUE.
        END DO
        IF (((nrrec(ng).eq.0).and.(nAVG(ng).gt.ntimes(ng))).or.         &
     &      (nAVG(ng).eq.0)) THEN
          LdefAVG(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nAVG2(ng).gt.ntimes(ng))).or.        &
     &      (nAVG2(ng).eq.0)) THEN
          LdefAVG2(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nDIA(ng).gt.ntimes(ng))).or.         &
     &      (nDIA(ng).eq.0)) THEN
          LdefDIA(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nFLT(ng).gt.ntimes(ng))).or.         &
     &      (nFLT(ng).eq.0)) THEN
          LdefFLT(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nHIS(ng).gt.ntimes(ng))).or.         &
     &      (nHIS(ng).eq.0)) THEN
          LdefHIS(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nRST(ng).gt.ntimes(ng))).or.         &
     &      (nRST(ng).eq.0)) THEN
          LdefRST(ng)=.FALSE.
        END  IF
        IF (((nrrec(ng).eq.0).and.(nSTA(ng).gt.ntimes(ng))).or.         &
     &      (nSTA(ng).eq.0)) THEN
          LdefSTA(ng)=.FALSE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        WRITE (out,60) TRIM(title), TRIM(my_os), TRIM(my_cpu),          &
     &                 TRIM(my_fort), TRIM(my_fc), TRIM(my_fflags),     &
     &                 TRIM(svn_url), TRIM(svn_rev), TRIM(Rdir),        &
     &                 TRIM(Hdir), TRIM(Hfile), TRIM(Adir)
        DO ng=1,Ngrids
!
!  Report grid size and domain decomposition.  Check for correct tile
!  decomposition.
!
          WRITE (out,90) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
     &                   NtileI(ng), NtileJ(ng)
          IF (NtileI(ng)*NtileJ(ng).le.0) THEN
            WRITE (out,100) ng
            exit_flag=6
            RETURN
          END IF
          IF (MOD(NtileI(ng)*NtileJ(ng),numthreads).ne.0) THEN
            WRITE (out,100) ng
            exit_flag=6
            RETURN
          END IF
!
!  Report physical parameters.
!
          WRITE (out,110) ng
          WRITE (out,120) ntimes(ng), 'ntimes',                         &
     &          'Number of timesteps for 3-D equations.'
          WRITE (out,140) dt(ng), 'dt',                                 &
     &          'Timestep size (s) for 3-D equations.'
          WRITE (out,130) ndtfast(ng), 'ndtfast',                       &
     &          'Number of timesteps for 2-D equations between',        &
     &          'each 3D timestep.'
          WRITE (out,120) ERstr, 'ERstr',                               &
     &          'Starting ensemble/perturbation run number.'
          WRITE (out,120) ERend, 'ERend',                               &
     &          'Ending ensemble/perturbation run number.'
          WRITE (out,120) nrrec(ng), 'nrrec',                           &
     &          'Number of restart records to read from disk.'
          WRITE (out,170) LcycleRST(ng), 'LcycleRST',                   &
     &          'Switch to recycle time-records in restart file.'
          WRITE (out,130) nRST(ng), 'nRST',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into restart fields.'
          WRITE (out,130) ninfo(ng), 'ninfo',                           &
     &          'Number of timesteps between print of information',     &
     &          'to standard output.'
          WRITE (out,170) ldefout(ng), 'ldefout',                       &
     &          'Switch to create a new output NetCDF file(s).'
          WRITE (out,130) nHIS(ng), 'nHIS',                             &
     &          'Number of timesteps between the writing fields',       &
     &          'into history file.'
          IF (ndefHIS(ng).gt.0) THEN
            WRITE (out,130) ndefHIS(ng), 'ndefHIS',                     &
     &            'Number of timesteps between creation of new',        &
     &            'history files.'
          END IF
          WRITE (out,130) ntsAVG(ng), 'ntsAVG',                         &
     &          'Starting timestep for the accumulation of output',     &
     &          'time-averaged data.'
          WRITE (out,130) nAVG(ng), 'nAVG',                             &
     &          'Number of timesteps between the writing of',           &
     &          'time-averaged data into averages file.'
          IF (ndefAVG(ng).gt.0) THEN
            WRITE (out,130) ndefAVG(ng), 'ndefAVG',                     &
     &            'Number of timesteps between creation of new',        &
     &            'time-averaged file.'
          END IF
          WRITE (out,130) ntsAVG2(ng), 'ntsAVG2',                       &
     &          'Starting timestep for the accumulation of output',     &
     &          'secondary time-averaged data.'
          WRITE (out,130) nAVG2(ng), 'nAVG2',                           &
     &          'Number of timesteps between the writing of',           &
     &          'secondary time-averaged data into averages file.'
          IF (ndefAVG2(ng).gt.0) THEN
            WRITE (out,130) ndefAVG2(ng), 'ndefAVG2',                   &
     &            'Number of timesteps between creation of new',        &
     &            'time-averaged file.'
          END IF
          DO itrc=1,NAT+NPT
            WRITE (out,190) tnu2(itrc,ng), 'tnu2', itrc,                &
     &            'Horizontal, harmonic mixing coefficient (m2/s)',     &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) visc2(ng), 'visc2',                           &
     &          'Horizontal, harmonic mixing coefficient (m2/s)',       &
     &          'for momentum.'
          DO itrc=1,NAT+NPT
            WRITE (out,190) Akt_bak(itrc,ng), 'Akt_bak', itrc,          &
     &            'Background vertical mixing coefficient (m2/s)',      &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Akv_bak(ng), 'Akv_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for momentum.'
          WRITE (out,200) rdrg(ng), 'rdrg',                             &
     &          'Linear bottom drag coefficient (m/s).'
          WRITE (out,200) rdrg2(ng), 'rdrg2',                           &
     &          'Quadratic bottom drag coefficient.'
          WRITE (out,200) Zob(ng), 'Zob',                               &
     &          'Bottom roughness (m).'
          WRITE (out,200) blk_ZQ(ng), 'blk_ZQ',                         &
     &          'Height (m) of surface air humidity measurement.'
          IF (blk_ZQ(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZQ.',                                  &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZT(ng), 'blk_ZT',                         &
     &          'Height (m) of surface air temperature measurement.'
          IF (blk_ZT(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZT.',                                  &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZW(ng), 'blk_ZW',                         &
     &          'Height (m) of surface winds measurement.'
          IF (blk_ZW(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZW.',                                  &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,120) lmd_Jwt(ng), 'lmd_Jwt',                       &
     &          'Jerlov water type.'
          IF ((lmd_Jwt(ng).lt.1).or.(lmd_Jwt(ng).gt.5)) THEN
            WRITE (out,260) 'lmd_Jwt.',                                 &
     &            'It must between one and five.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,120) Vtransform(ng), 'Vtransform',                 &
     &          'S-coordinate transformation equation.'
          WRITE (out,120) Vstretching(ng), 'Vstretching',               &
     &          'S-coordinate stretching function.'
          WRITE (out,200) theta_s(ng), 'theta_s',                       &
     &          'S-coordinate surface control parameter.'
          WRITE (out,200) theta_b(ng), 'theta_b',                       &
     &          'S-coordinate bottom  control parameter.'
          WRITE (out,160) Tcline(ng), 'Tcline',                         &
     &          'S-coordinate surface/bottom layer width (m) used',     &
     &          'in vertical coordinate stretching.'
          WRITE (out,140) rho0, 'rho0',                                 &
     &          'Mean density (kg/m3) for Boussinesq approximation.'
          WRITE (out,140) dstart, 'dstart',                             &
     &          'Time-stamp assigned to model initialization (days).'
          WRITE (out,150) time_ref, 'time_ref',                         &
     &          'Reference time for units attribute (yyyymmdd.dd)'
          DO itrc=1,NAT+NPT
            WRITE (out,190) Tnudg(itrc,ng), 'Tnudg', itrc,              &
     &            'Nudging/relaxation time scale (days)',               &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Znudg(ng), 'Znudg',                           &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for free-surface.'
          WRITE (out,210) M2nudg(ng), 'M2nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 2D momentum.'
          WRITE (out,210) M3nudg(ng), 'M3nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 3D momentum.'
          WRITE (out,210) obcfac(ng), 'obcfac',                         &
     &          'Factor between passive and active',                    &
     &          'open boundary conditions.'
          WRITE (out,140) T0(ng), 'T0',                                 &
     &          'Background potential temperature (C) constant.'
          WRITE (out,140) S0(ng), 'S0',                                 &
     &          'Background salinity (PSU) constant.'
          WRITE (out,160) gamma2(ng), 'gamma2',                         &
     &          'Slipperiness variable: free-slip (1.0) or ',           &
     &          '                     no-slip (-1.0).'
          IF (Hout(idFsur,ng)) WRITE (out,170) Hout(idFsur,ng),         &
     &       'Hout(idFsur)',                                            &
     &       'Write out free-surface.'
          IF (Hout(idUbar,ng)) WRITE (out,170) Hout(idUbar,ng),         &
     &       'Hout(idUbar)',                                            &
     &       'Write out 2D U-momentum component.'
          IF (Hout(idVbar,ng)) WRITE (out,170) Hout(idVbar,ng),         &
     &       'Hout(idVbar)',                                            &
     &       'Write out 2D V-momentum component.'
          IF (Hout(idUvel,ng)) WRITE (out,170) Hout(idUvel,ng),         &
     &       'Hout(idUvel)',                                            &
     &       'Write out 3D U-momentum component.'
          IF (Hout(idVvel,ng)) WRITE (out,170) Hout(idVvel,ng),         &
     &       'Hout(idVvel)',                                            &
     &       'Write out 3D V-momentum component.'
          IF (Hout(idWvel,ng)) WRITE (out,170) Hout(idWvel,ng),         &
     &       'Hout(idWvel)',                                            &
     &       'Write out W-momentum component.'
          IF (Hout(idOvel,ng)) WRITE (out,170) Hout(idOvel,ng),         &
     &       'Hout(idOvel)',                                            &
     &       'Write out omega vertical velocity.'
          DO itrc=1,NAT
            IF (Hout(idTvar(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idTvar(itrc),ng), 'Hout(idTvar)',                  &
     &          'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          IF (Hout(idUair,ng)) WRITE (out,170) Hout(idUair,ng),         &
     &       'Hout(idUair)',                                            &
     &       'Write out surface U-winds.'
          IF (Hout(idVair,ng)) WRITE (out,170) Hout(idVair,ng),         &
     &       'Hout(idVair)',                                            &
     &       'Write out surface V-winds.'
          IF (Hout(idUsms,ng)) WRITE (out,170) Hout(idUsms,ng),         &
     &       'Hout(idUsms)',                                            &
     &       'Write out surface U-momentum stress.'
          IF (Hout(idVsms,ng)) WRITE (out,170) Hout(idVsms,ng),         &
     &       'Hout(idVsms)',                                            &
     &       'Write out surface V-momentum stress.'
          IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idUbms,ng),         &
     &       'Hout(idUbms)',                                            &
     &       'Write out bottom U-momentum stress.'
          IF (Hout(idVbms,ng)) WRITE (out,170) Hout(idVbms,ng),         &
     &       'Hout(idVbms)',                                            &
     &       'Write out bottom V-momentum stress.'
          IF (Hout(idWamp,ng)) WRITE (out,170) Hout(idWamp,ng),         &
     &       'Hout(idWamp)',                                            &
     &       'Write out wave height.'
          IF (Hout(idWlen,ng)) WRITE (out,170) Hout(idWlen,ng),         &
     &       'Hout(idWlen)',                                            &
     &       'Write out wave length.'
          IF (Hout(idWdir,ng)) WRITE (out,170) Hout(idWdir,ng),         &
     &       'Hout(idWdir)',                                            &
     &       'Write out wave direction.'
          IF (Hout(idTsur(itemp),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(itemp),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net heat flux.'
          IF (Hout(idTsur(isalt),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(isalt),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net salt flux.'
          IF (Hout(idSrad,ng)) WRITE (out,170) Hout(idSrad,ng),         &
     &       'Hout(idSrad)',                                            &
     &       'Write out shortwave radiation flux.'
          IF (Hout(idLrad,ng)) WRITE (out,170) Hout(idLrad,ng),         &
     &       'Hout(idLrad)',                                            &
     &       'Write out longwave radiation flux.'
          IF (Hout(idLhea,ng)) WRITE (out,170) Hout(idLhea,ng),         &
     &       'Hout(idLhea)',                                            &
     &       'Write out latent heat flux.'
          IF (Hout(idShea,ng)) WRITE (out,170) Hout(idShea,ng),         &
     &       'Hout(idShea)',                                            &
     &       'Write out sensible heat flux.'
          IF (Hout(idevap,ng)) WRITE (out,170) Hout(idevap,ng),         &
     &       'Hout(idevap)',                                            &
     &       'Write out evaporation rate.'
          IF (Hout(idrain,ng)) WRITE (out,170) Hout(idrain,ng),         &
     &       'Hout(idrain)',                                            &
     &       'Write out rain rate.'
          IF (Hout(idDano,ng)) WRITE (out,170) Hout(idDano,ng),         &
     &       'Hout(idDano)',                                            &
     &       'Write out density anomaly.'
          IF (Hout(idVvis,ng)) WRITE (out,170) Hout(idVvis,ng),         &
     &       'Hout(idVvis)',                                            &
     &       'Write out vertical viscosity coefficient.'
          IF (Hout(idTdif,ng)) WRITE (out,170) Hout(idTdif,ng),         &
     &       'Hout(idTdif)',                                            &
     &       'Write out vertical T-diffusion coefficient.'
          IF (Hout(idSdif,ng)) WRITE (out,170) Hout(idSdif,ng),         &
     &       'Hout(idSdif)',                                            &
     &       'Write out vertical S-diffusion coefficient.'
          IF (Hout(idHsbl,ng)) WRITE (out,170) Hout(idHsbl,ng),         &
     &       'Hout(idHsbl)',                                            &
     &       'Write out depth of surface boundary layer.'
          IF (Hout2(idFsur,ng)) WRITE (out,170) Hout2(idFsur,ng),       &
     &       'Hout2(idFsur)',                                           &
     &       'Write out free-surface.'
          IF (Hout2(idUbar,ng)) WRITE (out,170) Hout2(idUbar,ng),       &
     &       'Hout2(idUbar)',                                           &
     &       'Write out 2D U-momentum component.'
          IF (Hout2(idVbar,ng)) WRITE (out,170) Hout2(idVbar,ng),       &
     &       'Hout2(idVbar)',                                           &
     &       'Write out 2D V-momentum component.'
          IF (Hout2(idUvel,ng)) WRITE (out,170) Hout2(idUvel,ng),       &
     &       'Hout2(idUvel)',                                           &
     &       'Write out 3D U-momentum component.'
          IF (Hout2(idVvel,ng)) WRITE (out,170) Hout2(idVvel,ng),       &
     &       'Hout2(idVvel)',                                           &
     &       'Write out 3D V-momentum component.'
          IF (Hout2(idWvel,ng)) WRITE (out,170) Hout2(idWvel,ng),       &
     &       'Hout2(idWvel)',                                           &
     &       'Write out W-momentum component.'
          IF (Hout2(idOvel,ng)) WRITE (out,170) Hout2(idOvel,ng),       &
     &       'Hout2(idOvel)',                                           &
     &       'Write out omega vertical velocity.'
          DO itrc=1,NAT
            IF (Hout2(idTvar(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(idTvar(itrc),ng), 'Hout2(idTvar)',                &
     &          'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          IF (Hout2(idUair,ng)) WRITE (out,170) Hout2(idUair,ng),       &
     &       'Hout2(idUair)',                                           &
     &       'Write out surface U-winds.'
          IF (Hout2(idVair,ng)) WRITE (out,170) Hout2(idVair,ng),       &
     &       'Hout2(idVair)',                                           &
     &       'Write out surface V-winds.'
          IF (Hout2(idUsms,ng)) WRITE (out,170) Hout2(idUsms,ng),       &
     &       'Hout2(idUsms)',                                           &
     &       'Write out surface U-momentum stress.'
          IF (Hout2(idVsms,ng)) WRITE (out,170) Hout2(idVsms,ng),       &
     &       'Hout2(idVsms)',                                           &
     &       'Write out surface V-momentum stress.'
          IF (Hout2(idUbms,ng)) WRITE (out,170) Hout2(idUbms,ng),       &
     &       'Hout2(idUbms)',                                           &
     &       'Write out bottom U-momentum stress.'
          IF (Hout2(idVbms,ng)) WRITE (out,170) Hout2(idVbms,ng),       &
     &       'Hout2(idVbms)',                                           &
     &       'Write out bottom V-momentum stress.'
          IF (Hout2(idWamp,ng)) WRITE (out,170) Hout2(idWamp,ng),       &
     &       'Hout2(idWamp)',                                           &
     &       'Write out wave height.'
          IF (Hout2(idWlen,ng)) WRITE (out,170) Hout2(idWlen,ng),       &
     &       'Hout2(idWlen)',                                           &
     &       'Write out wave length.'
          IF (Hout2(idWdir,ng)) WRITE (out,170) Hout2(idWdir,ng),       &
     &       'Hout2(idWdir)',                                           &
     &       'Write out wave direction.'
          IF (Hout2(idTsur(itemp),ng)) WRITE (out,170)                  &
     &        Hout2(idTsur(itemp),ng), 'Hout2(idTsur)',                 &
     &       'Write out surface net heat flux.'
          IF (Hout2(idTsur(isalt),ng)) WRITE (out,170)                  &
     &        Hout2(idTsur(isalt),ng), 'Hout2(idTsur)',                 &
     &       'Write out surface net salt flux.'
          IF (Hout2(idSrad,ng)) WRITE (out,170) Hout2(idSrad,ng),       &
     &       'Hout2(idSrad)',                                           &
     &       'Write out shortwave radiation flux.'
          IF (Hout2(idLrad,ng)) WRITE (out,170) Hout2(idLrad,ng),       &
     &       'Hout2(idLrad)',                                           &
     &       'Write out longwave radiation flux.'
          IF (Hout2(idLhea,ng)) WRITE (out,170) Hout2(idLhea,ng),       &
     &       'Hout2(idLhea)',                                           &
     &       'Write out latent heat flux.'
          IF (Hout2(idShea,ng)) WRITE (out,170) Hout2(idShea,ng),       &
     &       'Hout2(idShea)',                                           &
     &       'Write out sensible heat flux.'
          IF (Hout2(idevap,ng)) WRITE (out,170) Hout2(idevap,ng),       &
     &       'Hout2(idevap)',                                           &
     &       'Write out evaporation rate.'
          IF (Hout2(idrain,ng)) WRITE (out,170) Hout2(idrain,ng),       &
     &       'Hout2(idrain)',                                           &
     &       'Write out rain rate.'
          IF (Hout2(idDano,ng)) WRITE (out,170) Hout2(idDano,ng),       &
     &       'Hout2(idDano)',                                           &
     &       'Write out density anomaly.'
          IF (Hout2(idVvis,ng)) WRITE (out,170) Hout2(idVvis,ng),       &
     &       'Hout2(idVvis)',                                           &
     &       'Write out vertical viscosity coefficient.'
          IF (Hout2(idTdif,ng)) WRITE (out,170) Hout2(idTdif,ng),       &
     &       'Hout2(idTdif)',                                           &
     &       'Write out vertical T-diffusion coefficient.'
          IF (Hout2(idSdif,ng)) WRITE (out,170) Hout2(idSdif,ng),       &
     &       'Hout2(idSdif)',                                           &
     &       'Write out vertical S-diffusion coefficient.'
          IF (Hout2(idHsbl,ng)) WRITE (out,170) Hout2(idHsbl,ng),       &
     &       'Hout2(idHsbl)',                                           &
     &       'Write out depth of surface boundary layer.'
!
!-----------------------------------------------------------------------
!  Report output/input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,220)
          WRITE (out,230) '           Output Restart File:  ',          &
     &                    TRIM(RSTname(ng))
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).eq.0) THEN
              WRITE (out,230) '           Output History File:  ',      &
     &                        TRIM(HISname(ng))
            ELSE
              Lstr=LEN_TRIM(HISname(ng))
              WRITE (out,230) '      Prefix for History Files:  ',      &
     &                        HISname(ng)(1:Lstr-3)
            END IF
          END IF
          IF (ndefAVG(ng).eq.0) THEN
            WRITE (out,230) '          Output Averages File:  ',        &
     &                      TRIM(AVGname(ng))
          ELSE
            Lstr=LEN_TRIM(AVGname(ng))
            WRITE (out,230) '     Prefix for Averages Files:  ',        &
     &                      AVGname(ng)(1:Lstr-3)
          END IF
          IF (ndefAVG2(ng).eq.0) THEN
            WRITE (out,230) '          Output Secondary File:  ',       &
     &                      TRIM(AVG2name(ng))
          ELSE
            Lstr=LEN_TRIM(AVG2name(ng))
            WRITE (out,230) '     Prefix for Secondary Files:  ',       &
     &                      AVG2name(ng)(1:Lstr-3)
          END IF
          fname=GRDname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '               Input Grid File:  ',          &
     &                    TRIM(fname)
          fname=INIname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '  Input Nonlinear Initial File:  ',          &
     &                    TRIM(fname)
          DO i=1,nFfiles(ng)
            fname=FRCname(i,ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,310) '         Input Forcing File ', i, ':  ',   &
     &                      TRIM(fname)
          END DO
          fname=CLMname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '        Input Climatology File:  ',          &
     &                    TRIM(fname)
          fname=bparnam
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '       Biology Parameters File:  ',          &
     &                    TRIM(fname)
          fname=varname
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '  IO Variable Information File:  ',          &
     &                    TRIM(fname)
          GO TO 40
  30      IF (Master) WRITE (out,270) TRIM(fname)
          exit_flag=4
          RETURN
  40      CONTINUE
        END DO
        IF (Nuser.gt.0) THEN
          WRITE (out,230) '        Input/Output USER File:  ',          &
     &                    TRIM(USRname)
        END IF
!
!-----------------------------------------------------------------------
!  Report generic USER parameters.
!-----------------------------------------------------------------------
!
        IF (Nuser.gt.0) THEN
          WRITE (out,240)
          DO i=1,Nuser
            WRITE (out,250) user(i), i, i
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Rescale active tracer parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO itrc=1,NAT+NPT
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
          tnu4(itrc,ng)=SQRT(ABS(tnu4(itrc,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(itrc,ng).gt.0.0_r8) THEN
            Tnudg(itrc,ng)=1.0_r8/(Tnudg(itrc,ng)*86400.0_r8)
          ELSE
            Tnudg(itrc,ng)=0.0_r8
          END IF
        END DO
      END DO
  50  FORMAT (/,' READ_PhyPar - Error while processing line: ',/,a)
  60  FORMAT (/,1x,a,/,                                                 &
     &        /,1x,'Operating system : ',a,                             &
     &        /,1x,'CPU/hardware     : ',a,                             &
     &        /,1x,'Compiler system  : ',a,                             &
     &        /,1x,'Compiler command : ',a,                             &
     &        /,1x,'Compiler flags   : ',a,/,                           &
     &        /,1x,'SVN Root URL  : ',a,                                &
     &        /,1x,'SVN Revision  : ',a,/,                              &
     &        /,1x,'Local Root    : ',a,                                &
     &        /,1x,'Header Dir    : ',a,                                &
     &        /,1x,'Header file   : ',a,                                &
     &        /,1x,'Analytical Dir: ',a)
  70  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Nodes: ',i3,',',2x,'Tiling: ',i3.3,      &
     &        'x',i3.3)
  80  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i2.2,1x,          &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI * NtileJ  must be equal to the number of ', &
     &        'parallel nodes.',                                        &
     &        /,12x,'Change -np value to mpirun or',                    &
     &        /,12x,'change domain partition in input script.')
  90  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Threads: ',i2,',',2x,'Tiling: ',i3.3,    &
     &        'x',i3.3)
 100  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i3.3,1x,          &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI*NtileJ must be a positive multiple of the', &
     &        ' number of threads.',                                    &
     &        /,12x,'Change number of threads (environment variable) ', &
     &        'or',/,12x,'change domain partition in input script.')
 110  FORMAT (/,/,' Physical Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)  
 120  FORMAT (1x,i10,2x,a,t30,a)
 130  FORMAT (1x,i10,2x,a,t30,a,/,t32,a)
 140  FORMAT (f11.3,2x,a,t30,a)
 150  FORMAT (f11.2,2x,a,t30,a)
 160  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
 170  FORMAT (10x,l1,2x,a,t30,a)
 180  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 195  FORMAT (1p,e11.4,2x,a,t30,a,i2.2,':',1x,a)
 200  FORMAT (1p,e11.4,2x,a,t30,a)
 210  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 220  FORMAT (/,' Output/Input Files:',/)
 230  FORMAT (2x,a,a)
 240  FORMAT (/,' Generic User Parameters:',/)
 250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t30,                         &
     &        'User parameter ',i2.2,'.')
 260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,/,15x,a)
 270  FORMAT (/,' READ_PHYPAR - could not find input file:  ',a)
 280  FORMAT (/,' READ_PHYPAR - variable info not yet loaded, ', a)
 290  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,i4,    &
     &        /,15x,a)
 300  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,'(',   &
     &        i2.2,')',/,15x,a)
 310  FORMAT (2x,a,i2.2,a,a)
      RETURN
      END SUBROUTINE read_PhyPar
      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, is, itrc, ng, status
      integer :: decode_line, load_i, load_l, load_r
      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc
      real(r8), dimension(NBT,Ngrids) :: Rbio
      real(r8), dimension(300) :: Rval
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(300) :: Cval
! ==================================================================== !
! READ Bering Sea  PARAMS
! ==================================================================== !
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, Parfrac)
!----------------------------------
!  Vertical mixing tuning parameter
!----------------------------------
          ELSE IF (TRIM(KeyWord).eq.'VertMixIncr') THEN 
            Npts=load_r(Nval, Rval, 1, VertMixIncr)         
!------------------
!  Bio- conversions
!------------------
          ELSE IF (TRIM(KeyWord).eq.'xi') THEN
            Npts=load_r(Nval, Rval, 1, xi)
          ELSE IF (TRIM(KeyWord).eq.'ccr') THEN
            Npts=load_r(Nval, Rval, 1, ccr)
          ELSE IF (TRIM(KeyWord).eq.'ccrPhL') THEN 
            Npts=load_r(Nval, Rval, 1, ccrPhL)
!-------------------------
!  extinction coefficients
!-------------------------
          ELSE IF (TRIM(KeyWord).eq.'k_ext') THEN
            Npts=load_r(Nval, Rval, 1, k_ext)
          ELSE IF (TRIM(KeyWord).eq.'k_chl') THEN
            Npts=load_r(Nval, Rval, 1, k_chl)
!-------------------
!  PhS growth params
!-------------------
          ELSE IF (TRIM(KeyWord).eq.'DiS') THEN
            Npts=load_r(Nval, Rval, 1, DiS)
          ELSE IF (TRIM(KeyWord).eq.'DpS') THEN
            Npts=load_r(Nval, Rval, 1, DpS)
          ELSE IF (TRIM(KeyWord).eq.'alphaPhS') THEN
            Npts=load_r(Nval, Rval, 1, alphaPhS)
          ELSE IF (TRIM(KeyWord).eq.'psiPhS') THEN
            Npts=load_r(Nval, Rval, 1, psiPhS)
          ELSE IF (TRIM(KeyWord).eq.'k1PhS') THEN
            Npts=load_r(Nval, Rval, 1, k1PhS)
          ELSE IF (TRIM(KeyWord).eq.'k2PhS') THEN
            Npts=load_r(Nval, Rval, 1, k2PhS)
          ELSE IF (TRIM(KeyWord).eq.'aPS') THEN
            Npts=load_r(Nval, Rval, 1, aPS)
!-------------------
!  PhL growth params
!-------------------
          ELSE IF (TRIM(KeyWord).eq.'DiL') THEN
            Npts=load_r(Nval, Rval, 1, DiL)
          ELSE IF (TRIM(KeyWord).eq.'DpL') THEN
            Npts=load_r(Nval, Rval, 1, DpL)
          ELSE IF (TRIM(KeyWord).eq.'alphaPhL') THEN
            Npts=load_r(Nval, Rval, 1, alphaPhL)
          ELSE IF (TRIM(KeyWord).eq.'psiPhL') THEN
            Npts=load_r(Nval, Rval, 1, psiPhL)
          ELSE IF (TRIM(KeyWord).eq.'k1PhL') THEN
            Npts=load_r(Nval, Rval, 1, k1PhL)
          ELSE IF (TRIM(KeyWord).eq.'k2PhL') THEN
            Npts=load_r(Nval, Rval, 1, k2PhL)
          ELSE IF (TRIM(KeyWord).eq.'aPL') THEN
            Npts=load_r(Nval, Rval, 1, aPL)
!-----------------------
!  MZS preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSMZS') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSMZS)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLMZS') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLMZS)
!-------------------------------
!  MZS growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eMZS') THEN
            Npts=load_r(Nval, Rval, 1, eMZS)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZS') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZS)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZST') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZST)
          ELSE IF (TRIM(KeyWord).eq.'fMZS') THEN
            Npts=load_r(Nval, Rval, 1, fMZS)
          ELSE IF (TRIM(KeyWord).eq.'kMZS') THEN
            Npts=load_r(Nval, Rval, 1, kMZS)
          ELSE IF (TRIM(KeyWord).eq.'gammaMZS') THEN
            Npts=load_r(Nval, Rval, 1, gammaMZS)
!-----------------------
!  MZL preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSMZL') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSMZL)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLMZL') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLMZL)
!-------------------------------
!  MZL growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eMZL') THEN
            Npts=load_r(Nval, Rval, 1, eMZL)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZL') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZL)
          ELSE IF (TRIM(KeyWord).eq.'Q10MZLT') THEN
            Npts=load_r(Nval, Rval, 1, Q10MZLT)
          ELSE IF (TRIM(KeyWord).eq.'fMZL') THEN
            Npts=load_r(Nval, Rval, 1, fMZL)
          ELSE IF (TRIM(KeyWord).eq.'kMZL') THEN
            Npts=load_r(Nval, Rval, 1, kMZL)
          ELSE IF (TRIM(KeyWord).eq.'gammaMZL') THEN
            Npts=load_r(Nval, Rval, 1, gammaMZL)
!-----------------------
!  Cop preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSCop') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSCop)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLCop') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLCop)
          ELSE IF (TRIM(KeyWord).eq.'fpMZSCop') THEN
            Npts=load_r(Nval, Rval, 1, fpMZSCop)
          ELSE IF (TRIM(KeyWord).eq.'fpMZLCop') THEN
            Npts=load_r(Nval, Rval, 1, fpMZLCop)
!-------------------------------
!  Cop growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eCop') THEN
            Npts=load_r(Nval, Rval, 1, eCop)
          ELSE IF (TRIM(KeyWord).eq.'Q10Cop') THEN
            Npts=load_r(Nval, Rval, 1, Q10Cop)
          ELSE IF (TRIM(KeyWord).eq.'Q10CopT') THEN
            Npts=load_r(Nval, Rval, 1, Q10CopT)
          ELSE IF (TRIM(KeyWord).eq.'fCop') THEN
            Npts=load_r(Nval, Rval, 1, fCop)
          ELSE IF (TRIM(KeyWord).eq.'gammaCop') THEN
            Npts=load_r(Nval, Rval, 1, gammaCop)
          ELSE IF (TRIM(KeyWord).eq.'kCop') THEN
            Npts=load_r(Nval, Rval, 1, kCop)
!-----------------------
!  NCa preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSNCa)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLNCa)
          ELSE IF (TRIM(KeyWord).eq.'fpMZSNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpMZSNCa)
          ELSE IF (TRIM(KeyWord).eq.'fpMZLNCa') THEN
            Npts=load_r(Nval, Rval, 1, fpMZLNCa)
!-------------------------------
!  NCa growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eNCa') THEN
            Npts=load_r(Nval, Rval, 1, eNCa)
          ELSE IF (TRIM(KeyWord).eq.'Q10NCa') THEN
            Npts=load_r(Nval, Rval, 1, Q10NCa)
          ELSE IF (TRIM(KeyWord).eq.'Q10NCaT') THEN
            Npts=load_r(Nval, Rval, 1, Q10NCaT)
          ELSE IF (TRIM(KeyWord).eq.'fNCa') THEN
            Npts=load_r(Nval, Rval, 1, fNCa)
          ELSE IF (TRIM(KeyWord).eq.'gammaNCa') THEN
            Npts=load_r(Nval, Rval, 1, gammaNCa)
          ELSE IF (TRIM(KeyWord).eq.'kNCa') THEN
            Npts=load_r(Nval, Rval, 1, kNCa)
!-----------------------
!  Eup preference params
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'fpPhSEup') THEN
            Npts=load_r(Nval, Rval, 1, fpPhSEup)
          ELSE IF (TRIM(KeyWord).eq.'fpPhLEup') THEN
            Npts=load_r(Nval, Rval, 1, fpPhLEup)
          ELSE IF (TRIM(KeyWord).eq.'fpMZSEup') THEN
            Npts=load_r(Nval, Rval, 1, fpMZSEup)
          ELSE IF (TRIM(KeyWord).eq.'fpMZLEup') THEN
            Npts=load_r(Nval, Rval, 1, fpMZLEup)
          ELSE IF (TRIM(KeyWord).eq.'fpCopEup') THEN
            Npts=load_r(Nval, Rval, 1, fpCopEup)
!-------------------------------
!  Eup growth and feeding params
!-------------------------------
          ELSE IF (TRIM(KeyWord).eq.'eEup') THEN
            Npts=load_r(Nval, Rval, 1, eEup)
          ELSE IF (TRIM(KeyWord).eq.'Q10Eup') THEN
            Npts=load_r(Nval, Rval, 1, Q10Eup)
          ELSE IF (TRIM(KeyWord).eq.'Q10EupT') THEN
            Npts=load_r(Nval, Rval, 1, Q10EupT)
          ELSE IF (TRIM(KeyWord).eq.'fEup') THEN
            Npts=load_r(Nval, Rval, 1, fEup)
          ELSE IF (TRIM(KeyWord).eq.'gammaEup') THEN
            Npts=load_r(Nval, Rval, 1, gammaEup)
          ELSE IF (TRIM(KeyWord).eq.'kEup') THEN
            Npts=load_r(Nval, Rval, 1, kEup)
!--------------------------
!  Jellyfish  param
!--------------------------
          ELSE IF (TRIM(KeyWord).eq.'fpCopJel') THEN
            Npts=load_r(Nval, Rval, 1, fpCopJel)
          ELSE IF (TRIM(KeyWord).eq.'fpEupJel') THEN
            Npts=load_r(Nval, Rval, 1, fpEupJel)
          ELSE IF (TRIM(KeyWord).eq.'fpNCaJel') THEN
            Npts=load_r(Nval, Rval, 1, fpNCaJel)
          ELSE IF (TRIM(KeyWord).eq.'eJel') THEN
            Npts=load_r(Nval, Rval, 1, eJel)
          ELSE IF (TRIM(KeyWord).eq.'Q10Jelr') THEN
            Npts=load_r(Nval, Rval, 1, Q10Jelr) 
          ELSE IF (TRIM(KeyWord).eq.'Q10JelTr') THEN
            Npts=load_r(Nval, Rval, 1, Q10JelTr)
          ELSE IF (TRIM(KeyWord).eq.'Q10Jele') THEN
            Npts=load_r(Nval, Rval, 1, Q10Jele) 
          ELSE IF (TRIM(KeyWord).eq.'Q10JelTe') THEN
            Npts=load_r(Nval, Rval, 1, Q10JelTe)
          ELSE IF (TRIM(KeyWord).eq.'gammaJel') THEN
            Npts=load_r(Nval, Rval, 1, gammaJel) 
          ELSE IF (TRIM(KeyWord).eq.'mpredJel') THEN
            Npts=load_r(Nval, Rval, 1, mpredJel) 
          ELSE IF (TRIM(KeyWord).eq.'respJel') THEN
            Npts=load_r(Nval, Rval, 1,respJel)
          ELSE IF (TRIM(KeyWord).eq.'bmJ') THEN
            Npts=load_r(Nval, Rval, 1,bmJ)
          ELSE IF (TRIM(KeyWord).eq.'ktbmJ') THEN
            Npts=load_r(Nval, Rval, 1,ktbmJ)
          ELSE IF (TRIM(KeyWord).eq.'TrefJ') THEN
            Npts=load_r(Nval, Rval, 1,TrefJ)
          ELSE IF (TRIM(KeyWord).eq.'fJel') THEN
            Npts=load_r(Nval, Rval, 1,fJel)
!--------------------------
!  Phytoplankton senescence
!--------------------------
          ELSE IF (TRIM(KeyWord).eq.'mPhS') THEN
            Npts=load_r(Nval, Rval, 1, mPhS)
          ELSE IF (TRIM(KeyWord).eq.'mPhL') THEN
            Npts=load_r(Nval, Rval, 1, mPhL)
          ELSE IF (TRIM(KeyWord).eq.'NcritPhS') THEN
            Npts=load_r(Nval, Rval, 1, NcritPhS)
          ELSE IF (TRIM(KeyWord).eq.'minmPhL') THEN
            Npts=load_r(Nval, Rval, 1, minmPhL)
          ELSE IF (TRIM(KeyWord).eq.'maxmPhL') THEN
            Npts=load_r(Nval, Rval, 1, maxmPhL)
          ELSE IF (TRIM(KeyWord).eq.'NcritPhL') THEN
            Npts=load_r(Nval, Rval, 1, NcritPhL)
!-----------------------
!  Zoopkankton mortality
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'mMZS') THEN
            Npts=load_r(Nval, Rval, 1, mMZS)
          ELSE IF (TRIM(KeyWord).eq.'mMZL') THEN
            Npts=load_r(Nval, Rval, 1, mMZL)
          ELSE IF (TRIM(KeyWord).eq.'mCop') THEN
            Npts=load_r(Nval, Rval, 1, mCop)
          ELSE IF (TRIM(KeyWord).eq.'mNCa') THEN
            Npts=load_r(Nval, Rval, 1, mNCa)
          ELSE IF (TRIM(KeyWord).eq.'mEup') THEN
            Npts=load_r(Nval, Rval, 1, mEup)
!-------------------
!  predation closure
!-------------------
         ELSE IF (TRIM(KeyWord).eq.'mpredMZL') THEN
            Npts=load_r(Nval, Rval, 1, mpredMZL)
          ELSE IF (TRIM(KeyWord).eq.'mpredCop') THEN
            Npts=load_r(Nval, Rval, 1, mpredCop)
          ELSE IF (TRIM(KeyWord).eq.'mpredNCa') THEN
            Npts=load_r(Nval, Rval, 1, mpredNCa)
          ELSE IF (TRIM(KeyWord).eq.'mpredEup') THEN
            Npts=load_r(Nval, Rval, 1, mpredEup)
!--------------------------------
!  sinking 
!--------------------------------
          ELSE IF (TRIM(KeyWord).eq.'wPhS') THEN
            Npts=load_r(Nval, Rval, 1, wPhS)
          ELSE IF (TRIM(KeyWord).eq.'wPhL') THEN
            Npts=load_r(Nval, Rval, 1, wPhL)
          ELSE IF (TRIM(KeyWord).eq.'wDet') THEN
            Npts=load_r(Nval, Rval, 1, wDet)
	  ELSE IF (TRIM(KeyWord).eq.'wDetF') THEN
            Npts=load_r(Nval, Rval, 1, wDetF)
!------------------------
!  Respiration parameters
!------------------------
          ELSE IF (TRIM(KeyWord).eq.'respMZS') THEN
            Npts=load_r(Nval, Rval, 1, respMZS)
          ELSE IF (TRIM(KeyWord).eq.'respMZL') THEN
            Npts=load_r(Nval, Rval, 1, respMZL)
          ELSE IF (TRIM(KeyWord).eq.'respPhS') THEN
            Npts=load_r(Nval, Rval, 1, respPhS)
          ELSE IF (TRIM(KeyWord).eq.'respPhL') THEN
            Npts=load_r(Nval, Rval, 1, respPhL)
          ELSE IF (TRIM(KeyWord).eq.'respCop') THEN
            Npts=load_r(Nval, Rval, 1, respCop)
          ELSE IF (TRIM(KeyWord).eq.'respNCa') THEN
            Npts=load_r(Nval, Rval, 1, respNCa)
          ELSE IF (TRIM(KeyWord).eq.'respEup') THEN
            Npts=load_r(Nval, Rval, 1, respEup)
          ELSE IF (TRIM(KeyWord).eq.'ktbmC') THEN
            Npts=load_r(Nval, Rval, 1,ktbmC)
	  ELSE IF (TRIM(KeyWord).eq.'TrefC') THEN
            Npts=load_r(Nval, Rval, 1,TrefC)
	  ELSE IF (TRIM(KeyWord).eq.'ktbmN') THEN
            Npts=load_r(Nval, Rval, 1,ktbmN)
	  ELSE IF (TRIM(KeyWord).eq.'TrefN') THEN
            Npts=load_r(Nval, Rval, 1,TrefN)
	  ELSE IF (TRIM(KeyWord).eq.'ktbmE') THEN
            Npts=load_r(Nval, Rval, 1,ktbmE)
	  ELSE IF (TRIM(KeyWord).eq.'TrefE') THEN
            Npts=load_r(Nval, Rval, 1,TrefE)
          ELSE IF (TRIM(KeyWord).eq.'TmaxPhS') THEN
            Npts=load_r(Nval, Rval, 1, TmaxPhS)
          ELSE IF (TRIM(KeyWord).eq.'TminPhS') THEN
            Npts=load_r(Nval, Rval, 1,TminPhS)
          ELSE IF (TRIM(KeyWord).eq.'Topt_PhS') THEN
            Npts=load_r(Nval, Rval, 1, Topt_PhS)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_PhS') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_PhS)
          ELSE IF (TRIM(KeyWord).eq.'TmaxPhL') THEN
            Npts=load_r(Nval, Rval, 1, TmaxPhL)
          ELSE IF (TRIM(KeyWord).eq.'TminPhL') THEN
            Npts=load_r(Nval, Rval, 1,TminPhL)
          ELSE IF (TRIM(KeyWord).eq.'Topt_PhL') THEN
            Npts=load_r(Nval, Rval, 1, Topt_PhL)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_PhL') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_PhL)
          ELSE IF (TRIM(KeyWord).eq.'TmaxMZS') THEN
            Npts=load_r(Nval, Rval, 1, TmaxMZS)
          ELSE IF (TRIM(KeyWord).eq.'KtBm_MZS') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_MZS)
	  ELSE IF (TRIM(KeyWord).eq.'TmaxMZL') THEN
            Npts=load_r(Nval, Rval, 1, TmaxMZL)
	  ELSE IF (TRIM(KeyWord).eq.'KtBm_MZL') THEN
            Npts=load_r(Nval, Rval, 1,KtBm_MZL)
!------------------------
!  Iron climatology terms
!------------------------
          ELSE IF (TRIM(KeyWord).eq.'Feinlo') THEN
            Npts=load_r(Nval, Rval, 1, Feinlo)
          ELSE IF (TRIM(KeyWord).eq.'Feinhi') THEN
            Npts=load_r(Nval, Rval, 1, Feinhi)
          ELSE IF (TRIM(KeyWord).eq.'Feinh') THEN
            Npts=load_r(Nval, Rval, 1, Feinh)
          ELSE IF (TRIM(KeyWord).eq.'Feofflo') THEN
            Npts=load_r(Nval, Rval, 1, Feofflo)
          ELSE IF (TRIM(KeyWord).eq.'Feoffhi') THEN
            Npts=load_r(Nval, Rval, 1, Feoffhi)
          ELSE IF (TRIM(KeyWord).eq.'Feoffh') THEN
            Npts=load_r(Nval, Rval, 1, Feoffh)
!-----------------------
!  Iron limitation terms
!-----------------------
          ELSE IF (TRIM(KeyWord).eq.'kfePhS') THEN
            Npts=load_r(Nval, Rval, 1, kfePhS)
          ELSE IF (TRIM(KeyWord).eq.'kfePhL') THEN
            Npts=load_r(Nval, Rval, 1, kfePhL)
          ELSE IF (TRIM(KeyWord).eq.'FeC') THEN
            Npts=load_r(Nval, Rval, 1, FeC)
!----------
!  Diapause
!----------
          ELSE IF (TRIM(KeyWord).eq.'NCmaxz') THEN
            Npts=load_r(Nval, Rval, 1, NCmaxz)
          ELSE IF (TRIM(KeyWord).eq.'wNCrise') THEN
            Npts=load_r(Nval, Rval, 1, wNCrise)
          ELSE IF (TRIM(KeyWord).eq.'wNCsink') THEN
            Npts=load_r(Nval, Rval, 1, wNCsink)
          ELSE IF (TRIM(KeyWord).eq.'RiseStart') THEN
            Npts=load_r(Nval, Rval, 1, RiseStart)
          ELSE IF (TRIM(KeyWord).eq.'RiseEnd') THEN
            Npts=load_r(Nval, Rval, 1, RiseEnd)
          ELSE IF (TRIM(KeyWord).eq.'SinkStart') THEN
            Npts=load_r(Nval, Rval, 1, SinkStart)
          ELSE IF (TRIM(KeyWord).eq.'SinkEnd') THEN
            Npts=load_r(Nval, Rval, 1, SinkEnd)
!-----------------------------------
!Remineralization and Nitrification
!----------------------------------
          ELSE IF (TRIM(KeyWord).eq.'regen') THEN
            Npts=load_r(Nval, Rval, 1, regen)
          ELSE IF (TRIM(KeyWord).eq.'dgrad') THEN
            Npts=load_r(Nval, Rval, 1, dgrad)
          ELSE IF (TRIM(KeyWord).eq.'Pv0') THEN
            Npts=load_r(Nval, Rval, 1, Pv0)
	        ELSE IF (TRIM(KeyWord).eq.'PvT') THEN
            Npts=load_r(Nval, Rval, 1, PvT)
	        ELSE IF (TRIM(KeyWord).eq.'Nitr0') THEN
            Npts=load_r(Nval, Rval, 1,Nitr0)
	        ELSE IF (TRIM(KeyWord).eq.'KnT') THEN
            Npts=load_r(Nval, Rval, 1,KnT)
	        ELSE IF (TRIM(KeyWord).eq.'ToptNtr') THEN
            Npts=load_r(Nval, Rval, 1,ToptNtr)
	        ELSE IF (TRIM(KeyWord).eq.'ktntr') THEN
            Npts=load_r(Nval, Rval, 1,ktntr)
	        ELSE IF (TRIM(KeyWord).eq.'KNH4Nit') THEN
            Npts=load_r(Nval, Rval, 1,KNH4Nit)
	        ELSE IF (TRIM(KeyWord).eq.'tI0') THEN
            Npts=load_r(Nval, Rval, 1,tI0)
	        ELSE IF (TRIM(KeyWord).eq.'KI') THEN
            Npts=load_r(Nval, Rval, 1,KI)  
!-------------------
!Benthic Parameters
!-------------------
         ELSE IF (TRIM(KeyWord).eq.'Hout(idBvar)') THEN
           Npts=load_l(Nval, Cval,NBEN*Ngrids, Ltrc)
           DO ng=1,Ngrids
             DO itrc=1,NBEN
               i=idBvar(idben(itrc))
               Hout(i,ng)=Ltrc(itrc,ng)
             END DO
           END DO
         ELSE IF (TRIM(KeyWord).eq.'iremin') THEN
           Npts=load_r(Nval, Rval, 1, iremin)
         ELSE IF (TRIM(KeyWord).eq.'q10') THEN
           Npts=load_r(Nval, Rval, 1, q10)
         ELSE IF (TRIM(KeyWord).eq.'q10r') THEN
           Npts=load_r(Nval, Rval, 1, q10r)
         ELSE IF (TRIM(KeyWord).eq.'Rup') THEN
           Npts=load_r(Nval, Rval, 1, Rup)
         ELSE IF (TRIM(KeyWord).eq.'KupD') THEN
           Npts=load_r(Nval, Rval, 1, KupD)
         ELSE IF (TRIM(KeyWord).eq.'KupP') THEN
           Npts=load_r(Nval, Rval, 1, KupP)
         ELSE IF (TRIM(KeyWord).eq.'LupD') THEN
           Npts=load_r(Nval, Rval, 1, LupD)
         ELSE IF (TRIM(KeyWord).eq.'LupP') THEN
           Npts=load_r(Nval, Rval, 1, LupP)
         ELSE IF (TRIM(KeyWord).eq.'Qres') THEN
           Npts=load_r(Nval, Rval, 1, Qres)
         ELSE IF (TRIM(KeyWord).eq.'Rres') THEN
           Npts=load_r(Nval, Rval, 1, Rres)
         ELSE IF (TRIM(KeyWord).eq.'rmort') THEN
           Npts=load_r(Nval, Rval, 1, rmort)
         ELSE IF (TRIM(KeyWord).eq.'eex') THEN
           Npts=load_r(Nval, Rval, 1, eex)
         ELSE IF (TRIM(KeyWord).eq.'eexD') THEN
           Npts=load_r(Nval, Rval, 1, eexD)
         ELSE IF (TRIM(KeyWord).eq.'prefD') THEN
           Npts=load_r(Nval, Rval, 1, prefD)
         ELSE IF (TRIM(KeyWord).eq.'prefPL') THEN
           Npts=load_r(Nval, Rval, 1, prefPL)
         ELSE IF (TRIM(KeyWord).eq.'prefPS') THEN
           Npts=load_r(Nval, Rval, 1, prefPS)
         ELSE IF (TRIM(KeyWord).eq.'T0ben') THEN
           Npts=load_r(Nval, Rval, 1, T0ben)
         ELSE IF (TRIM(KeyWord).eq.'T0benr') THEN
           Npts=load_r(Nval, Rval, 1, T0benr)
         ELSE IF (TRIM(KeyWord).eq.'BenPred') THEN
           Npts=load_r(Nval, Rval, 1, BenPred)
         ELSE IF (TRIM(KeyWord).eq.'bmB') THEN
           Npts=load_r(Nval, Rval, 1,bmB)
         ELSE IF (TRIM(KeyWord).eq.'ktbmB') THEN
           Npts=load_r(Nval, Rval, 1,ktbmB)
         ELSE IF (TRIM(KeyWord).eq.'TrefB') THEN
           Npts=load_r(Nval, Rval, 1,TrefB)
!--------------
!ice bio params
!--------------
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
!	    Rbio(15,1)=360_r8
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
	    END DO
	END IF 
	print*,'KEYWORD=',Keyword
          IF ( Lwrite .and.                                             &
     &         TRIM(KeyWord).ne.'TNU2' .and.                            &
     &         TRIM(KeyWord).ne.'TNU4' .and.                            &
     &         TRIM(KeyWord).ne.'AKT_BAK' .and.                         &
     &         TRIM(KeyWord).ne.'TNUDG' .and.                           &
     &         TRIM(KeyWord).ne.'Hout(idTvar)' ) THEN
            write(6,'(a15,i3,20e12.5)') TRIM(KeyWord),Nval,Rval(1:Nval)
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO
  30  FORMAT (/,' READ_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Biology Parameters, Grid: ',i2.2,                   &
     &        /,  ' ============================',/)
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
!  70  FORMAT (f11.3,2x,a,t28,a)
!  80  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
! 100  FORMAT (1p,e11.4,2x,a,t28,a)
! 110  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)
      RETURN
      END SUBROUTINE read_BioPar
      FUNCTION decode_line (line_text, KeyWord, Nval, Cval, Rval)
!
!=======================================================================
!                                                                      !
!  This function decodes lines of text from input script files.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
! Imported variable declarations.
!
      character (len=*), intent(in) :: line_text
      character (len=40), intent(inout) :: KeyWord
      integer, intent(inout) :: Nval
      character (len=160), dimension(300), intent(inout) :: Cval
      real(r8), dimension(300), intent(inout) :: Rval
!
! Local variable declarations
!
      logical :: IsString, Kextract, decode, nested
      integer :: Iblank, Icmm, Kstr, Kend, Linp
      integer :: Lend, LenS, Lstr, Lval, Nmul, Schar
      integer :: copies, i, ic, ie, is, j, status
      integer, dimension(20) :: Imul
      integer :: decode_line
      character (len=1 ), parameter :: blank = ' '
      character (len=160) :: Vstring, line, string
!
!------------------------------------------------------------------------
!  Decode input line.
!------------------------------------------------------------------------
!
!  Initialize.
!
      DO i=1,LEN(line)
        line(i:i)=blank
        Vstring(i:i)=blank
        string(i:i)=blank
      END DO
!
!  Get length of "line", remove leading and trailing blanks.
!
      line=TRIM(ADJUSTL(line_text))
      Linp=LEN_TRIM(line)
!
!  If not a blank or comment line [char(33)=!], decode and extract input
!  values.  Find equal sign [char(61)].
!
      status=-1
      nested=.FALSE.
      IF ((Linp.gt.0).and.(line(1:1).ne.CHAR(33))) THEN
        status=1
        Kstr=1
        Kend=INDEX(line,CHAR(61),BACK=.FALSE.)-1
        Lstr=INDEX(line,CHAR(61),BACK=.TRUE.)+1
!
! Determine if KEYWORD is followed by double equal sign (==) indicating
! nested parameter.
!
        IF ((Lstr-Kend).eq.3) nested=.TRUE.
!
! Extract KEYWORD, trim leading and trailing blanks.
!
        Kextract=.FALSE.
        IF (Kend.gt.0) THEN
          Lend=Linp
          KeyWord=line(Kstr:Kend)
          Nval=0
          Kextract=.TRUE.
        ELSE
          Lstr=1
          Lend=Linp
          Kextract=.TRUE.
        END IF
!
! Extract parameter values string.  Remove comments [char(33)=!] or
! continuation symbol [char(92)=\], if any.  Trim leading trailing
! blanks.
!
        IF (Kextract) THEN
          Icmm=INDEX(line,CHAR(33),BACK=.FALSE.)
          IF (Icmm.gt.0) Lend=Icmm-1
          Icmm=INDEX(line,CHAR(92),BACK=.FALSE.)
          IF (Icmm.gt.0) Lend=Icmm-1
          Vstring=ADJUSTL(line(Lstr:Lend))
          Lval=LEN_TRIM(Vstring)
!
! The TITLE KEYWORD is a special one since it can include strings,
! numbers, spaces, and continuation symbol.
!
          IsString=.FALSE.
          IF (TRIM(KeyWord).eq.'TITLE') THEN
            Nval=Nval+1
            Cval(Nval)=Vstring(1:Lval)
            IsString=.TRUE.
          ELSE
!
! Check if there is a multiplication symbol [char(42)=*] in the variable
! string indicating repetition of input values.
!
            Nmul=0
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(42)) THEN
                Nmul=Nmul+1           
                Imul(Nmul)=i
              END IF
            END DO
            ic=1
!
! Check for blank spaces [char(32)=' '] between entries and decode.
!
            is=1
            ie=Lval
            Iblank=0
            decode=.FALSE.
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(32)) THEN
                IF (Vstring(i+1:i+1).ne.CHAR(32)) decode=.TRUE.
                Iblank=i
              ELSE
                ie=i
              ENDIF                                    
              IF (decode.or.(i.eq.Lval)) THEN
                Nval=Nval+1
!
! Processing numeric values.  Check starting character to determine
! if numeric or character values. It is possible to have both when
! processing repetitions via the multiplication symbol.              
!
                Schar=ICHAR(Vstring(is:is))
                IF (((48.le.Schar).and.(Schar.le.57)).or.               &
     &              (Schar.eq.43).or.(Schar.eq.45)) THEN
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Schar=ICHAR(Vstring(Imul(ic)+1:Imul(ic)+1))
                    IF ((43.le.Schar).and.(Schar.le.57)) THEN
                      READ (Vstring(Imul(ic)+1:ie),*) Rval(Nval)
                      DO j=1,copies-1
                        Rval(Nval+j)=Rval(Nval)
                      END DO
                    ELSE
                      string=Vstring(Imul(ic)+1:ie)
                      LenS=LEN_TRIM(string)
                      Cval(Nval)=string(1:LenS)
                      DO j=1,copies-1
                        Cval(Nval+j)=Cval(Nval)
                      END DO
                    END IF
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    LenS=LEN_TRIM(string)
                    READ (string(1:LenS),*) Rval(Nval)
                  END IF
                ELSE
!
! Processing character values (logicals and strings).
!
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Cval(Nval)=Vstring(Imul(ic)+1:ie)
                    DO j=1,copies-1
                      Cval(Nval+j)=Cval(Nval)
                    END DO
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    Cval(Nval)=TRIM(ADJUSTL(string))
                  END IF
                  IsString=.TRUE.
                END IF
                is=Iblank+1
                ie=Lval
                decode=.FALSE.
              END IF
            END DO
          END IF
        END IF             
        status=Nval
      END IF
      decode_line=status       
      RETURN
      END FUNCTION decode_line
      FUNCTION load_i (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model integer     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     load_i     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(r8), intent(in) :: Vinp(Ninp)
      integer, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_i
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(Ninp))
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
      END IF
      load_i=ic
      RETURN
      END FUNCTION load_i
      FUNCTION load_l (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model logical     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     load_l     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      character (len=*), intent(in) :: Vinp(Ninp)
      logical, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_l
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.(Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          IF ((Vinp(Ninp)(1:1).eq.'T').or.(Vinp(Ninp)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.(Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
      END IF
      load_l=ic
      RETURN
      END FUNCTION load_l
      FUNCTION load_r (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model real        !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output real variable.                                 !
!     load_r     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(r8), intent(in) :: Vinp(Ninp)
      real(r8), intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_r
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=Vinp(Ninp)
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
      END IF
      load_r=ic
      RETURN
      END FUNCTION load_r
