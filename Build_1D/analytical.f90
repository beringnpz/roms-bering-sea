      MODULE analytical_mod
!
!svn $Id: analytical.F 983 2009-05-23 01:07:05Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   PACKAGE:                                                 !
!                                                                      !
!  This package is used to provide various analytical fields to the    !
!  model when appropriate.                                             !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
      SUBROUTINE ana_biology (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_ocean
      USE mod_grid
      USE mod_grid
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL ana_biology_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % h,                              &
     &                       OCEAN(ng) % bt,                            &
     &                       OCEAN(ng) % t)
!
! Set analytical header file name used.
!
      IF (Lanafile.and.(tile.eq.0)) THEN
        ANANAME( 1)="ROMS/Functionals/ana_biology.h"
      END IF
      RETURN
      END SUBROUTINE ana_biology
!
!***********************************************************************
      SUBROUTINE ana_biology_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             z_r, h,                              &
     &                             bt,                                  &
     &                             t)
!***********************************************************************
!
      USE mod_param
      USE mod_biology
      USE mod_scalars
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: bt(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: i, is, itrc, j, k
      real(r8) :: var1, var2, var3, var4, var5, var6, var7
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: biod
      real(r8), dimension(NT(ng)) :: deepval
      real(r8), dimension(NT(ng)) :: loval
      real(r8), parameter :: eps = 1.0E-20_r8
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
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions for the BEST NPZ model .                                       !
!                                                                      !
!=======================================================================
!
!
! Make the curves.  First come parameters that are constant through 
! entire water column.
!
      deepval(iNO3) = 30.0_r8
      deepval(iNH4) = eps
      deepval(iPhS) = eps
      deepval(iPhL) = eps
      deepval(iMZS) = 0.0_r8
      deepval(iMZL) = 0.0_r8
      deepval(iCop) = 0.0_r8
      deepval(iNCaS) = 0.0_r8
      deepval(iEupS) = 0.0_r8
      deepval(iNCaO) = 0.0_r8
      deepval(iEupO) = 0.0_r8
      deepval(iDet) = eps
      deepval(iDetF) = eps
      deepval(iJel) = eps
      do k=1,N(ng)
         do j=JstrR,JendR
            do i=IstrR,IendR
               t(i,j,k,1,iNO3) = 18.0_r8 !  18.0_r8
               t(i,j,k,1,iNH4) =  0.0_r8 !2.0_r8
               t(i,j,k,1,iDet) =  eps
	        t(i,j,k,1,iDetF) =  eps
               biod(i,j,k) = -1 * ( z_r(i,j,k) + 2.5_r8 )
	       t(i,j,k,1,iJel) =  0.1_r8
            enddo
         enddo
      enddo
!
! PS - a combination of 2 curves: 2nd order polynomial curve to create
! a subsurface maximum at ~20 m and an exponentially decreasing curve to
! assure values aren't negative at depth.
!
      var1 = 25.58_r8
      var2 = -0.250_r8 / (5.0_r8**2)
      var3 =  0.008_r8 / (5.0_r8**3)
      var4 = 3.82_r8 * 5.0_r8
      var5 = 75.0_r8
      var6 = var5 - var4 
      var7 = var1 + var2*(var6**2) + var3*(var6**3)
      do i=IstrR,IendR
        do j=JstrR,JendR
          do k=1,N(ng)
            if ( biod(i,j,k) .le. var5 ) then
              var6 = biod(i,j,k) - var4
              t(i,j,k,1,iPhS) =                                         &
     &          var1 + var2*(var6**2) + var3*(var6**3)
            else
              t(i,j,k,1,iPhS) =                                         &
     &          var7 * exp( ( -1.0_r8*biod(i,j,k) + var5 ) / 5.0_r8 )
            endif
          enddo
        enddo
      enddo
!
! PL - exponentially decreasing with depth.
!
      var1 = 8.25_r8
      var2 = 0.322_r8 / 5.0_r8
      do i=IstrR,IendR
        do j=JstrR,JendR
          do k=1,N(ng)
            t(i,j,k,1,iPhL) = var1 * exp( -1.0_r8 * var2 * biod(i,j,k) )
          enddo
        enddo
      enddo
!
! Microzooplankton, from Howell-Kubler, 1996
! approximated with a straight line and an exponentially decreasing
! curve to assure values aren't negative at depth.  Curves meet ~60m
!
      var1 = 3.1714_r8
      var2 = -0.1865_r8 / 5.0_r8
      var3 = 60.0_r8
      var4 = 0.5_r8
      var5 = var1 + var2 * var3
      do i=IstrR,IendR
        do j=JstrR,JendR
          do k=1,N(ng)
            if ( biod(i,j,k) .le. var3 ) then
              t(i,j,k,1,iMZL) = var1 + var2 * biod(i,j,k)
            else
              t(i,j,k,1,iMZL) = var4 +                                  &
     &          ( var5 - var4 ) * exp( ( var3 - biod(i,j,k) ) / 5.0_r8 )
            endif
!            t(i,j,k,1,iMZS) = t(i,j,k,1,iMZL)
	    t(i,j,k,1,iMZS) = 0.0_r8
          enddo
        enddo
      enddo
!
! Pseudocalanus Copepods - from Shelikof data via Shelikof NPZ
! Step at ~ 50 m
!
      var1 = pi / 2.0_r8
      var2 = 2.876_r8 / pi
      var3 = 5.0_r8 / 5.0_r8
      var4 = 9.151_r8 * 5.0_r8
      do i=IstrR,IendR
        do j=JstrR,JendR
          do k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4 )
            t(i,j,k,1,iCop) = var1 - var2 * atan( var5 )
          enddo
        enddo
      enddo
!
! Neocalanus, from Shelikof NPZ. Step at ~ 30m
!
      var1 = 4.0_r8
      var2 = 1.3_r8
      var3 = 5.0_r8 / 5.0_r8
      var4 = 5.2_r8 * 5.0_r8
!On and off shelf neocalanus - at depth for january
      do i=IstrR,IendR
        do j=JstrR,JendR
          do k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4 )
	     if (z_r(i,j,1) .le. -400.0_r8) then 
	      if(z_r(i,j,k).le.-400.and.z_r(i,j,k).ge.-800)THEN
               t(i,j,k,1,iNCaO) = 3.0_r8
	      endif
	     elseif (z_r(i,j,1) .ge. -200.0_r8) then 
  	      if(k.le.5)THEN
               t(i,j,k,1,iNCaS) =  8.0_r8
	      endif
	     endif
          enddo
        enddo
      enddo
!
!
! Euphausiids, Wild guesses. Step at ~ 30m
!
      do i=IstrR,IendR
        do j=JstrR,JendR
	  var1 = grid(ng) % h(i,j)
          do k=1,N(ng)
	  if(var1.gt.100.0_r8.and.var1.lt.2000.0_r8) THEN
	   if (-z_r(i,j,k) .lt. 40.0_r8) then
	    t(i,j,k,1,iEupS) =  0.0_r8
	    t(i,j,k,1,iEupO) =  20.0_r8
	   end if
           else if (var1.gt.50_r8.and.var1.lt.100.0_r8) THEN 
	   if (-z_r(i,j,k) .lt. 25.0_r8) then
            t(i,j,k,1,iEupS) =  20.0_r8
	    t(i,j,k,1,iEupO) =  0_r8
            end if
           else 
            t(i,j,k,1,iEupS) =  0.0_r8
	    t(i,j,k,1,iEupO) =  0_r8
	   end if
          enddo
        enddo
      enddo
! Nitrate - depth dependent 
      do i=IstrR,IendR
        do j=JstrR,JendR
          var1 = MAX(0._r8,MIN(1._r8,                                   &
     &             (GRID(ng)%h(i,j)-20_r8)/(100_r8-20_r8)))
          var2 = 2_r8 + var1*(18_r8-2_r8)
          var3 = 2_r8 + var1*(30_r8-2_r8)
          var1 = (var3-var2) / 300._r8
          do k=1,N(ng)
	    if(h(i,j).gt.100_r8.and.z_r(i,j,k).ge.-50_r8)THEN
             t(i,j,k,1,iNO3) = var2  
	    else 
	    t(i,j,k,1,iNO3) = MIN(var3, var2 - z_r(i,j,k)*var1)
	    endif
          enddo
        enddo
      enddo
!
! Concentrations of everything below 100m - i.e. below
! depths where calculations are performed.  Have linear slope
! between values above and below.
! Iron deep values have already been determined.
!
      do i=IstrR,IendR
        do j=JstrR,JendR
          do k=N(ng),1,-1
            if ( biod(i,j,k) .gt. 300.0_r8 ) then       !was 120 for GOA
	      t(i,j,k,1,iNO3) = deepval(iNO3)
              t(i,j,k,1,iNH4) = deepval(iNH4)
              t(i,j,k,1,iPhS) = deepval(iPhS)
              t(i,j,k,1,iPhL) = deepval(iPhL)
              t(i,j,k,1,iMZS) = deepval(iMZS)
              t(i,j,k,1,iMZL) = deepval(iMZL)
              t(i,j,k,1,iCop) = deepval(iCop)
              t(i,j,k,1,iNCaS) = deepval(iNCaS)
              t(i,j,k,1,iEupS) = deepval(iEupS)
	      t(i,j,k,1,iNCaO) = deepval(iNCaO)
              t(i,j,k,1,iEupO) = deepval(iEupO)
              t(i,j,k,1,iDet) = deepval(iDet)
	      t(i,j,k,1,iDetF) = deepval(iDetF)
	      t(i,j,k,1,iJel) = deepval(iJel)
            else if ( biod(i,j,k) .gt. 100.0_r8 .and.                   &
     &                biod(i,j,k) .le. 300.0_r8) then
              var1 = ( 100.0_r8 - biod(i,j,k) ) / ( 100.0_r8-300.0_r8 )
              t(i,j,k,1,iNO3) = loval(iNO3) +                           &
     &                          ( deepval(iNO3) - loval(iNO3) ) * var1
              t(i,j,k,1,iNH4) = loval(iNH4) +                           &
     &                          ( deepval(iNH4) - loval(iNH4) ) * var1
              t(i,j,k,1,iPhS) = loval(iPhS) +                           &
     &                          ( deepval(iPhS) - loval(iPhS) ) * var1
              t(i,j,k,1,iPhL) = loval(iPhL) +                           &
     &                          ( deepval(iPhL) - loval(iPhL) ) * var1
              t(i,j,k,1,iMZS) = loval(iMZS) +                           &
     &                          ( deepval(iMZS) - loval(iMZS) ) * var1
              t(i,j,k,1,iMZL) = loval(iMZL) +                           &
     &                          ( deepval(iMZL) - loval(iMZL) ) * var1
              t(i,j,k,1,iCop) = loval(iCop) +                           &
     &                          ( deepval(iCop) - loval(iCop) ) * var1
              t(i,j,k,1,iNCaS) = loval(iNCaS) +                         &
     &                          ( deepval(iNCaS) - loval(iNCaS) ) * var1
              t(i,j,k,1,iEupS) = loval(iEupS) +                         &
     &                          ( deepval(iEupS) - loval(iEupS) ) * var1
              t(i,j,k,1,iNCaO) = loval(iNCaO) +                         &
     &                          ( deepval(iNCaO) - loval(iNCaO) ) * var1
              t(i,j,k,1,iEupO) = loval(iEupO) +                         &
     &                          ( deepval(iEupO) - loval(iEupO) ) * var1
              t(i,j,k,1,iDet) = loval(iDet) +                           &
     &                          ( deepval(iDet) - loval(iDet) ) * var1
              t(i,j,k,1,iDetF) = loval(iDetF) +                         &
     &                          ( deepval(iDetF) - loval(iDetF) ) * var1
              t(i,j,k,1,iJel) = loval(iJel) +                           &
                                ( deepval(iJel) - loval(iJel) ) * var1
            else
              loval(iNO3) = t(i,j,k,1,iNO3)
              loval(iNH4) = t(i,j,k,1,iNH4)
              loval(iPhS) = t(i,j,k,1,iPhS)
              loval(iPhL) = t(i,j,k,1,iPhL)
              loval(iMZS) = t(i,j,k,1,iMZS)
              loval(iMZL) = t(i,j,k,1,iMZL)
              loval(iCop) = t(i,j,k,1,iCop)
!              loval(iNCaS) = t(i,j,k,1,iNCaS)
              loval(iEupS) = t(i,j,k,1,iEupS)
!	      loval(iNCaO) = t(i,j,k,1,iNCaO)
              loval(iEupO) = t(i,j,k,1,iEupO)
              loval(iDet) = t(i,j,k,1,iDet)
	      loval(iDetF) = t(i,j,k,1,iDetF)
	      loval(iJel) = t(i,j,k,1,iJel) 
            endif
          enddo
        enddo
      enddo
!reduce the initial concs by 2 orders of magnitude 
!test for Bering Sea - initialize in Jan 
           do i=IstrR,IendR
              do j=JstrR,JendR
                do k=1,N(ng)
              t(i,j,k,1,iPhS) = 10_r8 !t(i,j,k,1,iPhS)*.01_r8
              t(i,j,k,1,iPhL) = 10_r8 !t(i,j,k,1,iPhL)*.01_r8
              t(i,j,k,1,iMZS) = 0.0_r8
              t(i,j,k,1,iMZL) = 10_r8
              t(i,j,k,1,iCop) = 3.0_r8
!              t(i,j,k,1,iNCaS) = t(i,j,k,1,iNCaS)*.0001_r8
!              t(i,j,k,1,iEupS) = t(i,j,k,1,iEupS)*.0001_r8
!	      t(i,j,k,1,iNCaO) = t(i,j,k,1,iNCaO)*.0001_r8
!              t(i,j,k,1,iEupO) = t(i,j,k,1,iEupO)*.0001_r8
           enddo
         enddo
       enddo
! Check for size, set other time index, and periodic BC's
!
      do i=IstrR,IendR
         do j=JstrR,JendR
            do k=1,N(ng)
               DO is=1,NBT
                  itrc=idbio(is)
                  t(i,j,k,1,itrc) = MAX(t(i,j,k,1,itrc),eps)
                  t(i,j,k,2,itrc) = t(i,j,k,1,itrc)
               enddo
            enddo
         enddo
      enddo
         do i=IstrR,IendR
           do j=JstrR,JendR
             do k=1,NBL(ng)
                 bt(i,j,k,1,iBen) = 8000.0_r8    !was eps 
		 bt(i,j,k,1,iBenDet) = 500.0_r8! 00.0_r8    !was eps  500
		 bt(i,j,k,2,iBen) = bt(i,j,k,1,iBen)
		 bt(i,j,k,2,iBenDet) = bt(i,j,k,1,iBenDet)
		  enddo
		 enddo  
            enddo
!g#ifdef 
!g         do i=IstrR,IendR
!g           do j=JstrR,JendR
!g             do k=1,NBL(ng)
!g                 bt(i,j,k,1,iBen) = MAX(bt(i,j,k,1,iBen),eps)
!g		 bt(i,j,k,1,iBenDet) = MAX(bt(i,j,k,1,iBenDet),eps)
!g		 bt(i,j,k,2,iBen) = bt(i,j,k,1,iBen)
!g		 bt(i,j,k,2,iBenDet) = bt(i,j,k,1,iBenDet)
!g		  enddo
!g		 enddo  
!g            enddo
!g#endif
      RETURN
      END SUBROUTINE ana_biology_tile
      SUBROUTINE ana_btflux (ng, tile, model, itrc)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic bottom flux of tracer type variables    !
!  (tracer units m/s).                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
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
      CALL ana_btflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      FORCES(ng) % btflx)
!
! Set analytical header file name used.
!
      IF (Lanafile.and.(tile.eq.0)) THEN
        ANANAME( 3)="ROMS/Functionals/ana_btflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_btflux
!
!***********************************************************************
      SUBROUTINE ana_btflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j
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
!  Set kinematic bottom heat flux (degC m/s) at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s) at horizontal RHO-points,
!  scaling by bottom salinity is done elsewhere.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_btflux_tile
      SUBROUTINE ana_hmixcoef (ng, tile, model)
!
!                                                                      !
!  This routine rescales horizontal mixing coefficients according      !
!  to the grid size.  Also,  if applicable,  increases horizontal      !
!  in sponge areas.                                                    !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL ana_hmixcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        MIXING(ng) % diff2,                       &
     &                        MIXING(ng) % visc2_p,                     &
     &                        MIXING(ng) % visc2_r,                     &
     &                        GRID(ng) % grdscl,                        &
     &                        GRID(ng) % xr,                            &
     &                        GRID(ng) % yr)
!
! Set analytical header file name used.
!
      IF (Lanafile.and.(tile.eq.0)) THEN
        ANANAME( 8)="ROMS/Functionals/ana_hmixcoef.h"
      END IF
      RETURN
      END SUBROUTINE ana_hmixcoef
!
!***********************************************************************
      SUBROUTINE ana_hmixcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              diff2,                              &
     &                              visc2_p,                            &
     &                              visc2_r,                            &
     &                              grdscl, xr, yr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
      real(r8), intent(inout) :: diff2(LBi:,LBj:,:)
      real(r8), intent(inout) :: visc2_p(LBi:,LBj:)
      real(r8), intent(inout) :: visc2_r(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: Iwrk, i, j, itrc
      real(r8) :: cff, cff1, cff2, fac
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
!  Scale horizontal viscosity according to the grid size.
!-----------------------------------------------------------------------
!
      cff=visc2(ng)/grdmax(ng)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          visc2_r(i,j)=cff*grdscl(i,j)
        END DO
      END DO
      cff=0.25_r8*cff
      DO j=Jstr,JendR
        DO i=Istr,IendR
          visc2_p(i,j)=cff*(grdscl(i,j  )+grdscl(i-1,j  )+              &
     &                      grdscl(i,j-1)+grdscl(i-1,j-1))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Scale horizontal diffusion according to the grid size.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        cff=tnu2(itrc,ng)/grdmax(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            diff2(i,j,itrc)=cff*grdscl(i,j)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        visc2_r)
      CALL exchange_p2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        visc2_p)
      DO itrc=1,NT(ng)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          diff2(:,:,itrc))
      END DO
      RETURN
      END SUBROUTINE ana_hmixcoef_tile
      SUBROUTINE ana_m2clima (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine sets analytical 2D momentum climatology fields.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL ana_m2clima_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng) % ubarclm,                       &
     &                       CLIMA(ng) % vbarclm)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(11)="Apps/1DBio/ana_m2clima.h"
      END IF
      RETURN
      END SUBROUTINE ana_m2clima
!
!***********************************************************************
      SUBROUTINE ana_m2clima_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ubarclm, vbarclm)
!***********************************************************************
!
      USE mod_param
!
      USE exchange_2d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(out) :: ubarclm(LBi:,LBj:)
      real(r8), intent(out) :: vbarclm(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
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
!  Set 2D momentum climatology.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ubarclm(i,j)=0.0_r8
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vbarclm(i,j)=0.0_r8
        END DO
      END DO
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ubarclm)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vbarclm)
      RETURN
      END SUBROUTINE ana_m2clima_tile
      SUBROUTINE ana_nudgcoef (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine set nudging coefficients time-scales (1/s).            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL ana_nudgcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
      IF (Lanafile.and.(tile.eq.0)) THEN
        ANANAME(16)="ROMS/Functionals/ana_nudgcoef.h"
      END IF
      RETURN
      END SUBROUTINE ana_nudgcoef
!
!***********************************************************************
      SUBROUTINE ana_nudgcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_boundary
      USE mod_clima
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Iwrk, i, itrc, j
      real(r8) :: cff1, cff2, cff3
      real(r8), parameter :: IniVal = 0.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
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
!  Set up nudging towards data time-scale coefficients (1/s).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          wrk(i,j)=0.0_r8
        END DO
      END DO
!
!  Default nudging coefficients.  Set nudging coefficients uniformly to
!  the values specified in the standard input file.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          CLIMA(ng)%M2nudgcof(i,j)=M2nudg(ng)
        END DO
      END DO
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            CLIMA(ng)%Tnudgcof(i,j,itrc)=Tnudg(itrc,ng)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set nudging coefficients (1/s) for passive/active (outflow/inflow)
!  open boundary conditions.  Weak nudging is expected in passive
!  outflow conditions and strong nudging is expected in active inflow
!  conditions.  Notice that interior nudging coefficient defined
!  above are zero out when boundary condition nudging.  The USER needs
!  to adapt this to his/her application!
!-----------------------------------------------------------------------
!
!  Free-surface nudging coefficients.
!
!
!  2D momentum nudging coefficients.
!
!
!  Tracers nudging coefficients.
!
!
!  3D momentum nudging coefficients.
!
      RETURN
      END SUBROUTINE ana_nudgcoef_tile
      SUBROUTINE ana_srflux (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This subroutine sets kinematic surface solar shortwave radiation    !
!  flux "srflx" (degC m/s) using an analytical expression.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
      CALL ana_srflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
     &                      FORCES(ng) % srflx)
!
! Set analytical header file name used.
!
      IF (Lanafile.and.(tile.eq.0)) THEN
        ANANAME(27)="ROMS/Functionals/ana_srflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_srflux
!
!***********************************************************************
      SUBROUTINE ana_srflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            lonr, latr,                           &
     &                            srflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: srflx(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: iday, month, year
      real(r8) :: Dangle, Hangle, LatRad
      real(r8) :: cff1, cff2, hour, yday
      real(r8) :: cff
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
!  Compute shortwave radiation (degC m/s):
!
!  ALBEDO option: Compute shortwave radiation flux using the Laevastu
!                 cloud correction to the Zillman equation for cloudless
!  radiation (Parkinson and Washington 1979, JGR, 84, 311-337).  Notice
!  that flux is scaled from W/m2 to degC m/s by dividing by (rho0*Cp).
!
!   option: Modulate shortwave radiation SRFLX (which
!                         read and interpolated elsewhere) by the local
!  diurnal cycle (a function of longitude, latitude and day-of-year).
!  This option is provided for cases where SRFLX computed by SET_DATA is
!  an average over >= 24 hours. For "diurnal_srflux" to work ana_srflux
!  must be undefined. If you want a strictly analytical diurnal cycle
!  enter it explicitly at the end of this subroutine or use the "albedo"
!  option.
!
!  For a review of shortwave radiation formulations check:
!
!    Niemela, S., P. Raisanen, and H. Savijarvi, 2001: Comparison of
!      surface radiative flux parameterizations, Part II, Shortwave
!      radiation, Atmos. Res., 58, 141-154.
!
!-----------------------------------------------------------------------
!
!  Assume time is in modified Julian day.  Get hour and year day.
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44_r8*COS((172.0_r8-yday)*2.0_r8*pi/365.25_r8)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
      Hangle=(12.0_r8-hour)*pi/12.0_r8
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!
!  Local daylight is a function of the declination (Dangle) and hour 
!  angle adjusted for the local meridian (Hangle-lonr(i,j)/15.0). 
!  The 15.0 factor is because the sun moves 15 degrees every hour.
!
          LatRad=latr(i,j)*deg2rad
          cff1=SIN(LatRad)*SIN(Dangle)
          cff2=COS(LatRad)*COS(Dangle)
!
!  SRFLX is reset on each time step in subroutine SET_DATA which 
!  interpolates values in the forcing file to the current date.
!  This  option is provided so that SRFLX values
!  corresponding to a greater or equal daily average can be modulated
!  by the local length of day to produce a diurnal cycle with the 
!  same daily average as the original data.  This approach assumes 
!  the net effect of clouds is incorporated into the SRFLX data. 
!
!  Normalization = (1/2*pi)*INTEGRAL{ABS(a+b*COS(t)) dt}  from 0 to 2*pi
!                = (a*ARCCOS(-a/b)+SQRT(b**2-a**2))/pi    for |a| < |b|
!  
          srflx(i,j)=MAX(0.0_r8, srflx(i,j))
          IF (ABS(cff1) > ABS(cff2)) THEN
            IF (cff1*cff2.gt.0.0_r8) THEN
              cff=cff1                                 ! All day case
              srflx(i,j)=MAX(0.0_r8,                                    &
     &                       srflx(i,j)/cff*                            &
     &                       (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
            ELSE
              srflx(i,j)=0.0_r8                        ! All night case
            END IF
          ELSE
            cff=(cff1*ACOS(-cff1/cff2)+SQRT(cff2*cff2-cff1*cff1))/pi
            srflx(i,j)=MAX(0.0_r8,                                      &
     &                     srflx(i,j)/cff*                              &
     &                     (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
          END IF
        END DO
      END DO
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        srflx)
      RETURN
      END SUBROUTINE ana_srflux_tile
      SUBROUTINE ana_stflux (ng, tile, model, itrc)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface flux of tracer type variables   !
!  "stflx" (tracer units m/s) using analytical expressions.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
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
      CALL ana_stflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % stflx)
!
! Set analytical header file name used.
!
      IF (Lanafile.and.(tile.eq.0)) THEN
        ANANAME(31)="ROMS/Functionals/ana_stflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_stflux
!
!***********************************************************************
      SUBROUTINE ana_stflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            srflx,                                &
     &                            stflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j
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
!  Set kinematic surface heat flux (degC m/s) at horizontal
!  RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater flux (m/s) at horizontal
!  RHO-points, scaling by surface salinity is done in STEP3D.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        stflx(:,:,itrc))
      RETURN
      END SUBROUTINE ana_stflux_tile
      END MODULE analytical_mod
