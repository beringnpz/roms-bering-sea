      MODULE mod_stepping
!
!svn $Id: mod_stepping.F 999 2009-06-09 23:48:31Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains time stepping indices.                         !
!                                                                      !
!  Lnew      New descent algorithm state solution index.               !
!  Lold      Previous descent algorithm state solution index.          !
!                                                                      !
!  knew      Barotropic (fast) time-step index corresponding to the    !
!              newest values for 2D primitive equation variables.      !
!  krhs      Barotropic (fast) time-step index used to compute the     !
!              right-hand-terms of 2D primitive equation variables.    !
!  kstp      Barotropic (fast) time-step index to which the current    !
!              changes are added to compute new 2D primitive equation  !
!              variables.                                              !
!                                                                      !
!  nfm3      Float index for time level "n-3".                         !
!  nfm2      Float index for time level "n-2".                         !
!  nfm1      Float index for time level "n-1".                         !
!  nf        Float index for time level "n".                           !
!  nfp1      Float index for time level "n+1".                         !
!                                                                      !
!  nnew      Baroclinic (slow) time-step index corresponding to the    !
!              newest values for 3D primitive equation variables.      !
!  nrhs      Baroclinic (slow) time-step index used to compute the     !
!              right-hand-terms of 3D primitive equation variables.    !
!  nstp      Baroclinic (slow) time-step index to which the current    !
!              changes are added to compute new 3D primitive equation  !
!              variables.                                              !
!                                                                      !
!=======================================================================
!
        USE mod_param
        implicit none
        integer, private :: ig
        integer, dimension(Ngrids) :: knew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: krhs = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: kstp = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nnew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nrhs = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nstp = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: Lnew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: Lold = (/ (1, ig=1,Ngrids) /)
      END MODULE mod_stepping
