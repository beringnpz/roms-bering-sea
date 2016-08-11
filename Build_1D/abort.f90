      SUBROUTINE abort (status)
!
!svn $Id: abort.F 895 2009-01-12 21:06:20Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This subroutine terminates execution after flushing all buffers and  !
! closing IO files.                                                    !
!                                                                      !
!=======================================================================
!
      USE ocean_control_mod, ONLY : ROMS_finalize
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: status
!
!-----------------------------------------------------------------------
!  Terminate execution due to fatal error.
!-----------------------------------------------------------------------
!
!  Finalize ROMS component.
!
      CALL ROMS_finalize
!
!  Stop execution.
!
      STOP
      END SUBROUTINE abort   
