      MODULE mod_parallel
!
!svn $Id: mod_parallel.F 911 2009-01-27 23:36:21Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
        USE mod_param
        USE mod_strings, ONLY: Nregion
!
        implicit none
!
!  Switch to identify master processor. In serial and shared-memory
!  applications it is always true.
!
        logical :: Master
!
!  Switch to identify which thread is processing input/output files.
!  In distributed-memory applications, this thread can be the master
!  thread or all threads in case of parallel output. In serial and
!  shared-memory applications it is always true.
!
        logical :: InpThread
        logical :: OutThread
!
!  Number of shared-memory parallel threads.  In distributed memory
!  configurations, its value must be equal to one.
!
        integer :: numthreads = 1
!
!  Number distributed memory nodes.
!
        integer :: numnodes = 0
!
!  Parallel nodes assined to the ocean model.
!
        integer :: peOCN_frst          ! first ocean parallel node
        integer :: peOCN_last          ! last  ocean parallel node
!
!  Parallel threads/nodes counters used in critical parallel regions.
!
        integer :: tile_count = 0
        integer :: block_count = 0
        integer :: thread_count = 0
!
!  Profiling variables as function of parallel thread:
!
!    proc          Parallel process ID.
!    Cstr          Starting time for program region.
!    Cend          Ending time for program region.
!    Csum          Accumulated time for progam region.
!
        integer  :: proc(0:1,4,Ngrids)
        real(r8) :: Cstr(0:Nregion,4,Ngrids)
        real(r8) :: Cend(0:Nregion,4,Ngrids)
        real(r8) :: Csum(0:Nregion,4,Ngrids)
!
!  Distributed-memory master process and rank of the local process.
!
        integer, parameter :: MyMaster = 0
        integer :: MyRank = 0
        CONTAINS
        SUBROUTINE initialize_parallel
!
!=======================================================================
!                                                                      !
!  This routine initializes and spawn distribute-memory nodes.         !
!                                                                      !
!=======================================================================
!
          USE mod_param
          USE mod_iounits
          USE mod_scalars
          USE mod_strings, ONLY: Nregion
!
!  Local variable declarations.
!
          integer :: i
          integer :: my_numthreads
!
!-----------------------------------------------------------------------
!  Initialize shared-memory (OpenMP) or serial configuration.
!-----------------------------------------------------------------------
!
!  Inquire number of threads in parallel region.
!
          numthreads=my_numthreads()
          Master=.TRUE.
          InpThread=.TRUE.
          OutThread=.TRUE.
!
!-----------------------------------------------------------------------
!  Initialize profiling variables.
!-----------------------------------------------------------------------
!
          proc(0:1,1:4,1:Ngrids)=0
          Cstr(0:Nregion,1:4,1:Ngrids)=0.0_r8
          Cend(0:Nregion,1:4,1:Ngrids)=0.0_r8
          Csum(0:Nregion,1:4,1:Ngrids)=0.0_r8
          RETURN
        END SUBROUTINE initialize_parallel
      END MODULE mod_parallel
