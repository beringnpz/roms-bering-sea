# svn $Id: Linux-ifort.mk 1020 2009-07-10 23:10:30Z kate $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2009 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for Intel IFORT (version 10.1) compiler on Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fotran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := ifort
#            FFLAGS := -heap-arrays -fp-model precise
#            FFLAGS := -heap-arrays
#            FFLAGS := -fp-model precise
            FFLAGS := 
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
          LDFLAGS := -Vaxlib
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_NETCDF4
    NETCDF_INCDIR ?= /home/aydink/include
    NETCDF_LIBDIR ?= /home/aydink/lib
      HDF5_LIBDIR ?= /home/aydink/lib
else
    NETCDF_INCDIR ?= /home/aydink/include
    NETCDF_LIBDIR ?= /home/aydink/lib
endif
# KAK 2017/10/11 netcdf library may be split across two files 
# (libnetcdf.a and libnetcdff.a)... this part is very 
# computer-specific at the moment, designed so code will work on mox with 
# two files (HOSTNAME is eithermox1.hyak.local or mox2.hyak.local) or 
# beast or cluster1 with a single file.
#ifeq ($(shell echo $(HOSTNAME) | head -c 3), mox)
#             LIBS := -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf
#else
#             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf
#endif
LIBS := -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf
ifdef USE_NETCDF4
             LIBS += -L$(HDF5_LIBDIR) -lhdf5_hl -lhdf5 -lz
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/intelsoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/intelsoft/PARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -lfmpi-pgi -lmpi-pgi 
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -openmp
endif

ifdef USE_DEBUG
#          FFLAGS += -g -check bounds -traceback
#           FFLAGS += -g -check uninit -ftrapuv -traceback -CB
           FFLAGS += -g -zero -CB -debug all -debug inline-debug-info -traceback
else
           FFLAGS += -ip -O3
 ifeq ($(CPU),i686)
           FFLAGS += -pc80 -xW
 endif
 ifeq ($(CPU),x86_64)
           FFLAGS += -xW
 endif
endif

ifdef USE_MCT
       MCT_INCDIR ?= /opt/intelsoft/mct/include
       MCT_LIBDIR ?= /opt/intelsoft/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

ifdef USE_ESMF
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) -lesmf -lC
endif

       clean_list += ifc* work.pc*

#
# Use full path of compiler.
#
               FC := $(shell which ${FC})
               LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -free
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += -free
$(SCRATCH_DIR)/analytical.o: FFLAGS += -free
$(SCRATCH_DIR)/biology.o: FFLAGS += -free
ifdef USE_ADJOINT
$(SCRATCH_DIR)/ad_biology.o: FFLAGS += -free
endif
ifdef USE_REPRESENTER
$(SCRATCH_DIR)/rp_biology.o: FFLAGS += -free
endif
ifdef USE_TANGENT
$(SCRATCH_DIR)/tl_biology.o: FFLAGS += -free
endif

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -nofree
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -nofree
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swanser.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -nofree
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -nofree
$(SCRATCH_DIR)/m_constants.o: FFLAGS += -free
$(SCRATCH_DIR)/m_fileio.o: FFLAGS += -free
$(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += -free
$(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -free

endif
