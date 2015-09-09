# This script calls the ROMS makefile with the appropriate options for 
# the NEP5 application (Bering Sea 10km Grid), compiled for parallel runs 
# on either cluster1 or beast.  It compiles 3 separate variants: 
# physics-only, bestnpz, and feast.  
#


# Things that are the same on all computers

export   ROMS_APPLICATION=NEP5 
export        NestedGrids=1
export        MY_ROOT_DIR=${PWD}                  # Location of roms-bering-sea clone
export     MY_PROJECT_DIR=${MY_ROOT_DIR}/Apps/NEP # Location of analyticals, headers, etc
export        MY_ROMS_SRC=${MY_ROOT_DIR}          # Main ROMS source code

export           USE_MPI=on
export        USE_MPIF90=on
export              FORT=ifort
export         USE_LARGE=on


# Things that vary by computer

case $HOSTNAME in
  afsc-s45.afsc.noaa.gov) # beast
    export           PATH=/usr/mpi/intel/openmpi-1.4.1/bin:$PATH  # mpif90 location
    export  NETCDF_INCDIR=/home/aydink/include                    # netcdf include
    export  NETCDF_LIBDIR=/home/aydink/lib                        # netcdf lib
	export      NC_CONFIG=home/aydink/bin/nc-config
    export           PATH=/home/aydink/bin:$PATH
    ;; 
  afsc-s29.afsc.noaa.gov) # cluster1
    export           PATH=/opt/intel/openmpi/163/bin:$PATH
    ;;
  *)   
    echo "Not set up to compile on this computer"
    exit 1
    ;;
esac

# These variables are based on those set above

export     MY_HEADER_DIR=${MY_PROJECT_DIR}  # Where header files are
export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}  # Where analytical files are
export            BINDIR=${MY_ROOT_DIR}     # Where the compiled program goes

# Make

echo "***********************************"
echo "Compiling Bering 10K: physics-only"
echo "***********************************"
echo "  Location of ifort:"
which ifort
echo "  Location of mpif90"
which mpif90
echo "  Location of netcdf"
which nc-config
echo "  Environment variables"
printenv
echo "***********************************"


cp ${MY_PROJECT_DIR}/GK_nep5.h-phys-floats ${MY_PROJECT_DIR}/nep5.h
export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build1

make clean
make -j

mv oceanM oceanM_phys


