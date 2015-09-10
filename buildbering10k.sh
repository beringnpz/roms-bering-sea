# This script calls the ROMS makefile with the appropriate options for 
# the NEP5 application (Bering Sea 10km Grid), compiled for parallel runs 
# on either cluster1 or beast.  It compiles 3 separate variants: 
# physics-only, bestnpz, and feast.  
#

#--------------
# Setup
#--------------

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
    export            PATH=/usr/mpi/intel/openmpi-1.4.1/bin:$PATH  # mpif90 location
    export   NETCDF_INCDIR=/home/aydink/includeM                   # netcdf include
    export   NETCDF_LIBDIR=/home/aydink/libM                       # netcdf lib
	;; 
  afsc-s29.afsc.noaa.gov) # cluster1
    export            PATH=/opt/intel/openmpi/163/bin:$PATH
    export   NETCDF_INCDIR=/home/aydink/include                    # netcdf include
    export   NETCDF_LIBDIR=/home/aydink/lib                        # netcdf lib
    ;;
  *)   
    echo "Not set up to compile on this computer"
    exit 1
    ;;
esac

source /opt/intel/Compiler/11.1/069/bin/ifortvars.sh intel64 # TODO: What does this do?

# These variables are based on those set above

export     MY_HEADER_DIR=${MY_PROJECT_DIR}  # Where header files are
export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}  # Where analytical files are
export            BINDIR=${MY_ROOT_DIR}     # Where the compiled program goes

#--------------
# Make
#--------------

# Compile physics

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build1

make clean &>/dev/null
echo "Compiling physics-only variant"
make -j &> buildouterr.txt
if [ $? -ne 0 ]; then
	mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
else
	mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanM oceanM_phys
	echo "  Success: oceanM_phys created"
fi

# Compile bestnpz

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build2
export      MY_CPP_FLAGS="-DBIOLOGY -DBEST_NPZ"

make clean &>/dev/null
echo "Compiling bestnpz variant"
make -j &> buildouterr.txt
if [ $? -ne 0 ]; then
	mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
else
	mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanM oceanM_npz
	echo "  Success: oceanM_npz created"
fi

# Compile Feast

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build3
export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DFEAST"

make clean &>/dev/null
echo "Compiling feast variant"
make -j &> buildouterr.txt
if [ $? -ne 0 ]; then
	mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
else
	mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanM oceanM_feast
	echo "  Success: oceanM_feast created"
fi





