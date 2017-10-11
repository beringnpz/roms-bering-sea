# This script calls the ROMS makefile with the appropriate options for 
# the NEP5 application (Bering Sea 10km Grid), compiled for parallel runs 
# on either cluster1 or beast.  It compiles 3 separate variants: 
# physics-only, bestnpz, and feast.  
#
# See useage statement for syntax.


USEAGE="Usage: buildbering10k [-s <suffix>] [-e <epath>] [-pPnNfF] [-h]

where:

  -s suffix:  add suffix string to the end of the ocean[M/G] 
              executables. Useful if compiling a version based on a 
              branch without wanting to overwrite the master-compiled 
              version.  If not included, the default names (oceanM_phys, 
              oceanM_npz, and oceanM_feast) will be used.

  -e epath:   path to folder where compiled executables should be 
              placed.  Default is ../romsexecs/

  -p:         compile physics-only version

  -P:         compile physics-only version in debug mode

  -n:         compile BEST_NPZ version

  -N:         compile BEST_NPZ version in debug mode

  -f:         compile FEAST version

  -F:         compile FEAST version in debug mode
  
  -h:         show this help text
"

if [[ $# -eq 0 ]] ; then
    echo "This script now requires arguments:

$USEAGE"
    exit 0
fi

pfile="oceanM_phys"
nfile="oceanM_npz"
ffile="oceanM_feast"
Pfile="oceanG_phys"
Nfile="oceanG_npz"
Ffile="oceanG_feast"

epath="../romsexecs/"

pflag=false
nflag=false
fflag=false
Pflag=false
Nflag=false
Fflag=false

while getopts ":s:e:pPnNfFh" opt; do
  case $opt in
    s) pfile="oceanM_phys_${OPTARG}"
       nfile="oceanM_npz_${OPTARG}"
       ffile="oceanM_feast_${OPTARG}"
       Pfile="oceanG_phys_${OPTARG}"
       Nfile="oceanG_npz_${OPTARG}"
       Ffile="oceanG_feast_${OPTARG}"
      ;;
    e) epath=$OPTARG;;
    p) pflag=true;;
    P) Pflag=true;;
    n) nflag=true;;
    N) Nflag=true;;
    f) fflag=true;;
    F) Fflag=true;;
    h) echo "$USEAGE"
       exit
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$USEAGE" >&2
       exit 1
       ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "$USEAGE" >&2
      exit 1
      ;;
  esac
done

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
    
    source /opt/intel/Compiler/11.1/069/bin/ifortvars.sh intel64 # TODO: What does this do?
    
  ;;
  afsc-s29.afsc.noaa.gov) # cluster1
    export            PATH=/opt/intel/openmpi/163/bin:$PATH
    export   NETCDF_INCDIR=/home/aydink/include                    # netcdf include
    export   NETCDF_LIBDIR=/home/aydink/lib                        # netcdf lib
    
    source /opt/intel/Compiler/11.1/069/bin/ifortvars.sh intel64 # TODO: What does this do?
    
    ;;
  mox[12].hyak.local) # hyak-mox
    export   PATH=/gscratch/sw/intel-201703/compilers_and_libraries_2017.2.174/linux/mpi/intel64/bin/mpif90:$PATH
    export   NETCDF_INCDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/include # netcdf include
    export   NETCDF_LIBDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/lib     # netcdf lib
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

#--------------
# Make
#--------------

# Compile physics

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_phys
export      USE_DEBUG=

if [ "$pflag" = true ]; then
  make clean &>/dev/null
  echo "Compiling physics-only variant"
  make -j &> buildouterr.txt
  if [ $? -ne 0 ]; then
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
  else
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanM ${epath%%/}/${pfile}
    echo "  Success: $pfile created"
  fi
fi

# Compile physics (debugging)

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_physdb
export      USE_DEBUG=on

if [ "$Pflag" = true ]; then
  make clean &>/dev/null
  echo "Compiling physics-only (debugging) variant"
  make -j &> buildouterr.txt
  if [ $? -ne 0 ]; then
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
  else
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanG ${epath%%/}/${Pfile}
    echo "  Success: $Pfile created"
  fi
fi

# Compile bestnpz

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_npz
export      USE_DEBUG=
export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBEST_NPZ"

if [ "$nflag" = true ]; then
  make clean &>/dev/null
  echo "Compiling bestnpz variant"
  make -j &> buildouterr.txt
  if [ $? -ne 0 ]; then
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
      echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
  else
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
      mv oceanM ${epath%%/}/${nfile}
    echo "  Success: $nfile created"
  fi
fi

# Compile bestnpz (debugging)

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_npzdb
export      USE_DEBUG=on

if [ "$Nflag" = true ]; then
  make clean &>/dev/null
  echo "Compiling bestnpz (debugging) variant"
  make -j &> buildouterr.txt
  if [ $? -ne 0 ]; then
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
      echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
  else
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanG ${epath%%/}/${Nfile}
    echo "  Success: $Nfile created"
  fi
fi

# Compile Feast

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_feast
export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DFEAST"
export      USE_DEBUG=

if [ "$fflag" = true ]; then
  make clean &>/dev/null
  echo "Compiling feast variant"
  make -j &> buildouterr.txt
  if [ $? -ne 0 ]; then
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
      echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
  else
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
      mv oceanM ${epath%%/}/${ffile}
    echo "  Success: $ffile created"
  fi
fi

# Compile Feast (debugging)

export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_feastdb
export      USE_DEBUG=on

if [ "$Fflag" = true ]; then
  make clean &>/dev/null
  echo "Compiling feast (debugging) variant"
  make -j &> buildouterr.txt
  if [ $? -ne 0 ]; then
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
      echo "  Compilation failed: see ${SCRATCH_DIR}/buildouterr.txt for details"
  else
    mv buildouterr.txt ${SCRATCH_DIR}/buildouterr.txt
    mv oceanG ${epath%%/}/${Ffile}
    echo "  Success: $Ffile created"
  fi
fi




