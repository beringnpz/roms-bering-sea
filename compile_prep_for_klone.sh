# Computer-specific path setup for compiling roms-bering-sea
# Note: this sets enviromental variables so needs to be sourced in the existing bash environment, not just executed, e.g.:
# source compile_prep_for_klone.sh

# ifort compile and mpif90

module load intel/oneAPI/2023.2.1
module load goaclim/netcdf-c-4.7.4/intel-classic_2023.2.1
module load goaclim/netcdf-fortran-4.5.3/intel-classic_2023.2.1
module load goaclim/openmpi-4.1.6/intel-classic_2023.2.1

export I_MPI_F90=ifort

# netcdf-fortran

export NETCDF_INCDIR=/sw/contrib/goaclim-src/netcdf-fortran-4.5.3/intel-classic_2023.2.1/include
export NETCDF_LIBDIR=/sw/contrib/goaclim-src/netcdf-fortran-4.5.3/intel-classic_2023.2.1/lib

# build default K20 version of roms-bering-sea code

export MY_CPP_FLAGS="-DPI_CONSTANT -DGPPMID"




