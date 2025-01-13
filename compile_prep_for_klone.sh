# Computer-specific path setup for compiling roms-bering-sea

# ifort compile and mpif90

module load intel/oneAPI

# netcdf-fortran

export   NETCDF_INCDIR=/gscratch/goaclim/GR042989_mCDR/Software/Programming/netcdf-fortran-4.5.3/intel-classic_2023.2.1/include
export   NETCDF_LIBDIR=/gscratch/goaclim/GR042989_mCDR/Software/Programming/netcdf-fortran-4.5.3/intel-classic_2023.2.1/lib

# build default K20 version

export CPP_FLAGS="-DPI_CONSTANT -DGPPMID"




