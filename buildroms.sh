#!/bin/sh
#
# This script compiles ROMS using the NEP5 options.



# Check that we're on one of the appropriate computers, and set some 
# compiler stuff

case $HOSTNAME in
  afsc-s29.afsc.noaa.gov) # beast
  	echo "Compiling"
	export PATH=/usr/mpi/intel/openmpi-1.4.1/bin:$PATH
	;; 
  afsc-s45.afsc.noaa.gov) # cluster1
  	echo "Compiling"
	export PATH=/opt/intel/openmpi/163/bin:$PATH
	;;
  *)   
  	echo "Not set up to compile on this computer"
    exit 1
	;;
esac

source /opt/intel/Compiler/11.1/069/bin/ifortvars.sh intel64
export LD_LIBRARY_PATH=/opt/intel/Compiler/11.1/069/lib/intel64/:LD_LIBRARY_PATH

# Build

make clean
make -j
