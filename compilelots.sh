# ./buildbering10k.sh -n -s PIvavg

export MY_CPP_FLAGS="-DGPPMID"
./buildbering10k.sh -n -s PIvmid

export MY_CPP_FLAGS="-DPI_CONSTANT"
./buildbering10k.sh -n -s PIcavg
