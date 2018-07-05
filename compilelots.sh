# ./buildbering10k.sh -n -s PIvavg

./buildbering10k.sh -n -s PIvGPPavg
export MY_CPP_FLAGS="-DPI_CONSTANT"
./buildbering10k.sh -n -s PIcGPPavg
export MY_CPP_FLAGS="-DPI_CONSTANT -DGPPMID"
./buildbering10k.sh -n -s PIcGPPmid
export MY_CPP_FLAGS="-DGPPMID"
./buildbering10k.sh -n -s PIvGPPmid
