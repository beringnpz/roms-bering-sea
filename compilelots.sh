export MY_CPP_FLAGS="-DNEWSHADE"
./buildbering10k.sh -n -s PIvAtten1avg
export MY_CPP_FLAGS="-DNEWSHADESHALLOW"
./buildbering10k.sh -n -s PIvAtten2avg
export MY_CPP_FLAGS="-DCOKELET"
./buildbering10k.sh -n -s PIvAtten3avg

export MY_CPP_FLAGS="-DNEWSHADE -DGPPMID"
./buildbering10k.sh -n -s PIvAtten1mid
export MY_CPP_FLAGS="-DNEWSHADESHALLOW -DGPPMID"
./buildbering10k.sh -n -s PIvAtten2mid
export MY_CPP_FLAGS="-DCOKELET -DGPPMID"
./buildbering10k.sh -n -s PIvAtten3mid

export MY_CPP_FLAGS="-DNEWSHADE -DPI_CONSTANT"
./buildbering10k.sh -n -s PIcAtten1avg
export MY_CPP_FLAGS="-DNEWSHADESHALLOW -DPI_CONSTANT"
./buildbering10k.sh -n -s PIcAtten2avg
export MY_CPP_FLAGS="-DCOKELET -DPI_CONSTANT"
./buildbering10k.sh -n -s PIcAtten3avg

