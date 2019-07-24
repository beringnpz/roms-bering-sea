---
title: "Compiling the Bering10K code"
tags:
  - documentation
---


In the top level of the roms-bering-sea code, you'll find a script named `buildbering.sh`.  The goal of this particular build script is to allow users to compile the Bering10K model with a variety of different options without cluttering the code repository with temporary build files or continuously-changing header files.  It tells ROMS the location of the MPI and NetCDF libraries on a computer, and then sets a number of ROMS-specific environmental variables.

The MPI and NetCDF library part of the script is currently hard-coded for the specific computers on which this code has been compiled to date.  This includes the AFSC machines cluster1 and beast, as well as UW's hyak-mox cluster.  To compile on a new computer, you'll need to modify this portion of the build script appropriately.

The remainder of the script deals with setting ROMS environmental variables in such a way that one doesn't need to mess with any of the header files to compile.  The [ROMS wiki](https://www.myroms.org/wiki/build_Script) provides some detail on these various options.  The key things to know about this hard-coded parts of the build script:

 - $ROMS_APPLICATION = NEP5
 - $MY_ROOT_DIR = $MY_ROMS_SRC = top level of this repo
 - $MY_PROJECT_DIR = `./Apps/NEP`
 - $SCRATCH_DIR = `./Apps/NEP/Build_[phys/npz/feast/physdb/npzdb/feastdb]` (the part in brackets depends on the -npfNPF options)

Usage:
```
buildbering10k [-s <suffix>] [-e <epath>] [-pPnNfF] [-h]
```
Options

 - `-s suffix`:  add suffix string to the end of the ocean[M/G] 
              executables. Useful if compiling a version based on a 
              branch without wanting to overwrite the master-compiled 
              version.  If not included, the default names (oceanM_phys, 
              oceanM_npz, and oceanM_feast) will be used.
 - `-e epath`:   path to folder where compiled executables should be 
              placed.  Default is `../romsexecs/`
 - `-p`:         compile physics-only version
 - `-P`:         compile physics-only version in debug mode
 - `-n`:         compile BEST_NPZ version
 - `-N`:         compile BEST_NPZ version in debug mode
 - `-f`:         compile FEAST version
 - `-F`:         compile FEAST version in debug mode
 - `-h`:         show this help text

## Examples

A straightforward compilation of the physics-only code:

```shell
./buildbering.sh -p
```

This creates an executable, `../romsexecs/oceanM_phys`.  All the intermediate build files will be found in `./Apps/NEP/Build_phys/`, and you may need to examine them if the compilation fails.

The next example demonstrates how to add addition definitions without changing the primary header file.  The current "best" version of the BESTNPZ model (as of this post) requires a couple biology-specific compilation flags to be defined, in addition to those in the nep.h file.  It was compiled using:

```shell
export MY_CPP_FLAGS="-DGPPMID -DPI_CONSTANT"
./buildbering.sh -n -s 20190412 
```

The result is an executable, `../romsexecs/oceanM_npz_20190412`.  Note that I didn't need to define the BEST_NPZ option as a compilation flag; this is automatically added by the -n option.



