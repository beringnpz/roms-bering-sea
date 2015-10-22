## ROMS: Bering Sea domain

This repository holds the source code for the Bering Sea ROMS domain, including biological code for the BESTNPZ and FEAST models.  The code originated in 2009 as a clone of Kate Hedstrom's branch.  This repository was initiated on 08/17/15 in an attempt to bring several slightly different versions of the code together into one cohesive set.

The ROMS folder holds the main source code.  The Apps/NEP folder holds most of the Bering 10k-related files.

### Compiling the NEP5 ROMS code

A script, `buildbering10k.sh`, is located in the top level folder of this repository, and can be used to compile code for all variants of the Bering 10k application, and is currently set up to compile 3 variants of the source code (physics-only, with BEST_NPZ, and with FEAST) on either beast or cluster1.  Additional compilation options (for example, adding FLOATS) can be modeled on this example.

With the exception of this build script, please try not to modify any of the files in this repository (including the makefile and nep5.h header file) unless the changes are intended to apply to all future use of the code.  Other changes, e.g. adding CPP flags for single set of experiments, should be done via the build script by adding to the `MY_CPP_FLAGS` exported variable. 


### Not in this repo

The compiled executable is not tracked by this repository.  Each user should compile the code themselves whenever updates to the source code are downloaded.

I have also not added any input parameter files (ocean.in, biology.in, varinfo.dat, etc.) files or any input forcing files to this repo.  To keep things as clean as possible, I believe any files that will be changed regularly for individual experiments should be kept outside the repo.  

### Tutorial

If you're new to git, read [the tutorial html file](https://cdn.rawgit.com/kakearney/roms-bering-sea/master/tutorial/roms_git.html) in the tutorial folder of this repository before you get started.

