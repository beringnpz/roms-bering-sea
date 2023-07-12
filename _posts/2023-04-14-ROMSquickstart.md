---
title: "Quick-start guide to our ROMS setup"
tags:
  - documentation
---

This guide will get new users up and running with our current setup of ROMS on the UW computer cluster, [mox-hyak](https://wiki.cac.washington.edu/display/hyakusers/WIKI+for+Hyak+users).  

For those familiar with ROMS, we will set up the usual pieces (ROMS source code, application data, and input datasets) and pair that with my [ROMS Communication Toolbox](https://github.com/beringnpz/romscom) python tools.  This organization system allows us to mix and match the pieces needed for the many, many variations of ROMS that we are currently supporting across different projects while minimizing duplication and one-off edits to files (because those one-off edits make it _so_ much harder to keep a good record of each experiment!)

## Hyak basics

This tutorial assumes you have received a UW Hyak user account and have been added to the `bumblereem` group.  For the basics of using Hyak, I'll refer you to the [wiki]((https://wiki.cac.washington.edu/display/hyakusers/WIKI+for+Hyak+user).

### Folder structure

Here's the basic folder structure on mox, including folders we will reference in this tutorial:

    /usr/lusers
      /user1                 <- user folder
      /user2
    /gscratch
        /jisao
        /bumblereem
            /ROMS_Datasets   <- Shared data folder
            /groupmember1    <- personal folder
                /roms        <- ROMS source code
                /beringApps  <- ROMS application data
                /Project1    <- ROMS project working folder
            /groupmember2
            /...
        /other_groups...


The majority of our work will be placed in the scratch folder: `/gscratch/bumblereem`.  This folder is shared by all members of the group.  File and folder permissions vary based on who created them, but in general, most files here are fully accessible to all group members.

Under /gscratch/bumblereem, each group member can create their own personal folder.  Note that this is separate from your user folder, i.e. the folder you land in when you first log in.  Your user folder is assigned when you get an account, and has a very low quota; you should try not to store much there.  The personal folder under `/gscratch/bumblereem` is where most of your stuff should go instead (although note that this drive is not backed up).

### Modules

Most software is available on HYAK through preinstalled modules that can be loaded as needed.  To run and compile ROMS, we will need MPI, NetCDF, and python3.  Load these with the following commands:

```
module load icc_17-impi_2017
module load netcdf_fortran+c_4.4.1.1-icc_17
module load anaconda3_4.3.1
```

(These versions are out of date... but we tend to follow an "if it ain't broke, don't fix it" policy regarding finicky ROMS dependencies.  This works, so we'll stick with it for now.)

I've added these lines to my `~/.bash_profile` file so they're loaded on login.

## ROMS code, data, and toolboxes

Our ROMS workflow relies on several parts:

 - The source code (2 variants)
 - Application data (for several domains and compilations variants)
 - romscom python package
 - Input data (many, many variants)
 - Project-specific files
 
The first 3 of these are files shared by the group, version controlled by git, and hosted on GitHub.  The input forcing files are also shared (but not version-controlled).

In the past, we tried sharing a single copy of the code repos, but that led to  a mess of mixed file and folder permission issues when trying to compile.  So our new guidelines recommend everyone keep their own local copies of the repos.  Here are the steps to do that.  All repos should be cloned into your personal folder.

### Clone the "new ROMS" source code

"New ROMS" is a fork of Kate Hedstrom's version of ROMS.  Until 2021, Kate regularly synced this copy with the main myroms.org trunk version of ROMS.  However, she is no longer maintaining it, and so I am no longer actively pulling in any updates from her branch.

From your personal folder in `/gscratch/bumblereem`, clone the code:

```
$ git clone git@github.com:beringnpz/roms.git
```

Because I’m balancing so many different versions of ROMS in my own collection, I’ve chosen to rename this one to a less generic `roms-kate-ice` (since Kate’s version included many updates to the ice code) whenever I clone it:

```
$ git clone git@github.com:beringnpz/roms.git roms-kate-ice
```

But that’s just a personal preference; rename or not as you wish.

**branches**

When you first clone this repository, you will find yourself on the default branch, `kate\_svn`.  Somewhat confusingly, this is _not_ the branch you want to use to compile and develop code; you want to work with the `main` branch.  (`kate_svn` instead mirrors the myroms.org code -- or at least it did until Kate stopped updating -- so it doesn't contain the most up to date ice and biology options, or any of our custom Alaksa-region stuff).  You should switch to the `main` branch to compile the primary stable version of our ROMS domains:

```
$ git checkout main
```

When developing, please use a feature branch workflow: create a new branch off of `main` to do your coding and testing, and only merge the changes back into `main` when you're confident it compiles cleanly, works as intended, and ideally doesn't introduce any conflicts with existing compiled versions.  If you're not familiar with git and git workflows, I refer you to the many online tutorials on the topic.


### Clone older ROMS source code

The roms-bering-sea code holds the “old” version of ROMS (i.e. the K20P19) and BESTNPZ.  This is the version used for the operational hindcast, ACLIM, and FEAST.  If you need this one, clone as above:

```
$ git clone git@github.com:beringnpz/roms-bering-sea.git
```

### Clone ROMS application data

A ROMS application refers to a specific implementation of ROMS: the grid for the domain, related input files, variable info, header files, etc.  In early versions of ROMS this was included within the code repository (under the Apps folder), but now the standard recommendation is to keep that stuff fully separate.  The beringApps repo holds our application data. 

```
$ git clone git@github.com:beringnpz/beringApps.git
```

The repo currently includes two apps:

 - B10KH16\_with\_FEAST: The example App I put together for FEAST runs.  These files should be paired with a ROMS executable compiled off the roms-bering-sea aclimfeast branch.
 - Bering\_BGC\_variants: The app to run new ROMS in with either the BIO\_BANAS, BEST\_NPZ, or BIO\_COBALT biology (or with no bgc, ocean and ice only) and for either the Bering10K, NEP, or EBS3K domains.  Still a bit of a work in progress.

The subApps folder includes all the pieces that go into these apps.  The details of the subapps are outside the scope of this tutorial; talk to me if interested!

In many ROMS workflows, the application data folder doubles as a working directory or project folder.  However, in our workflow, this application data is considered a default template only.  Please do not be make any changes to the files in this repo unless they are intended to be permanent updates for all users and to all model configurations going forward.

### Clone python toolboxes

There are two different versions of this toolbox.  Both are designed to programmatically manipulate ROMS input files.  The original (romsascii) module is the one I used up until last 2021ish, and it is still used by some of my older scripts (like the one that updates the operational hindcast).  In prep for the BCG project, I overhauled things to be lot more flexible, robust, and intuitive based on lessons learned over the past few years, and the result is the newer ROMS Communication Toolbox (romscom).  I now consider romsascii deprecated, and I’m not making any updates to it.  But you may need it if you need back-compatibility with old scripts.  You’ll also need romsascii if you are running FEAST (because I haven’t yet ported some of the more niche FEAST-specific functions).  You can clone these as follows:

```
$ git clone https://bitbucket.org/kakearney/romsascii.git
$ git clone git@github.com:beringnpz/romscom.git
```

## Set up python environments

We use the conda package manager to download a few packages beyond those in the default Anaconda3 package loaded from hyak's pre-installed options.  In an ideal world, we'd just install those two packages in our root environment and go on with our lives.  But it turns out there are several conflicts involved.  The netCDF library on mox doesn't interface properly with the default version of the netCDF4 python interface that comes with Anaconda.  We can get a version from the conda-forge channel (a channel is just a server that stores python modules and other code) that does interface properly, but it also overloads some details of the netCDF library, making it incompatible with the ROMS compiler.  Likewise, the NCO package brings along with it a version of mpi that ends up overshadowing mox's mpirun command, and that doesn't work to run our ROMS batch jobs.  Luckily, conda allows us to set up environments that support the necessary combinations.  We end up needing 3 isolated environments: 

 1. root: used to compile ROMS
 2. pythonnc: includes the conda-forge netCDF4.  Must be active to interact with netCDF files in python (also Matlab, ncdump, etc... possibly R?).  Can't be active when *compiling* ROMS, but doesn't affect *running* of ROMS jobs.
 3. mynco: includes the conda-forge nco and also the conda-forge netCDF4.  Must be active to use any NCO tools.  Cannot be active when compiling *or* running ROMS.

First, before setting up these environments, you need to tell conda where to install its packages.  By default it will put these in your home directory (/usr/lusers/[username]), but mox sets a very low quota on that folder, so you want to redirect it elsewhere.

So, open `~/.condarc` in your favorite text editor.  Then edit the file (which is probably empty if you haven’t used conda before) with this text:

    auto_update_conda: true
    envs_dirs:
      - /gscratch/bumblereem/[personalfolder]/anaconda/dirs
    pkgs_dirs:
      - /gscratch/bumblereem/[personalfolder]/anaconda/pkgs
    conda-build:
      root-dir: /gscratch/bumblereem/[personalfolder]/anaconda/conda-builds
    channels:
      - defaults
      - conda-forge

Save and close.  Also, go ahead and create the indicated folders under your personal folder if they don’t exist.

Now we need to create the new environments:

```
 $ conda create --prefix [envs_dirs path]/pythonnc
 $ conda create --prefix [envs_dirs path]/mynco
 $ source activate pythonnc
 $ conda install -c conda-forge netcdf4
 $ source activate mynco
 $ conda install -c conda-forge netcdf4
 $ conda install -c conda-forge nco
```

(You can choose whatever names you want in place of “pythonnc” and “mynco”)

### Install the python toolboxes

Because they're still in development, I haven't made romsascii or romscom available via common package managers (yet).  So they need to be installed from your cloned source code.  You'll need to do this for all conda environments where you might use them.  To install, first enter an interactive build node:

```
$ srun -p build --time=4:00:00 --mem=100G --pty /bin/bash
```

Then build in each environment:

```
$ cd /gscratch/bumblereem/[personalfolder]/romscom/
$ source deactivate
$ python setup.py develop
$ source activate pythonnc
$ python setup.py develop
$ source activate mynco
$ python setup.py develop
```

Repeat for romsascii if necessary.  You can then exit the build node.

```
$ exit
```

## Compile ROMS

For compiling, I follow the recommended method of using a build script.  Please do not edit the makefiles within the ROMS source code repository; that code is shared by several users and we do not want to mess with others' compilation setups or add unecessesary changes to the repositories.

### Compiling new ROMS

Compilation for ROMS should always be performed from a working folder.  You can place this folder anywhere you like, probably under your personal folder somewhere, but *not* in any of the version-controlled folders (i.e. not in the source code or application data folders).  We'll refer to this as the project folder.

Before compiling, you first need to copy a few files from the Application Data folder to your project folder:

```
$ cd /gscratch/bumblereem/[personalfolder]/[projectfolder]
$ cp /gscratch/bumblereem/beringApps/Apps/Bering_BGC_variants/build_roms.bash .
$ cp /gscratch/bumblereem/beringApps/Apps/Bering_BGC_variants/ana_psource.h .
$ cp /gscratch/bumblereem/beringApps/Apps/Bering_BGC_variants/bering_10k.h .
```

You can then modify your *local* copy of the build script (build_roms.bash) to point to your locations:

```
line 147 MY_PROJECT_DIR=/gscratch/bumblereem/[personalfolder]/[projectfolder]
line 159 MY_ROMS_SRC=/gscratch/bumblereem/[personalfolder]/
```

Alternatively, you can keep the build script as is (with my paths as the default) and override those defaults with exported shell variables:

```
$ export MY_PROJECT_DIR=/gscratch/bumblereem/[personalfolder]/[projectfolder]
$ export MY_ROMS_SRC=/gscratch/bumblereem/[personalfolder]/
```

(I generally prefer the latter method, since I switch regularly between projects and prefer minimizing any editing of files, but both methods are fine).

This version of the build script offers a few extra options beyond those usually found there in standard ROMS versions:

    #    -j [N]       Compile in parallel using N CPUs                      :::
    #                  omit argument for all available CPUs                 :::
    #                                                                       :::
    #    -p macro     Prints any Makefile macro value. For example,         :::
    #                                                                       :::
    #                  build.bash -p FFLAGS                                 :::
    #                                                                       :::
    #    -noclean     Do not clean already compiled objects                 :::
    #                                                                       :::
    #    -db          Compile in debug mode (USE_DEBUG)                     :::
    #                                                                       :::
    #    -variant [x] CPP-combo variant to compile.  Can be:                :::
    #                 phys:    no bio, default ice [default]                :::
    #                 cobalt:  BIO_COBALT bio, default ice                  :::
    #                 bestnpz: BEST_NPZ bio, default ice                    :::
    #                 banas:   BIO_BANAS bio, default ice                   :::
    #                                                                       :::
    #    -serial      Compile serially with only one node                   :::
    #                                                                       :::
    #    -datestr     Use the specified 12-digit YYYYMMDDhhmm string rather :::
    #                 than the current time to label the simulation.  Needs :::
    #                 to be passed if you want to do a -noclean recompile   :::


The extra options were added to make it easier to switch between the different compilation variants without always having to manually add to the `MY\_CPP\_FLAGS` shell variable or modify the script.

To compile the physics version, first you need to check out a node from mox's allocation:

```
$ srun -p build --time=4:00:00 --mem=100G --pty /bin/bash
```

Export any ROMS-related variables as needed, e.g.:

```
$ export MY_PROJECT_DIR=/gscratch/bumblereem/[personalfolder]/[projectfolder]
$ export MY_ROMS_SRC=/gscratch/bumblereem/[personalfolder]/
```
(If you've modified your copy of the build script, you don't need to do anything here.  This is also where you can specify dynamic CPP flags by setting the MY_CPP_FLAGS variable.)

Then build:

```
./build_roms.bash -variant phys
```

The script should create an executable in your project folder with the name romsM_phys_YYYYMMDDhhmm, where the 12-digit string reflects the time of compilation.  Feel free to rename to something more meaningful.

Hopefully, the code compiled happily at that point.  If it fails, it will tell you so and point you to the log file where you can find all the detailed compile-time standard output.

You can then exit the interactive mode:

```
exit
```

### Compiling older ROMS

The workflow for older ROMS is a little different, since it was developed before we decided to stop trying to share code folders and executables, and when the application data folder was still under the umbrella of the source code.

Rather than compiling from the project folder and pointing the script to your source code, the older build script uses the opposite convention: you compile from the source code folder and point to where you want to place the executable.

Move to the source code folder:

```
$ cd /gscratch/bumblereem/[personalfolder]/roms-bering-sea
```

Enter interactive mode:

```
$ srun -p build --time=4:00:00 --mem=100G --pty /bin/bash
```

If compiling NPZ, activate a few compilation flags (these started as optional but are part of the "operational" setup):

```
$ export MY_CPP_FLAGS=”-DPI_CONSTANT -DGPPMID”
```

And build:

```
$ ./buildbering.sh -n -s [mytestname] -e [execpath]
```

The default path to place the executable is `../romsexecs`; you'll want to change that to point to your project folder.  The [mytestname] suffix will be appended to the romsM name to identify this specific compilation; set it to whatever you want (or omit it to just stick with romsM).

The exit interactive mode:

```
exit
```

## Running ROMS

At this point, you can run ROMS the "standard" way, manually creating and modifying lots of .in files.  But most of the files in the Application Data repository is geared toward running things "my way", i.e. with the romscom toolbox.  

To get familiar with the toolbox, I refer you to the [romscom](https://github.com/beringnpz/romscom) documentation and examples.  One of the example notebooks is tailored to this particular ROMS setup, and should walk you through the basics of writing a run script.  I won't cover the specific details of these scripts here since they can be very experiment-dependent.  

Assuming you've written an appropriate run script, you can now submit the job.  See the HYAK wiki for full details of using the queues.  But in general, my submission scripts usually look like this:


    #!/bin/bash
    
    ## Job Name
    #SBATCH --job-name=myjob
    
    ## Allocation Definition
    
    ## On mox and ikt, the account and partition options should be the same.
    #SBATCH --account=bumblereem
    #SBATCH --partition=bumblereem
    
    ## Resources
    ## Nodes
    #SBATCH --nodes=5
    ## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
    #SBATCH --ntasks-per-node=28
    
    ## Walltime 
    #SBATCH --time=10:00:00
    
    ## Memory per node
    #SBATCH --mem=10G
    
    module load icc_17-impi_2017
    module load netcdf_fortran+c_4.4.1.1-icc_17
    source activate pythonnc
    
    python myrunscript.py

Save this to a file (here, we’ll call it roms.sub).  The submit the job:

```
$ sbatch roms.sub
```

And your ROMS simulation should be off and running!




















