---
title: "The Bering10K dataset"
tags:
  - documentation
---

Numerous Bering 10K ROMS model simulations have been run to date, including hindcasts of the past few decades, long-term forecasts under CMIP5 and CMIP6 emissions scenarios, and seasonal retropective forecasts.  Data and metadata related to these simulations are held in a number of locations.  This page serves as a centralized hub for this data and metadata.

## The model

Model source code is available on GitHub: [beringnpz/roms-bering-sea](https://github.com/beringnpz/roms-bering-sea)

## The documentation

A few guides for working with the Bering10K output dataset can be found

- [The Bering10K Dataset documentation](https://zenodo.org/record/4586950/files/Bering10K_dataset_documentation.pdf): A pdf describing the dataset, including:
    1. A description of the various simulations (base model versions, parent model forcing datasets, and biological modules) and the output naming scheme for each.
    2. A tutorial on working with the curvilinear, sigma-coordinate ROMS spatial grid, aimed at those who will be working with data on the native grid.
    3. An overview of the ACLIM index dataset; this is a set of time series derived from the Bering10K output and intended for Alaska Climate Integrated Modeling (ACLIM) collaborators.
- [Bering10K Simulaton Variables ](https://zenodo.org/record/4586950/files/Bering10K_simulation_variables.xlsx?download=1): A spreadsheet listing all simulations and the archived output variables associated with each, updated periodically as new simulations are run or new variables are made available.  Note that this includes both data available on both public and private servers (see below).

Please also see the [Literature page](../literature) for a list scientific publications related to the model, including model description and skill assessment.

## The output data

The Bering10K simulation output data is available in a few different locations, depending on the status of the project they are associated with.

**The PMEL server**

This server hosts all data from completed projects, as well as for ongoing projects where the underlying model setup has been used in published papers.  As of this time, this includes:

- Hindcast simulations: K20 operational variant and H16 variant from ACLIM Phase 1
- ACLIM Phase 1 long-term forecasts: CMIP5-forced simulations using the H16 model variant

The data can be accessed through three different front-ends:

- [THREDDS](https://data.pmel.noaa.gov/aclim/thredds/): Catalog listing where data and metadata can be accessed and/or downloaded
- [Live Action Server](https://data.pmel.noaa.gov/aclim/las/): An interactive web interface with plotting and mapping capabilities, primarily for data exploration (though some limited download can be achieved from here)
- [ERDDAP](https://data.pmel.noaa.gov/aclim/erddap/): Web interface to access and download tabular data.  Note that only a small subset of the model output (primarily Level 3 indices) is able to be formatted for access through this interface.

Note that these data portals are currently in a testing phase, and documentation related to data access is being actively added to the [the Bering10K Dataset documentation](https://zenodo.org/record/4586950/files/Bering10K_dataset_documentation.pdf) on an ongoing basis.

**Google Drive (via NOAA's G Suite)**

Some data associated with the ongoing ACLIM 2.0 project has been uploaded to the ACLIM project folder on Google Drive.  This is accessible to participants in the ACLIM study, who have/will receive access to the relevant folder (`00_ACLIM_shared`) directly.  Data is in the following subfolder:

`00_ACLIM_shared/02_Data/Newest/ACLIM_ROMS_indices/roms_for_aclim/`

Due to size restrictions, only the Level 3 regional and survey-replicated indices data (Level 3) are included in this location.  For Level 1-2, see the hyak-mox option. 

**UW hyak mox computer**

All simulation output data, including the data mirrored on the PMEL server and Google Drive, is stored on the University of Washington hyak-mox computer cluster, and is accessible to collaborators only. If you are a current collaborator who needs access to this data, or are interested in beginning a collaboration, please [contact](mailto:kelly.kearney@noaa.gov) the Bering10K team.


## The tools

As part of the [ACLIM project](https://www.fisheries.noaa.gov/alaska/ecosystems/alaska-climate-integrated-modeling-project), a suite of R tools have been developed to assist with accessing, subsetting, and manipulating the Bering10K model data.  You can find these [on GitHub in the kholsman/ACLIM2 repository](https://github.com/kholsman/ACLIM2).


