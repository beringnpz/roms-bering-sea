---
title: "The Bering10K dataset"
tags:
  - documentation
---

Numerous Bering 10K ROMS model simulations have been run to date, including hindcasts of the past few decades, long-term forecasts under CMIP5 and CMIP6 emissions scenarios, and seasonal retropective forecasts.  Data and metadata related to these simulations are held in a number of locations.  This page serves as a centralized hub for this data and metadata.

## The model

Model source code is available on GitHub: [beringnpz/roms-bering-sea](https://github.com/beringnpz/roms-bering-sea)

## The documentation

A few guides for working with the Bering10K output dataset:

- [The Bering10K Dataset documentation](https://drive.google.com/file/d/1GlITTIvbs2gRBMNIxdDI15cZU6mH4ckg/view?usp=sharing): A pdf describing the dataset, including:
    1. A description of the various simulations (base model versions, parent model forcing datasets, and biological modules) and the output naming scheme for each.
    2. A tutorial on working with the curvilinear, sigma-coordinate ROMS spatial grid, aimed at those who will be working with data on the native grid.
    3. An overview of the ACLIM index dataset; this is a set of time series derived from the Bering10K output and intended for Alaska Climate Integrated Modeling (ACLIM) collaborators.
- [Bering10K Simulaton Variables ](https://drive.google.com/file/d/1C1FCxRMBm0uBv2wEKwrGfHmLnjt_gFvG/view?usp=sharing): A spreadsheet listing all simulations and the archived output variables associated with each, updated periodically as new simulations are run or new variables are made available.  Note that this includes both data available on both public and private servers (see below).

Please also see the [Literature page](../literature) for a list scientific publications related to the model, including model description and skill assessment.

## The output data

The following simulations are available for public download on the PMEL server via either [THREDDS](https://data.pmel.noaa.gov/aclim/thredds/), [ERDDAP](https://data.pmel.noaa.gov/aclim/erddap/), or the [Live Action Server](https://data.pmel.noaa.gov/aclim/las/):

- Hindcast simulations: K20 operational variant and H16 variant from ACLIM Phase 1
- ACLIM Phase 1 long-term forecasts: CMIP5-forced simulations using the H16 model variant

These data portals are currently in a testing phase, and the data files may not match the above documentation.  Updated documentation encompassing these data portals will be made available soon.

The remaining output data is stored on the University of Washington hyak-mox computer cluster, and is accessible to collaborators only.  If interested in beginning a collaboration, please [contact](mailto:kelly.kearney@noaa.gov) the Bering10K team.

