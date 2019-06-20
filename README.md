## ROMS: Bering Sea domain

This repository holds the source code used for the Bering Sea Regional Ocean Modeling System (ROMS) domain, including biological code for the BESTNPZ (lower trophic level) and FEAST (fisheries food web) models.  

This code is not necessarily domain-specific (input forcing files, including the Bering Sea grid file, are not included here).  However, the biological code is tailored to the Bering Sea ecosystem, and all application data (Apps/ folder) reflect the model's use with the Bering 10K grid.  Compilation may be possible for other ROMS grids, but the portability has not been tested.

This code originated in 2009 as a clone of Kate Hedstrom's branch of the ROMS code, a variant that includes a sea ice module. This git repository was initiated on 08/17/15 in an attempt to bring several slightly different versions of Bering Sea ROMS code together into one cohesive set.  It remains under active development, particularly with respect to the biological models.  It is not in sync with more recent versions of ROMS (myroms.org), though we have attempted to incorporate some key bug fixes from that version when applicable.

For a detailed description of the Bering Sea ROMS model, please see the following publications:

- Hermann AJ, Gibson GA, Bond NA, Curchitser EN, Hedstrom K, Cheng W, Wang M, Stabeno PJ, Eisner L, Cieciel KD (2013) A multivariate analysis of observed and modeled biophysical variability on the Bering Sea shelf: Multidecadal hindcasts (1970-2009) and forecasts (2010-2040). Deep Res Part II Top Stud Oceanogr 94:121–139
- Hermann AJ, Curchitser EN, Hedstrom K, Cheng W, Bond NA, Wang M, Aydin K, Stabeno PJ, Cokelet ED, Gibson GA (2016) Projected future biophysical states of the Bering Sea. Deep Sea Res Part II Top Stud Oceanogr 134:30–47

The Bering Sea ROMS model was derived from the larger Northeast Pacific (NEP5) ROMS model.  Much of the preliminary validation of that model, including sea ice and tidal dynamics, are equally applicable to the Bering Sea sub-domain, and are documented in:

- Danielson S, Curchitser E, Hedstrom K, Weingartner T, Stabeno P (2011) On ocean and sea ice modes of variability in the Bering Sea. J Geophys Res Ocean 116:1–24 

The lower trophic level BESTNPZ model is described here:

- Gibson GA, Spitz YH (2011) Impacts of biological parameterization, initial conditions, and environmental forcing on parameter sensitivity and uncertainty in a marine ecosystem model for the Bering Sea. J Mar Syst 88:214–231
- Kearney KA, Hermann A, Cheng W, Ortiz I, Aydin K (in prep) A coupled pelagic-benthic-sympagic biogeochemical model for the Bering Sea: documentation and validation of the BESTNPZ model within a high-resolution regional ocean model.


<!-- Please follow [this link](http://kakearney.github.io/roms-bering-sea/) for a description of the code, as well as a tutorial for new git users. -->