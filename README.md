# Code and data for "The Berkeley High Resolution Tropospheric NO2 Product"

This repository contains the Matlab code used to carry out the analysis
in the listed paper as well as netCDF file versions of the average NO2 
VCDs and other quantities for each change increment described in the paper.

The top-level `behr3_dataset_plots.m` has the driver code used to generate
the figures and tables in the paper. The netCDF files of the increment averages
can be found in Dependencies/BEHR-v3-update-analysis/Workspaces/IncrementTests.

If you want this code to function, you will have to run `BEHR_initial_setup.m` 
in the Dependencies/BEHR-core-utils folder, in order to generate the `behr_paths`
class and compile the PSM/CVM gridding code. Some modfication would also be necessary
as it expects to load Matlab files for the averages rather than the netCDF files,
plus it will not have access to the BEHR .mat files or the WRF output files.
