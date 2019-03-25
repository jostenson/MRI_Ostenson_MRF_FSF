# SETUP #
MRI, Ostenson, Damon, Welch, MR Fingerprinting with Simultaneous T1, T2, and Fat Signal Fraction Estimation with Integrated B0 Correction Reduces Bias in Water T1 and T2 Estimates
=============================

Requirements:
-------------
* MATLAB  - for reconstruction, processing, and figure generation
* C compiler - for generating MEX files from contributing code
* Berkeley Advanced Reconstruction Toolbox (BART) - for coil combination and image reconstruction (https://mrirecon.github.io/bart/)
* Michigan Image Reconstruction Toolbox - for image simulations (https://web.eecs.umich.edu/~fessler/code/)


Tested Configuration:
---------------------
* Windows 10 Enterprise (v1803) Subsystem for Linux (Ubuntu 16.04)
* MATLAB R2017b (v9.3.0)
* (MATLAB R2018b in Windows 10 used for figure generation script to improve graphics quality)
* Berkeley Advanced Reconstruction Toolbox (v0.4.01)

* 12-core Intel Xeon CPU E5-2687W v4
* 128 GB RAM

Installation Options:
---------------------
* Download the zip file of this repository to your local machine
* OR clone the git repository to your local machine
* OR fork to your own repository and then clone the git repository to your local machine


* Download contributing sample density compensation code (https://zenodo.org/record/401057/files/sdc3_nrz_11aug.zip) into `./contrib/`
* Download contributing gridding code (for acquisition time map generation) (https://zenodo.org/record/401036/files/grid3_dct_11aug.zip) into `./contrib/`
* Compile necessary mex files in respective contributing code directories:
    -(e.g. `mex -compatibleArrayDims grid3_MAT.c`)
    -(e.g. `mex -compatibleArrayDims sdc3_MAT.c`)
* Download contributing graph cut fat-water separation code (https://github.com/welcheb/FattyRiot/tree/master/FattyRiot_toolbox/hernando) into `./contrib/`
* Download contributing code for figure annotation and capture:
    -annotation https://www.mathworks.com/matlabcentral/fileexchange/278-arrow
    -capture https://github.com/altmany/export_fig/archive/f0af704d84608f5a69e3de82581869e7b6161d4f.zip	

Usage:
------
* Download input data (see Releases)
* Go to `./code/` and modify the setup parameters in the first section of `run_all.m`:
    -modify `bart_path` to point to BART path
    -modify other paths to other contributing code (if different from default)
* Go to `./code` and modify the setup parameters in the first section of `conventional_T1T2fatsep.m`:
    -modify `BASEPATH` to point to Hernando et al.'s graph cut fat-water separation code
* Run `run_all.m`
* (Figure generation: Go to `./code/` and modify the first section of `figure_generation.m` to correctly point to contributing annotation and figure export code)


* Code is not optimized for speed, the total processing time maybe ~days if all processing is performed
* The MRF fat-water separation algorithm is memory intensive, requiring ~100 GB of memory in some instances
* Figure generation quality is display dependent

Folder Structure:
--------

* `./code/` - contains all code (exluding contributing code noted above) necessary to reconstruct, process, and generate figures
* `./contrib/` - the downloaded contributing code
* `./data_in/` - the data input directory
* `./data_out/` - the reconstruction and processing output directory
* `./figures/` - the figure output directory

