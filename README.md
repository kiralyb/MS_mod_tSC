# MATLAB code for 'The medial septum controls hippocampal supra-theta oscillations'

## Description

This repository contains MATLAB code for data analysis on the connection between theta-nested spectral components (tSCs), single unit activities and optogenetic stimulus events. Single unit clusters and hippocampal tSCs were extracted with Kilosort (https://github.com/MouseLand/Kilosort) and the tSC extraction package (https://data.mrc.ox.ac.uk/data-set/tsc).

The figures and supplementary figures of 'The medial septum modulates hippocampal oscillations beyond the theta rhythm' by Király et al. can be produced with the MS_mod_tSC_main function relying on four datasets:
- freely moving mouse (hippocampal LFP - septal units) 
- anesthetized rat (hippocampal LFP - septal units)
- anesthetized mouse (hippocampal LFP - septal units)
- optogenetic stimulation of PV experessing septal neurons in mice (hippocampal LFP, hippocampal units - septal stimulus events)

Demo session data, showing the required data format and hierarchy, is avaibale at the following link: https://figshare.com/articles/dataset/Demo_data_for_the_MATLAB_code_for_analyzing_the_connection_between_the_medial_septum_and_hippocampal_oscillations_beyond_the_theta_rhythm_/22060964
(content: raw downsampled lfp; position data of the animal (optional); stimulation times (optional); preprocessed data after spike sorting and tSC extraction)
(e.g. tSC_run_session_analyses ('L:\Balint\tsc\awake_mouse\','20161989','107108','10',1,1,'awake_mouse') for the demo session, see code documentation for further details on the inputs)

Session data can be preprocessed with the tSC_run_session_analyses function. Preprocessed data is stored in the Matrix.mat and ses_Matrix.mat files in the root directory.

## Content

- .m files for data analysis
- license file

## Installation

Move the .m files on your MATLAB path (should take around a few seconds). No further installation is needed. 

## Dependencies

This code uses functions from the hangya-matlab-code package for spike triggered average analyses (https://github.com/hangyabalazs/Hangya-Matlab-code).

Classification of MS neurons based on their rhythmicity was performed with the ms_sync_analysis package (https://github.com/hangyabalazs/ms_sync_analysis). 

Optogenetically stimulated hippocampal neurons were examined using Cellbase (https://github.com/hangyabalazs/CellBase). Note that when running the first time, Cellbase needs to be initialized first.

Not ordinary data formats are loaded with the (https://github.com/kwikteam/npy-matlab/blob/master/npy-matlab/readNPY.m) and the (https://github.com/open-ephys/analysis-tools/blob/master/load_open_ephys_data.m) functions.

## System requirements

- Windows 10 64bits
- Intel i7
- 32 GB RAM
- MatlabR2016a (Signal Processing Toolbox, Statistics and Machine Learning Toolbox, Wavelet Toolbox, Bioinformatics Toolbox)
- Tested on MatlabR2016a and MatlabR2018
