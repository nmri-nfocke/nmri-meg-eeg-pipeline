# nmri-meg-eeg-pipeline

This is a collection of Matlab scripts deviced my member of AG Focke (https://neurologie.umg.eu/forschung/arbeitsgruppen/epilepsie-und-bildgebungsforschung/) wrapping around Fieldtrip and other common tool for scripted MEG and EEG analyis.

The provided repo should contain all relevant Matlab scripts to run, but requieres some external tools that need to be either placed in the Matlab path, or be included in the scripts/utilities directory.


## External Tools

-Fieldtrip (required, tested with version 20191127, 20220713, other version are likely to work as well), https://www.fieldtriptoolbox.org/

-SPM12 (required, tested with version 7487, other version are likely to work as well), https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

-Freesurfer (required, tested with version 6.0.0, 7.X versions are likely to work), https://surfer.nmr.mgh.harvard.edu/

-PALM (optional, for group comparisions, not needed by core scripts), https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM

## Installation / Setup

Download / clone the repo

Make sure the external tools are available or in the Matlab path

Change to the rootpath / cloned directory

`addpath(genpath('scripts'))`


## Analysis Steps

TBD

## Reference
Please cite either one of the below references if you use these scripts (or parts/concepts thereof) in your work:

_EEG+MEG in Epilepsy_
Stier, et al., Combined electrophysiological and morphological phenotypes in patients with genetic generalized epilepsy and their healthy siblings, Epilepsia 2022; 63(7):1643-57, doi:10.1111/epi.17258

_MEG in Epilepsy_
Stier, et al., Heritability of Magnetoencephalography Phenotypes Among Patients With Genetic Generalized Epilepsy and Their Siblings, Neurology 2021; 97(2):e166-e77, doi:10.1212/wnl.0000000000012144

_Technical Description and Reliability_
Marquetand, et al., Reliability of Magnetoencephalography and High-Density Electroencephalography Resting-State Functional Connectivity Metrics, Brain Connectivity 2019; 9(7):539-53, doi:10.1089/brain.2019.0662




