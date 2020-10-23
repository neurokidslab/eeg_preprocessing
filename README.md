# EEG preprocessing for infant data

> Authors: Ana Flo (anaflom@gmail.com)

This is the code repository for pre-preprocessing developmental EEG data of the UNICOG BabyLab, NeuroSpin.
All routines are implemented in MATLAB.

#### Requirements
* <a href="https://mathworks.com/" target_="blank">MATLAB</a>
* <a href="https://sccn.ucsd.edu/eeglab/" target_="blank">EEGLAB toolbox</a> (recommend >=14.0)
* <a href="http://audition.ens.fr/adc/NoiseTools/" target_="blank">NoiseTools</a>, only necessary to perform robust detrending
* <a href="https://github.com/irenne/MARA" target_="blank">MARA EEGLAB plugin</a>, only necessary to automatically select IC after ICA
* <a href="https://www.fieldtriptoolbox.org/" target_="blank">FieldTrip toolbox</a>, only necessary to perform a cluster based permutation analysis
