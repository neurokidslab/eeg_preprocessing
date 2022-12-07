These scripts exemplify how to pre-process data using APICE.

The data to analyze is located in the folder 'DATA'. For each subject, there is one folder where all further pre-processing steps are saved (BIDS format).  

An example data file already imported to the EELAB format is located in \examples\DATA\S01

=====================================
ex2_details_IMPORT
=====================================

This code is an example of importing .raw data to EEGLAB (.set/.fdt). 
Other datafiles can be imported using the EEGLAB functions.

=====================================
ex2_details_PP.m
=====================================

This code is an example of how to pre-process continuous data using APICE. 

It takes input data already imported to the EEGLAB format (.det/.fdt) and loops over all the subject folders. 

All the parameters are set at the beginning based on the scripts in the parameters folder. 
The parameters for the artifacts' detection are similar to those used in (Fló et al. 2021). The only differences are:
- Artifacts are not masked after each algorithm is applied to avoid excessive rejection. Thus, rejection is restricted to samples identified as bad by the algorithms.
- A mask 0.5 s is applied when defining bad times (BT). 
- The threshold is set to 2 for all algorithms besides eega_tRejAmpElecVar, which uses 3.
- Include rejected segments shorter than 40 ms (instead of 20 ms). Thus, the segments on which target PCA is applied are at least 40 ms long. 
- Reject good segments between bad segments shorter than 0.5 s (instead of 2 s)
- The last artifact rejection does not overwrite the previous rejection.


NOTE for APICE+W-ICA
-------------------------------------
We found some compatibility problems in iMARA regarding the electrodes layouts 
To solve this 
1) We indicate how to change the names of the electrodes to be compatible with the standard system using the function eega_changechlable.
See the example script to see how to indicate it in the structure input
2) In the function iMARA from the iMARA toolbox we commented the lines 
clab(7) = {'T7'};
clab(24) = {'T8'};
i = contains(loc_labels, clab(7)); loc_labels(i) = {EEG.chanlocs(7).labels};
j = contains(loc_labels, clab(24)); loc_labels(j) = {EEG.chanlocs(24).labels};
3) In the function getchanloc from the iMARA toolbox, we replaced 
EEGchanslabels = {'Fp1';'AF3';'F7';'F3';'FC1';'FC5';'T7(T3)';'C3';'CP1';'CP5';'P7';'P3';'Pz';'PO3';'O1';'Oz';'O2';'PO4';'P4';'P8';'CP6';'CP2';'C4';'T8(T4)';'FC6';'FC2';'F4';'F8';'AF4';'Fp2';'Fz';'Cz'};
by 
EEGchanslabels = {'Fp1';'AF3';'F7';'F3';'FC1';'FC5';'T7';'C3';'CP1';'CP5';'P7';'P3';'Pz';'PPO1';'O1';'Oz';'O2';'PPO2';'P4';'P8';'CP6';'CP2';'C4';'T8';'FC6';'FC2';'F4';'F8';'AF4';'Fp2';'Fz';'Cz'};

ex_ch_iMARA
This script shows how to construct a cell to use for W-ICA to convert the labels of the channels to the 10-10 standard system used by MARA/iMARA  

=====================================
ex2_details_ERPs.m
=====================================

This code is an example of how to epoch data and obtain ERPs. 

It takes input continuous data pre-process using APICE and loops over all the subject folders. 


=====================================
ex2_compact_PP
=====================================

This code shows how to pre-process the data using the function eega_RunAll_BIDS.

eega_RunAll_BIDS runs a set of functions on the data located in each folder within the data folder. 
The functions to run and their parameters are provided as inputs. 
It runs on 
- all the input files 
- only ththe input files not yet analyzed (for which an output file is not found in the subject data folder). 
- on the input files for which the user decides to apply it 

