We provide different example codes.

An example data file already imported to the EELAB format is located in \examples\DATA\set
the corresponding event file is located in \examples\DATA\evt

=====================================
A. Detailed example code
=====================================

-------------------------------------
A.1 example_details
-------------------------------------
This code exemplifies how to run APICE more transparently and with a detailed explanation of each step.
The parameters a the begin allow avoiding steps that are omitted in the different versions of the pipeline described in the manuscript (Fló et al. 2021) 
The code first exemplified how to obtained continuous preprocess data and then how to obtain ERPs.

Notice that the code runs on data already imported to the EEGLAB format. 
To import data to the EEGLAB format, use the pertinent EEGLAB toolbox or see B.1 to import data acquired using EGI


=====================================
B. Codes to run different version 
of the pipeline for multiple subjects 
using the function eega_RunAll
=====================================

These codes (named example_mlt_*) show how to use the function eega_RunAll to apply a set of functions to a dataset.

eega_RunAll runs the specified functions (with their inputs) on the files in the inputs folder matching the input data name
The advantage of the function is that the process can be applied:
- on all the input files found
- on the input files for which the output file does not exist (new files)
- on the files for which the user decides to apply it 

Different codes are provided to run the different versions of the pipeline described in the manuscript (Fló et al. 2021)

All the codes take as input the raw data imported to the EEGLAB format

Additionally, we provide a code showing how to import the data recorded using the EGI system.
Data recorded using other systems can be imported using the pertinent EEGLAB toolbox.  


-------------------------------------
B.1 example_mlt_impordata
-------------------------------------
This code shows how to import data recorded using the EGI system to EEGLAB.


-------------------------------------
B.2 example_mlt_APICE_preprocessing
-------------------------------------
This code shows how to run the APICE pipeline to obtained continuos preprocess data


-------------------------------------
B.3 example_mlt_APICEa_preprocessing
-------------------------------------
This code shows how to run the APICEa pipeline to obtained continuos preprocess data


-------------------------------------
B.4 example_mlt_APICE_WICA_preprocessing
-------------------------------------
This code shows how to run the APICE+W-ICA pipeline to obtained continuos preprocess data

NOTE: We found some compatibility problems in iMARA regarding the electrodes layouts 
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


-------------------------------------
B.5 example_mlt_APICE_ERPs
-------------------------------------
This code takes continuous preprocessed data and performs all the steps to obtain the ERPs


-------------------------------------
B.6 example_mlt_APICE_ERPs_DSS
-------------------------------------
This code takes continuous preprocessed data and performs all the steps to obtain the ERPs.
Before averaging applies the DSS to clean the data

=====================================
C. Codes to run different parts 
of the pipeline for multiple subjects 
using the function eega_RunAll
=====================================

These codes (named example_step_*) show how to use the function eega_RunAll to apply a set of functions dedicated to a specific step to a dataset.

Different codes are provided to run the different processes.

All the codes take as input data in the EEGLAB format.

-------------------------------------
C.1 example_step_artifactsdetection
-------------------------------------
This code shows how to perform the artifacts' detection on a dataset imported to EEGLAB.

-------------------------------------
C.2 example_step_artifactscorrection
-------------------------------------
This code shows how to perform the artifacts' correction on a dataset imported to EEGLAB.  
Artifacts need to have already been detected.

=====================================
D. Extras
=====================================

-------------------------------------
D.1 example_ch_iMARA
-------------------------------------
This script shows how to construct a cell to use for W-ICA to convert the labels of the channels to the 10-10 standard system used by MARA/iMARA  