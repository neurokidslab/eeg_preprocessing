% .........................................................................
% This code shows how to perform the artifacts' correction on a dataset imported to EEGLAB.  
% Artifacts need to have already been detected.
%
% Ana Fló, August 2021
% .........................................................................

clear, close all
restoredefaultpath
clc

%% ========================================================================
%% PATHS
% Specify the paths to the different relevant folders

% Folder where EEGLAB is:
Path2EEGLAB = 'xxx';
% Folder where the function for APICE are:
Path2APICE = 'xxx';
% Current path:
Path0 = '\examples';
% Folder where the scripts defining the parameters are:
Path2Parameters = fullfile(Path0,'parameters');
% Folder where the data imported to EEGLAB is: 
Path2InputData = 'xxx';
% Folder where the continuos preprocess data will be saved
Path2OutputData = 'xxx';


%% ========================================================================
%% ADD TOOLBOXES AND PATHS TO USEFUL FOLDERS

% Add EEGLAB to the path
cd(Path2EEGLAB)
eeglab
close all
% Add the functions for APICE (https://github.com/neurokidslab/eeg_preprocessing)
addpath(genpath(Path2APICE))
% Add the path to the folder with the parameters
addpath(Path2Parameters)
% Got to the original path
cd(Path0)


%% ========================================================================
%% PARAMETERS
% Define the parameters

% -------------------------------------------------------------------------
% Parameters for importing the data
% -------------------------------------------------------------------------

% Name of the files to read from the Path2InputData folder
Ppp.Files2Read = '*.set';


% -------------------------------------------------------------------------
% Parameters for transient artifacts interpolation
% -------------------------------------------------------------------------

% Generate the structures with the parameters for the artifact correction
% The functions to generate the structures should be in the path (in Path2Parameters)
% This function should be modified to change the algorithms applied for
% artifacts interpolation
Ppp.Int = example_APICE_Interpolation;

% Create a cell array with all the input to use later
inputsInterPCA = {Ppp.Int.PCA.nSV, Ppp.Int.PCA.vSV,...
    'maxTime', Ppp.Int.PCA.maxTime,'maskTime', Ppp.Int.PCA.maskTime,'splicemethod', Ppp.Int.PCA.splicemethod};

inputsInterSegmentsSpline = {Ppp.Int.Spl.p, 'pneigh', Ppp.Int.Spl.pneigh, 'splicemethod', Ppp.Int.Spl.splicemethod,...
    'mingoodtime', Ppp.Int.Spl.minGoodTime, 'minintertime', Ppp.Int.Spl.minInterTime, 'masktime', Ppp.Int.Spl.maskTime};


% -------------------------------------------------------------------------
% Parameters to define Bad Times (BT) and Bad Channels (BC) 
% -------------------------------------------------------------------------

% Limits for the proportion of BT to define a BC (the last value is the final/effective one)
Ppp.BCall.nbt           = [0.70 0.50 0.30];
% Limits for the proportion of BC to define a BT (the last value is the final/effective one)
Ppp.BTall.nbc           = [0.70 0.50 0.30];
% Shorter intervals between bad segments will be marked as bad
Ppp.BTall.minGoodTime   = 1.000;
% Shorter periods will not be considered as bad
Ppp.BTall.minBadTime    = 0.100;   
% Also mark as bad surronding samples within this value 
Ppp.BTall.maskTime      = 0.500;        

% Create a cell array with all the input to use later
inputsDefBTBC = {Ppp.BTall.nbc, Ppp.BCall.nbt, Ppp.BCall.nbt, 'keeppre', 0,...
    'minBadTime', Ppp.BTall.minBadTime, 'minGoodTime', Ppp.BTall.minGoodTime, 'maskTime', Ppp.BTall.maskTime};


% -------------------------------------------------------------------------
% Parameters to generate the report
% -------------------------------------------------------------------------

% Pattern to find in the names to create the table (subjetcs identifier). 
% If empty, the whole files name are used
Ppp.report.patternname = 'data';  
% Number of caracters to use from the begiging of the patter to create the table. 
% If empty, the name is extracted from the begign of the patter till the
% end of the name
Ppp.report.patternn = [];

%% ========================================================================
%% RUN THE SET OF FUNCTIONS FOR ALL THE FILES

% Name of the files to run the processes
FilesIn = 'art_*.set';

% String to add at the begign of the output filename
FilesOut = 'int_';

% Define how to run eega_RunAll
% - If 1, run to all the files in the input folder
% - If 2, run to the files in the input folder do not present in the ouput folder (files still not analyzed)
runto = 1;

% Import all the files in the input folder

eega_RunAll(FilesIn,...
            Path2InputData,...
            Path2OutputData,...
            'eega_tDefBTBC',                inputsDefBTBC,...
            'eega_tTargetPCAxElEEG',        inputsInterPCA,...
            'eega_demean',                  {},...
            'pop_eegfiltnew',               {Ppp.filt_highpass, [], [], 0, [], [], 0},...
            'eega_tDefBTBC',                inputsDefBTBC,...
            'eega_tInterpSpatialSegmentEEG',inputsInterSegmentsSpline,...
            'eega_demean',                  {},...
            'pop_eegfiltnew',               {Ppp.filt_highpass, [], [], 0, [], [], 0},...
            'eega_tInterpSpatialEEG',       {Ppp.Int.Spl.p, 'pneigh', Ppp.Int.Spl.pneigh},...
            'prename',FilesOut,'runto', runto);

% -------------------------------------------------------------------------
% Print a preprocessing report for all the subjetcs
eega_printsummary('int_*.set', Path2OutputData, Ppp.report.patternname, Ppp.report.patternn, 0);

