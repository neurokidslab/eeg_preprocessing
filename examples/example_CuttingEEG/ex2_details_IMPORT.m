%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code exemplifies how to preprocess continuos EEG data using APICE
% 
% Data is organized in one folder per subject where all the pre-processing
% steps are saved
%
% Ana Fló, December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, close all
restoredefaultpath
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PATHS
%
% Specify the paths to the different relevant folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder where EEGLAB is:
Path2EEGLAB = 'C:\Users\anafl\Documents\MATLAB\eeglab2022.1';
% Folder where the function for APICE are:
Path2APICE = 'C:\Users\anafl\Documents\GitHub\eeg_preprocessing';

% Folder where the main scripts are:
Path0 = 'C:\Users\anafl\Documents\GitHub\eeg_preprocessing\examples\example_CuttingEEG';
% Folder where the scripts with the parameters are:
Path2Parameters = fullfile(Path0, 'parameters');
% Folder where the data is:
Path2Data = fullfile(Path0, 'DATA');
% File with the hannels layout 
filechanloc = fullfile(Path0,'ElectrodesLayout','GSN-HydroCel-128.sfp');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD TOOLBOXES AND PATHS TO USEFUL FOLDERS
%
% Add to the path the required toolboxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add EEGLAB to the path
cd(Path2EEGLAB)
eeglab
close all

% Add the functions for APICE (https://github.com/neurokidslab/eeg_preprocessing)
addpath(genpath(Path2APICE))

% Add the path to the folder with the parameters
addpath(Path2Parameters)

% Got to the path where the main scripts are
cd(Path0)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPORT 
%
% Import the raw data to EEGLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% .........................................................................
% Raw files names
% .........................................................................
rawfiles = '*.raw';

% Generate a list of files to analyze 
SBJs = eega_getfilesinfolders(Path2Data, rawfiles);

% Loop over all the files
for sbj=1:length(SBJs)
    
    % ---------------------------------------------------------------------
    % Import
    EEG = eega_importdata(SBJs{sbj});

    % ---------------------------------------------------------------------
    % Save the continuos pre-process data
    [sbjfolder,name,~] = fileparts(SBJs{sbj});
    newname = ['raw_' name '.set'];
    pop_saveset(EEG, newname, sbjfolder);
end
