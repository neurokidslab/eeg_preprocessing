% .........................................................................
% This code shows how to import data recorded using the EGI system to EEGLAB.
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
Path2EEGLAB = fullfile('xxx','eeglab2020_0');
% Folder where the function for APICE are:
Path2APICE = fullfile('xxx','eeg_preprocessing');
% Current path:
Path0 = 'xxx\eeg_preprocessing\examples';
% Folder where the raw data is: 
Path2DataRaw = fullfile(Path0,'DATA','raw');
% Folder where the data imported to EEGLAB will be saved
Path2DataSet = fullfile(Path0,'DATA','set');


%% ========================================================================
%% ADD TOOLBOXES AND PATHS TO USEFUL FOLDERS

% Add EEGLAB to the path
cd(Path2EEGLAB)
eeglab
close all
% Add the functions for APICE (https://github.com/neurokidslab/eeg_preprocessing)
addpath(genpath(Path2APICE))
% Got to the original path
cd(Path0)

%% ========================================================================
%% IMPORT THE RAW DATA TO EEGLAB FOR ALL THE SUBJECTS

% Name of the files to import from the folder FolderDataRaw
FilesIn = '*.raw';

% Define how to run eega_RunAll
% - If 1, run to all the files in the input folder
% - If 2, run to the files in the input folder do not present in the ouput folder (files still not analyzed)
runto = 1;

% Import all the files in the input folder
eega_RunAll(FilesIn,...
            Path2DataRaw,...
            Path2DataSet,...
            'eega_importdata', {},...
            'prename','','runto', runto);
