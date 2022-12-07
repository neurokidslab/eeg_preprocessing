% .........................................................................
% This code shows how to run the APICE pipeline to obtained continuos preprocess data
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
Path2EEGLAB = fullfile('xxx\toolbox','eeglab2020_0');
% Folder where the function for APICE are:
Path2APICE = fullfile('xxx','eeg_preprocessing');
% Current path:
Path0 = 'xxx\eeg_preprocessing\examples';
% Folder where the scripts defining the parameters are:
Path2Parameters = fullfile(Path0,'parameters');
% Folder where the continuos preprocess data will be saved
Path2DataPrp = fullfile(Path0,'DATA','prp');
% Folder where the ERPs will be saved
Path2DataERP = fullfile(Path0,'DATA','erp');


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
% Parameters for filtering
% -------------------------------------------------------------------------

% High pass filter
Perp.filt_highpass = 0.2;

% Low pass filter
Perp.filt_lowpass = 20;


% -------------------------------------------------------------------------
% Parameters for epoching
% -------------------------------------------------------------------------

% Time window to epoch in seconds
Perp.Epoch.tw = [-1.600, 2.200];  

% Events relative to which the data epoched 
Perp.Epoch.ev = {'Iout'};  

% -------------------------------------------------------------------------
% Parameters for artifacts interpolation
% -------------------------------------------------------------------------

% Generate the structures with the parameters for the artifact correction
% The functions to generate the structures should be in the path (in Path2Parameters)
% This function should be modified to change the algorithms applied for
% artifacts interpolation
Perp.Int = example_APICE_Interpolation;

% -------------------------------------------------------------------------
% Parameters to define Bad Times (BT) and Bad Channels (BC) 
% -------------------------------------------------------------------------

% Limits for the proportion of BT to define a BC during the whole recording (the last value is the final/effective one, here 30 %)
Perp.BCall.nbt           = [0.70 0.50 0.30 0.30];
% Limits for the proportion of BT to define a BC during each epoch (the last value is the final/effective one, here 100 ms)
Perp.BCep.nbt            = [0.70 0.50 0.30 0.10/diff(Perp.Epoch.tw)];   
% Limits for the proportion of BC to define a BT (the last value is the final/effective one, here 30 %)
Perp.BTep.nbc            = [0.70 0.50 0.30 0.30]; 
% Shorter intervals between bad segments will be marked as bad
Perp.BTep.minGoodTime    = 1.00;            
% Shorter periods will not be considered as bad
Perp.BTep.minBadTime     = 0.100;            
% Also mark as bad surronding samples within this value 
Perp.BTep.maskTime       = 0;                

% Create a cell array with all the input to use later
inputsDefBTBC = {Perp.BTep.nbc, Perp.BCep.nbt, Perp.BCall.nbt, 'keeppre', 0,...
    'minBadTime', Perp.BTep.minBadTime, 'minGoodTime', Perp.BTep.minGoodTime, 'maskTime', Perp.BTep.maskTime}; 


% -------------------------------------------------------------------------
% Parameters to define Bad Epochs (BE) based on the amount of bad data
% -------------------------------------------------------------------------

% Maximun proportion of bad data per epoch  
Perp.DefBEa.limBCTa      = 1.00;   
% Maximun proportion of bad times per epoch  
Perp.DefBEa.limBTa       = 0.00;   
% Maximun proportion of bad channels per epoch  
Perp.DefBEa.limBCa       = 0.30;   
% Maximun proportion of interpolated data per epoch  
Perp.DefBEa.limCCTa      = 0.50;   

% -------------------------------------------------------------------------
% Parameters for defining experimental factors and averaging
% -------------------------------------------------------------------------

% Define based on the events the factors that will be used to determine conditions
% They should be events' properties (e.g., EEG.epoch(i).eventCong) 
% This factors can be the used for averaging across one or multiple of them to obtain the different experimental conditions
% The factors are defined in EEG.F For each factor i it contains:
% - EEG.F{i}.name: name of the factor (e.g., 'eventCong')
% - EEG.F{i}.val: possible values of the factor (e.g., {'x0'  'x1'})
% - EEG.F{i}.g: vector with length equal to the number of epochs with the indexes
% to EEG.F{i}.val indicating the value of the factor (EEG.F.val(EEG.F{i}.g(k)) is the the value of the facto i for epoch k)
Perp.factors = {'eventCong' 'eventProb'};

% Parameters for averaging across some factors to obtain the ERPs for different conditions
% After averaging, EEG.data has size (channles) x (samples) x (possible conditions) 
Perp.avgcnd = {'eventCong' 'eventProb'};

% -------------------------------------------------------------------------
% Parameters to perform DSS
% -------------------------------------------------------------------------

% Run or not DSS
% (1) apply | (0) do not apply
Perp.DSS.apply      = 0;  
% components to keep in the first PCA
Perp.DSS.k          = 50; 
% components to keep in the second PCA
Perp.DSS.n          = 15; 
% Define the trials to use to bias the filter. 
% Different filters can be used for different trials sets
% If empty, one filter is created using all trials
Perp.DSS.fbias      = [];
% Define the trials to which the DSS is applied 
% Different filters can be applied for different trials sets
% If empty, the single filter is applied to all trials
Perp.DSS.fapply     = [];


% -------------------------------------------------------------------------
% Parameters for baseline correction
% -------------------------------------------------------------------------

% Time windos (in ms) to compute the basiline
Perp.BL.tw = [-100 100];  


% -------------------------------------------------------------------------
% Parameters to generate the report
% -------------------------------------------------------------------------

% Pattern to find in the names to create the table (subjetcs identifier). 
% If empty, the whole files name are used
Perp.report.patternname = 'data';  
% Number of caracters to use from the begiging of the patter to create the table. 
% If empty, the name is extracted from the begign of the patter till the
% end of the name
Perp.report.patternn = [];

%% ========================================================================
%% RUN THE SET OF FUNCTIONS FOR ALL THE FILES

% Name of the files to run the processes
FilesIn = '*.set';

% String to add at the begign of the output filename
FilesOut = 'erp_';

% Define how to run eega_RunAll
% - If 1, run to all the files in the input folder
% - If 2, run to the files in the input folder do not present in the ouput folder (files still not analyzed)
runto = 1;

% Import all the files in the input folder

eega_RunAll(FilesIn,...
            Path2DataPrp,...
            Path2DataERP,...
            'pop_eegfiltnew',           {[], Perp.filt_lowpass,  [], 0, [], [], 0},...
            'pop_eegfiltnew',           {Perp.filt_highpass, [], [], 0, [], [], 0},...
            'eega_epoch',               {Perp.Epoch.ev, Perp.Epoch.tw},...
            'eega_definefactors',       {Perp.factors},... 
            'eega_tDefBTBC',            inputsDefBTBC,...
            'eega_tInterpSpatialEEG',   {Perp.Int.Spl.p,'pneigh', Perp.Int.Spl.pneigh},...
            'eega_tDefBTBC',            inputsDefBTBC,...
            'eega_tDefBEbaddata',       {Perp.DefBEa, 'keeppre', 0, 'plot', 0},...
            'eega_rmvbadepochs',        {},...
            'pop_reref',                {[]},...
            'dss_denoise_EEG',          {Perp.DSS.k, Perp.DSS.n, Perp.DSS.fbias, Perp.DSS.fapply},...
            'eega_rmbaseline',          {Perp.BL.tw},...
            'eega_avgdatabyfactors',    {Perp.avgcnd, 'dim2avg', 3},...
            'prename',FilesOut,'runto', runto);

% -------------------------------------------------------------------------
% Print a preprocessing report for all the subjetcs
eega_printsummary('erp_*.set', Path2DataERP, Perp.report.patternname, Perp.report.patternn, 0);

