%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code exemplifies how to obtained ERPs from the continuos preprocess
% data obtained using APICE
% 
% Data is organized in one folder per subject where all the pre-processing
% steps and further steps are saved there
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
%% PARAMETERS
%
% Define all the required parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% .........................................................................
% Pre-processed files names
% .........................................................................
prpfiles = 'prp_raw_*.set';
epochfiles = 'epch_';
erpfiles = 'erp_';

% .........................................................................
% Parameters for filtering
% .........................................................................
Perp.filt = ex2_Perp_Filtering;

% .........................................................................
% Parameters for epoching
% .........................................................................
% Time window to epoch in seconds
Perp.Epoch.tw = [-1.600, 2.200];  
% Events relative to which the data epoched 
Perp.Epoch.ev = {'Iout'};  

% .........................................................................
% Parameters for artifacts interpolation using spherical spline
% .........................................................................
Perp.Int.spl = ex2_Perp_Int_Spl;

% .........................................................................
% Parameters to define Bad Times (BT) and Bad Channels (BC) 
% .........................................................................
Perp.BTBC = ex2_Perp_DefBTBC;

% .........................................................................
% Parameters to define Bad Epochs (BE) based on the amount of bad data
% .........................................................................
Perp.BE = ex2_Perp_DefBE;

% .........................................................................
% Parameters for defining experimental factors and averaging
% .........................................................................

% Define based on the events the factors that will be used to determine conditions
% They should be events' properties (e.g., EEG.epoch(i).eventCong) 
% This factors can be the used for averaging across one or multiple of them to obtain the different experimental conditions
% The factors are defined in EEG.F For each factor i it contains:
% - EEG.F{i}.name: name of the factor (e.g., 'eventCong')
% - EEG.F{i}.val: possible values of the factor (e.g., {'x0'  'x1'})
% - EEG.F{i}.g: vector with length equal to the number of epochs with the indexes
% to EEG.F{i}.val indicating the value of the factor (EEG.F.val(EEG.F{i}.g(k)) is the the value of the facto i for epoch k)
Perp.Factors = {'eventCong' 'eventProb'};

% Parameters for averaging across some factors to obtain the ERPs for different conditions
% After averaging, EEG.data has size (channles) x (samples) x (possible conditions) 
Perp.Avg = {'eventCong' 'eventProb'};

% .........................................................................
% Parameters to perform DSS
% .........................................................................
Perp.dss = ex2_Perp_DSS;
do_applyDSS = 0; % run (1) or not (0)

% .........................................................................
% Parameters for baseline correction
% .........................................................................
% Time windos (in ms) to compute the basiline
Perp.BL.tw = [-100 0];  

% .........................................................................
% Plotting
% .........................................................................
do_plotrejection = 0; % run (1) or not (0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EPOCHS
%
% Epoch the pre-processed data for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================================================================
%% EPOCH THE DATA FOR ALL THE SUBJECTS

% Generate a list of files to analyze 
SBJs = eega_getfilesinfolders(Path2Data, prpfiles);

% Loop over all the files
SBJs_epch = cell(size(SBJs));
for sbj=1:length(SBJs)
    
    % ---------------------------------------------------------------------
    % Load the data
    EEG = pop_loadset(SBJs{sbj});
    
    % ---------------------------------------------------------------------
    % Filter
    
    % low-pass
    EEG = pop_eegfiltnew(EEG, [], Perp.filt.lowpass,  [], 0, [], [], 0);
    
    % high-pass
    EEG = pop_eegfiltnew(EEG, Perp.filt.highpass, [], [], 0, [], [], 0);
    
    % ---------------------------------------------------------------------
    % Epoch
    EEG = eega_epoch(EEG, Perp.Epoch.ev, Perp.Epoch.tw);
    
    % ---------------------------------------------------------------------
    % Define the factors that will enable to define the conditions
    EEG = eega_definefactors(EEG, Perp.Factors);
    
    % ---------------------------------------------------------------------
    % Define Bad Times (BT) and Bad Channels (BC)
    EEG = eega_tDefBTBC(EEG, Perp.BTBC.BT.nbc, Perp.BTBC.BCep.nbt, Perp.BTBC.BCall.nbt,...
        'minBadTime', Perp.BTBC.BT.minBadTime, 'minGoodTime', Perp.BTBC.BT.minGoodTime,...
        'maskTime', Perp.BTBC.BT.maskTime,'keeppre', 0);
    
    % ---------------------------------------------------------------------
    % Apply spherical spline interpolation on each each epoch
    EEG = eega_tInterpSpatialEEG(EEG, Perp.Int.spl.p,'pneigh', Perp.Int.spl.pneigh);
    
    % ---------------------------------------------------------------------
    % Define Bad Times (BT) and Bad Channels (BC)
    EEG = eega_tDefBTBC(EEG, Perp.BTBC.BT.nbc, Perp.BTBC.BCep.nbt, Perp.BTBC.BCall.nbt,...
        'minBadTime', Perp.BTBC.BT.minBadTime, 'minGoodTime', Perp.BTBC.BT.minGoodTime,...
        'maskTime', Perp.BTBC.BT.maskTime,'keeppre', 0);
    
    % ---------------------------------------------------------------------
    % Define bad epochs
    EEG = eega_tDefBEbaddata(EEG, Perp.BE, 'keeppre', 0, 'plot', 0);
    
    % ---------------------------------------------------------------------
    % Plot the rejection matrix
    if do_plotrejection
        % Plot the rejection matrix
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts rejection 2')
    end
    
    % ---------------------------------------------------------------------
    % Reference average
    EEG = pop_reref(EEG,[]);
    
    % ---------------------------------------------------------------------
    % DSS
    if do_applyDSS
        EEG = dss_denoise_EEG(EEG, Perp.dss.k, Perp.dss.n, Perp.dss.fbias, Perp.dss.fapply);
    end
    
    % ---------------------------------------------------------------------
    % Save the epochs
    [sbjfolder,name,~] = fileparts(SBJs{sbj});
    newname = [epochfiles name '.set'];
    pop_saveset(EEG, newname, sbjfolder);
    SBJs_epch{sbj} = fullfile(sbjfolder,newname);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ERPs
%
% Average trials across factors of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a list of files to analyze 
SBJs = eega_getfilesinfolders(Path2Data, [epochfiles prpfiles]);

% Loop over all the files
SBJs_erp = cell(size(SBJs));
for sbj=1:length(SBJs)
    
    % ---------------------------------------------------------------------
    % Load the data
    EEG = pop_loadset(SBJs{sbj});
    
    % ---------------------------------------------------------------------
    % Baseline correction
    EEG = eega_rmbaseline(EEG, Perp.BL.tw);
    
    % ---------------------------------------------------------------------
    % Remove the bad epochs
    EEG = eega_rmvbadepochs(EEG);
    
    % ---------------------------------------------------------------------
    % Average per conditions
    EEG = eega_avgdatabyfactors(EEG, Perp.Avg);

    % ---------------------------------------------------------------------
    % Save the ERPs
    [sbjfolder,name,~] = fileparts(SBJs{sbj});
    newname = [erpfiles name '.set'];
    pop_saveset(EEG, newname, sbjfolder);
    SBJs_erp{sbj} = fullfile(sbjfolder,newname);

end





