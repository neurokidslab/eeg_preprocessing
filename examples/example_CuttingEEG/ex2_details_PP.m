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
% Folder where iMARA is:
Path2iMARA = 'C:\Users\anafl\Documents\MATLAB\iMARA-main';

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

% Add iMARA (https://github.com/Ira-marriott/iMARA/tree/main)
addpath(genpath(Path2iMARA))

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
% Raw files names
% .........................................................................
rawfiles = 'raw_*.set';

% .........................................................................
% Parameters dealing with events 
% .........................................................................
do_arreangeevents = 0; % run (1) or not (0)
Ppp.events = ex2_Ppp_Events;

% .........................................................................
% Parameters for filtering
% .........................................................................
Ppp.filt = ex2_Ppp_Filtering;

% .........................................................................
% Parameters for artifacts detection
% .........................................................................
% Algorithms to detect bad electrodes
Ppp.Art.BadEl = ex2_Ppp_Art_BadEl;
% Algorithms to detect jumps in the signal
Ppp.Art.Jump = ex2_Ppp_Art_Jump;
% Algorithms to detect motion artifacts
Ppp.Art.Mot1 = ex2_Ppp_Art_Mot1;
% Algorithms to detect motion artifacts
Ppp.Art.Mot2 = ex2_Ppp_Art_Mot2;

% .........................................................................
% Parameters for transient artifacts interpolation using target PCA
% .........................................................................
Ppp.Int.pca = ex2_Ppp_Int_tPCA;

% .........................................................................
% Parameters for artifacts interpolation using spherical spline
% .........................................................................
Ppp.Int.spl = ex2_Ppp_Int_Spl;

% .........................................................................
% Parameters to define Bad Times (BT) and Bad Channels (BC) 
% .........................................................................
Ppp.BTBC = ex2_Ppp_DefBTBC;

% .........................................................................
% Parameters for W-ICA
% .........................................................................
Ppp.ica = ex2_Ppp_wICA;
do_applyICA = 1; % run (1) or not (0)

% .........................................................................
% Plotting
% .........................................................................
do_plotrejection = 0; % run (1) or not (0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRE-PROCESSING
%
% Run the pre-processing for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a list of files to analyze 
SBJs = eega_getfilesinfolders(Path2Data, rawfiles);

% Loop over all the files
SBJs_prp = cell(size(SBJs));
for sbj=1:length(SBJs)
   
    [sbjfolder,sbjname,~] = fileparts(SBJs{sbj});
    
    % ---------------------------------------------------------------------
    % Load the data
    EEG = pop_loadset(SBJs{sbj});
    

    % ---------------------------------------------------------------------
    % Import the channels location file
    EEG = pop_chanedit(EEG, 'load', {filechanloc 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );
    

    % ---------------------------------------------------------------------
    % Deal with events

    if do_arreangeevents

        % Add extra information from the event file
        if Ppp.events.import.apply
            EEG = eega_importinfoevents(EEG, sbjfolder, Ppp.events.import.ev0);
        end

        % Correct the latency of the events by the DINs
        if ~isempty(Ppp.events.corrlatencies.ev)
            for i=1:length(Ppp.events.corrlatencies.ev)
                EEG = eega_latencyevent(EEG, Ppp.events.corrlatencies.ev(i).event2, Ppp.events.corrlatencies.ev(i).event1);
            end
        end

        % Remove unuseful events
        if ~isempty(Ppp.events.rmv)
            EEG = eega_removeevent(EEG, Ppp.events.rmv.ev);
        end
    end

    % ---------------------------------------------------------------------
    % Filter
    
    % remove the mean
    EEG = eega_demean(EEG);
    
    % low-pass
    EEG = pop_eegfiltnew(EEG, [], Ppp.filt.lowpass,  [], 0, [], [], 0);
    
    % high-pass
    EEG = pop_eegfiltnew(EEG, Ppp.filt.highpass, [], [], 0, [], [], 0);
    
    % notch 
    EEG = pop_eegfiltnew(EEG, Ppp.filt.notch(1), Ppp.filt.notch(2), [], 1, [], [], 0);
    
    % ---------------------------------------------------------------------
    % Detect artifacts
    
    % Detect bad channels
    EEG = eega_tArtifacts(EEG, Ppp.Art.BadEl, 'KeepRejPre', 1);
    
    % Detect motion artifacts
    EEG = eega_tArtifacts(EEG, Ppp.Art.Mot1, 'KeepRejPre', 1);
    
    % Detect jumps in the signal
    EEG = eega_tArtifacts(EEG, Ppp.Art.Jump, 'KeepRejPre', 1);
    
    % Define Bad Times (BT) and Bad Channels (BC)
    EEG = eega_tDefBTBC(EEG, Ppp.BTBC.BT.nbc, Ppp.BTBC.BC.nbt, Ppp.BTBC.BC.nbt,...
        'minBadTime', Ppp.BTBC.BT.minBadTime,'minGoodTime', Ppp.BTBC.BT.minGoodTime,...
        'maskTime', Ppp.BTBC.BT.maskTime, 'keeppre', 0);
    
    % ---------------------------------------------------------------------
    % Plot the rejection matrix   
    if do_plotrejection      
        % Plot the rejection matrix 
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts rejection 1')
        % plot the data and bad segments and channels
        eega_plot_rejection(EEG, 1, 1, 1, 120)
        set(gcf, 'Name', 'artifacts rejection 1')
    end
    
    % ---------------------------------------------------------------------
    % Correct brief jumps in the signal using target PCA

    % Apply target PCA
    EEG = eega_tTargetPCAxElEEG(EEG, Ppp.Int.pca.nSV, Ppp.Int.pca.vSV,...
        'maxTime', Ppp.Int.pca.maxTime,'maskTime',...
        Ppp.Int.pca.maskTime,'splicemethod', Ppp.Int.pca.splicemethod);
    
    % High-pass filter the data to remove drifts
    EEG = eega_demean(EEG);
    EEG = pop_eegfiltnew(EEG, Ppp.filt.highpass, [], [], 0, [], [], 0);

    % Define Bad Times (BT) and Bad Channels (BC)
    EEG = eega_tDefBTBC(EEG, Ppp.BTBC.BT.nbc, Ppp.BTBC.BC.nbt, Ppp.BTBC.BC.nbt,...
        'minBadTime', Ppp.BTBC.BT.minBadTime,'minGoodTime', Ppp.BTBC.BT.minGoodTime,...
        'maskTime', Ppp.BTBC.BT.maskTime, 'keeppre', 0);
    
    % ---------------------------------------------------------------------
    % Correct channels not working during some time using spherical spline

    % Apply spherical spline interpolation
    EEG = eega_tInterpSpatialSegmentEEG(EEG, Ppp.Int.spl.p, 'pneigh',...
        Ppp.Int.spl.pneigh, 'splicemethod', Ppp.Int.spl.splicemethod,...
        'mingoodtime', Ppp.Int.spl.minGoodTime,...
        'minintertime', Ppp.Int.spl.minInterTime, 'masktime', Ppp.Int.spl.maskTime);

    % High-pass filter the data to remove drifts 
    EEG = eega_demean(EEG);
    EEG = pop_eegfiltnew(EEG, Ppp.filt.highpass, [], [], 0, [], [], 0);

    % ---------------------------------------------------------------------
    % Perform ICA

    if do_applyICA
        EEG = eega_pcawtica(EEG,'filthighpass',Ppp.ica.filthighpass, 'filtlowpass', Ppp.ica.filtlowpass,...
            'npc', Ppp.ica.npc, 'classifyIC', Ppp.ica.classifyIC, 'classifyICfun', Ppp.ica.classifyICfun,...
            'changelabelch',Ppp.ica.changelabelch,'labelch',Ppp.ica.labelch,...
            'saveica',Ppp.ica.saveica,'icaname',Ppp.ica.icaname,'icapath',Ppp.ica.icapath);
    end
    
    % ---------------------------------------------------------------------
    % Spatially interpolate channels not working during the whole recording

    % Apply spherical spline interpolation
    EEG = eega_tInterpSpatialEEG(EEG, Ppp.Int.spl.p, 'pneigh', Ppp.Int.spl.pneigh);

    % ---------------------------------------------------------------------
    % Detect artifacts

    % Detect motion artifacts
    EEG = eega_tArtifacts(EEG, Ppp.Art.Mot2, 'KeepRejPre', 1);
    
    % Define Bad Times (BT) and Bad Channels (BC)
    EEG = eega_tDefBTBC(EEG, Ppp.BTBC.BT.nbc, Ppp.BTBC.BC.nbt, Ppp.BTBC.BC.nbt,...
        'minBadTime', Ppp.BTBC.BT.minBadTime,'minGoodTime', Ppp.BTBC.BT.minGoodTime,...
        'maskTime', Ppp.BTBC.BT.maskTime, 'keeppre', 0);
    
    % ---------------------------------------------------------------------
    % Plot the rejection matrix
    if do_plotrejection
        % Plot the rejection matrix 
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts rejection 2')    
        % plot the data and bad segments and channels
        eega_plot_rejection(EEG, 1, 1, 1, 120)
        set(gcf, 'Name', 'artifacts rejection 2')
    end
    
    % ---------------------------------------------------------------------
    % Save the continuos pre-process data
    newsbjname = ['prp_' sbjname '.set'];
    pop_saveset(EEG, newsbjname, sbjfolder);
    SBJs_prp{sbj} = fullfile(sbjfolder,newsbjname);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REPORT
%
% Print a preprocessing report for all the subjetcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = eega_printsummary(SBJs_prp, Path2Data, [], [], 0);

