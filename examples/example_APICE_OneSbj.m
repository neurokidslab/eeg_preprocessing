% This is an example of a preprocessing pipeline for one subject 
% It can be extended to multiple subjetc loopinng over subjects files or
% using the the script for multiple subjects 
%
% --> to save data use: pop_saveset(EEG, FILE.set, FOLDER);
% --> to load data use: EEG = pop_loadset('FILE.set', FOLDER);


clear, close all
restoredefaultpath
clc

%% ========================================================================
%% PATHS
% Folder where EEGLAB is
Path2EEGLAB = 'C:\Users\an251296\Documents\MATLAB\toolbox\eeglab2020_0';
% Folder where the function for APICE are
Path2APICE = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eeg_preprocessing';
% Folder where iMARA is
Path2iMARA = 'C:\Users\an251296\Documents\MATLAB\toolbox\iMARA-main';
% Current path
Path0 = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eeg_preprocessing\examples';
% Folder where the data is and will be saved
Path2Data = fullfile(Path0,'DATA');
% Folder where the event files are
Path2DataEvent = fullfile(Path0,'DATA','evt');
% Path to the parameters
Path2Parameters = fullfile(Path0,'parameters');
% Channels location file
filechanloc128 = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eeg_preprocessing\electrodeslayout\GSN-HydroCel-128.sfp';


%% ========================================================================
%% ADD TOOLBOXES

% Add EEGLAB to the path
cd(Path2EEGLAB)
eeglab
close all
% Add the functions for APICE
addpath(genpath(Path2APICE))
% Add iMARA 
addpath(genpath(Path2iMARA))
% Add the path to the folder with the parameters
addpath(Path2Parameters)
% Got to the original path
cd(Path0)


%% ========================================================================
%% PRE-PROCESSING

%% ------------------------------------------------------------------------
%% Parameters for the preprocessing 
%% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parameters for filtering
% -------------------------------------------------------------------------
filt_highpass   = 0.1;
filt_lowpass    = 40;

% -------------------------------------------------------------------------
% Generate the structures with the parameters for the artifact detection
% -------------------------------------------------------------------------
% for the algorithms to detect bad electrodes
ArtBadEl = example_APICE_ArtPP_BadEl(3);
% for the algorithms to detect jumps in the signal
ArtJump = example_APICE_ArtPP_Jump(3);
% for the algorithms to detect motion artifacts
ArtMot1 = example_APICE_ArtPP_Mot1(3);
% for the algorithms to detect motion artifacts
ArtMot2 = example_APICE_ArtPP_Mot2(3);

% -------------------------------------------------------------------------
% Generate the structures with the parameters for the artifact correction
% -------------------------------------------------------------------------
Int = example_APICE_Interpolation;

% -------------------------------------------------------------------------
% Parameters to define Bad Times (BT) and Bad Channels (BC) 
% -------------------------------------------------------------------------
BCall.nbt           = [0.70 0.50 0.30];
BTall.nbc           = [0.70 0.50 0.30];
BTall.minGoodTime   = 1.000;        % shorter intervals between bad segments will be marked as bad
BTall.minBadTime    = 0.100;        % too short periods will not be considered as bad
BTall.maskTime      = 0.500;        % also mark as bad surronding samples

% -------------------------------------------------------------------------
% Parameters for ICA
% -------------------------------------------------------------------------
ICA.do          = 0;    % run or not ICA
ICA.npc         = 50;   % number of components to keep in the PCA
ICA.flthpass    = 2;    % high pass filter before ICA

%% ------------------------------------------------------------------------
%% Run the prerpocessing
%% ------------------------------------------------------------------------

% In order to properly used the pipeline files in the EEGLAB fromat (.set) are requeired
% Importing the data and eventually extra information for the events can be
% done using the functions provided here or any other method. 

% -------------------------------------------------------------------------
% Import data, filter
% -------------------------------------------------------------------------
Path2DataRaw = fullfile(Path2Data, 'raw');
FilesnameIn = 'example_dataset_01.raw';
EEG = eega_importdata(fullfile(Path2DataRaw, FilesnameIn) );
EEG = eega_importinfoevents(EEG, Path2DataEvent, 'STRT');
EEG = eega_addchinfo(EEG, filechanloc128, 'filetype', 'sfp');
    
% -------------------------------------------------------------------------
% Filter
% -------------------------------------------------------------------------
EEG = eega_demean(EEG);
EEG = pop_eegfiltnew(EEG, [], filt_lowpass,  [], 0, [], [], 0);
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], 0);

% -------------------------------------------------------------------------
% Detect artifacts
% -------------------------------------------------------------------------
% detect bad channels
EEG = eega_tArtifacts(EEG, ArtBadEl, 'KeepRejPre', 1);
% detect motion artifacts
EEG = eega_tArtifacts(EEG, ArtMot1, 'KeepRejPre', 1);
% detect jumps in the signal
EEG = eega_tArtifacts(EEG, ArtJump, 'KeepRejPre', 1);
% define Bad Times and Bad Channels
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);
    
% -------------------------------------------------------------------------
% Correct artifacts 
% -------------------------------------------------------------------------
% correct brief jumps in the signal using target PCA
EEG = eega_tTargetPCAxElEEG(EEG, Int.PCA.nSV, Int.PCA.vSV, 'maxTime', Int.PCA.maxTime,'maskTime', Int.PCA.maskTime,'splicemethod', Int.PCA.splicemethod);
EEG = eega_demean(EEG );
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], 0);
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);
% spatially interpolate channels not working during some time
EEG = eega_tInterpSpatialSegmentEEG(EEG, Int.Spl.p,'pneigh',Int.Spl.pneigh,'splicemethod',Int.Spl.splicemethod,'mingoodtime',Int.Spl.minGoodTime,'minintertime',Int.Spl.minInterTime,'masktime',Int.Spl.maskTime);
EEG = eega_demean(EEG );
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], 0);
% perform ICA
if ICA.do
    EEG = eega_pcawtica(EEG, 'icaname','wica_','filthighpass',ICA.flthpass,'filtlowpass',[],'npc',ICA.npc,'funrejica','iMARA');
end
% spatially interpolate channels not working during the whole recording
EEG = eega_tInterpSpatialEEG(EEG, Int.Spl.p,'pneigh',Int.Spl.pneigh);

% -------------------------------------------------------------------------
% Detect artifacts again
% -------------------------------------------------------------------------
EEG = eega_tArtifacts(EEG, ArtMot2, 'KeepRejPre', 0);
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);


%% ========================================================================
%% EPOCH AND OBTAIN THE ERP

%% ------------------------------------------------------------------------
%% Parameters for the ERP
%% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parameters for filtering
% -------------------------------------------------------------------------
filt_highpass = 0.2;
filt_lowpass = 20;

% -------------------------------------------------------------------------
% Parameters for epoching
% -------------------------------------------------------------------------
Epoch.tw = [-1.600, 2.200];
Epoch.ev = {'Iout'};

% -------------------------------------------------------------------------
% Generate the structures with the parameters for the artifact correction
% -------------------------------------------------------------------------
Int = example_APICE_Interpolation;

% -------------------------------------------------------------------------
% Parameters to define Bad Times (BT) and Bad Channels (BC) 
% -------------------------------------------------------------------------
BCall.nbt           = [0.70 0.50 0.30 0.30];
BCep.nbt            = [0.70 0.50 0.30 0.10/diff(Epoch.tw)];   
BTep.nbc            = [0.70 0.50 0.30 0.30]; 
BTep.minGoodTime    = 1.00;             % shorter intervals between bad segments will be marked as bad
BTep.minBadTime     = 0.10;             % too short periods will not be considered as bad
BTep.maskTime       = 0;                % also mark as bad surronding samples

% -------------------------------------------------------------------------
% Parameters to define Bad Epochs (BE) based on the amount of bad data
% -------------------------------------------------------------------------
DefBEa.maxloops     = 1;
DefBEa.limBCTa      = 1.00;   % maximun bbad data per epoch
DefBEa.limBTa       = 0.00;   % maximun bad times per epoch
DefBEa.limBCa       = 0.30;   % maximun channles with some bad data
DefBEa.limCCTa      = 0.50;   % maximun interpolated data

% -------------------------------------------------------------------------
% Parameters to perform DSS
% -------------------------------------------------------------------------
DSS.do      = 0;    % apply or not the DSS
DSS.k   	= 50;   % components to keep in the first PCA
DSS.n   	= 15;   % components to keep in the second PCA
DSS.fbias   = [];
DSS.fapply  = [];
% DSS.fbias   = {{'eventCnd' 'eventItm'} {'STDx4'};{'eventCnd' 'eventItm'} {'STDx5'};{'eventCnd' 'eventItm'} {'DEVx4'};{'eventCnd' 'eventItm'} {'DEVx5'}};
% DSS.fapply  = DSS.fbias;

% -------------------------------------------------------------------------
% Parameters for baseline correction
% -------------------------------------------------------------------------
BL.tw          = [-100 100]; 

%% ------------------------------------------------------------------------
%% Run the ERP
%% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Arreange events
% -------------------------------------------------------------------------
% correct the latency of the events by the DINs
% The latency of the third input (e.g., 'Icue') is corrected by the second input (e.g., 'DIN6')
EEG = eega_latencyevent(EEG, 'DIN6', {'Icue'});
EEG = eega_latencyevent(EEG, 'DIN6', {'Iout'});
EEG = eega_latencyevent(EEG, 'DIN6', {'Ieye'});
% remove unuseful events
EEG = eega_removeevent(EEG, 'DIN6',  [], 'all');

% -------------------------------------------------------------------------
% Filter
% -------------------------------------------------------------------------
EEG = pop_eegfiltnew(EEG, [], filt_lowpass,  [], 0, [], [], 0);
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], 0);

% -------------------------------------------------------------------------
% Epoch, define conditions, and BT and BC on the epoched data
% -------------------------------------------------------------------------

% Epoch
% time for the epoching (Epoch.tw) and referece events to use as time zero (Epoch.ev) should be modified
EEG = eega_epoch(EEG, Epoch.ev, Epoch.tw);

% define based on the events the factors that will be used to determine conditions
% the second input indicates events properties (e.g., in EEG.epoch(i).eventCong) to use for defining factors of analysis 
% this factors can be the used for averaging
EEG = eega_definefactors(EEG, {'eventCong' 'eventProb'});

% define Bad Times and Bad Channels
EEG = eega_tDefBTBC(EEG, BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime);

% -------------------------------------------------------------------------
% Interpolate the bad channels for each epoch 
% -------------------------------------------------------------------------
EEG = eega_tInterpSpatialEEG(EEG, Int.Spl.p,'pneigh', Int.Spl.pneigh);
EEG = eega_tDefBTBC(EEG, BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime);

% -------------------------------------------------------------------------
% Define bad epochs and remove them
% -------------------------------------------------------------------------
% Define bad epochs
EEG = eega_tDefBEbaddata(EEG, DefBEa,'keeppre',0,'plot', 0);
% Remove bad epochs
EEG = eega_rmvbadepochs(EEG);

% -------------------------------------------------------------------------
% Further process
% -------------------------------------------------------------------------
% reference average
EEG = pop_reref(EEG,[]);
% DSS
if DSS.do
    EEG = dss_denoise_EEG(EEG, DSS.k, DSS.n, DSS.fbias, DSS.fapply);
end
% baseline correction
EEG = eega_rmbaseline(EEG, BL.tw);
% average per conditions
EEG = eega_avgdatabyfactors(EEG, {'eventCong' 'eventProb'}, 'dim2avg', 3);


