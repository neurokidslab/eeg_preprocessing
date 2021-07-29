% This is an example of a preprocessing pipeline for one dataset

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
ICA.npc         = 50;   % number of components to keep in the PCA
ICA.flthpass    = 2;    % high pass filter before ICA

%% ------------------------------------------------------------------------
%% Run the prerpocessing
%% ------------------------------------------------------------------------

% The data provided as example starts from step 1, a dataset in the EEGLAB 
% format saved as .set
% Importing the data and eventually extra information for the events can be
% done using the functions provided here or any other method

% -------------------------------------------------------------------------
% 0)Import data 
% -------------------------------------------------------------------------
FilesnameIn = '*.raw';
FilesnameOut = '';

eega_RunAll( FilesnameIn,...
    fullfile(Path2Data, 'raw'),...
    fullfile(Path2Data, 'set'),...
    'eega_importdata',              {},...
    'eega_importinfoevents',        {Path2DataEvent, 'STRT'},...
    'eega_addchinfo',               {filechanloc128, 'filetype', 'sfp'},...
    'eega_resample',                {250},...
    'prename',FilesnameOut,'runto', 1);

% -------------------------------------------------------------------------
% 1)Filter
% -------------------------------------------------------------------------
FilesnameIn = '*.set';
FilesnameOut = sprintf('fh%1.0e_fl%d_',filt_highpass,filt_lowpass);

eega_RunAll( FilesnameIn,...
    fullfile(Path2Data, 'set'),...
    fullfile(Path2Data, 'flt'),...
    'eega_demean',                  {},...
    'pop_eegfiltnew',               {[], filt_lowpass,  [], 0, [], [], 0},...
    'pop_eegfiltnew',               {filt_highpass, [], [], 0, [], [], 0},...
    'prename',FilesnameOut,'runto', 1);

% -------------------------------------------------------------------------
% 2) Detect artefacts
% -------------------------------------------------------------------------
FilesnameIn = [sprintf('fh%1.0e_fl%d_',filt_highpass,filt_lowpass) '*.set'];
FilesnameOut = 'a_';

eega_RunAll( FilesnameIn,...
    fullfile(Path2Data, 'flt'),...
    fullfile(Path2Data, 'art'),...
    'eega_tArtifacts',              {ArtBadEl, 'KeepRejPre', 1},...
    'eega_tArtifacts',              {ArtMot1, 'KeepRejPre', 1},...
    'eega_tArtifacts',              {ArtJump, 'KeepRejPre', 1},...
    'eega_tDefBTBC',                {BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime},...
    'prename',FilesnameOut,'runto', 2);

% -------------------------------------------------------------------------
% 3a) Correct localized artifacts, and detect artifacts again 
% APICE
% -------------------------------------------------------------------------
FilesnameIn = ['a_' sprintf('fh%1.0e_fl%d_',filt_highpass,filt_lowpass) '*.set'];
FilesnameOut = 'a_i_is_';

eega_RunAll( FilesnameIn,...
    fullfile(Path2Data, 'art'),...
    fullfile(Path2Data, 'art'),...
    'eega_tTargetPCAxElEEG',        {Int.PCA.nSV, Int.PCA.vSV, 'maxTime', Int.PCA.maxTime,'maskTime', Int.PCA.maskTime,'splicemethod', Int.PCA.splicemethod},...
    'eega_demean',                  {},...
    'pop_eegfiltnew',               {filt_highpass, [], [], 0, [], [], 0},...
    'eega_tDefBTBC',                {BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime},...
    'eega_tInterpSpatialSegmentEEG',{Int.Spl.p,'pneigh',Int.Spl.pneigh,'splicemethod',Int.Spl.splicemethod,'mingoodtime',Int.Spl.minGoodTime,'minintertime',Int.Spl.minInterTime,'masktime',Int.Spl.maskTime},...
    'eega_demean',                  {},...
    'pop_eegfiltnew',               {filt_highpass, [], [], 0, [], [], 0},...
    'eega_tInterpSpatialEEG',       {Int.Spl.p,'pneigh',Int.Spl.pneigh},...
    'eega_tArtifacts',              {ArtMot2, 'KeepRejPre', 0},...
    'eega_tDefBTBC',                {BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime},...
    'prename',FilesnameOut,'runto', 2);


% -------------------------------------------------------------------------
% 3b) Correct localized artifacts, apply ICA, and detect artifacts again
% APICE + ICA
% -------------------------------------------------------------------------
FilesnameIn = ['a_' sprintf('fh%1.0e_fl%d_',filt_highpass,filt_lowpass) '*.set'];
FilesnameOut = 'a_i_wica_is_';

eega_RunAll( FilesnameIn ,...
    fullfile(Path2Data, 'art'),...
    fullfile(Path2Data, 'art'),...
    'eega_tTargetPCAxElEEG',        {Int.PCA.nSV, Int.PCA.vSV, 'maxTime', Int.PCA.maxTime,'maskTime', Int.PCA.maskTime,'splicemethod', Int.PCA.splicemethod},...
    'eega_demean',                  {},...
    'pop_eegfiltnew',               {filt_highpass, [], [], 0, [], [], 0},...
    'eega_tDefBTBC',                {BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime},...
    'eega_tInterpSpatialSegmentEEG',{Int.Spl.p,'pneigh',Int.Spl.pneigh,'splicemethod',Int.Spl.splicemethod,'mingoodtime',Int.Spl.minGoodTime,'minintertime',Int.Spl.minInterTime,'masktime',Int.Spl.maskTime},...
    'eega_demean',                  {},...
    'pop_eegfiltnew',               {filt_highpass, [], [], 0, [], [], 0},...
    'eega_tDefBTBC',                {BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime},...
    'eega_pcawtica',                {'icaname','wica_','filthighpass',ICA.flthpass,'filtlowpass',[],'npc',ICA.npc,'funrejica','iMARA'},...
    'eega_tInterpSpatialEEG',       {Int.Spl.p,'pneigh',Int.Spl.pneigh},...
    'eega_tArtifacts',              {ArtMot2, 'KeepRejPre', 0},...
    'eega_tDefBTBC',                {BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime},...
    'prename',FilesnameOut,'runto', 2);

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
% 4b) Epoch, define experimental factors, filter, interpolate bad channels, 
%       remove bad epochs, average reference, baseline correction, average  
% APICE
% -------------------------------------------------------------------------
FilesnameIn = ['a_i_is_a_' sprintf('fh%1.0e_fl%d_',0.1,40) '*.set'];
FilesnameOut = sprintf('avg_b%dto%d_r_i_e%s%dto%d_fl%d_fh%1.0e_',BL.tw(1),BL.tw(2),Epoch.ev{1},1000*Epoch.tw(1),1000*Epoch.tw(2),filt_lowpass,filt_highpass);

eega_RunAll( FilesnameIn,...
    fullfile(Path2Data, 'art'),...
    fullfile(Path2Data, 'erp'),...
    'eega_latencyevent',        {'DIN6', {'Icue'}},...
    'eega_latencyevent',        {'DIN6', {'Iout'}},...
    'eega_latencyevent',        {'DIN6', {'Ieye'}},...
    'pop_eegfiltnew',           {filt_highpass, [], [], 0, [], [], 0},...
    'pop_eegfiltnew',           {[], filt_lowpass, [], 0, [], [], 0},...
    'eega_epoch',               {Epoch.ev, Epoch.tw},...
    'eega_definefactors',       {{'eventCong' 'eventProb'}},...
    'eega_tDefBTBC',            {BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime},...
    'eega_tInterpSpatialEEG',   {Int.Spl.p,'pneigh', Int.Spl.pneigh},...
    'eega_tDefBTBC',            {BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime},...
    'eega_tDefBEbaddata',       {DefBEa,'keeppre',0,'plot', 0},...
    'eega_rmvbadepochs',        {},...
    'pop_reref',                {[]},...
	'pop_rmbase',               {BL.tw},...
    'eega_avgdatabyfactors',    {{'eventCong' 'eventProb'}, 'dim2avg', 3},...
    'prename', FilesnameOut, 'saveformat', 'set', 'runto', 2);


% -------------------------------------------------------------------------
% 4b) Epoch, define experimental factors, filter, interpolate bad channels, 
%       remove bad epochs, average reference, DSS, baseline correction, average  
% APICE + DSS
% -------------------------------------------------------------------------
FilesnameIn = ['a_i_is_a_' sprintf('fh%1.0e_fl%d_',0.1,40) '*.set'];
FilesnameOut = sprintf('avg_b%dto%d_dss_r_i_e%s%dto%d_fl%d_fh%1.0e_',BL.tw(1),BL.tw(2),Epoch.ev{1},1000*Epoch.tw(1),1000*Epoch.tw(2),filt_lowpass,filt_highpass);

eega_RunAll( FilesnameIn,...
    fullfile(Path2Data, 'art'),...
    fullfile(Path2Data, 'erp'),...
    'eega_latencyevent',        {'DIN6', {'Icue'}},...
    'eega_latencyevent',        {'DIN6', {'Iout'}},...
    'eega_latencyevent',        {'DIN6', {'Ieye'}},...
    'pop_eegfiltnew',           {filt_highpass, [], [], 0, [], [], 0},...
    'pop_eegfiltnew',           {[], filt_lowpass, [], 0, [], [], 0},...
    'eega_epoch',               {Epoch.ev, Epoch.tw},...
    'eega_definefactors',       {{'eventCong' 'eventProb'}},...
    'eega_tDefBTBC',            {BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime},...
    'eega_tInterpSpatialEEG',   {Int.Spl.p,'pneigh', Int.Spl.pneigh},...
    'eega_tDefBTBC',            {BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime},...
    'eega_tDefBEbaddata',       {DefBEa,'keeppre',0,'plot', 0},...
    'eega_rmvbadepochs',        {},...
    'pop_reref',                {[]},...
	'dss_denoise_EEG',          {DSS.k, DSS.n, DSS.fbias, DSS.fapply},...
    'pop_rmbase',               {BL.tw},...
    'eega_avgdatabyfactors',    {{'eventCong' 'eventProb'}, 'dim2avg', 3},...
    'prename', FilesnameOut, 'saveformat', 'set', 'runto', 2);












