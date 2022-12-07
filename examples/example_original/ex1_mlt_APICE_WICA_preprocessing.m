% .........................................................................
% This code shows how to run the APICE+W-ICA pipeline to obtained continuos preprocess data
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
% Folder where iMARA is:
Path2iMARA = fullfile('xxx','iMARA-main');
% Current path:
Path0 = 'xxx\eeg_preprocessing\examples';
% Folder where the scripts defining the parameters are:
Path2Parameters = fullfile(Path0,'parameters');
% Channels location file
filechanloc = fullfile(Path0,'DATA','ElectrodesLayout','GSN-HydroCel-128.sfp');
% Folder where the event files are:
Path2DataEvent = fullfile(Path0,'DATA','evt');
% Folder where the data imported to EEGLAB is: 
Path2DataSet = fullfile(Path0,'DATA','set');
% Folder where the continuos preprocess data will be saved
Path2DataPrp = fullfile(Path0,'DATA','prp');


%% ========================================================================
%% ADD TOOLBOXES AND PATHS TO USEFUL FOLDERS

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
% Got to the original path
cd(Path0)


%% ========================================================================
%% PARAMETERS
% Define the parameters

% -------------------------------------------------------------------------
% Parameters for importing the data
% -------------------------------------------------------------------------

% Name of the files to read from the Path2DataSet folder
Ppp.Files2Read = '*.set';

% Add extra information in an even file (.evt) located in Path2DataEventimport 
% .........................................................................
% To remove this step, delete the line with eega_importinfoevents!!!
% .........................................................................
% This step should be avoided unless extra information in the event 
% file is not imported from the raw data to the EEGLAB structure. 
% The function adding events has been created to add events sent by 
% Psychtoolbox to EGI. The addition of new events will not work with other
% acquisition systems. Customized functions might be required

% Event to use to align differences in time offsets between the event file 
% and the EEG structure
Ppp.importevnt0 = 'STRT';   
	

% -------------------------------------------------------------------------
% Parameters for correcting the events
% -------------------------------------------------------------------------

% Correct the latency of a given event using Digital Input events (DINs)
% .........................................................................
% To remove this step, delete the line with eega_latencyevent!!!
% .........................................................................
% The latency of the events 'Icue' 'Iout' 'Ieye' are corrected by the latency of 'DIN6'
Ppp.eventcorrlat.event1 = {'Icue' 'Iout' 'Ieye'};  % event which latencies will be corrected
Ppp.eventcorrlat.event2 = 'DIN6';    % event to use to correct the latency

% Delete unuseful events
% .........................................................................
% To remove this step, delete the line with eega_removeevent!!!
% .........................................................................
% Events to delete
Ppp.eventrmv = {'DIN6'};


% -------------------------------------------------------------------------
% Parameters for filtering
% -------------------------------------------------------------------------

% High pass filter
Ppp.filt_highpass   = 0.1;

% Low pass filter
Ppp.filt_lowpass    = 40;


% -------------------------------------------------------------------------
% Parameters for artifacts detection
% -------------------------------------------------------------------------

% Generate the structures with the parameters for the artefact detection
% The functions to generate the structures should be in the path (Path2Parameters)
% These functions should be modified to change the algorithms applied for
% artifacts detection

% In this example, we use a relative threshold of 3 for all types of
% artifacts (APICE(3))

% Algorithms to detect bad electrodes
Ppp.ArtBadEl = example_APICE_ArtPP_BadEl(3);
% Algorithms to detect jumps in the signal
Ppp.ArtJump = example_APICE_ArtPP_Jump(3);
% Algorithms to detect motion artifacts
Ppp.ArtMot1 = example_APICE_ArtPP_Mot1(3);
% Algorithms to detect motion artifacts
Ppp.ArtMot2 = example_APICE_ArtPP_Mot2(3);


% -------------------------------------------------------------------------
% Parameters for transient artifacts interpolation
% -------------------------------------------------------------------------

% Generate the structures with the parameters for the artifact correction
% The functions to generate the structures should be in the path (in Path2Parameters)
% This function should be modified to change the algorithms applied for
% artifacts interpolation
Ppp.Int = example_APICE_Interpolation;

% Create a cell array with all the input to use later
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
% Parameters for W-ICA
% -------------------------------------------------------------------------

% Run or not ICA
% (1) apply | (0) do not apply
Ppp.ICA.apply           = 1;    
% Number of components to keep in the PCA beofre runing ICA. If 0 it is not applied. Default 0
Ppp.ICA.npc             = 50;   
% High pass filter applied before ICA. Default 2
Ppp.ICA.filthighpass    = 2;    
% Low pass filter applied before ICA. Default 40
Ppp.ICA.filtlowpass     = [];    
% Apply an authomatic classification of IC components. Default (1)
Ppp.ICA.classifyIC      = 1;    
% Function use to classify IC. Default 'iMARA'
Ppp.ICA.classifyICfun   = 'iMARA';
% Change the labels of the channels to be consistent with the classification algorithm. Default False.
% Consider that MARA and iMARA use the international 10-10 electrodes position. 
% If channels labels correspond to another system, they need to be adapted.
% The function performing the ICA provides the possibility to change the labels. 
% Use this parameter to apply this step. 
% (1) apply | (0) do not apply
Ppp.ICA.changelabelch   = 1; 
% It can be:
% - a cell of size n x 2 with the labels in chanlocs in the EEGLAB structure and the new names in the 10-10 electrodes position (see example_ch_iMARA.m). 
% - the name (with path) of a text file with the old and new names (see example_ch_iMARA.txt). 
% - empty, then the default is the conversion from a EGI 129 layout to the 10-10 electrodes position
Ppp.ICA.labelch         = [];   
% Save a file with the ICA weights. 
% The file is saved in the folder specified in EEG.filename
% (1) save | (0) do not save. Default True. 
Ppp.ICA.saveica         = 1;    
% Name added to save the ICA decomposition. Default 'ica_'
Ppp.ICA.icaname         = 'wica_'; 
% Folder were the ICA decomposition is seved. If empty, in EEG.filepath
Ppp.ICA.icapath         = Path2DataPrp;

% Create a cell array with all the input to use later
inputsICA = {'filthighpass',Ppp.ICA.filthighpass, 'filtlowpass', Ppp.ICA.filtlowpass,...
            'npc', Ppp.ICA.npc, 'classifyIC', Ppp.ICA.classifyIC, 'classifyICfun', Ppp.ICA.classifyICfun,...
            'changelabelch',Ppp.ICA.changelabelch,'labelch',Ppp.ICA.labelch,...
            'saveica',Ppp.ICA.saveica,'icaname',Ppp.ICA.icaname,'icapath',Ppp.ICA.icapath};
        
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
FilesIn = '*.set';

% String to add at the begign of the output filename
FilesOut = 'prp_';

% Define how to run eega_RunAll
% - If 1, run to all the files in the input folder
% - If 2, run to the files in the input folder do not present in the ouput folder (files still not analyzed)
runto = 1;

% Import all the files in the input folder

eega_RunAll(FilesIn,...
            Path2DataSet,...
            Path2DataPrp,...
            'pop_chanedit',                 {'load', {filechanloc 'filetype' 'autodetect'}},...
            'eeg_checkset',                 {},...
            'eega_importinfoevents',        {Path2DataEvent, Ppp.importevnt0},...
            'eega_latencyevent',            {Ppp.eventcorrlat.event2, Ppp.eventcorrlat.event1},...
            'eega_removeevent',             {Ppp.eventrmv,  [], 'all'},...
            'eega_demean',                  {},...
            'pop_eegfiltnew',               {[], Ppp.filt_lowpass,  [], 0, [], [], 0},...
            'pop_eegfiltnew',               {Ppp.filt_highpass, [], [], 0, [], [], 0},...
            'eega_tArtifacts',              {Ppp.ArtBadEl, 'KeepRejPre', 1},...
            'eega_tArtifacts',              {Ppp.ArtMot1, 'KeepRejPre', 1},...
            'eega_tArtifacts',              {Ppp.ArtJump, 'KeepRejPre', 1},...
            'eega_tDefBTBC',                inputsDefBTBC,...
            'eega_tInterpSpatialSegmentEEG',inputsInterSegmentsSpline,...
            'eega_demean',                  {},...
            'pop_eegfiltnew',               {Ppp.filt_highpass, [], [], 0, [], [], 0},...
            'eega_pcawtica',                inputsICA,...
            'eega_tInterpSpatialEEG',       {Ppp.Int.Spl.p, 'pneigh', Ppp.Int.Spl.pneigh},...
            'eega_tArtifacts',              {Ppp.ArtMot2, 'KeepRejPre', 0},...
            'eega_tDefBTBC',                inputsDefBTBC,...
            'prename',FilesOut,'runto', runto);

% -------------------------------------------------------------------------
% Print a preprocessing report for all the subjetcs
eega_printsummary('prp_*.set', Path2DataPrp, Ppp.report.patternname, Ppp.report.patternn, 0);

