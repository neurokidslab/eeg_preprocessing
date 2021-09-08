% .........................................................................
% This code exemplifies how to run APICE transparently and with a detailed 
% explanation of each step.
% The parameters a the begin allow avoiding steps that are omitted in the 
% different versions of the pipeline described in the manuscript 
% (Fló et al. 2021) 
% The code first exemplified how to obtained continuous preprocess data and
% then how to obtain ERPs.
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
Path2EEGLAB = fullfile('C:\Users\an251296\Documents\MATLAB\toolbox','eeglab2020_0');
% Folder where the function for APICE are:
Path2APICE = fullfile('C:\Users\an251296\nextCloud\MATLAB\mytoolbox','eeg_preprocessing');
% Folder where iMARA is:
Path2iMARA = fullfile('C:\Users\an251296\Documents\MATLAB\toolbox','iMARA-main');
% Current path:
Path0 = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eeg_preprocessing\examples';
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
% Add iMARA (https://github.com/Ira-marriott/iMARA/tree/main)
addpath(genpath(Path2iMARA))
% Add the path to the folder with the parameters
addpath(Path2Parameters)
% Got to the original path
cd(Path0)


%% ========================================================================
%% PARAMETERS
% Define the parameters

%% ------------------------------------------------------------------------
%% Parameters for plotting

% Plot rejection matrixes at different stages
% (1) plot | (0) do not plot
Pplot.rejection = 0;

%% ------------------------------------------------------------------------
%% Parameters for the preprocessing (continuos data)

% -------------------------------------------------------------------------
% Parameters for importing the data
% -------------------------------------------------------------------------

% Name of the files to read from the Path2DataSet folder
Ppp.Files2Read = '*.set';

% Add extra information in an even file (.evt) located in Path2DataEventimport 
% This step should be avoided unless extra information in the event 
% file is not imported from the raw data to the EEGLAB structure. 
% The function adding events has been created to add events sent by 
% Psychtoolbox to EGI. The addition of new events will not work with other
% acquisition systems. Customized functions might be required
% (1) add | (0) do not add
Ppp.importevntapply = 1;  

% Event to use to align differences in time offsets between the event file 
% and the EEG structure
% Only necessary if events are added (Ppp.importevntapply = 1)
Ppp.importevnt0 = 'STRT';   
	

% -------------------------------------------------------------------------
% Parameters for correcting the events
% -------------------------------------------------------------------------

% Correct the latency of a given event using Digital Input events (DINs)
% (1) apply | (0) do not apply
Ppp.eventcorrlatapply = 1;  

% The latency of the events 'Icue' 'Iout' 'Ieye' are corrected by the latency of 'DIN6'
% Only necessary if events are added (Ppp.eventcorrlatapply = 1)
Ppp.eventcorrlat = [];
Ppp.eventcorrlat(1).event1 = {'Icue'};  % event which latencies will be corrected
Ppp.eventcorrlat(1).event2 = 'DIN6';    % event to use to correct the latency
Ppp.eventcorrlat(2).event1 = {'Iout'};  % event which latencies will be corrected
Ppp.eventcorrlat(2).event2 = 'DIN6';    % event to use to correct the latency
Ppp.eventcorrlat(3).event1 = {'Ieye'};  % event whdich latencies will be corrected
Ppp.eventcorrlat(3).event2 = 'DIN6';    % event to use to correct the latency

% Delete unuseful events
% (1) apply | (0) do not apply
Ppp.eventrmvapply = 1;
% Events to delete
% Only necessary if events are added (Ppp.eventrmvapply = 1)
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

% Apply the interpolation of transient artifacts using target PCA
% (1) apply (APICE) | (0) do not apply (APICEa, APICE+W-ICA)
Ppp.IntTransientArtPCA = 1;

% Apply the interpolation of transient artifacts using spherical spline
% (1) apply (APICE, APICE+W-ICA) | (0) do not apply (APICEa)
Ppp.IntTransientArtSpline = 1;

% Generate the structures with the parameters for the artifact correction
% The functions to generate the structures should be in the path (in Path2Parameters)
% This function should be modified to change the algorithms applied for
% artifacts interpolation
Ppp.Int = example_APICE_Interpolation;


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

% -------------------------------------------------------------------------
% Parameters for W-ICA
% -------------------------------------------------------------------------

% Run or not ICA
% (1) apply | (0) do not apply
Ppp.ICA.apply           = 0;    
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


%% ------------------------------------------------------------------------
%% Parameters for ERP analysis
%% ------------------------------------------------------------------------

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
%% RUN THE PREPROCESSING FOR ALL THE SUBJECTS

files2pp = dir(fullfile(Path2DataSet,Ppp.Files2Read));
files2pp = files2pp(~ismember({files2pp.name},{'.', '..'}));
for sbj=1:length(files2pp)
    
    FilesnameIn = files2pp(sbj).name;
    
    % ---------------------------------------------------------------------
    % Load the data
    EEG = pop_loadset(fullfile(Path2DataSet, FilesnameIn) );
    
    % ---------------------------------------------------------------------
    % Import the channels location file
    EEG = pop_chanedit(EEG, 'load', {filechanloc 'filetype' 'autodetect'});
    EEG = eeg_checkset( EEG );
    
    % ---------------------------------------------------------------------
    % Add extra information from the evnet file
    if Ppp.importevntapply
        EEG = eega_importinfoevents(EEG, Path2DataEvent, Ppp.importevnt0);
    end
    
    % ---------------------------------------------------------------------
    % Arreange events
    
    % Correct the latency of the events by the DINs
    if Ppp.eventcorrlatapply
        for i=1:length(Ppp.eventcorrlat)
            EEG = eega_latencyevent(EEG, Ppp.eventcorrlat(i).event2, Ppp.eventcorrlat(i).event1);
        end
    end
    
    % Remove unuseful events
    if Ppp.eventrmvapply
        for i=1:length(Ppp.eventrmv)
            EEG = eega_removeevent(EEG, Ppp.eventrmv(i),  [], 'all');
        end
    end
    
    % ---------------------------------------------------------------------
    % Filter
    
    EEG = eega_demean(EEG);
    EEG = pop_eegfiltnew(EEG, [], Ppp.filt_lowpass,  [], 0, [], [], 0);
    EEG = pop_eegfiltnew(EEG, Ppp.filt_highpass, [], [], 0, [], [], 0);
    
    % ---------------------------------------------------------------------
    % Detect artifacts
    
    % detect bad channels
    EEG = eega_tArtifacts(EEG, Ppp.ArtBadEl, 'KeepRejPre', 1);
    % detect motion artifacts
    EEG = eega_tArtifacts(EEG, Ppp.ArtMot1, 'KeepRejPre', 1);
    % detect jumps in the signal
    EEG = eega_tArtifacts(EEG, Ppp.ArtJump, 'KeepRejPre', 1);
    % define Bad Times and Bad Channels
    EEG = eega_tDefBTBC(EEG, Ppp.BTall.nbc, Ppp.BCall.nbt, Ppp.BCall.nbt, 'keeppre', 0, 'minBadTime', Ppp.BTall.minBadTime, 'minGoodTime', Ppp.BTall.minGoodTime, 'maskTime', Ppp.BTall.maskTime);
    
    % Plot the rejection matrix
    if Pplot.rejection
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts rejection 1')
        
        eega_plot_rejection(EEG, 1, 1, 1, 120)
        set(gcf, 'Name', 'artifacts rejection 1')
    end
    
    % ---------------------------------------------------------------------
    % Correct transient artifacts artifacts
    
    % correct brief jumps in the signal using target PCA
    if Ppp.IntTransientArtPCA
        EEG = eega_tTargetPCAxElEEG(EEG, Ppp.Int.PCA.nSV, Ppp.Int.PCA.vSV, 'maxTime', Ppp.Int.PCA.maxTime,'maskTime', Ppp.Int.PCA.maskTime,'splicemethod', Ppp.Int.PCA.splicemethod);
        EEG = eega_demean(EEG);
        EEG = pop_eegfiltnew(EEG, Ppp.filt_highpass, [], [], 0, [], [], 0);
        EEG = eega_tDefBTBC(EEG, Ppp.BTall.nbc, Ppp.BCall.nbt, Ppp.BCall.nbt, 'keeppre', 0, 'minBadTime', Ppp.BTall.minBadTime, 'minGoodTime', Ppp.BTall.minGoodTime, 'maskTime', Ppp.BTall.maskTime);
    end
    
    % spatially interpolate channels not working during some time
    if Ppp.IntTransientArtSpline
        EEG = eega_tInterpSpatialSegmentEEG(EEG, Ppp.Int.Spl.p, 'pneigh', Ppp.Int.Spl.pneigh, 'splicemethod', Ppp.Int.Spl.splicemethod, 'mingoodtime', Ppp.Int.Spl.minGoodTime, 'minintertime', Ppp.Int.Spl.minInterTime, 'masktime', Ppp.Int.Spl.maskTime);
        EEG = eega_demean(EEG);
        EEG = pop_eegfiltnew(EEG, Ppp.filt_highpass, [], [], 0, [], [], 0);
    end
    
    % ---------------------------------------------------------------------
    % Perform ICA
    
    if Ppp.ICA.apply
        EEG = eega_pcawtica(EEG,'filthighpass',Ppp.ICA.filthighpass, 'filtlowpass', Ppp.ICA.filtlowpass,...
            'npc', Ppp.ICA.npc, 'classifyIC', Ppp.ICA.classifyIC, 'classifyICfun', Ppp.ICA.classifyICfun,...
            'changelabelch',Ppp.ICA.changelabelch,'labelch',Ppp.ICA.labelch,...
            'saveica',Ppp.ICA.saveica,'icaname',Ppp.ICA.icaname,'icapath',Ppp.ICA.icapath);
    end
    
    % ---------------------------------------------------------------------
    % Spatially interpolate channels not working during the whole recording
    EEG = eega_tInterpSpatialEEG(EEG, Ppp.Int.Spl.p, 'pneigh', Ppp.Int.Spl.pneigh);
    
    % Plot the rejection matrix
    if Pplot.rejection
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts correction')
        
        
        eega_plot_rejection(EEG, 1, 1, 1, 120)
        set(gcf, 'Name', 'artifacts correction')
    end
    
    % ---------------------------------------------------------------------
    % Detect artifacts again
    EEG = eega_tArtifacts(EEG, Ppp.ArtMot2, 'KeepRejPre', 0);
    EEG = eega_tDefBTBC(EEG, Ppp.BTall.nbc, Ppp.BCall.nbt, Ppp.BCall.nbt, 'keeppre', 0, 'minBadTime', Ppp.BTall.minBadTime, 'minGoodTime', Ppp.BTall.minGoodTime, 'maskTime', Ppp.BTall.maskTime);
    
    % Plot the rejection matrix
    if Pplot.rejection
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts rejection 2')
        
        eega_plot_rejection(EEG, 1, 1, 1, 120)
        set(gcf, 'Name', 'artifacts rejection 2')
    end
    
    % ---------------------------------------------------------------------
    % Save the continuos pre-process data
    [~,name,~] = fileparts(EEG.filename);
    newfilename = ['prp_' name];
    pop_saveset(EEG, newfilename,Path2DataPrp);
    
end

% -------------------------------------------------------------------------
% Print a preprocessing report for all the subjetcs
eega_printsummary('prp_*.set', Path2DataPrp, Ppp.report.patternname, Ppp.report.patternn, 0);

%% ========================================================================
%% OBTAIN ERP FOR ALL THE SUBJECTS

files2pp = dir(fullfile(Path2DataPrp,'prp_*.set'));
files2pp = files2pp(~ismember({files2pp.name},{'.', '..'}));
for sbj=1:length(files2pp)
    
    FilesnameIn = files2pp(sbj).name;
    
    % ---------------------------------------------------------------------
    % Load the data
    EEG = pop_loadset(fullfile(Path2DataPrp, FilesnameIn) );
    
    % ---------------------------------------------------------------------
    % Filter
    EEG = pop_eegfiltnew(EEG, [], Perp.filt_lowpass,  [], 0, [], [], 0);
    EEG = pop_eegfiltnew(EEG, Perp.filt_highpass, [], [], 0, [], [], 0);
    
    % ---------------------------------------------------------------------
    % Epoch
    EEG = eega_epoch(EEG, Perp.Epoch.ev, Perp.Epoch.tw);
    
    % ---------------------------------------------------------------------
    % Define the factors that will enable to define the conditions
    EEG = eega_definefactors(EEG, Perp.factors);
    
    % ---------------------------------------------------------------------
    % Define Bad Times and Bad Channels on the epoched data
    EEG = eega_tDefBTBC(EEG, Perp.BTep.nbc, Perp.BCep.nbt, Perp.BCall.nbt, 'keeppre', 0, 'minBadTime', Perp.BTep.minBadTime, 'minGoodTime', Perp.BTep.minGoodTime, 'maskTime', Perp.BTep.maskTime);
    
    % ---------------------------------------------------------------------
    % Interpolate the bad channels for each epoch
    EEG = eega_tInterpSpatialEEG(EEG, Perp.Int.Spl.p,'pneigh', Perp.Int.Spl.pneigh);
    EEG = eega_tDefBTBC(EEG, Perp.BTep.nbc, Perp.BCep.nbt, Perp.BCall.nbt, 'keeppre', 0, 'minBadTime', Perp.BTep.minBadTime, 'minGoodTime', Perp.BTep.minGoodTime, 'maskTime', Perp.BTep.maskTime);
    
    % ---------------------------------------------------------------------
    % Define bad epochs
    EEG = eega_tDefBEbaddata(EEG, Perp.DefBEa, 'keeppre', 0, 'plot', 0);
    
    % Plot the rejection matrix
    if Pplot.rejection
        eega_plot_artifacts(EEG)
        set(gcf, 'Name', 'artifacts rejection 2')
    end
    
    % Remove the bad epochs
    EEG = eega_rmvbadepochs(EEG);
    
    % ---------------------------------------------------------------------
    % Reference average
    EEG = pop_reref(EEG,[]);
    
    % ---------------------------------------------------------------------
    % DSS
    if Perp.DSS.apply
        EEG = dss_denoise_EEG(EEG, Perp.DSS.k, Perp.DSS.n, Perp.DSS.fbias, Perp.DSS.fapply);
    end
    
    % ---------------------------------------------------------------------
    % Baseline correction
    EEG = eega_rmbaseline(EEG, Perp.BL.tw);
    
    % ---------------------------------------------------------------------
    % Average per conditions
    EEG = eega_avgdatabyfactors(EEG, Perp.avgcnd, 'dim2avg', 3);
    
    % ---------------------------------------------------------------------
    % Save the ERPs
    [~,name,~] = fileparts(EEG.filename);
    newfilename = ['erp_' name];
    pop_saveset(EEG, newfilename,Path2DataERP);
end

% -------------------------------------------------------------------------
% Print a preprocessing report for all the subjetcs
eega_printsummary('erp_*.set', Path2DataERP, Perp.report.patternname, Perp.report.patternn, 0);
