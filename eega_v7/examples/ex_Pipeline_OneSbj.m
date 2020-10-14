% This is an example of a pre-processing pipeline for one subject 
% It can be extended to multiple subjetc loopinng over subjects files or
% using the the script for multiple subjects 
%
% Once the data id in the EEGLAB structure 
% --> save data use: pop_saveset(EEG, FILE.set, FOLDER);
% --> to load data use: EEG = pop_loadset('FILE.set', FOLDER);

clear, close all
restoredefaultpath
clc

%% ------------------------------------------------------------------------
%% PATHS
% Folder where EEGLAB is
Path2Mat1 = 'C:\Users\an251296\Documents\MATLAB\toolbox\eeglab14_0_0b';
% Folder where the function eega_ are
Path2Mat2 = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eega_v6_used';
% Folder where the data is and will be saved
Path2Data = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eega_v6_used\examples\DATA';
% Folder where the event files are
Path2DataEvent = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eega_v6_used\examples\DATA\evt';
% Current path
Path0 = 'C:\Users\an251296\nextCloud\MATLAB\mytoolbox\eega_v6_used\examples';
% Channels location file
filechanloc128 = fullfile(Path0, 'ElectrodesLayout','GSN-HydroCel-128.sfp');

%% ------------------------------------------------------------------------
%% ADD TOOLBOXES
% add EEGLAB
cd(Path2Mat1)
eeglab
close all
% Add the functions
addpath(genpath(Path2Mat2))
% Add the path to the data
addpath(genpath(fullfile(Path2Data)))
% Got to the original path
cd(Path0)
% Add the path to the folder with the parameters
addpath(fullfile(Path0, 'parameters'))

%% ------------------------------------------------------------------------
%% IMPORT DATA AND EVENTS

Path2DataRaw = fullfile(Path2Data, 'raw');
Path2DataMFF = fullfile(Path2Data, 'mff');

% OPTION 1: import .raw files and .evt event files
EEG = eega_importdata( fullfile(Path2DataRaw,'Session 20170602 BB30.raw') );
EEG = eega_importinfoevents(EEG, Path2DataEvent);

% % OPTION 2: import .raw files and .evt event files
% EEG = eega_importdata( fullfile(Path2Data,'Session 20170602 BB30.mff') );
% EEG = eega_addeventsfromevt(EEG, Path2DataEvent, 4);
% EEG = eega_eventextra2fields(EEG, [],'%d');
  
% Add the channel info
EEG = eega_addchinfo(EEG, filechanloc128, 'filetype', 'sfp');

%% ------------------------------------------------------------------------
%% DEAL WITH ISSUES WITH THE DINs and EVENTS 
% Remove DINs that are too close
EEG = eega_rmvcloseeventsB(EEG, 'DIN1', 0.500, 0);
% Show all the events
eega_displayevent(EEG);
% Show the events specified 
eega_displayevent(EEG, {'rs' 'Fr' 'Ftl' 'Fts' 'Tw'});


%% ------------------------------------------------------------------------
%% ARTIFACTS DETECTION AND CORRECTION ON CONTINUOS DATA

ex_parameters_preprocessing     % the parameters 
filt_highpass = 0.2;
filt_lowpass = 40;
minphase = 0;

% %down sample if necessary
% EEG = pop_resample(EEG, 250);
%demean, low pass-filter, and high pass- filter
EEG = eega_demean(EEG);
EEG = pop_eegfiltnew(EEG, [], filt_lowpass,  [], 0, [], [], minphase);
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], minphase);
%detect artifacts
EEG = eega_tArtifacts(EEG, ArtPwrS, 'FilterDo', 0, 'KeepRejPre', 0);
EEG = eega_tArtifacts(EEG, ArtMot1, 'FilterDo', 1, 'KeepRejPre', 1);
EEG = eega_tArtifacts(EEG, ArtJump, 'FilterDo', 1, 'KeepRejPre', 1);
%define bad channels and times
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);
%correct short bad segments of data (+ deman and high pass filter)
EEG = eega_tTargetPCAxElEEG(EEG, Int.PCA.nSV, Int.PCA.vSV, 'maxTime', Int.PCA.maxTime,'maskTime', Int.PCA.maskTime,'splicemethod', Int.PCA.splicemethod);
EEG = eega_demean(EEG);
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], minphase);
%define bad channels and times
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);
%spatially interpolate segments of bad data (+ deman and high pass filter)
EEG = eega_tInterpSpatialSegmentEEG(EEG, Int.Spl.p,'pneigh',Int.Spl.pneigh,'splicemethod',Int.Spl.splicemethod,'mingoodtime',Int.Spl.minGoodTime,'minintertime',Int.Spl.minInterTime,'masktime',Int.Spl.maskTime);
EEG = eega_demean(EEG);
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], minphase);
EEGforica = EEG;  % --> this data may be useful for ICA
%spatially interpolate channels not working during the whole recording
EEG = eega_tInterpSpatialEEG(EEG, Int.Spl.p,'pneigh',Int.Spl.pneigh);        
%detect remaning artifacts
EEG = eega_tArtifacts(EEG, ArtMot2, 'FilterDo', 1, 'KeepRejPre', 1);
%define bad channels and times
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);

% % ICA
% % ICA has to be run before interpolating channels not working during the
% % whole experimetn and before averea referencing 
% rmvart = 1;
% npc = 30;  % number of component to keep in the PCA (how many can be kept it depends of the lenght of the data)
% nameica = 'ica_';
% EEGica = eega_pcawtica(EEGforica, 'icaname',nameica,'filthighpass', 1, 'filtlowpass', 40, 'npc', npc, 'level', 6, 'threshtype', 'xlevel', 'applymara',1, 'usealphaband', 0,'rmvart',rmvart);
% %spatially interpolate channels not working during the whole recording
% EEGica = eega_tInterpSpatialEEG(EEGica, Int.Spl.p,'pneigh',Int.Spl.pneigh);        
% %remove components
% EEGicarm = eega_rmvic(EEGica);         
           
%% ------------------------------------------------------------------------
%% EPOCH AND OBTAIN THE ERP

ex_parameters_erp   % the parameters 
filt_highpass_ERP = 1;
filt_lowpass_ERP  = 15;
minphase = 0;

%correct the latency of the events by the DINs and remove the DINs
EEG = eega_latencyevent(EEG, 'DIN1', {'Tw'}, [1 1]);
%remove unuseful events
EEG = eega_removeevent(EEG, 'DIN1',  [], 'all');
%high pass filter 
EEG = pop_eegfiltnew(EEG, filt_highpass_ERP, [], [], 0, [], [], minphase);
%epoch
EEG = eega_epoch(EEG, {'Tw'}, OPT.epoch.tw);
%define bad channels and times
EEG = eega_tDefBTBC(EEG, BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',1,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime );
%interpolate the bad channels for each epoch
EEG = eega_tInterpSpatialEEG(EEG, Int.Spl.p,'pneigh', Int.Spl.pneigh);
%check there are not big artifacts in the epoch data
EEG = eega_tArtifacts(EEG, ArtMot2, 'FilterDo', 1, 'KeepRejPre', 1);
%define bad channels and times
EEG = eega_tDefBTBC(EEG, BTep.nbc,BCep.nbt,BCall.nbt,'keeppre',1,'minBadTime',BTep.minBadTime,'minGoodTime',BTep.minGoodTime,'maskTime',BTep.maskTime );
%define bad epochs based on the amount of bad data, corrected data
EEG = eega_tDefBEbaddata(EEG, DefBEa,'where',DefBEa.where,'keeppre',0,'plot', 0);
%define bad epochs based on the distance of the activity for each epoch across time and electrodes
EEG = eega_tDefBEdistT(EEG, DefBEdT.limMean, DefBEdT.limMax,'where',DefBEdT.where, 'maxloops', DefBEdT.maxloops, 'keeppre',1, 'plot', 0, 'savefigure', 0);
EEG = eega_tDefBEdistE(EEG, DefBEdE.limMean, DefBEdE.limMax,'where',DefBEdE.where, 'maxloops', DefBEdE.maxloops, 'keeppre',1, 'plot', 0, 'savefigure', 0);
%low pass filter 
EEG = pop_eegfiltnew(EEG, [], filt_lowpass_ERP, [], 0, [], [], minphase);
%reference avereage
EEG = eega_refavg(EEG);
%normalize
EEG = eega_normalization(EEG, 'epochs',OPT.norm.epochs,'electrodes',OPT.norm.electrodes,'latency',OPT.norm.latency,'mean',OPT.norm.mean,'sd',OPT.norm.sd,'BadData',OPT.norm.baddata);
%remove a baseline
EEG = eega_rmbaseline(EEG, OPT.bl.tw);
%define based on the events the factors that will be used to determine conditions
EEG = eega_definefactors(EEG, {'eventCnd'});
%remove bad epochs
EEG = eega_rmvbadepochs(EEG);
%average
EEG = eega_avgdatabyfactors(EEG, {'eventCnd'}, 'dim2avg', 3);

