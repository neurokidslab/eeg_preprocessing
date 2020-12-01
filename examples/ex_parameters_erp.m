%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Parameters requiered for the analysis %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Epoch
OPT.epoch.tw = [-1.750, 3.250];

%% Artifacts Detection
ArtMot2 = ArtPP_Mot2;

%% Artifacts Correction
Int = Interpolation;

%% Define Bad Segments and Electrodes 

BCall.nbt           = [0.70 0.40 0.40];

BCep.nbt            = [0.50 0.25 0]; %[0.15 0.05];

BTep.nbc            = [0.70 0.40 0.30];
BTep.minGoodTime    = 1.000;         % shorter intervals between bad segments will be marked as bad
BTep.minBadTime     = 0.080;         % too short periods will not be considered as bad
BTep.maskTime       = 0.500;         % also mark as bad surronding samples


%% Define bad epoch

% based on the amount of bad data absolute threshold
DefBEa.maxloops   = 1;
DefBEa.where      = [];
DefBEa.limBCTa    = 0.10;  % 10% of bad data maximun
DefBEa.limBTa     = 0.00; 
DefBEa.limBCa     = 0.10;  
DefBEa.limCCTa    = 0.50;
% based on the distance to the mean across samples
DefBEdT.maxloops    = 2;
DefBEdT.where       = [];
DefBEdT.limMean     = 1.50;     % if the time average distance from the mean is larger than 2 SD
DefBEdT.limMax      = 1.50;     % if the time maximun distance from the mean is larger than 2 SD
% based on the distance to the mean acoss electrodes
DefBEdE.maxloops    = 2;
DefBEdE.where       = [];
DefBEdE.limMean     = 1.50;     % if the time average distance from the mean is larger than 2 SD
DefBEdE.limMax      = 1.50;     % if the time maximun distance from the mean is larger than 2 SD

%% Normalization
OPT.norm.electrodes   = 'all';  % 'all' / 'single'
OPT.norm.epochs       = 'single';  % 'all' / 'single' / 'timewindow'
OPT.norm.latency      = 'all'; % 1000*OPT.epoch.tw;  % [] / [opt.epoch.tw(1) 0]
OPT.norm.mean         = 0;
OPT.norm.sd           = [];  % only divide by the standar desviation of the epoch
OPT.norm.baddata      = 'none'; % none / replacebynan

%% Baseline
OPT.bl.tw          = [-500 0]; 
OPT.bl.type        = 'xtrial'; % xtrial / xcond / alltrials
OPT.bl.baddata     = 'none'; % none / replacebynan
    
                
                
