%% -----------------------------------------------------------------------
% ARTIFACTS DETECTION: MOTION ARTIFACTS - first time
%
% This script sets the input structure to provide to eega_tArtifacts to run
% algorithms that detect motion artifacts in the signal 
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Art_Mot1
i=1;

% NON average reference data, ABSOLUTE threshold for ALL electrodes
% -------------------------------------------------------------------------

% REJECT: if the amplitud is bigger than a threshold 
P(i).algorithm        = 'eega_tRejAmp';
P(i).loops            = 1;
P(i).P.refdata        = 0;
P(i).P.dozscore       = 0;
P(i).P.thresh         = 500; 
P(i).P.relative       = 0;
P(i).P.xelectrode     = 0;
P(i).P.mask           = 0.500;
i = i+1;

% NON average reference data, RELATIVE threshold for EACH electrode
% -------------------------------------------------------------------------

% REJECT: if the amplitud is bigger than a threshold 
P(i).algorithm        = 'eega_tRejAmp';
P(i).loops            = [2 3];
P(i).P.refdata        = 0;
P(i).P.dozscore       = 0;
P(i).P.thresh         = 2;  
P(i).P.relative       = 1;
P(i).P.xelectrode     = 1;
P(i).P.mask           = 0.00;
i = i+1;

% REJECT: electrode based on time variance 
P(i).algorithm        = 'eega_tRejTimeVar';
P(i).loops            = [2 3];
P(i).P.refdata        = 0;
P(i).P.dozscore       = 0;
P(i).P.twdur          = 0.500;
P(i).P.twstep         = 0.100;
P(i).P.thresh         = 2*[-1 1];
P(i).P.relative       = 1; 
P(i).P.xelectrode     = 1;
i = i+1;

% REJECT: using the weighted running average Net Station algorithm 
P(i).algorithm        = 'eega_tRejRunningAvg';
P(i).loops            = [2 3];
P(i).P.refdata        = 0;
P(i).P.dozscore       = 0;
P(i).P.thresh_fa      = 2;  
P(i).P.thresh_da      = 2;
P(i).P.relative       = 1;
P(i).P.xelectrode     = 1;
P(i).P.mask           = 0.00;
i = i+1;

% Average reference data, RELATIVE threshold for ALL electrodes
% -------------------------------------------------------------------------

% REJECT: if the amplitud is bigger than a threshold 
P(i).algorithm        = 'eega_tRejAmp';
P(i).loops            = [4 5];
P(i).P.refdata        = 1;
P(i).P.refbaddata     = 'replacebynan';
P(i).P.dozscore       = 0;
P(i).P.thresh         = 2;  
P(i).P.relative       = 1;
P(i).P.xelectrode     = 0;
P(i).P.mask           = 0.00;
i = i+1;

% REJECT: electrode based on time variance 
P(i).algorithm        = 'eega_tRejTimeVar';
P(i).loops            = [4 5];
P(i).P.refdata        = 1;
P(i).P.refbaddata     = 'replacebynan';
P(i).P.dozscore       = 0;
P(i).P.twdur          = 0.500;
P(i).P.twstep         = 0.100;
P(i).P.thresh         = 2*[-1 1];
P(i).P.relative       = 1; 
P(i).P.xelectrode     = 0;
i = i+1;

% REJECT: using the weighted running average Net Station algorithm 
P(i).algorithm        = 'eega_tRejRunningAvg';
P(i).loops            = [4 5];
P(i).P.refdata        = 1;
P(i).P.refbaddata     = 'replacebynan';
P(i).P.dozscore       = 0;
P(i).P.thresh_fa      = 2;  
P(i).P.thresh_da      = 2;
P(i).P.relative       = 1;
P(i).P.xelectrode     = 0;
P(i).P.mask           = 0.00;
i = i+1;

% REJECT: using the variance across electrodes
P(i).algorithm         = 'eega_tRejAmpElecVar';
P(i).loops             = 5;
P(i).P.refdata         = 1;
P(i).P.refbaddata      = 'replacebynan';
P(i).P.dozscore        = 0;
P(i).P.thresh          = 3;
P(i).P.mask            = 0.00;
i = i+1;

% Include / reject data based on rejected data
% -------------------------------------------------------------------------

% REJECT: too short rejected segments are marked as good
P(i).algorithm        = 'eega_tIncShortBad';
P(i).loops            = 5;
P(i).P.timelim        = 0.040;
i = i+1;

% REJECT: too short not included segments
P(i).algorithm        = 'eega_tRejShortGood';
P(i).loops            = 5;
P(i).P.timelim        = 0.500;
i = i+1;

end