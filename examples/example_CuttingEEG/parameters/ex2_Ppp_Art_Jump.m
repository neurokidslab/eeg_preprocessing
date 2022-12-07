%% -----------------------------------------------------------------------
% ARTIFACTS DETECTION: FAST CHANGES
%
% This script sets the input structure to provide to eega_tArtifacts to run
% algorithms that detect fast jumps in the signal 
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Art_Jump
i=1;

% NON average reference data, RELATIVE threshold for EACH electrode
% -------------------------------------------------------------------------

% REJECT: if a big change happens in a small interval (RELATIVE THRESHOLD)
P(i).algorithm         = 'eega_tRejFastChange';
P(i).loops             = 1;
P(i).P.refdata         = 0;
P(i).P.dozscore        = 0;
P(i).P.tmotion         = 0.020;
P(i).P.thresh          = 2;  
P(i).P.relative        = 1;
P(i).P.xelectrode      = 1;
i = i+1;


% Average reference data, RELATIVE threshold for ALL electrodes
% -------------------------------------------------------------------------

% REJECT: if a big change happens in a small interval (RELATIVE THRESHOLD)
P(i).algorithm         = 'eega_tRejFastChange';
P(i).loops             = 2;
P(i).P.refdata         = 1;
P(i).P.refbaddata      = 'replacebynan';
P(i).P.dozscore        = 0;
P(i).P.tmotion         = 0.020;
P(i).P.thresh          = 2;  
P(i).P.relative        = 1;
P(i).P.xelectrode      = 0;
i = i+1;

% Include / reject data based on rejected data
% -------------------------------------------------------------------------

% REJECT: too short rejected segments are marked as good
P(i).algorithm         = 'eega_tIncShortBad';
P(i).loops             = 2;
P(i).P.timelim         = 0.040;
i = i+1;

end