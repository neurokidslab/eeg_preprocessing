
%% -----------------------------------------------------------------------
% ARTIFACTS DETECTION: BAD ELECTRODES
%
% This script sets the input structure to provide to eega_tArtifacts to run
% algorithms that detect non-functional channels 
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Art_BadEl
i=1;

% NON average reference data, RELATIVE threshold
% -------------------------------------------------------------------------

% REJECT: electrodes by epochs based on the correlation between channels
P(i).algorithm      = 'eega_tRejCorrCh';
P(i).loops          = 1;
P(i).P.refdata      = 0;
P(i).P.dozscore     = 0;
P(i).P.twdur        = 4;
P(i).P.twstep       = 1;
P(i).P.topcorrch    = 5;  
P(i).P.thresh       = 0.4; % threshold for the correlation
P(i).P.relative     = 0; 
i = i+1;

% REJECT: electrodes by epochs based on the  power spectrum 
P(i).algorithm      = 'eega_tRejPwr';
P(i).loops          = 2;
P(i).P.refdata      = 0;
P(i).P.dozscore     = 1;
P(i).P.twdur        = 4;
P(i).P.twstep       = 2;
P(i).P.frqband      = [1 10; 20 40];    %(Hz) frequency bands
P(i).P.thresh       = [-3 Inf; -Inf 3]; % threshold for the power 
P(i).P.relative     = [1; 1]; 
i = i+1;

% Include / reject data based on rejected data
% ------------------------------------------------------------------------

% REJECT: too short included segments
P(i).algorithm      = 'eega_tRejShortGood';
P(i).loops          = 2;
P(i).P.timelim      = 0.5;
i = i+1;

end