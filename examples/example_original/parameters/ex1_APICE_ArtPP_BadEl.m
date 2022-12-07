%% ARTIFACTS REJECTION PARAMETERS
% Identify bad electrodes

function Art = example_APICE_ArtPP_BadEl(thresh)

i=1;

%% NON average reference data, RELATIVE threshold

% REJECT: electrodes by epochs based on the correlation between channels
Art(i).algorithm        = 'eega_tRejCorrCh';
Art(i).loops            = [1];
Art(i).P.refdata        = 0;
Art(i).P.dozscore       = 0;
Art(i).P.twdur          = 4;
Art(i).P.twstep         = 2;
Art(i).P.topcorrch      = 5;  
Art(i).P.thresh         = 0.4; % threshold for the correlation
Art(i).P.relative       = 0; 
i = i+1;

% REJECT: electrodes by epochs based on the  power spectrum 
Art(i).algorithm        = 'eega_tRejPwr';
Art(i).loops            = [2];
Art(i).P.refdata        = 0;
Art(i).P.dozscore       = 1;
Art(i).P.twdur          = 4;
Art(i).P.twstep         = 2;
Art(i).P.frqband        = [1 10; 20 40];    %(Hz) frequency bands
Art(i).P.thresh         = [-thresh Inf; -Inf thresh]; % threshold for the power 
Art(i).P.relative       = [1; 1]; 
i = i+1;

%% Include / reject data based on rejected data

% REJECT: too short included segments
Art(i).algorithm            = 'eega_tRejShortGood';
Art(i).loops                = 2;
Art(i).P.timelim            = 2.000;
i = i+1;

end
