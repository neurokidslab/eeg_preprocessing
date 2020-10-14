%% ARTIFACTS REJECTION PARAMETERS

function Art = ArtPP_PwrS

i=1;

%% NON average reference data, RELATIVE threshold for EACH electrode

% REJECT: electrodes by epochs based on the  power spectrum 
Art(i).algorithm        = 'eega_tRejPwr';
Art(i).loops            = [1 2];
Art(i).P.refdata        = 0;
Art(i).P.dozscore       = 1;
Art(i).P.twdur          = 5;
Art(i).P.twstep         = 2;
Art(i).P.frqband        = [1 10; 20 40];    %(Hz) frequency bands
Art(i).P.thresh         = [-3 Inf; -Inf 3]; % threshold for the power 
Art(i).P.relative       = [1; 1]; 
i = i+1;

end
