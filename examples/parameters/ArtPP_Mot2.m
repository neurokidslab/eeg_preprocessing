%% ARTIFACTS REJECTION PARAMETERS
%  First round

function Art = ArtPP_Mot2

i=1;

%% NON average reference data, RELATIVE threshold for EACH electrode

%% Average reference data, RELATIVE threshold for ALL electrodes

% REJECT: if the amplitud is bigger than a threshold (RELATIVE THRESHOLD)
Art(i).algorithm         = 'eega_tRejAmp';
Art(i).loops             = [1 2];
Art(i).P.refdata         = 1;
Art(i).P.refbaddata      = 'none';
Art(i).P.dozscore        = 0;
Art(i).P.thresh          = 2.0;  
Art(i).P.relative        = 1;
Art(i).P.xelectrode      = 0;
Art(i).P.mask            = 0;
i = i+1;

% REJECT: using the weighted running average Net Station algorithm 
Art(i).algorithm         = 'eega_tRejRunningAvg';
Art(i).loops             = [1 2];
Art(i).P.refdata         = 1;
Art(i).P.refbaddata      = 'none';
Art(i).P.dozscore        = 0;
Art(i).P.thresh_fa       = 2;  
Art(i).P.thresh_da       = 2;
Art(i).P.relative        = 1;
Art(i).P.xelectrode      = 0;
Art(i).P.mask            = 0;
i = i+1;

% REJECT: using the weighted running average Net Station algorithm 
Art(i).algorithm         = 'eega_tRejAmpElecVar';
Art(i).loops             = [1 2];
Art(i).P.refdata         = 1;
Art(i).P.refbaddata      = 'none';
Art(i).P.dozscore        = 0;
Art(i).P.thresh          = 3;
Art(i).P.mask            = 0;
i = i+1;

%% Include / reject data based on rejected data

% REJECT: too short rejected segments are marked as good
Art(i).algorithm            = 'eega_tIncShortBad';
Art(i).loops                = 2;
Art(i).P.timelim            = 0.080;
i = i+1;

% REJECT: too short not included segments
Art(i).algorithm            = 'eega_tRejShortGood';
Art(i).loops                = 2;
Art(i).P.timelim            = 0.500; %0.200;
i = i+1;


end