function Art = example_APICE_ArtPP_Jump(thresh)

i=1;

%% NON average reference data, RELATIVE threshold for EACH electrode
% REJECT: if a big change happens in a small interval (RELATIVE THRESHOLD)
Art(i).algorithm         = 'eega_tRejFastChange';
Art(i).loops             = [1];
Art(i).P.refdata         = 0;
Art(i).P.dozscore        = 0;
Art(i).P.tmotion         = 0.020;
Art(i).P.thresh          = thresh;  
Art(i).P.relative        = 1;
Art(i).P.xelectrode      = 1;
i = i+1;


%% Average reference data, RELATIVE threshold for ALL electrodes
% REJECT: if a big change happens in a small interval (RELATIVE THRESHOLD)
Art(i).algorithm         = 'eega_tRejFastChange';
Art(i).loops             = [2];
Art(i).P.refdata         = 1;
Art(i).P.refbaddata      = 'replacebynan';
Art(i).P.dozscore        = 0;
Art(i).P.tmotion         = 0.020;
Art(i).P.thresh          = thresh;  
Art(i).P.relative        = 1;
Art(i).P.xelectrode      = 0;
i = i+1;

%% Include / reject data based on rejected data

% REJECT: too short rejected segments are marked as good
Art(i).algorithm         = 'eega_tIncShortBad';
Art(i).loops             = 2;
Art(i).P.timelim         = 0.020;
i = i+1;

% REJECT: too short not included segments
Art(i).algorithm        = 'eega_tRejShortGood';
Art(i).loops            = 2;
Art(i).P.timelim        = 0.100;
i = i+1;

end
