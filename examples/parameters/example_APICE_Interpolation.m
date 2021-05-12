function Int = example_APICE_Interpolation

%% Artifacts Correction by PCA
Int.PCA.maxTime         = 0.100;
Int.PCA.nSV             = [];
Int.PCA.vSV             = 0.90;
Int.PCA.splicemethod    = 1;   % 0 / 1 / 2 / 3
Int.PCA.maskTime        = 0.050;

%% Artifacts Correction by spline interpolation
Int.Spl.p               = 0.50;    % a bit bigger than the limit defining BT to make the interpolation faster
Int.Spl.pneigh          = 1;
Int.Spl.splicemethod    = 1; %1;   % 0 / 1 / 2 / 3
Int.Spl.minGoodTime     = 2.000;   % segments that are too short and between bad segments are marked as bad
Int.Spl.minInterTime    = 0.100;   % segments shorter than this are not interpolated
Int.Spl.maskTime        = 1.000;   % mask surronding samples for interpolation (mask before was 0.500)

end
