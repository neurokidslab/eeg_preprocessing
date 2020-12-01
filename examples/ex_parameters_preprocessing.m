%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Parameters requiered for the analysis %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Artifacts Detection
ArtPwrS = ArtPP_PwrS;
ArtJump = ArtPP_Jump;
ArtMot1 = ArtPP_Mot1;
ArtMot2 = ArtPP_Mot2;

%% Artifacts Correction
Int = Interpolation;

%% Define Bad Segments and Electrodes 

BCall.nbt           = [0.70 0.50 0.30];

BTall.nbc           = [0.70 0.50 0.30];
BTall.minGoodTime   = 1.000;         % shorter intervals between bad segments will be marked as bad
BTall.minBadTime    = 0.080;         % too short periods will not be considered as bad
BTall.maskTime      = 0.500;         % also mark as bad surronding samples









