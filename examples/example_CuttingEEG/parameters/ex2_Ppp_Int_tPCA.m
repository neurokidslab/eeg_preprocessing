%% -----------------------------------------------------------------------
% ARTIFACTS CORRECTION: target PCA
%
% This script sets the parameters to run target PCA to correct brief 
% artifacts in the signal (e.i., jumps) using the function eega_tTargetPCAxElEEG  
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Int_tPCA

P.maxTime         = 0.100;
P.nSV             = [];
P.vSV             = 0.90;
P.splicemethod    = 1;   
P.maskTime        = 0.050;

end
