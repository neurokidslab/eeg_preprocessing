%% -----------------------------------------------------------------------
% ARTIFACTS CORRECTION: target PCA
%
% This script sets the parameters to run spherical spline interpolation of 
% channels not working during the whole recording or epoch using the 
% function eega_tInterpSpatialEEG 
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Perp_Int_Spl

P.p               = 0.50;    % a bit bigger than the limit defining BT to make the interpolation faster
P.pneigh          = 1;

end