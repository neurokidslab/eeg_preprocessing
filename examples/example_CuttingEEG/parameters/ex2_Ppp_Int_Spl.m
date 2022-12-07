%% -----------------------------------------------------------------------
% ARTIFACTS CORRECTION: target PCA
%
% This script sets the parameters to run spherical spline interpolation of 
% bad segments of data or channels
% It can be done using the functions:
% - eega_tInterpSpatialSegmentEEG (interpolates segments on continuos data)
% - eega_tInterpSpatialEEG (interpolates a channel during the whole
% recording or during one epoch)
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Int_Spl

P.p               = 0.50;    % a bit bigger than the limit defining BT to make the interpolation faster
P.pneigh          = 1;
P.splicemethod    = 1; 
P.minGoodTime     = 2.000;   % segments that are too short and between bad segments are marked as bad
P.minInterTime    = 0.100;   % segments shorter than this are not interpolated
P.maskTime        = 1.000;   % mask surronding samples for interpolation 

end
