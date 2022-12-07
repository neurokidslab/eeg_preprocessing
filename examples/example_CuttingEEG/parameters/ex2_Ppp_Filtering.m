%% -----------------------------------------------------------------------
% Filtering
%
% This script sets the parameters for filtering during the pre-processing
% of continuos data
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Filtering

% High pass filter
% ------------------------------------------------------------------------
P.highpass   = 0.1;

% Low pass filter
% ------------------------------------------------------------------------
P.lowpass    = 98;

% Notch
% ------------------------------------------------------------------------
P.notch      = [48 52];

end
