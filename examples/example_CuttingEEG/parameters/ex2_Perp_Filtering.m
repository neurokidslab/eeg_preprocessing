%% -----------------------------------------------------------------------
% Filtering
%
% This script sets the parameters for filtering before epoching the data
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Perp_Filtering

% High pass filter
% ------------------------------------------------------------------------
P.highpass   = 0.2;

% Low pass filter
% ------------------------------------------------------------------------
P.lowpass    = 20;

end
