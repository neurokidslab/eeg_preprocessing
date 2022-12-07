%% -----------------------------------------------------------------------
% W-ICA
%
% This script sets the parameters to run ICA to remove physiological
% artifacts using the function eega_pcawtica
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_wICA

% Number of components to keep in the PCA beofre runing ICA. If 0 it is not applied. Default 0
P.npc             = 50;

% High pass filter applied before ICA. Default 2
P.filthighpass    = 2;

% Low pass filter applied before ICA. Default 40
P.filtlowpass     = [];

% Apply an authomatic classification of IC components. Default (1)
P.classifyIC      = 1;

% Function use to classify IC. Default 'iMARA'
P.classifyICfun   = 'iMARA';

% Change the labels of the channels to be consistent with the classification algorithm. Default False.
% Consider that MARA and iMARA use the international 10-10 electrodes position.
% If channels labels correspond to another system, they need to be adapted.
% The function performing the ICA provides the possibility to change the labels.
% Use this parameter to apply this step.
% (1) apply | (0) do not apply
P.changelabelch   = 1;

% It can be:
% - a cell of size n x 2 with the labels in chanlocs in the EEGLAB structure and the new names in the 10-10 electrodes position (see example_ch_iMARA.m).
% - the name (with path) of a text file with the old and new names (see example_ch_iMARA.txt).
% - empty, then the default is the conversion from a EGI 129 layout to the 10-10 electrodes position
P.labelch         = [];

% Save a file with the ICA weights.
% The file is saved in the folder specified in EEG.filename
% (1) save | (0) do not save. Default True.
P.saveica         = 1;

% Name added to save the ICA decomposition. Default 'ica_'
P.icaname         = 'wica_';

% Folder were the ICA decomposition is seved. If empty, in EEG.filepath
P.icapath         = [];

end
