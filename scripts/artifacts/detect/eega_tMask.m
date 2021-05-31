% -------------------------------------------------------------------------
% Function that reject data before and after a mad segment
%
% INPUTS
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - tmask   : duration of the mask in seconds. Default 0.20
%
% OUTPUTS
%   EEG     output data
%   BCT     bad data 
%
% -------------------------------------------------------------------------

function [ EEG, BCT ] = eega_tMask( EEG, varargin )

fprintf('### Mask around artifacts ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.tmask = 0.050;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tMask: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Get data and check that the artifact structure exists 
[nEl, nS, nEp] = size(EEG.data);
EEG = eeg_checkart(EEG);

%% ------------------------------------------------------------------------
%% Algorithm
BCT = eega_maskmatrix(EEG.artifacts.BCT, P.tmask, EEG.srate);

%% ------------------------------------------------------------------------
%% Display rejected data
n = nEl*nS*nEp;
new = BCT & ~EEG.artifacts.BCT;
new = sum(new(:));
fprintf('Total data rejected %3.2f %%\n', new/n*100 )

%% ------------------------------------------------------------------------
%% Update the rejection matrix
if P.updateBCT
    EEG.artifacts.BCT = EEG.artifacts.BCT | BCT;
end
if P.updatesummary
    EEG.artifacts.summary = eega_summaryartifacts(EEG);
end
if P.updatealgorithm
    EEG.artifacts.algorithm.parameters = cat(1,EEG.artifacts.algorithm.parameters(:),{P});
    f = dbstack;
    EEG.artifacts.algorithm.stepname = cat(1,EEG.artifacts.algorithm.stepname(:),{f(1).name});
    EEG.artifacts.algorithm.rejxstep = cat(1,EEG.artifacts.algorithm.rejxstep(:),new);
end

fprintf('\n' )
end





