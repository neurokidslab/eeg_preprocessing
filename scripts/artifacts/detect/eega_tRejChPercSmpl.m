% -------------------------------------------------------------------------
% This functions rejects all the samples for an electrode in an epoch if 
% too muach data was rejected
%
% INPUT
% EEG   EEG structure
%
% OPTIONAL INPUTS
%    - maxBadSmpl   maximun porportion of bad time for each channel 
%                   accepted (dfautult value 0.50)
%
% OUTPUTS
%   EEG     output data
%   BCT     bad data 
%   isbadch
%
% -------------------------------------------------------------------------

function [ EEG, BCT, isbadch ] = eega_tRejChPercSmpl( EEG, varargin )

fprintf('### Rejecting Channels based on the rejected Samples ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.maxBadSmpl = 0.50;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejChPercSmpl: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Get data and check that the artifact structure exists 
[nEl, nS, nEp] = size(EEG.data);
EEG = eeg_checkart(EEG);
BCTin = EEG.artifacts.BCT;

%% ------------------------------------------------------------------------
%% Define artifacts
BCT = false(nEl,nS,nEp);
nbad = sum(EEG.artifacts.BCT,2);
isbadch = nbad/nS > P.maxBadSmpl;
for e=1:nEp
    BCT(isbadch(:,1,e),:,e) = 1;
end

%% ------------------------------------------------------------------------
%% Display rejected data
n = nEl*nS*nEp;
new = BCT & ~BCTin;
new = sum(new(:));
fprintf('Total new data rejected %3.2f %%\n', new/n*100 )

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
    EEG.artifacts.algorithm.rejxstep = cat(1,EEG.artifacts.algorithm.rejxstep(:),sum(BCT(:)));
end

fprintf('\n' )
end
