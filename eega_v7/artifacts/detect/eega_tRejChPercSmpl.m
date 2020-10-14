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

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejChPercSmpl: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Get the data from the EEG structure
BCTin = EEG.artifacts.BCT;
nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

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
EEG.artifacts.BCT = EEG.artifacts.BCT | BCT;
EEG.artifacts.summary = eega_summaryartifacts(EEG);

fprintf('\n' )
end
