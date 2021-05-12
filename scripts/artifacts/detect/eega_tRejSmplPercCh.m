% -------------------------------------------------------------------------
% This functions defines which samples are bad based on the amount of 
% rejected channels. 
%
% INPUT
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - maxBadCh     maximun porportion of bad channels in order to not reject
%                  a sample. Default 0.35
%
% OUTPUTS
%   EEG         output data
%   BCT         bad data 
%   isbadsmpl   
%
% -------------------------------------------------------------------------


function [ EEG, BCT, isbadsmpl ] = eega_tRejSmplPercCh( EEG, varargin )

fprintf('### Rejecting Samples based on the rejected Channels ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.maxBadCh = 0.30;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejSmplPercCh: Non recognized inputs')
end


%% ------------------------------------------------------------------------
%% Get data and check that the artifact structure exists 
[nEl, nS, nEp] = size(EEG.data);
EEG = eeg_checkart(EEG);
BCTin = EEG.artifacts.BCT;

%% ------------------------------------------------------------------------
%% Define artifacts
BCT = false(nEl,nS,nEp);
nbad = sum(EEG.artifacts.BCT,1);
isbadsmpl = nbad/nEl > P.maxBadCh;
for e=1:nEp
    BCT(:,isbadsmpl(1,:,e),e) = 1;
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

