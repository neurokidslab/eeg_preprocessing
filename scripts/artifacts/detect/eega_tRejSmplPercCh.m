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
P.maxBadCh = 0.35;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejSmplPercCh: Non recognized inputs')
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
EEG.artifacts.BCT = EEG.artifacts.BCT | BCT;
EEG.artifacts.summary = eega_summaryartifacts(EEG);

fprintf('\n' )

end

