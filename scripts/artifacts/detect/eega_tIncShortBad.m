% -------------------------------------------------------------------------
% This function includes back segments that are rejected but are too short
% 
% INPUTS
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - timelim       limit to reinclude, in seconds (default 0.04)
%
% OUTPUTS
%   EEG     output data
%   Inc     
%
% -------------------------------------------------------------------------

function [ EEG, Inc ] = eega_tIncShortBad( EEG, varargin )

fprintf('### Do not reject too short segments ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.timelim = 0.020;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tIncShortBad: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Get data and check that the artifact structure exists 
[nEl, nS, nEp] = size(EEG.data);
EEG = eeg_checkart(EEG);

timelim = round(P.timelim*EEG.srate);

restore = 0;
Inc = false(nEl,nS,nEp);
for ep=1:nEp
    for el=1:nEl
        badt = EEG.artifacts.BCT(el,:,ep);
        isbad_i = [badt(1) ((badt(2:end) - badt(1:end-1)) == 1) ];
        isbad_i = find(isbad_i);
        isbad_f = [ ((badt(1:end-1) - badt(2:end)) == 1) badt(end) ];
        isbad_f = find(isbad_f);
        if ~isempty(isbad_i)
            durgood = isbad_f - isbad_i;
            isshort = find(durgood<timelim);
            for i=1:numel(isshort)
                Inc(el,isbad_i(isshort(i)):isbad_f(isshort(i)),ep) = 1;
                restore = restore + (isbad_f(isshort(i)) - isbad_i(isshort(i)) + 1);
            end
        end
    end 
end

%% ------------------------------------------------------------------------
%% Display rejected data
n = nEl*nS*nEp;
fprintf('Total data re-included %3.2f %%\n', restore/n*100 )

%% ------------------------------------------------------------------------
%% Update the rejection matrix
if P.updateBCT
    EEG.artifacts.BCT(Inc) = 0;
end
if P.updatesummary
    EEG.artifacts.summary = eega_summaryartifacts(EEG);
end
if P.updatealgorithm
    EEG.artifacts.algorithm.parameters = cat(1,EEG.artifacts.algorithm.parameters(:),{P});
    f = dbstack;
    EEG.artifacts.algorithm.stepname = cat(1,EEG.artifacts.algorithm.stepname(:),{f(1).name});
    EEG.artifacts.algorithm.rejxstep = cat(1,EEG.artifacts.algorithm.rejxstep(:),-sum(Inc(:)));
end

fprintf('\n' )
end    


