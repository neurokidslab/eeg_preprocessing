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
P.timelim = 0.040;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tIncShortBad: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Get data
timelim = round(P.timelim*EEG.srate);

nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

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
EEG.artifacts.BCT(Inc) = 0;
EEG.artifacts.summary = eega_summaryartifacts(EEG);

fprintf('\n' )
end    


