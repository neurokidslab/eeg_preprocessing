% -------------------------------------------------------------------------
% This function rejects good segments sourronded by bad data if they are 
% shorte than a limit. Default value 0.2 s
%
% INPUT
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - timelim   time limit for the rejection (seconds) 
%
% OUTPUTS
%   EEG     output data
%   BCT     bad data 
%
% -------------------------------------------------------------------------

function [ EEG, BCT ] = eega_tRejShortGood( EEG, varargin )


fprintf('### Rejecting short good segments ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.timelim = 2.000;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejShortGood: Non recognized inputs')
end


%% ------------------------------------------------------------------------
%% Get data and check that the artifact structure exists 
[nEl, nS, nEp] = size(EEG.data);
EEG = eeg_checkart(EEG);
BCTin = EEG.artifacts.BCT;
BCT = false(nEl,nS,nEp);

timelim = round(P.timelim*EEG.srate);
if nS<=timelim
    timelim = nS-1;
end

%% ------------------------------------------------------------------------
%% Algorithm

for ep=1:nEp
    bct = BCTin(:,:,ep);
    for el=1:nEl
        goodt = ~bct(el,:);
        isgood_i = [goodt(1) ((goodt(2:end) - goodt(1:end-1)) == 1) ];
        isgood_i = find(isgood_i);
        isgood_f = [ ((goodt(1:end-1) - goodt(2:end)) == 1) goodt(end) ];
        isgood_f = find(isgood_f);
        if ~isempty(isgood_i)
            durgood = isgood_f - isgood_i + 1;
            isshort = find(durgood<=timelim);
            for i=1:numel(isshort)
                bct(el,isgood_i(isshort(i)):isgood_f(isshort(i))) = 1;
            end
        end
    end
    BCT(:,:,ep) = bct & ~BCTin(:,:,ep);
end

%% ------------------------------------------------------------------------
%% Display rejected data
n = nEl*nS*nEp;
new = BCT & ~BCTin;
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
    EEG.artifacts.algorithm.rejxstep = cat(1,EEG.artifacts.algorithm.rejxstep(:),sum(BCT(:)));
end

fprintf('\n' )

end
