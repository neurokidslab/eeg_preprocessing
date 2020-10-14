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
P.timelim = 0.200;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejShortGood: Non recognized inputs')
end


%% ------------------------------------------------------------------------
%% Get data
nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

BCTin = EEG.artifacts.BCT;
BCT = false(nEl,nS,nEp);

timelim = round(P.timelim*EEG.srate);

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
            durgood = isgood_f - isgood_i;
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
EEG.artifacts.BCT = EEG.artifacts.BCT | BCT;
EEG.artifacts.summary = eega_summaryartifacts(EEG);

fprintf('\n' )

end
