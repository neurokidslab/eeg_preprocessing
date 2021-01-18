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
P.tmask = 0.200;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tMask: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Get data
nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

BCTin = EEG.artifacts.BCT;
BCT = false(nEl,nS,nEp);

%% ------------------------------------------------------------------------
%% Algorithm

art_buffer = round(P.tmask*EEG.srate); % time in seconds times sample rate

for ep=1:nEp
    for el = 1:nEl
        bad_i = diff([0; BCTin(el,:,ep)'])==1;  % begining of bad segments
        bad_f = diff([BCTin(el,:,ep)'; 0])==-1; % end of bad segments
        bad_i = find(bad_i);
        bad_f = find(bad_f);
        for j=1:length(bad_i)
            bad_idx_i = ((bad_i(j)-art_buffer):(bad_i(j)-1));
            bad_idx_i(bad_idx_i<=0) = [];
            bad_idx_f = ((bad_f(j)+1):(bad_f(j)+art_buffer));
            bad_idx_f(bad_idx_f>nS) = [];
            BCT(el,bad_idx_i,ep) = 1;
            BCT(el,bad_idx_f,ep) = 1;
        end
    end
    
end


%% ------------------------------------------------------------------------
%% Display rejected data
n = nEl*nS*nEp;
new = BCT & ~BCTin;
new = sum(new(:));
fprintf('Total data rejected %3.2f %%\n', new/n*100 )

%% ------------------------------------------------------------------------
%% Update the rejection matrix
EEG.artifacts.BCT = BCTin | BCT;
EEG.artifacts.summary = eega_summaryartifacts(EEG);

fprintf('\n' )
end





