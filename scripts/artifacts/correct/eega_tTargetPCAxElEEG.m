% -------------------------------------------------------------------------
% This function applies eega_tTargetPCAxEl on the EEG structure.
%
% INPUTS
% EEG   : EEG structure. Needs to have
%           * artifacts.BCT (defining bad electrodes x samples x epochs
%           * artifacts.BC (defining bad electrodes x 1 x epochs)
% nSV   : number of principal components to remove from the data
% vSV   : proportion of variance to remove from the data
%
% OPTIONAL INPUTS
% Plus all the optional inputs of eega_tTargetPCAxEl
%
% OUTPUT
% EEG   : EEG structure after target PCA correction
%
% USAGE
% [EEG] = eega_tTargetPCAxElEEG(EEG, [], 0.98, 'maxTime',0.08)
%
% -------------------------------------------------------------------------


function [EEG] = eega_tTargetPCAxElEEG( EEG, nSV, vSV, varargin)

fprintf('### Target PCA ###\n')

%% ------------------------------------------------------------------------
%% Parameters
savecorrected = 1;
alltime = 0;
allchannels = 0;
idxrmv = [];
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'maxtime') || strcmpi(varargin{i},'masktime') || strcmpi(varargin{i},'wsize')
        varargin{i+1} = round(varargin{i+1}*EEG.srate);
    end
    if strcmp(varargin{i},'savecorrected')
        savecorrected = varargin{i+1};
        idxrmv=[idxrmv i i+1];
    end
    if strcmp(varargin{i},'alltime')
        alltime = varargin{i+1};
        idxrmv=[idxrmv i i+1];
    end
    if strcmp(varargin{i},'allchannels')
        allchannels = varargin{i+1};
        idxrmv=[idxrmv i i+1];
    end
end
varargin(idxrmv) = [];

%% ------------------------------------------------------------------------
%% Interpolate
if alltime
    intertime=[];
else
    intertime=~EEG.artifacts.BT;
end
if allchannels
    interch=[];
else
    interch=~EEG.artifacts.BC;
end
if ~isempty(EEG.data)
    [ EEG.data, Interp ] = eega_tTargetPCAxEl( EEG.data,...
        EEG.artifacts.BCT,...
        intertime,...
        interch,...
        nSV,...
        vSV,...
        varargin{:});
    
    if savecorrected
        if ~isfield(EEG.artifacts,'CCT')
            EEG.artifacts.CCT = false(size(EEG.data));
        end
        EEG.artifacts.BCT(logical(Interp) & EEG.artifacts.BCT) = 0;
        EEG.artifacts.CCT = logical(EEG.artifacts.CCT) | logical(Interp);
        EEG.artifacts.summary = eega_summaryartifacts(EEG);
    end
    
else
    fprintf('No data, nothing will be done\n')
end
fprintf('\n')

end