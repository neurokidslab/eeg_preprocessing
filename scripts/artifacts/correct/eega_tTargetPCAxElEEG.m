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
P.savecorrected = 1;
P.alltime       = 'nobadtime';
P.allchannels   = 'nobadch';
P.idxrmv        = [];
P.maxTime       = 0.100;
P.maskTime      = 0.050;
P.splicemethod  = 1; % 
P.order         = 3; % order of polynomial to fit the trend if P.splicemethod = 4
P.wsize         = 4; % time to mask the bad segment to detrend if P.splicemethod = 4
P.silent        = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tTargetPCAxElEEG: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Interpolate
switch P.alltime
    case 'all'
        intertime=[];
    case 'nobadtime'
        intertime=~EEG.artifacts.BT;
    case 'badtime'
        intertime=EEG.artifacts.BT;
end
switch P.allchannels
    case 'all'
        interch=[];
    case 'nobadch'
        interch=~EEG.artifacts.BC;
    case 'badch'
        interch=EEG.artifacts.BC;
end
if ~isempty(EEG.data)
    
    % target PCA
    [ EEG.data, Interp ] = eega_tTargetPCAxEl( EEG.data,...
        EEG.artifacts.BCT,...
        intertime,...
        interch,...
        nSV,...
        vSV,...
        'maxTime',round(P.maxTime*EEG.srate),...
        'maskTime',round(P.maskTime*EEG.srate),...
        'order',P.order,...
        'splicemethod',P.splicemethod,...
        'wsize',round(P.wsize*EEG.srate),...
        'silent',P.silent);

    % mark the interpolated data                        
    if P.savecorrected
        if ~isfield(EEG.artifacts,'CCT')
            EEG.artifacts.CCT = false(size(EEG.data));
        end
        EEG.artifacts.BCT(logical(Interp) & EEG.artifacts.BCT) = 0;
        EEG.artifacts.CCT = logical(EEG.artifacts.CCT) | logical(Interp);
    end
    if exist('eega_summarypp','file')==2
        EEG = eega_summarypp(EEG);
    end
else
    fprintf('No data, nothing will be done\n')
end
fprintf('\n')

end