% -------------------------------------------------------------------------
% The function takes as input the EEG structure and interpolates for each
% epoch the channels marked as bad in EEG.artifacts.BC
%
% INPUTS:
% EEG   : EEG structure. Needs to have
%           * artifacts.BCT (defining bad electrodes x samples x epochs
%           * artifacts.BC (defining bad electrodes x 1 x epochs)
%           * chanlocs (defining the localization of the channels)
% p_int : maximun proportion of bad channels in a segment to interpolate
%
% OPTIONAL INPUTS
% savecorrected     : save the corrected data in the EEG structure and
%                     update the artifacts field
% Plus all the optional inputs of eega_tInterpSpatial
%
% OUTPUT:
% EEG   : EEG structure after interpolation
%
% USAGE
% [EEG] = eega_tInterpSpatialEEG(EEG, 0.4, 'distmethod', 'triangulation')
%
% -------------------------------------------------------------------------

function [ EEG ] = eega_tInterpSpatialEEG(EEG, p_int, varargin )

%% ------------------------------------------------------------------------
%% Parameters
savecorrected = 1;

for i=1:2:length(varargin)
    if strcmp(varargin{i},'savecorrected')
        savecorrected = varargin{i+1};
        varargin(i:i+1) = [];
    end
end

%% ------------------------------------------------------------------------
%% Interpolate
if ~isempty(EEG.data)
    
    % interpolate
    [ EEG.data, Interp ] = eega_tInterpSpatial( EEG.data, ~EEG.artifacts.BC, EEG.chanlocs, p_int, varargin{:});
    
    if savecorrected
        
        if ~isfield(EEG.artifacts,'CCT')
            EEG.artifacts.CCT = false(size(EEG.data));
        end
        
        EEG.artifacts.BCT(logical(Interp) & EEG.artifacts.BCT) = 0;
        EEG.artifacts.CCT = logical(EEG.artifacts.CCT) | logical(Interp);
        
        % mark samples for the interpolated channels as bad during BT
        idx = Interp & repmat(EEG.artifacts.BT,[size(EEG.data,1) 1 1]);
        EEG.artifacts.BCT(idx) = 1;
        
        % mark samples for the interpolated channels as bad if bad data was used for the interpolation
        bct = EEG.artifacts.BCT;
        bct(repmat(~EEG.artifacts.BC,[1 size(bct,2) 1]))=0;
        idx = Interp & repmat(~EEG.artifacts.BT,[size(EEG.data,1) 1 1]) & repmat(any(bct,1),[size(EEG.data,1) 1 1]);
        EEG.artifacts.BCT(idx) = 1;
        
        % mark interpolated channels as good in BC
        EEG.artifacts.BC(all(Interp,2)) = 0;
        
    end
    if exist('eega_summarypp','file')==2
        EEG = eega_summarypp(EEG);
    end
else
    fprintf('No data, nothing will be done\n')
end

fprintf('\n')
end

