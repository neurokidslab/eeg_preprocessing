% -------------------------------------------------------------------------
% The function takes as input the EEG structure and interpolates for each
% channels the bad segments if they are comprised in the InterTime
%
% INPUTS
% EEG   : EEG structure. Needs to have
%           * artifacts.BCT (defining bad electrodes x samples x epochs
%           * artifacts.BC (defining bad electrodes x 1 x epochs)
%           * chanlocs (defining the localization of the channels)
% p_int : maximun proportion of bad channels in a segment to interpolate
%
% OPTIONAL INPUTS
% savecorrected     : save the corrected data in the EEG structure and
%                     update the artifacts field
% Plus all the optional inputs of eega_tInterpSpatialSegment
%
% OUTPUT:
% EEG   : EEG structure after interpolation
%
% USAGE
% [EEG] = eega_tInterpSpatialSegmentEEG(EEG, 0.4, 'distmethod', 'triangulation')
%
% -------------------------------------------------------------------------

function [ EEG ] = eega_tInterpSpatialSegmentEEG(EEG, p_int, varargin )

fprintf('### Spatial interpolatiom of bad segments ###\n')

%% ------------------------------------------------------------------------
%% Parameters
savecorrected = 1;
mingoodtime = 0;
idxrmv = false(1,length(varargin));
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'minintertime') || strcmp(varargin{i},'masktime')
        varargin{i+1} = round(varargin{i+1}*EEG.srate);
    end
    if strcmpi(varargin{i},'savecorrected')
        savecorrected = varargin{i+1};
        idxrmv(i:i+1) = 1;
    end
    if strcmpi(varargin{i},'mingoodtime')
        mingoodtime = varargin{i+1};
        idxrmv(i:i+1) = 1;
    end
end
varargin(idxrmv) = [];

%% ------------------------------------------------------------------------
%% Interpolate
if ~isempty(EEG.data)
    
    % mark as bad too short bad segments
    if mingoodtime~=0
        [ EEG, ~ ] = eega_tRejShortGood( EEG, 'timelim', mingoodtime );
    end
    
    % interpolate
    [ EEG.data, Interp ] = eega_tInterpSpatialSegment( EEG.data,...
                                EEG.artifacts.BCT,...
                                ~EEG.artifacts.BC,...
                                ~EEG.artifacts.BT,...
                                EEG.chanlocs,...
                                p_int,...
                                varargin{:});
    
    % mark the interpolated data                        
    if savecorrected
        if ~isfield(EEG.artifacts,'CCT')
            EEG.artifacts.CCT = false(size(EEG.data));
        end
        
        EEG.artifacts.BCT(Interp & EEG.artifacts.BCT) = 0;
        EEG.artifacts.CCT = EEG.artifacts.CCT | Interp;
    end
    if exist('eega_summarypp','file')==2
        EEG = eega_summarypp(EEG);
    end
else
    fprintf('No data, nothing will be done\n')
end

end

