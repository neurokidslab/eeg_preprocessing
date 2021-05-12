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
P.savecorrected     = 1;
P.mingoodtime       = 2.00;
P.minintertime      = 0.10;
P.masktime          = 1.00;
P.distneighbour     = [];
P.distmethod        = 'triangulation'; % 'distance', 'triangulation' or 'template'
P.pneigh            = 1;
P.maxloop           = 10;
P.splicemethod      = 1; % 0 / 1 / 2 / 3
P.silent            = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tInterpSpatialSegmentEEG: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Interpolate
if ~isempty(EEG.data)
    
    % mark as bad too short bad segments
    if P.mingoodtime~=0
        [ EEG, ~ ] = eega_tRejShortGood( EEG, 'timelim', P.mingoodtime );
    end
    
    % interpolate
    [ EEG.data, Interp ] = eega_tInterpSpatialSegment( EEG.data,...
                                EEG.artifacts.BCT,...
                                ~EEG.artifacts.BC,...
                                ~EEG.artifacts.BT,...
                                EEG.chanlocs,...
                                p_int,...
                                'minintertime',round(P.minintertime*EEG.srate),...
                                'masktime',round(P.masktime*EEG.srate),...
                                'distneighbour',P.distneighbour,...
                                'distmethod',P.distmethod,...
                                'pneigh',P.pneigh,...
                                'maxloop',P.maxloop,...
                                'splicemethod',P.splicemethod,...
                                'silent',P.silent);
    
    % mark the interpolated data                        
    if P.savecorrected
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

