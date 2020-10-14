% -------------------------------------------------------------------------
% Function that applies a baseline correction
%
% INPUTS
% EEG           EEG structure
% timerange     time window to compute the baseline
%
% OPTIONAL INPUTS
% DataField     fields in EEG containig the data. It can be a cell with 
%               multiple fields. In that case the baseline correction is 
%               applied to all fileds. Default {'data'}
% FactorField   field in EEG containg the factors characterizing epochs. 
%               Default 'F'
% BadData       how to treat bad data 
%               Possible values: 'replacebynan' | 'replacebymean' | 'none'
%               Default 'none'
% Type          How to compute the baseline 
%               Possible values: 'xtrial | 'alltrials' | 'xcond'
%               Default 'xtrial'
% TableCNDs     cell array with the factors used to define the conditions.
%               Only requiered if Type = 'xcond'
% timereference time point from where the time range is measured. Default 0
%
% OUTPUTS
% EEG           EEG structure
%
% -------------------------------------------------------------------------

function EEG = eega_rmbaseline(EEG, timerange, varargin)

%% ========================================================================
%% Parameters

P.DataField 	= {'data'};
P.FactorField 	= 'F';
P.BadData       = 'none';   % 'replacebynan' / 'replacebymean' / 'none'
P.Type          = 'xtrial';  % xtrial / alltrials / xcond
P.TableCNDs 	= [];
P.timereference = 0; % time point from where the time range is measured

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_MergeData: Non recognized inputs')
end

if ischar(P.DataField)
    P.DataField = {P.DataField};
end

%% ========================================================================
%% Remove bad data
[ eeg ] = eega_rmvbaddata(EEG, 'BadData', P.BadData, 'DataField', P.DataField{1});

%% ========================================================================
%% Remove baseline

if ischar(P.timereference)
    for iep=1:size(EEG.data,3)
        ev = find(strcmp(P.timereference,EEG.epoch(iep).eventtype),1);
        timeref = EEG.epoch(iep).eventlatency{ev};
        idxt = false(size(EEG.times));
        for i=1:size(timerange,1)
            idxt = idxt | ((EEG.times>=timerange(i,1)+timeref) & (EEG.times<=timerange(i,2)+timeref));
        end
        bl = nanmean(eeg.(P.DataField{1})(:,idxt,iep),2);
        EEG.(P.DataField{1})(:,:,iep) = bsxfun(@minus,EEG.(P.DataField{1})(:,:,iep),bl);
    end
elseif ismatrix(timerange)
    idxt = false(size(EEG.times));
    for i=1:size(timerange,1)
        idxt = idxt | ((EEG.times>=timerange(i,1)) & (EEG.times<=timerange(i,2)));
    end
    
    if strcmp(P.Type,'xtrial')
        bl = nanmean(eeg.(P.DataField{1})(:,idxt,:),2);
        EEG.(P.DataField{1})= bsxfun(@minus,EEG.(P.DataField{1}),bl);
    elseif strcmp(P.Type,'alltrials')
        bl = nanmean(nanmean(eeg.(P.DataField{1})(:,idxt,:),2),3);
        EEG.(P.DataField{1})= bsxfun(@minus,EEG.(P.DataField{1}),bl);
    elseif strcmp(P.Type,'xcond') && ~isempty(P.TableCNDs)
        F=EEG.(P.FactorField);
        TableCNDs = cnd_buildtablecond(P.TableCNDs, F);
        [condition, ~, ~, ~] = cnd_findtrialsxcond(TableCNDs, F);
        CNDs=unique(condition);
        ntheCND = length(CNDs);
        for i=1:ntheCND
            bl = nanmean(nanmean(eeg.(P.DataField{1})(:,idxt,condition==CNDs(i)),2),3);
            EEG.(P.DataField{1})(:,:,condition==CNDs(i)) = bsxfun(@minus,EEG.(P.DataField{1})(:,:,condition==CNDs(i)),bl);
        end
    elseif strcmp(P.Type,'xcond') && isempty(P.TableCNDs)
        error('conditions have to be provided')
    end
end

end

