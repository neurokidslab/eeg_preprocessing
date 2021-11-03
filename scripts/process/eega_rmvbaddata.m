% -------------------------------------------------------------------------
% Deal with bad data
%
% INPUTS
% EEG           EEG structure
%
% OPTIONAL INPUTS
% DataField     fields in EEG containig the data. It can be a cell with
%               multiple fields. In that case the baseline correction is
%               applied to all fileds. Default {'data'}
% FactorField   field in EEG containg the factors characterizing epochs.
%               Default 'F'
% BadData       how to treat bad data 'none' | 'replacebynan' | 'replacebyzero' | 'replacebymean' | 'replacebymeanxcnd'
%
% OUTPUTS
% EEG           EEG structure
%
% -------------------------------------------------------------------------

function [EEG, varargout ] = eega_rmvbaddata(EEG, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.DataField     = {'data'};
P.FactorField   = {'F'};
P.BadData       = 'none';   % 'none' | 'replacebynan' | 'replacebyzero' | 'replacebymean' | 'replacebymeanxcnd'
P.TableCNDs 	= [];
P.what          = 'all';  % 'all' | 'BCT' | 'BTBC' | 'BT' | 'BC'
P.Silent        = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_rmvbaddata: Non recognized inputs')
end

% Check the inputs
if ~any(strcmp(P.BadData,{'none', 'replacebynan' , 'replacebyzero', 'replacebymean','replacebymeanxcnd'}))
    error('eega_rmvbaddata: The option ''BadData'' can have values ''replacebynan'' / ''none''')
end
if ~iscell(P.DataField)
    P.DataField = {P.DataField};
end


%% ------------------------------------------------------------------------
%% Obtain some usefull data

sz = size(EEG.(P.DataField{1}));
data2rmv = false(sz);

if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BCT') && any(strcmp(P.what,{'all' 'BCT'}))
    if length(sz)<4
        data2rmv(EEG.artifacts.BCT) = 1;
    else
        idx = permute(EEG.artifacts.BCT,[1 4 2 3]);
        data2rmv(repmat(idx,[1 sz(2) 1 1])) = 1;
    end
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BT') && any(strcmp(P.what,{'all' 'BTBC' 'BT'}))
    if length(sz)<4
        data2rmv(repmat(EEG.artifacts.BT,[sz(1) 1 1])) = 1;
    else
        idx = permute(EEG.artifacts.BT,[1 4 2 3]);
        data2rmv(repmat(idx,[sz(1) sz(2) 1 1])) = 1;
    end
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BC') && any(strcmp(P.what,{'all' 'BTBC' 'BC'}))
    if length(sz)<4
        data2rmv(repmat(EEG.artifacts.BC,[1 sz(2) 1])) = 1;
    else
        idx = permute(EEG.artifacts.BC,[1 4 2 3]);
        data2rmv(repmat(idx,[1 sz(2) sz(3) 1])) = 1;
    end
end

%% ------------------------------------------------------------------------
%% Remove bad data
nBad=sum(data2rmv(:));
goodsbj = 1;
if nBad>0
    if ~P.Silent; fprintf('Percentage of bad data in the data : %d samples out of %d (%5.1f%%)\n',nBad,prod(sz),nBad/prod(sz)*100); end
    switch P.BadData
        case 'none'
            if ~P.Silent; fprintf(' ---> Bad data will be left \n'); end
            
        case 'replacebynan'
            if ~P.Silent; fprintf('---> Bad data will be replaced by NaNs \n'); end
            for i=1:length(P.DataField)
                EEG.(P.DataField{i})(data2rmv)=NaN;
            end
            
        case'replacebyzero'
            if ~P.Silent; fprintf('---> Bad data will be replaced by zeros \n'); end
            for i=1:length(P.DataField)
                EEG.(P.DataField{i})(data2rmv)=0;
            end
            
        case'replacebymean'
            if ~P.Silent; fprintf('---> Bad data will be replaced by the mean over all epochs \n'); end
            for i=1:length(P.DataField)
                M = nanmean(EEG.(P.DataField{i}),3);
                % scale the mean based on the standard deviations (otherwise the mean and the single trial data have different amplitudes)
                sdM = std(M,[],2);
                sdD = nanstd(reshape(EEG.(P.DataField{i}),[size(EEG.(P.DataField{i}),1) size(EEG.(P.DataField{i}),2)*size(EEG.(P.DataField{i}),3)]),[],2);
                M = M .* sdD ./ sdM;
                % replace the nans by the average
                M = repmat(M,[1 1 size(EEG.(P.DataField{i}),3)]);
                EEG.(P.DataField{i})(data2rmv)=M(data2rmv);
            end
            
        case'replacebymeanxcnd'
            if ~P.Silent; fprintf('---> Bad data will be replaced by the mean per condition \n'); end
            if isempty(P.TableCNDs)
                error('the optional input TableCNDs needs to be provided')
            end
            % Build the table defying the conditions
            if ~isempty(P.TableCNDs) && ~istable(P.TableCNDs)
                P.TableCNDs = cnd_buildtablecond(P.TableCNDs, EEG.(P.FactorField{1}));
            end
            if isempty(P.TableCNDs)
                goodsbj=0;
            else
                F = EEG.(P.FactorField{1});
                [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(P.TableCNDs, F);
                ntheCND = length(theCND);
            end
            if goodsbj
                for i=1:length(P.DataField)
                    for c=1:ntheCND
                        if sum(condition==c)>1
                            dcnd = EEG.(P.DataField{i})(:,:,condition==c);
                            data2rmvcnd = data2rmv(:,:,condition==c);
                            dcnd(data2rmvcnd) = nan;
                            M = nanmean(dcnd,3);
                            % scale the mean based on the standard deviations (otherwise the mean and the single trial data have different amplitudes)
                            sdM = std(M,[],2);
                            sdD = nanstd(reshape(dcnd,[size(dcnd,1) size(dcnd,2)*size(dcnd,3)]),[],2);
                            M = M .* sdD ./ sdM;
                            % replace the nans by the average
                            M = repmat(M,[1 1 size(dcnd,3)]);
                            dcnd(data2rmvcnd) = M(data2rmvcnd);
                            EEG.(P.DataField{i})(:,:,condition==c) = dcnd;
                        end
                    end
                    
                end
            end
            
    end
    
    varargout=cell(1,length(P.DataField));
    for i=1:length(P.DataField)
        varargout{i}=EEG.(P.DataField{i});
    end
    
    
end
end