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
% BadData       how to treat bad data 'replacebynan' | 'replacebymean' | 'none' | 'zero'
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
P.BadData       = 'none';   % 'replacebynan' / 'replacebymean' / 'none' / 'zero'
P.Silent        = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_rmvbaddata: Non recognized inputs')
end

% Check the inputs
if ~any(strcmp(P.BadData,{'replacebynan' , 'none', 'zero'}))
    error('eega_rmvbaddata: The option ''BadData'' can have values ''replacebynan'' / ''none''')
end
if ~iscell(P.DataField)
    P.DataField = {P.DataField};
end


%% ------------------------------------------------------------------------
%% Obtain some usefull data

sz = size(EEG.(P.DataField{1}));
data2rmv = false(sz);

if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BCT')
    if length(sz)<4
        data2rmv(EEG.artifacts.BCT) = 1;
    else
        idx = permute(EEG.artifacts.BCT,[1 4 2 3]);
        data2rmv(repmat(idx,[1 sz(2) 1 1])) = 1;
    end
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BT')
    if length(sz)<4
        data2rmv(repmat(EEG.artifacts.BT,[sz(1) 1 1])) = 1;
    else
        idx = permute(EEG.artifacts.BT,[1 4 2 3]);
        data2rmv(repmat(idx,[sz(1) sz(2) 1 1])) = 1;
    end
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BC')
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
if nBad>0
    if ~P.Silent; fprintf('Percentage of bad data in the data : %d samples out of %d (%5.1f%%)\n',nBad,prod(sz),nBad/prod(sz)*100); end
    if strcmp(P.BadData,'replacebynan')
        if ~P.Silent; fprintf('---> Bad data will be replaced by NaNs \n'); end
        for i=1:length(P.DataField)
            EEG.(P.DataField{i})(data2rmv)=NaN;
        end
    elseif strcmp(P.BadData,'zero')
        if ~P.Silent; fprintf('---> Bad data will be replaced by zeros \n'); end
        for i=1:length(P.DataField)
            EEG.(P.DataField{i})(data2rmv)=0;
        end
    elseif strcmp(P.BadData,'none')
        if ~P.Silent; fprintf(' ---> Bad data will be left \n'); end
    end
end
varargout=cell(1,length(P.DataField));
for i=1:length(P.DataField)
    varargout{i}=EEG.(P.DataField{i});
end


end