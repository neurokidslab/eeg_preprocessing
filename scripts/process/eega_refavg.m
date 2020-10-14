%Function that performs an average reference

function [ EEG, ref ] = eega_refavg( EEG, varargin )

fprintf('### Reference to the mean ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.BadData = 'none';   % 'replacebynan' / 'none' / 'zero'
P.SaveRef = 0;   % keep the reference as an electrode 

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_refavg: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Reference
nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

[~, data_good] = eega_rmvbaddata(EEG, 'BadData', P.BadData);

% Reference to the mean
ref = nanmean(data_good,1);
if any(isnan(data_good(:)))
    ref(isnan(ref)) = 0;
end

EEG.data = EEG.data - repmat(ref,[nEl 1 1]);
EEG.ref = 'averef';
EEG.refval = ref;

if P.SaveRef
    % Save the reference and the data
    % EEG.reference = ref;
    EEG.data = cat(1,EEG.data,ref);
    if isfield(EEG, 'artifacts')
        if isfield(EEG.artifacts, 'BCT')
            EEG.artifacts.BCT=cat(1,EEG.artifacts.BCT,false(1,nS,nEp));
        end
        if isfield(EEG.artifacts, 'CCT')
            EEG.artifacts.CCT=cat(1,EEG.artifacts.CCT,false(1,nS,nEp));
        end
        if isfield(EEG.artifacts, 'BC')
            EEG.artifacts.BC=cat(1,EEG.artifacts.BC,false(1,1,nEp));
        end
    end
end

fprintf('\n' )

end

