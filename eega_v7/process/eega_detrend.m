% implementation of the de Cheveigne detrending function to EEGlab data

function EEG = eega_detrend(EEG, order, weights , varargin)

fprintf('### Detrending ###\n')

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.basis     = [];
P.thresh    = [];
P.niter     = [];

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_detrend: Non recognized inputs')
end


%% ------------------------------------------------------------------------
%% Detrend

if isempty(weights)
    weights = ones(size(EEG.data));
end
if size(weights,2)~=size(EEG.data,2)
    error('weigths time dimentiona has to be equal o data time dimension')
end
if size(weights,1)~=size(EEG.data,1)
    weights=repmat(weights,[size(EEG.data,1) 1 1]);
end
if size(weights,3)~=size(EEG.data,3)
    weights=repmat(weights,[1 1 size(EEG.data,3)]);
end

for ep=1:size(EEG.data,3)
    [y,~,~] = nt_detrend(EEG.data(:,:,ep)', order, weights(:,:,ep)',...
        P.basis, P.thresh, P.niter);
    EEG.data(:,:,ep) = y';
end

fprintf('\n')

end