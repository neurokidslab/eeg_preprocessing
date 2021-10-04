% -------------------------------------------------------------------------
% This function applies the robust detrending by de Chevaigne (NoiseTools)
% to data in EEGLAB format
%
% INPUTS
% EEG           EEG structure
% order         order of polynomial or number of sin/cosine pairs
% w             weights to compute the robus detrending
%               - [] --> none (default)
%               - string --> the matrix EEG.(string) is used
%               - [ta1 ta2; tb1 tb2] time intervals --> w=1 in [t1 t2], w=0 out [t1 t2]
%               - logical matrix of the size of the data
% basis         'polynomials' [default] or 'sinusoids', or user-provided matrix
% thresh        threshold for outliers [default: 3 sd]
% niter         number of iterations [default: 3]
% wsize         window size for local detrending [default: all] in seconds
%
% OPTIONAL INPUTS
% all possible inputs of the nt_detrend function
%
% OUTPUTS
% EEG           EEG structure
%
% -------------------------------------------------------------------------


function EEG = eega_ntdetrend(EEG, order, varargin)

fprintf('### Robust detrending ###\n')

% check that the NoiseTools are in the path
if exist('nt_detrend','file')~=2
    error('The NoiseTools is not in the path. Add them to the path (http://audition.ens.fr/adc/NoiseTools/)')
end

%% ------------------------------------------------------------------------
%% Parameters

P.weights   = [];  %[] | 'art' | [t1 t2; t1 t2; ...] | logical matrix of the size of the data
P.basis     = 'polynomials'; % 'polynomials' | 'sinusoids'
P.thresh    = 3;
P.niter     = 3;
P.wsize     = [];

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_normalization: Non recognized inputs')
end

% check the inputs
if isempty(P.weights)
    fprintf('No initial weights were provided\n')
else
    if ischar(P.weights)
        if strcmp(P.weights,'arttifacts.BCT')
            if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BCT')
                fprintf('Artifacted samples have an initial weight equal to 0\n')
                P.weights = ~EEG.artifacts.BCT;
                weightsmat = 1;
            else
                warning('%s field not found. No initial weights will be used',P.weights)
                P.weights = [];
                weightsmat = 0;
            end
        else
            if isfield(EEG,P.weights)
                fprintf('Initial weights based on the matrix in field %s\n', P.weights)
                P.weights = EEG.(P.weights);
                if size(P.weights,1)==1
                    P.weights = repmat(P.weights,[size(EEG.data,1) 1 1]);
                end
                if size(P.weights,3)==1
                    P.weights = repmat(P.weights,[1 1 size(EEG.data,3)]);
                end
                if any(size(P.weights)~=size(EEG.data))
                    error('The weights matrix has a wrong size')
                end
                weightsmat = 1;
            else
                warning('%s field not found. No initial weights will be used',P.weights)
                P.weights = [];
                weightsmat = 0;
            end
        end
    elseif isnumeric(P.weights)
        if size(P.weights,2)==2
            fprintf('Initial weights equal to 1 for samples in the time windows\n')
            tlimits = P.weights;
            P.weights = false(size(EEG.data,2),1);
            for i=1:size(tlimits,1)
                idx = (EEG.times>=tlimits(i,1)) & (EEG.times<=tlimits(i,2));
                P.weights(idx) = 1;
            end
            weightsmat = 0;
        elseif size(P.weights,2)==size(EEG.data,2)
            fprintf('Initial weights based on the provided matrix\n')
            if size(P.weights,1)==1
                P.weights = repmat(P.weights,[size(EEG.data,1) 1 1]);
            end
            if size(P.weights,3)==1
                P.weights = repmat(P.weights,[1 1 size(EEG.data,3)]);
            end
            if any(size(P.weights)~=size(EEG.data))
                error('The weights matrix has a wrong size')
            end
            weightsmat = 1;
        else
            error('Problem with weights')
        end
    else
        error('Problem with weights')
    end
end
if ~isempty(P.wsize)
    P.wsize = P.wsize*EEG.srate;
end

%% ------------------------------------------------------------------------
%% Robust detrending
fprintf('Epoch % 4.0d', 0)
for epi=1:size(EEG.data,3)
    fprintf(repmat('\b',[1 4]))
    fprintf('% 4.0d', epi)
    for chi=1:size(EEG.data,1)
        y = EEG.data(chi,:,epi)';
        if ~isempty(P.weights) && weightsmat == 1
            w = P.weights(chi,:,epi)';
        else
            w = P.weights;
        end
        [yout, wout] = nt_detrend(y, order, w, P.basis, P.thresh, P.niter, P.wsize);
        yout = yout';
        wout = wout';
        
        EEG.data(chi,:,epi) = yout;
    end
    
end

fprintf('\n')

end