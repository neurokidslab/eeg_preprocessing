% INPUTS
% Y         data channels x samples x trials
% srate     sampling rate
% Ncmp      components to keep
% tbias     trials to compute the bias filter (trials x number of filters)
% tnoise    trials to compute the "noise" (trials x number of filters)
% tapply    trials to apply the filter (trials x number of filters)

function [freq, Xdft, W, E] = ems_freq(Y, srate, Ncmp, tbias, tnoise, tapply, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.fmin          = 0;
P.fmax          = 100;
P.fneig         = 0.25;

% Optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dft_ems: Non recognized inputs')
end
if size(tbias,1)~=size(Y,3)
    error('dft_ems: Different number of trials in the data and biased vector')
end
if isempty(tnoise)
    tnoise = true(size(tbias));
end
if isempty(tapply)
    tapply = tbias;
end
if size(tnoise,1)~=size(Y,3)
    error('dft_ems: Different number of trials in the data and noise vector')
end
if size(tapply,1)~=size(Y,3)
    error('dft_ems: Different number of trials in the data and apply vector')
end
if size(tnoise,2)~=size(tbias,2)
    error('dft_ems: Different number of filters in the bias and noise vector')
end
if size(tapply,2)~=size(tbias,2)
    error('dft_ems: Different number of filters in the bias and apply vector')
end
if any(sum(tapply,2)>1)
    warning('dft_ems: Two differen filters apply to the same trial. Only one will be kept')
end

%% ------------------------------------------------------------------------
%% compute DFT
[Ne, Ns, Nt] = size(Y);
f_res = srate / Ns;             % resolution (Hz)
freq  = (0:floor(Ns/2))*f_res;  % frequencies
idx_freq = logical(freq>=P.fmin & freq<=P.fmax);
freq = freq(idx_freq);
Ydft = fft(Y,Ns,2);
Ydft = Ydft(:,1:floor(Ns/2)+1,:);
Ydft = Ydft(:,idx_freq,:);
Nfrq = length(freq);

%% ------------------------------------------------------------------------
%% spatial filter
Nfilt = size(tbias,2);
Xdft = nan(Ncmp, Nfrq, Nt);
W = nan(Ne, Ncmp, Nfrq, Nt);
E = nan(Ncmp, Nfrq, Nt);
for iflt=1:Nfilt
    fprintf('EMS %d: bias on %d trials, noise based on %d trials, apply to %d trials...\n',iflt,sum(tbias(:,iflt)), sum(tnoise(:,iflt)), sum(tapply(:,iflt)))
    
    if (sum(tbias(:,iflt))>2) && (sum(tnoise(:,iflt))>2) && (sum(tapply(:,iflt))>1)
        
        idnoise = find(tnoise(:,iflt));
        idbias = find(tbias(:,iflt));
        idapply = find(tapply(:,iflt));
        
        expvar = nan(1,Nfrq);
        for ifrq=1:Nfrq
            
                     
            % compute the noise covariance
            C0 = zeros(Ne,Ne);
            for i=1:length(idnoise)
                ydft = squeeze(Ydft(:,ifrq,idnoise(i)));
                C0 = C0 + (ydft*ydft');
            end
            C0 = real(C0);
            C0 = C0 *(1/length(idnoise));
                        
            % compute the covariance of the mean (bias)
            Ydfm = mean(Ydft(:,ifrq,tbias(:,iflt)),3);
            C1 = (Ydfm*Ydfm');
            C1 = real(C1);
            
            % eigenvalues decomposition
            [V, S] = eig(C1,C0);
            L = abs(real(S));
            [L, idx] = sort(diag(L), 'descend') ;
            V = V(:,idx);
            
            if any(isinf(L))
                idxinf = isinf(L);
                L(idxinf) = 1;
                L(~idxinf) = 0;
            end
            expvarxcmp = L/sum(L)*100;
            expvar(:,ifrq) = sum(L(1:Ncmp))/sum(L)*100;
            E(:,ifrq,tapply(:,iflt)) = repmat(expvarxcmp(1:Ncmp), [1 1 sum(tapply(:,iflt))]);
        
            % keep the first components as spatial filters
            w = V;
%             w = C1 * V / (V' * C1 * V);
            w = bsxfun(@rdivide,w,sqrt(sum(w.^2,1))); % normalize vectors (not really necessary, but OK)
            w = w(:,1:Ncmp);
            W(:,:,ifrq,tapply(:,iflt)) = repmat(w,[1 1 sum(tapply(:,iflt))]);
            
            % project each trial
            y = squeeze(Ydft(:,ifrq,idapply));
            x = y' * w;
            Xdft(:,ifrq,idapply) = x';
                      
        end
        fprintf('- the explained variance by the %d first EMS at each frequency bin is in the range [%4.2f %4.2f]\n',Ncmp,min(expvar), max(expvar))
        
    else
        warning('Not enougth good data')
    end
    
end
end