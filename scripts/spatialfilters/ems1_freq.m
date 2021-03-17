% INPUTS
% Y         data channels x samples x trials
% srate     sampling rate
% tbias     trials to compute the bias filter (trials x number of filters)
% tnoise    trials to compute the "noise" (trials x number of filters)
% tapply    trials to apply the filter (trials x number of filters)

function [freq, Xdft, W] = ems1_freq(Y, srate, tbias, tnoise, tapply, varargin)

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
Xdft = nan(1, Nfrq, Nt);
W = nan(Ne, 1, Nfrq, Nt);
for iflt=1:Nfilt
    fprintf('EMS %d: bias on %d trials, noise based on %d trials, apply to %d trials...\n',iflt,sum(tbias(:,iflt)), sum(tnoise(:,iflt)), sum(tapply(:,iflt)))
    
    if (sum(tbias(:,iflt))>2) && (sum(tnoise(:,iflt))>2) && (sum(tapply(:,iflt))>1)
        
        idnoise = find(tnoise(:,iflt));
        idbias = find(tbias(:,iflt));
        idapply = find(tapply(:,iflt));
        

        for ifrq=1:Nfrq
            
            % compute the noise covariance
            C0 = zeros(Ne,Ne);
            for i=1:length(idnoise)
                ydft = squeeze(Ydft(:,ifrq,idnoise(i)));
                C0 = C0 + (ydft*ydft');
            end
            C0 = real(C0);
            C0 = C0 *(1/length(idnoise));
                        
            % compute the mean for the bias trials
            Ydfm = mean(Ydft(:,ifrq,tbias(:,iflt)),3);
            
            % get the filters
%             w = inv(C0)*Ydfm;
            w = C0\Ydfm;
%             w = pinv(C0)*Ydfm;
            w = real(w);
            w = bsxfun(@rdivide,w,sqrt(sum(w.^2,1))); % normalize vectors (not really necessary, but OK)
            W(:,:,ifrq,tapply(:,iflt)) = repmat(w,[1 1 sum(tapply(:,iflt))]);
            
            % project each trial
            y = squeeze(Ydft(:,ifrq,idapply));
            x = y' * w;
            Xdft(:,ifrq,idapply) = x';
            
            
                      
        end
            
    else
        warning('Not enougth good data')
    end
    
end
end