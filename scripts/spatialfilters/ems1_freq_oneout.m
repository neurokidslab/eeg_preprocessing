% INPUTS
% Y         data channels x samples x trials
% srate     sampling rate
% Ncmp      components to keep
% tbias     trials to compute the bias filter (trials x number of filters)
% tnoise    trials to compute the "noise" (trials x number of filters)
% tapply    trials to apply the filter (trials x number of filters)

function [freq, Xdft, W] = ems1_freq_oneout(Y, srate, tbias, tnoise, tapply, varargin)

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
    tnoise = tbias;
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
        for ifrq=1:Nfrq
            
            
            idnoise = find(tnoise(:,iflt));
            idbias = find(tbias(:,iflt));
            idapply = find(tapply(:,iflt));
            
            % compute the noise covariance
            % -------------------------------------------------------------
            C0 = zeros(Ne,Ne);
            for i=1:length(idnoise)
                ydft = squeeze(Ydft(:,ifrq,idnoise(i)));
                C0 = C0 + (ydft*ydft');
            end
            C0 = real(C0);
            C0 = C0 *(1/length(idnoise));
            
            % compute per condition to align the trials
            % -------------------------------------------------------------
            
            % compute the covariance of the mean (bias)
            Ydfm = mean(Ydft(:,ifrq,tbias(:,iflt)),3);
            
            % get the filters
%             w = inv(C0)*Ydfm;
            wcnd = C0\Ydfm;
            wcnd = real(wcnd);
            wcnd = bsxfun(@rdivide,wcnd,sqrt(sum(wcnd.^2,1))); % normalize vectors (not really necessary, but OK)
            
            
            % compute per trials
            % -------------------------------------------------------------
            for ti=1:length(idapply)
                
                t = idapply(ti);
                
                % find the trials to bias
                t0t = tbias(:,iflt);
                t0t(t) = 0;
                
                % compute the covariance of the mean (bias)
                Ydfm = mean(Ydft(:,ifrq,t0t),3);
                
                % get the filters
                w = C0\Ydfm;
                w = real(w);
                w = bsxfun(@rdivide,w,sqrt(sum(w.^2,1))); % normalize vectors (not really necessary, but OK)
                
                % align with the spacial filters of other trials
                for k=1:Ncmp
                    r = corrcoef(wcnd(:,k),w(:,k));
                    if r(1,2)<0
                        w_sign = -1;
                    else
                        w_sign = 1;
                    end
                    w(:,k) = w_sign * w(:,k); 
                end
                 W(:,:,ifrq,t) = w;
                
                % project each trial 
                y = Ydft(:,ifrq,t);
                x = y' * w;
                Xdft(:,ifrq,t) = x;
                
               
            end
                
            
        end
    else
        warning('Not enougth good data')
    end
    
end

end