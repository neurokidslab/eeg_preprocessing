% -------------------------------------------------------------------------
% This function identifies segments of bad data based on the amount of
% power in different frequency bands
% It is a aplied in a slidding time windows.
% In the case of a relative threshold it is determined across all time
% windows and across all electrodes
%
% INPUTS
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - thresh        upper and lower thresholds. Size n x 2 
%                   n = number of frequncy bands
%                   first column = lower limit
%                   second column = upper limit
%                   Default [-Inf 3]
%   - refdata       referenced average the data before (1) or not (0) (default 0)
%   - refbaddata    how to teat bad data when reference average ('replacebynan' / 'none' / 'zero', default 'none')
%   - dozscore      1/0 (default 1)
%   - relative      appply relative (1) or absolute (0) thresholds (default 1)
%   - twdur         time window length in seconds (default 5s)
%   - twstep        time window step in seconds (default 2s)
%   - frqband       frequncy band. Size n x 2 (default [20 40])
%                   n = number of frequency bands 
%   - mask          time to mask bad segments (default 0)
%
% OUTPUTS
%   EEG     output data
%   BCTS  rejection matrix
%   Ti    thresholds
%
% -------------------------------------------------------------------------


function [ EEG, BCT, Ti ] = eega_tRejPwr( EEG, varargin )

fprintf('### Rejecting based on the spectrum ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.refdata = 0;
P.refbaddata = 'none'; % 'replacebynan' / 'none' / 'zero'
P.twdur = 4;
P.twstep = 2;
P.dozscore = 1;
P.frqband = [1 10; 20 40];
P.relative = [1; 1];
P.thresh = [-3 Inf; -Inf 3];
P.mask = 0;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejSpectrum: Non recognized inputs')
end
if length(P.relative)<size(P.frqband,1)
    n = size(P.frqband,1)-length(P.relative);
    P.relative = [P.relative; P.relative(end)*ones(n,1) ];
end

% Check the inputs
if length(P.dozscore)>1 || any(P.dozscore(:)~=0 & P.dozscore(:)~=1)
    error('eega_tRejPwr: dozscore has to be 0 / 1')
end
if any(P.relative(:)~=0 & P.relative(:)~=1)
    error('eega_tRejPwr: relative has to have values 0 / 1')
end
if size(P.thresh,2)~=2
    error('eega_tRejPwr: The threshold has to be a matrix of size n x 2')
end
if size(P.thresh,1)~=size(P.frqband,1)
    error('eega_tRejPwr: The number of raws of the threshold matrix has to be equeal to number of raws of the frequency bands matrix')
end
if size(P.thresh,1)~=length(P.relative)
    warning('eega_tRejPwr: relative = %d for all frequency bands',P.relative)
    P.relative = repmat(P.relative(1),[size(P.thresh,1) 1]);
end

fprintf('- referenced data: %d\n',P.refdata)
fprintf('- z-score data: %d\n',P.dozscore)
fprintf('- relative threshold: %d\n',P.relative)
fprintf('\n')

%% ------------------------------------------------------------------------
%% Get data and check that the artifact structure exists 
[nEl, nS, nEp] = size(EEG.data);
EEG = eeg_checkart(EEG);

%% ------------------------------------------------------------------------
%% Reference data
if P.refdata
    [ EEG, reference ] = eega_refavg( EEG ,'BadData',P.refbaddata,'SaveRef',0);
end

%% ------------------------------------------------------------------------
%% Z-score
if P.dozscore
    [EEG, mu, sd] = eega_ZscoreForArt(EEG);
end

%% ------------------------------------------------------------------------
%% Define the time windows
if ~isempty(P.twdur) && ~(P.twdur==Inf)
    twdur = round(P.twdur*EEG.srate);
    twstep = round(P.twstep*EEG.srate);
    ntw = round((nS-twdur+1)/twstep)+1;
else
    ntw = 1;
    twdur = nS;
end
if any(ntw<=0)
    warning('The time window is too long')
    return
end
if ntw==1
    i_t = (1:nS)';
else
    i_t = linspace(1,nS-twdur+1,ntw);
    i_t = round(i_t);
    i_t = repmat(i_t,[twdur 1]) + repmat((0:twdur-1)',[1 size(i_t,2)]);
end
bct_tw = false(nEl,ntw,nEp);

%% ------------------------------------------------------------------------
%% Compute the spectrum

% Power in each frequency bin
L       = twdur;                % length of signal
Fs      = EEG.srate;            % sampling frequency
freq    = Fs*(0:(L/2))/L;       % frequencies
nfreq   = length(freq);
P1 = nan(nEl, nfreq, ntw, nEp);
for ep = 1:nEp
    for el = 1:nEl
        % data in each time window
        d = permute(EEG.data(el,:,ep),[2 1 3]);
        dtw = d(i_t);
        % Fourier transform
        fD = fft(dtw - repmat(mean(dtw,1),[size(dtw,1) 1]), L, 1); % fast fourier trasnform
        P2 = abs(fD).^2 / (L*Fs);  % two sides spectrum
        P1(el,:,:,ep) = P2(1:floor(L/2+1),:); % one side spectrum
    end
end

% Power in the different frequency bands
nfbands = size(P.frqband,1);
RRR = nan(nEl, nfbands, ntw, nEp);
for i=1:nfbands
    i_band = (freq > P.frqband(i,1)) & (freq <= P.frqband(i,2));
    
    p_band = log10( mean(P1(:,i_band,:,:),2) );
    p_base = nanmedian(p_band(:));
    p_band = p_band - p_base; 
    
%     p_band = mean( log10(P1(:,i_band,:,:)), 2);
%     p_base = mean(p_band,1);
%     p_band = p_band - repmat(p_base,[nEl 1 1 1]); 
    
    p_band = permute( p_band, [1 3 4 2] );
    RRR(:,i,:,:) = p_band;
end
% bad time windows
for i=1:ntw
    idx = i_t(:,i);
    bct_tw(:,i,:) = sum(EEG.artifacts.BCT(:,idx,:),2) >= (1.0*length(idx));
end


%% ------------------------------------------------------------------------
%% Reject bad electrodes per epoch and time window
Ru = false(size(RRR));
Rl = false(size(RRR));
for iband=1:nfbands
    if P.relative(iband)
        RRRi = RRR(:,iband,:,:);
        RRRi(permute(bct_tw,[1 4 2 3])) = nan;
        RRRi(isinf(RRRi)) = nan;
        RRRi(isnan(RRRi)) = [];
        perc = prctile(RRRi(:),[25 50 75]);
        IQ = perc(3) - perc(1);
        t_l = perc(1) + P.thresh(iband,1)*IQ;
        t_u = perc(3) + P.thresh(iband,2)*IQ;
        Ru(:,iband,:,:) = RRR(:,iband,:,:)>t_u;
        Rl(:,iband,:,:) = RRR(:,iband,:,:)<t_l;
        
%         for iep=1:nEp
%             for itw=1:ntw
%                 RRRi = RRR(:,iband,itw,iep);
%                 perc = prctile(RRRi,[25 50 75]);
%                 IQ   = perc(3) - perc(1);
%                 t_l  = perc(1) + P.thresh(iband,1)*IQ;
%                 t_u  = perc(3) + P.thresh(iband,2)*IQ;
%                 Ru(:,iband,itw,iep) = RRR(:,iband,:,:)>t_u;
%                 Rl(:,iband,itw,iep) = RRR(:,iband,:,:)<t_l;
%             end
%         end
    else
        t_u = P.thresh(iband,2);
        t_l = P.thresh(iband,1);
        Ru(:,iband,:,:) = RRR(:,iband,:,:)>t_u;
        Rl(:,iband,:,:) = RRR(:,iband,:,:)<t_l;
    end
    ru_sum = sum(sum(sum(Ru(:,iband,:,:),1),3),4);
    rl_sum = sum(sum(sum(Rl(:,iband,:,:),1),3),4);
    
    %% --------------------------------------------------------------------
    %% Display rejected data
    n = nEl*ntw*nEp;
    fprintf('Data rejected frequncy band [%2.0f %2.0f]: thresh_u %3.2f %%, thresh_l %3.2f %%\n',...
        P.frqband(iband,1), P.frqband(iband,2),ru_sum/n*100, rl_sum/n*100)
    
end
clear RRR

%% ------------------------------------------------------------------------
%% Go back to each sample
bcts_u = false(nEl,nS,nEp,nfbands);
bcts_l = false(nEl,nS,nEp,nfbands);
Ti = i_t(1,:);
Tf = i_t(end,:);
for i=1:nfbands
    for ep = 1:nEp
        for el = 1:nEl
            ru = Ru(el,i,:,ep);
            rl = Rl(el,i,:,ep);
            tiu = Ti(ru);
            tfu = Tf(ru);
            til = Ti(rl);
            tfl = Tf(rl);
            if length(tfu)>1
                iu = (tfu(1:end-1) - tiu(2:end)) >= 0 ;
                tiu(logical([0 iu])) = [];
                tfu(logical([iu 0])) = [];
            end
            if length(tfl)>1
                il = (tfl(1:end-1) - til(2:end)) >= 0 ;
                til(logical([0 il])) = [];
                tfl(logical([il 0])) = [];
            end
            for j=1:length(tiu)
                bcts_u(el,tiu(j):tfu(j),ep,i) = true;
            end
            for j=1:length(til)
                bcts_l(el,til(j):tfl(j),ep,i) = true;
            end
        end
    end
end
clear Ru Rl
BCT = any(bcts_u,4)  | any(bcts_l,4);

%% ------------------------------------------------------------------------
%% Update the rejection matrix
if P.updateBCT
    EEG.artifacts.BCT = EEG.artifacts.BCT | BCT;
end
if P.updatesummary
    EEG.artifacts.summary = eega_summaryartifacts(EEG);
end
if P.updatealgorithm && isfield(EEG.artifacts, 'algorithm')
    EEG.artifacts.algorithm.parameters = cat(1,EEG.artifacts.algorithm.parameters(:),{P});
    f = dbstack;
    EEG.artifacts.algorithm.stepname = cat(1,EEG.artifacts.algorithm.stepname(:),{f(1).name});
    EEG.artifacts.algorithm.rejxstep = cat(1,EEG.artifacts.algorithm.rejxstep(:),sum(BCT(:)));
end

%% ------------------------------------------------------------------------
%% Data back
if P.dozscore
    EEG.data = EEG.data.*repmat( sd, [1 nS nEp]) + repmat( mu, [1 nS nEp]);
end
if P.refdata
    EEG.data = EEG.data + repmat(reference,[size(EEG.data,1) 1 1]);
end

%% ------------------------------------------------------------------------
%% Mask around 
if ~isempty(P.mask) && P.mask~=0
    [ EEG, bctmask ] = eega_tMask( EEG, 'tmask', P.mask);
    BCT = BCT | bctmask;
    clear bctmask
end

%% ------------------------------------------------------------------------
%% Display the total rejected data
n = nEl*nS*nEp;
fprintf('Total data rejected %3.2f \n', sum(BCT(:))/n*100)
fprintf('\n' )


end
