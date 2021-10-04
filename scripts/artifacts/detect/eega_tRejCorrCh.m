function [ EEG, BCT, Ti ] = eega_tRejCorrCh( EEG, varargin )

fprintf('### Rejecting based on the channels correlations ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.thresh = 0.4;
P.refdata = 0;
P.refbaddata = 'none'; % 'replacebynan' / 'none' / 'zero'
P.twdur = 4;
P.twstep = 2;
P.dozscore = 0;
P.relative = 0;
P.topcorrch = 5;  % percentage of top correlations to compute the mean
P.mask = 0;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;


[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejCorrCh: Non recognized inputs')
end

% Check the inputs
if length(P.dozscore)>1 || ~any(P.dozscore==[0 1])
    error('eega_tRejCorrCh: dozscore has to be 0 / 1')
end
if any(P.relative(:)~=0 & P.relative(:)~=1)
    error('eega_tRejCorrCh: relative has to have values 0 / 1')
end
if length(P.thresh)>1
    error('eega_tRejCorrCh: The threshold has to be a single value')
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
% bad time windows
for i=1:ntw
    idx = i_t(:,i);
    bct_tw(:,i,:) = sum(EEG.artifacts.BCT(:,idx,:),2) >= (1.0*length(idx));
end

%% ------------------------------------------------------------------------
%% Compute the correlation
CC = nan(nEl, ntw, nEp);
for ep = 1:nEp
    for itw=1:ntw
        % take the data
        d = EEG.data(:,i_t(:,itw),ep);
        % compute the correlation with all the other channels
        corrs = abs(corrcoef(d'));
        % remove correlation with itself
        corrs(logical(eye(size(corrs,1)))) = NaN;
        % compute the average correlation for the top channels
        ptop = prctile(corrs,100-P.topcorrch,1);
        corrs(corrs<=repmat(ptop,[size(corrs,1) 1])) = NaN;
        avgcorrs = nanmean(corrs);
        % store the data
        CC(:,itw,ep) = avgcorrs;
    end
end

%% ------------------------------------------------------------------------
%% Reject bad electrodes per epoch and time window
if P.relative    
    CCi = CC;
    CCi(bct_tw) = nan;
    CCi(isnan(CCi)) = [];
    perc = prctile(CCi(:),[25 50 75]);
    IQ = perc(3) - perc(1);
    t_l = perc(1) - P.thresh*IQ;
    R = CC<t_l;
else
    t_l = P.thresh;
    R = CC<t_l;
end
rl_sum = sum(sum(sum(R,1),2),3);

%% --------------------------------------------------------------------
%% Display rejected data
n = nEl*ntw*nEp;
fprintf('Data rejected due to a low correlation: thresh_l %3.2f %%\n', rl_sum/n*100)

%% ------------------------------------------------------------------------
%% Go back to each sample
BCT = false(nEl,nS,nEp);
Ti = i_t(1,:);
Tf = i_t(end,:);
for ep = 1:nEp
    for el = 1:nEl
        rl = R(el,:,ep);
        til = Ti(rl);
        tfl = Tf(rl);
        if length(tfl)>1
            il = (tfl(1:end-1) - til(2:end)) >= 0 ;
            til(logical([0 il])) = [];
            tfl(logical([il 0])) = [];
        end
        for j=1:length(til)
            BCT(el,til(j):tfl(j),ep) = true;
        end
    end
end
clear R

%% ------------------------------------------------------------------------
%% Update the rejection matrix
if P.updateBCT
    EEG.artifacts.BCT = EEG.artifacts.BCT | BCT;
end
if P.updatesummary
    EEG.artifacts.summary = eega_summaryartifacts(EEG);
end
if P.updatealgorithm
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