% -------------------------------------------------------------------------
% Function that performs an artefact rejection algorithm based on the
% variace of the signal. The time windows showing an abnormal variance
% are rejected
%
% INPUT
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - thresh        lower and upper threshold [thresh_low, thresh_up] (default [-2 2]))
%   - refdata       referenced average the data before (1) or not (0) (default 0)
%   - refbaddata    how to treat bad data when reference average ('replacebynan' / 'none' / 'zero', default 'none')
%   - dozscore      z-score the data per electrodes before (1) or not (0) (default 0)
%   - relative      appply relative (1) or absolute (0) thresholds (default 1)
%   - xelectrode    appply the threhold per electrode (1) or over all electrodes (0) (default 1)
%   - twdur         time window length in seconds (default 0.5s)
%   - twstep        time window step in seconds (default 0.1s)
%   - mask          time to mask bad segments (default 0)
%
% OUTPUTS
%   EEG     output data
%   BCT     bad data 
%   T       threshold
%
% -------------------------------------------------------------------------

function [ EEG, BCT, T ] = eega_tRejTimeVar( EEG, varargin )

fprintf('### Rejecting based on the time variance ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.thresh = [-3 3];
P.refdata = 0;
P.refbaddata = 'none'; % 'replacebynan' / 'none' / 'zero'
P.dozscore = 0;
P.twdur = 0.5;
P.twstep = 0.1;
P.relative = 1;
P.xelectrode = 1;
P.mask = 0;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;


[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tDefBTBC: Non recognized inputs')
end

% Check the inputs
if length(P.dozscore)>1 || any(P.dozscore(:)~=0 & P.dozscore(:)~=1)
    error('eega_tDefBTBC: dozscore has to be 0 / 1')
end
if any(P.relative(:)~=0 & P.relative(:)~=1)
    error('eega_tDefBTBC: relative has to have values 0 / 1')
end
if any(P.xelectrode(:)~=0 & P.xelectrode(:)~=1)
    error('eega_tDefBTBC: xelectrode has to have values 0 / 1')
end
if length(P.twdur)>1 || length(P.twstep)>1
    error('eega_tDefBTBC: twdur and twstep have to be a number')
end
if size(P.thresh,2)~=2
    error('eega_tDefBTBC: The threshold has to be a matrix of size n x 2')
end
if size(P.thresh,1)~=length(P.relative)
    warning('eega_tDefBTBC: relative = %d for all thresholds',P.relative)
    P.relative = repmat(P.relative(1),[size(P.thresh,1) 1]);
end
if size(P.thresh,1)~=length(P.xelectrode)
    warning('eega_tDefBTBC: xelectrode = %d for all thresholds',P.relative)
    P.xelectrode = repmat(P.xelectrode(1),[size(P.thresh,1) 1]);
end

fprintf('- referenced data: %d\n',P.refdata)
fprintf('- z-score data: %d\n',P.dozscore)
fprintf('- relative threshold: %d\n',P.relative)
fprintf('- threshold per electrode: %d\n',P.xelectrode)
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
%% Time variability computation
STD = nan(nEl, ntw, nEp);
for ep = 1:nEp
    for el = 1:nEl
        d = permute(EEG.data(el,:,ep),[2 1 3]);
        dtw = d(i_t);
        if sum(~isnan(dtw))>2
            STD(el,:,ep) = nanstd(dtw,[],1);
        end
    end
end
for i=1:ntw
    idx = i_t(:,i);
    bct_tw(:,i,:) = sum(EEG.artifacts.BCT(:,idx,:),2) >= (1.0*length(idx));
end

%% ------------------------------------------------------------------------
%% Reject
nR = length(P.relative);
T = nan(nEl,2,nR);
Ru = false(size(STD));
Rl = false(size(STD));
for j=1:nR
    if P.relative(j)
        STD = log(STD/median(STD(:))); % siganl stability
        if P.xelectrode(j)
            ru_sum = 0;
            rl_sum = 0;
            for el = 1:nEl
                STD_el = STD(el,:,:);
                STD_el(bct_tw(el,:,:)) = nan;
                STD_el = STD_el(:);
                
                perc = prctile(STD_el,[25 50 75]);
                IQ   = perc(3) - perc(1);
                t_l_el  = perc(1) + P.thresh(j,1)*IQ;
                t_u_el  = perc(3) + P.thresh(j,2)*IQ;
                
                Ru(el,:,:) = Ru(el,:,:) | STD(el,:,:)>t_u_el;
                Rl(el,:,:) = Rl(el,:,:) | STD(el,:,:)<t_l_el;
                T(el,1,j) = t_u_el;
                T(el,2,j) = t_l_el;
                ru_sum = ru_sum + sum(sum(STD(el,:,:)>t_u_el,2),3);
                rl_sum = rl_sum + sum(sum(STD(el,:,:)<t_l_el,2),3);
        
                clear STD_el;
            end
        else
            perc = prctile(STD(~bct_tw),[25 50 75]);
            IQ   = perc(3) - perc(1);
            t_l  = perc(1) + P.thresh(j,1)*IQ;
            t_u  = perc(3) + P.thresh(j,2)*IQ;
            
            Ru = Ru | STD>t_u;
            Rl = Rl | STD<t_l;
            T(:,1,j) = t_l;
            T(:,2,j) = t_u;
            ru_sum = sum(STD(:)>t_u);
            rl_sum = sum(STD(:)<t_l);

        end
    else
        t_u = P.thresh(j,2);
        t_l = P.thresh(j,1);
        Ru = Ru | STD>t_u;
        Rl = Rl | STD<t_l;
        T(:,1,j) = t_l;
        T(:,2,j) = t_u;
        ru_sum = sum(STD(:)>t_u);
        rl_sum = sum(STD(:)<t_l);
    end
    
   
end
clear STD

%% ------------------------------------------------------------------------
%% Go back to each sample
bcts_u = false(nEl,nS,nEp);
bcts_l = false(nEl,nS,nEp);
Ti = i_t(1,:);
Tf = i_t(end,:);
for ep = 1:nEp
    for el = 1:nEl        
        ru = Ru(el,:,ep);
        rl = Rl(el,:,ep);
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
        for i=1:length(tiu)
            bcts_u(el,tiu(i):tfu(i),ep) = true;
        end
        for i=1:length(til)
            bcts_l(el,til(i):tfl(i),ep) = true;
        end
    end
end
clear Ru Rl
BCT = bcts_u  | bcts_l;

%% ------------------------------------------------------------------------
%% Mask around 
if ~isempty(P.mask) && P.mask~=0
    fprintf('- Mask around %4.3f s\n', P.mask)
    BCT = eega_maskmatrix(BCT, P.mask, EEG.srate);
end

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
%% Display the total rejected data
n = nEl*nS*nEp;
fprintf('Data rejected thresh_l %3.2f %%\n', sum(bcts_l(:))/n*100)
fprintf('Data rejected thresh_u %3.2f %%\n', sum(bcts_u(:))/n*100)
fprintf('Total data rejected %3.2f \n', sum(BCT(:))/n*100)
fprintf('\n' )



end

