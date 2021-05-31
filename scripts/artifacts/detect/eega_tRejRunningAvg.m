% -------------------------------------------------------------------------
% Function that performs a weighted running average algorithm to detect
% artifacts. Periods with a too big fast running average or with a too big
% difference between a fast and slow runnig average are discarted
%
% INPUTS
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - thresh_fa     fast running average threshold (default 2)
%   - thresh_da     difference running average threshold (default 2)
%   - refdata       referenced average the data before (1) or not (0) (default 0)
%   - refbaddata    how to teat bad data when reference average ('replacebynan' / 'none' / 'zero', default 'none')
%   - dozscore      z-score the data per electrodes before (1) or not (0) (default 0)
%   - relative      appply relative (1) or absolute (0) thresholds (default 1)
%   - xelectrode    appply the threhold per electrode (1) or over all electrodes (0) (default 1)
%   - mask          time to mask bad segments (default 0)
%
% OUTPUTS
%   EEG     output data
%   BCT     bad data 
%   T       threshold
%
% -------------------------------------------------------------------------

function [ EEG, BCT, T ] = eega_tRejRunningAvg( EEG, varargin )

fprintf('### Rejecting based on the running average algorithm ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.thresh_fa = 3;
P.thresh_da = 3;
P.refdata = 0;
P.refbaddata = 'none'; % 'replacebynan' / 'none' / 'zero'
P.dozscore = 0;
P.relative = 1;
P.xelectrode = 1;
P.mask = 0.05;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejRunningAvg: Non recognized inputs')
end

% Check the inputs
if length(P.dozscore)>1 || any(P.dozscore(:)~=0 & P.dozscore(:)~=1)
    error('eega_tRejRunningAvg: dozscore has to be 0 / 1')
end
if any(P.relative(:)~=0 & P.relative(:)~=1)
    error('eega_tRejRunningAvg: relative has to have values 0 / 1')
end
if any(P.xelectrode(:)~=0 & P.xelectrode(:)~=1)
    error('eega_tRejRunningAvg: xelectrode has to have values 0 / 1')
end
if size(P.thresh_fa,2)~=1
    error('eega_tRejRunningAvg: thresh_fa has to be a matrix of size n x 1')
end
if size(P.thresh_da,2)~=1
    error('eega_tRejRunningAvg: thresh_da has to be a matrix of size n x 1')
end
if size(P.thresh_da,2)~=size(P.thresh_fa,2)
    error('eega_tRejRunningAvg: thresh_fa and thresh_da have to have the same size')
end
if size(P.thresh_fa,1)~=length(P.relative)
    warning('eega_tRejRunningAvg: relative = %d for all thresholds',P.relative)
    P.relative = repmat(P.relative(1),[size(P.thresh_fa,1) 1]);
end
if size(P.thresh_fa,1)~=length(P.xelectrode)
    warning('eega_tRejRunningAvg: xelectrode = %d for all thresholds',P.relative)
    P.xelectrode = repmat(P.xelectrode(1),[size(P.thresh_fa,1) 1]);
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
%% Running average algorithm
AvgFast = nan(nEl,nS,nEp);
AvgSlow = nan(nEl,nS,nEp);

%remove a baseline to the data
baseline = nanmean(EEG.data,2);
EEG.data = EEG.data - baseline;

%compute the running average
AvgFast(:,1,:) = 0.800*zeros(nEl,1,nEp)+0.200*EEG.data(:,1,:);
AvgSlow(:,1,:) = 0.975*zeros(nEl,1,nEp)+0.025*EEG.data(:,1,:);
for j=2:nS
    AvgFast(:,j,:) = 0.800*AvgFast(:,j-1,:)+0.200*EEG.data(:,j,:);
    AvgSlow(:,j,:) = 0.975*AvgSlow(:,j-1,:)+0.025*EEG.data(:,j,:);
end
AvgDiff = AvgFast - AvgSlow;

%re-add the baseline to the data
EEG.data = EEG.data + baseline;

clear AvgSlow

%% ------------------------------------------------------------------------
%% Reject
nR = length(P.relative);
T = nan(nEl,2,nR);
BCT_fast = zeros(nEl,nS,nEp);
BCT_diff = zeros(nEl,nS,nEp);
for j=1:nR
    if P.relative(j)
        
        if P.xelectrode(j)
              
            %remove periods already signaled as bad
            AvgFast_el = abs(AvgFast);
            AvgFast_el(EEG.artifacts.BCT) = nan;
            AvgFast_el = reshape(AvgFast_el,[nEl nS*nEp]);
            
            AvgDiff_el = abs(AvgDiff);
            AvgDiff_el(EEG.artifacts.BCT) = nan;
            AvgDiff_el = reshape(AvgDiff_el,[nEl nS*nEp]);
            
            %obtain the thresholds
            perc_fast = prctile(AvgFast_el,75,2);
            IQ   = 2*perc_fast;  % becasue I'm taking half of the distribution, which is centered in zero
            t_fa  = perc_fast + P.thresh_fa(j)*IQ;
            
            perc_diff = prctile(AvgDiff_el,75,2);
            IQ   = 2*perc_diff;  % becasue I'm taking half of the distribution, which is centered in zero
            t_da  = perc_diff + P.thresh_da(j)*IQ;
            
            
        else
            
            %remove periods already signaled as bad
            AvgFast_all = abs(AvgFast);
            AvgFast_all(EEG.artifacts.BCT) = nan;
            AvgFast_all = AvgFast_all(:);
            
            AvgDiff_all = abs(AvgDiff);
            AvgDiff_all(EEG.artifacts.BCT) = nan;
            AvgDiff_all = AvgDiff_all(:);
            
            %obtain the thresholds
            perc_fast = prctile(AvgFast_all,75);
            IQ = 2*perc_fast;
            t_fa  = perc_fast + P.thresh_fa(j)*IQ;
            t_fa = repmat(t_fa,[nEl 1 1]);
            
            perc_diff = prctile(AvgDiff_all,75);
            IQ = 2*perc_diff;
            t_da  = perc_diff + P.thresh_da(j)*IQ;
            t_da = repmat(t_da,[nEl 1 1]);
            
        end
        
        %reject
        bctfast = AvgFast > repmat(t_fa,[1 nS nEp]);
        bctdiff = AvgDiff > repmat(t_da,[1 nS nEp]);
        BCT_fast = BCT_fast | bctfast;
        BCT_diff = BCT_diff | bctdiff;
        
        %store the thresholds
        T(:,1,j) = t_fa;
        T(:,2,j) = t_da;
        
    else
        
        %reject
        bctfast = AvgFast > repmat(P.thresh_fa(j),[nEl nS nEp]);
        bctdiff = AvgDiff > repmat(P.thresh_da(j),[nEl nS nEp]);
        BCT_fast = BCT_fast | bctfast;
        BCT_diff = BCT_diff | bctdiff;
        
        %store the thresholds
        T(:,1,j) = repmat(P.thresh_fa(j),[nEl 1 1]);
        T(:,2,j) = repmat(P.thresh_da(j),[nEl 1 1]);
        
    end
    
end
clear AvgFast AvgDiff
BCT = BCT_fast | BCT_diff;

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
%% Display rejected data
n = nEl*nS*nEp;
fprintf('Data rejected thresh_fa %3.2f %%\n', sum(BCT_fast(:))/n*100 )
fprintf('Data rejected thresh_da %3.2f %%\n', sum(BCT_diff(:))/n*100 )
fprintf('Total data rejected %3.2f %%\n', sum(BCT(:))/n*100 )
fprintf('\n' )


end

