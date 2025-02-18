% -------------------------------------------------------------------------
% Function that reject data when the maximun change in a given time window
% is bigger than a threshold 
%
% INPUTS
% EEG   EEG structure
%
% OPTIONAL INPUTS
%   - thresh        upper threshold (default 2)
%   - tmotion       time window where fast chenges are detected (default 0.020s)
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

function [ EEG, BCT, T ] = eega_tRejMaxChangeEpoch( EEG, varargin )

fprintf('### Rejecting maximun changes ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.thresh = 3;
P.refdata = 0;
P.refbaddata = 'none'; % 'replacebynan' / 'none' / 'zero'
P.dozscore = 0;
P.relative = 1;
P.xelectrode = 1;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejMaxChange: Non recognized inputs')
end

% Check the inputs
if length(P.dozscore)>1 || any(P.dozscore(:)~=0 & P.dozscore(:)~=1)
    error('eega_tRejMaxChange: dozscore has to be 0 / 1')
end
if any(P.relative(:)~=0 & P.relative(:)~=1)
    error('eega_tRejMaxChange: relative has to have values 0 / 1')
end
if any(P.xelectrode(:)~=0 & P.xelectrode(:)~=1)
    error('eega_tRejMaxChange: xelectrode has to have values 0 / 1')
end
if length(P.relative)~=length(P.thresh)
    error('eega_tRejMaxChange: relative has to have the same length than size(thresh,1)')
end
if length(P.xelectrode)~=length(P.thresh)
    error('eega_tRejMaxChange: xelectrode has to have the same length than size(thresh,1)')
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
%% Calculate the max change

change = max(EEG.data,[],2)-min(EEG.data,[],2);
change = abs(change);
change = log(change);

%% ------------------------------------------------------------------------
%% Reject
nR = length(P.relative);
T = nan(nEl,nR);
Ru = false([nEl 1 nEp]);

for j=1:nR
    if P.relative(j)
        if P.xelectrode(j)
            ru_sum = 0;
            for el = 1:nEl
                dd = change(el,:,:);
                dd(EEG.artifacts.BC(el,:,:)) = nan;
                dd = dd(:);
                perc = prctile(dd,[25 50 75]);
                IQ   = perc(3)-perc(1);  
                t_u_el  = perc(3) + P.thresh(j)*IQ;
               
                Ru(el,:,:) = Ru(el,:,:) | change(el,:,:)>t_u_el;
                T(el,j) = t_u_el;
                ru_sum = ru_sum + sum(sum(change(el,:,:)>t_u_el,2),3);
            end            
        else
            dd = change;    
            dd(EEG.artifacts.BC) = nan;
            dd = dd(:);
            perc = prctile(dd,[25 50 75]);
            IQ   = perc(3)-perc(1);
            t_u  = perc(3) + P.thresh(j)*IQ;
            
            Ru = Ru | change>t_u;
            T(:,j) = t_u;
            ru_sum = sum(change(:)>t_u);
        end
    else
        dd = change;        
        t_u = log(P.thresh(j));
        Ru = Ru | dd>t_u;
        T(:,j) = t_u;
        ru_sum = sum(change(:)>t_u);
    end
    
    %% --------------------------------------------------------------------
    %% Display rejected data
    n = nEl*nEp;
    fprintf('Data rejected thresh_u %3.2f %%\n', ru_sum/n*100)
    
end
BCT = repmat(Ru,[1 nS 1]);
clear Ru


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
fprintf('Total data rejected %3.2f %%\n', sum(BCT(:))/n*100 )
fprintf('\n' )


end



