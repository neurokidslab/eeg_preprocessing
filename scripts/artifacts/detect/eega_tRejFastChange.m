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

function [ EEG, BCT, T ] = eega_tRejFastChange( EEG, varargin )

fprintf('### Rejecting fast changes ###\n' )

%% ------------------------------------------------------------------------
%% Parameters
P.thresh = 3;
P.tmotion = 0.020;  
P.refdata = 0;
P.refbaddata = 'none'; % 'replacebynan' / 'none' / 'zero'
P.dozscore = 0;
P.relative = 1;
P.xelectrode = 1;
P.mask = 0;

P.updateBCT = 1;
P.updatesummary = 1;
P.updatealgorithm = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tRejFastChange: Non recognized inputs')
end

% Check the inputs
if length(P.dozscore)>1 || any(P.dozscore(:)~=0 & P.dozscore(:)~=1)
    error('eega_tRejFastChange: dozscore has to be 0 / 1')
end
if any(P.relative(:)~=0 & P.relative(:)~=1)
    error('eega_tRejFastChange: relative has to have values 0 / 1')
end
if any(P.xelectrode(:)~=0 & P.xelectrode(:)~=1)
    error('eega_tRejNetStation: xelectrode has to have values 0 / 1')
end
if length(P.relative)~=length(P.thresh)
    error('eega_tRejFastChange: relative has to have the same length than size(thresh,1)')
end
if length(P.xelectrode)~=length(P.thresh)
    error('eega_tRejFastChange: xelectrode has to have the same length than size(thresh,1)')
end
if length(P.tmotion)>1
    error('eega_tRejFastChange: tmotion has to be a number')
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

% number of samples in the time window
smpls_amp = round(P.tmotion*EEG.srate); 
if mod(smpls_amp,2)==1
    smpls_amp = smpls_amp-1;
end

% first derivate
dd1 = diff(EEG.data,1,2); 
dd1 = cat(2,dd1,EEG.data(:,end,:));

% mask to take the data
id = (-smpls_amp/2:smpls_amp/2-1);
id = repmat(id,[nS 1]) + repmat((1:nS)',[1 smpls_amp]);
idmask = ones(size(id));
idmask(id<1) = 0;
idmask(id>(nS-1)) = 0;
id(id<1) = 1;
id(id>(nS-1)) = nS-1;

% calculate the maximun change in each window
change = nan(nEl,nS,nEp);
for ep=1:nEp 
    for el = 1:nEl    
        dd1e = dd1(el,:,ep)';
        dd1e_ = dd1e(id) .* idmask;
        change_el = sum(dd1e_,2);
        change(el,:,ep) = change_el;
    end
end
change = abs(change);

%% ------------------------------------------------------------------------
%% Reject
nR = length(P.relative);
T = nan(nEl,nR);
Ru = false([nEl nS nEp]);

for j=1:nR
    if P.relative(j)
        if P.xelectrode(j)
            ru_sum = 0;
            for el = 1:nEl
                dd = change(el,:,:);
                dd(EEG.artifacts.BCT(el,:,:)) = nan;
                dd = dd(:);
                perc = prctile(dd,75);
                IQ   = 2*perc;  % becasue I'm taking half of the distribution, which is centered in zero
                t_u_el  = perc + P.thresh(j)*IQ;
               
                Ru(el,:,:) = Ru(el,:,:) | change(el,:,:)>t_u_el;
                T(el,j) = t_u_el;
                ru_sum = ru_sum + sum(sum(change(el,:,:)>t_u_el,2),3);
            end            
        else
            dd = change(~EEG.artifacts.BCT);          
            perc = prctile(dd,75);
            IQ   = 2*perc;
            t_u  = perc + P.thresh(j)*IQ;
            
            Ru = Ru | change>t_u;
            T(:,j) = t_u;
            ru_sum = sum(change(:)>t_u);
        end
    else
        dd = change;        
        t_u = P.thresh(j);
        Ru = Ru | dd>t_u;
        T(:,j) = t_u;
        ru_sum = sum(change(:)>t_u);
    end
    
    %% --------------------------------------------------------------------
    %% Display rejected data
    n = nEl*nS*nEp;
    fprintf('Data rejected thresh_u %3.2f %%\n', ru_sum/n*100)
    
end
BCT = Ru;
clear Ru

%% ------------------------------------------------------------------------
%% Reject around (in the length where the change was detected)

BCTin = BCT;
BCT = false(nEl,nS,nEp);

% set artifact buffer for tmsk seconds on each side of spike
art_buffer = smpls_amp/2; % time in seconds times sample rate
for ep=1:nEp
    bct = BCTin(:,:,ep);
    for el = 1:nEl
        bad = bct(el,:)';
        bad_idx = find(bad);
        % Eliminate time points before or after motion artifacts
        if ~isempty(bad_idx)
            bad_idx = repmat(bad_idx, 1, 2*art_buffer+1)+repmat(-art_buffer:art_buffer,length(bad_idx), 1);
            bad_idx = unique(bad_idx(:));
            bad_idx = bad_idx( (bad_idx>0) & (bad_idx<=nS) );
        end
        
        % Reject
        bct(el,bad_idx) = 1;
    end
    BCT(:,:,ep) = bct & ~BCTin(:,:,ep);
end
BCT = BCT | BCTin; 

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
fprintf('Total data rejected %3.2f %%\n', sum(BCT(:))/n*100 )
fprintf('\n' )


end



