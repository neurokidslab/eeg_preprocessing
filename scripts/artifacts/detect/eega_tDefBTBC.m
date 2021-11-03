% -------------------------------------------------------------------------
% This function defines where there are artifacts affecting all the
% electrodes and which are bad channels
%
% INPUTS
% EEG   EEG structure
% DefBT
% DefBC
%
% OPTIONAL INPUTS
%   - maskTime      time to mask bad times (s)
%   - minGoodTime   Periods with non artifact shorter than this (s) are
%                   marked as having artifacts
%   - minBadTime    Periods with artifacts shorter than this (s) are not
%                   included
%   - keeppre       keep previuos values
%   - plot          plot the rejection
%
% OUTPUTS
% EEG     output data
% BTnew     bad times, logical indexes (1 x samples x epochs)
% BCnew     bad electrodes, logial indexes (electrodes x 1 x epochs)
%
% -------------------------------------------------------------------------

function [ EEG, BTnew, BCnew ] = eega_tDefBTBC( EEG, DefBT, DefBC, DefBCAll, varargin )

fprintf('### Identifying bad samples and channels ###\n' )


%% ------------------------------------------------------------------------
%% Parameters
P.keeppre       = 0;
P.plot          = 0;
P.minBadTime    = 0;
P.maskTime      = 0;
P.minGoodTime   = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tDefBTBC: Non recognized inputs')
end

ncycle = unique([length(DefBT) length(DefBC) length(DefBC)]);
if length(ncycle)~=1
    error('eega_tDefBTBC: all the thresholds have to have the same length')
end

%% ------------------------------------------------------------------------
%% Get the data from the EEG structure
nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

%% ------------------------------------------------------------------------
%% Define artifacts

% If requiered keep the previous information
if isnumeric( P.keeppre )
    if (P.keeppre ==1)
        if isfield(EEG.artifacts,'BC')
            BCold = EEG.artifacts.BC;
        else
            BCold = false(nEl,1,nEp);
        end
        if isfield(EEG.artifacts,'BT')
            BTold = EEG.artifacts.BT;
        else
            BTold = false(1,nS,nEp);
        end
    else
        BCold = false(nEl,1,nEp);
        BTold = false(1,nS,nEp);
    end
elseif ischar( P.keeppre )
    switch P.keeppre
        case 'BT'
            if isfield(EEG.artifacts,'BT')
                BTold = EEG.artifacts.BT;
                BCold = false(nEl,1,nEp);
            else
                BTold = false(1,nS,nEp);
                BCold = false(nEl,1,nEp);
            end
        case 'BC'
            if isfield(EEG.artifacts,'BC')
                BCold = EEG.artifacts.BC;
                BTold = false(1,nS,nEp);
            else
                BCold = false(nEl,1,nEp);
                BTold = false(1,nS,nEp);
            end
    end
end

BCnew = BCold;
BTnew = BTold;
BCBTnew = zeros(nEl,nS,nEp);
BCBTnew(repmat(BCnew,[1 nS 1])) = 1;
BCBTnew(repmat(BTnew,[nEl 1 1])) = 1;
BCall = all(BCnew,3);
if isfield(EEG.artifacts,'BCmanual')
    BCall(EEG.artifacts.BCmanual) = true;
    BCnew(EEG.artifacts.BCmanual,:,:) = true;
end


for icyc=1:ncycle
    
    % DEFINE BAD SAMPLES BASED ON ABSOLUTE THRESHOLD
    %number of bad channels per sample
    thrsBadCh = DefBT(icyc);
    bct = EEG.artifacts.BCT;
    bct(repmat(BCnew,[1 nS 1])) = 0;
    nbadCh = sum(bct,1);
    pbadCh = nbadCh ./ repmat(sum(~BCnew,1),[1 nS 1]);
    
    % DEFINE BAD CHANNELS DURING THE WHOLE RECORDING ON ABSOLUTE THRESHOLD
    %number of bad samples per channel
    thrsBadSAll = DefBCAll(icyc);
    bct = EEG.artifacts.BCT;
    bct(repmat(BTnew,[nEl 1 1])) = 0;
    bct = reshape(bct,[nEl nS*nEp]);
    nbadSAll = sum(bct,2);
    pbadSAll = nbadSAll/ sum(~reshape(BTnew,[1 nS*nEp]),2);
    
    % DEFINE BAD CHANNELS PER EPOCH ON ABSOLUTE THRESHOLD
    % number of bad samples per channel
    thrsBadS = DefBC(icyc);
    bct = EEG.artifacts.BCT;
    bct(repmat(BTnew,[nEl 1 1])) = 0;
    nbadS = sum(bct,2);
    pbadS = nbadS ./ repmat(sum(~BTnew,2), [nEl 1 1]);
    
    % REJECT
    BTnew = BTnew | (pbadCh>thrsBadCh);
    BCall = BCall | (pbadSAll>thrsBadSAll);
    BCnew = BCnew | repmat(BCall, [1 1 nEp]);
    BCnew = BCnew | (pbadS>thrsBadS);
    
    % TEST IF THE DEFINITION HAS CHANGED
    bcbt = zeros(nEl,nS,nEp);
    bcbt(repmat(BCnew,[1 nS 1])) = 1;
    bcbt(repmat(BTnew,[nEl 1 1])) = 1;
    T = bcbt~=BCBTnew;
    BCBTnew = bcbt;
    fprintf('Cycle %d, new rejected data %3.2f %%\n',icyc, sum(T(:))/length(T(:))*100)
end

%% ------------------------------------------------------------------------
%% Update
BT = BTold | BTnew;
BC = BCold | BCnew;

%% ------------------------------------------------------------------------
%% Samples rejection or not based on the rejected data

% Remove too short artifacts
if ~isempty(P.minBadTime)
    P.minBadTime = P.minBadTime * EEG.srate;
    if P.minBadTime~=0
        for e=1:nEp
            isbadsmpl = BT(1,:,e);
            tbad = isbadsmpl(1,:);
            badi = [tbad(1) ( ( tbad(2:end) - tbad(1:end-1) ) == 1 ) ];
            badf = [ ( ( tbad(1:end-1) - tbad(2:end) ) == 1 ) tbad(end) ];
            badi = find( badi );
            badf = find( badf );
            
            lim = [badi' badf'];
            dur = lim(:,2) - lim(:,1);
            idx = find(dur < P.minBadTime);
            for i=1:numel(idx)
                isbadsmpl(1,badi(idx(i)):badf(idx(i))) = 0;
            end
            BT(1,:,e) = logical(isbadsmpl);
        end
    end
end

% Mask around
if ~isempty(P.maskTime) && P.maskTime~=0
    P.maskTime = P.maskTime * EEG.srate;
    art_buffer=round(P.maskTime);
    for ep=1:nEp
        bt = BT(1,:,ep);
        bad = bt(:);
        bad_idx = find(bad);
        
        % Eliminate time points before or after motion artifacts
        if ~isempty(bad_idx)
            bad_idx = repmat(bad_idx, 1, 2*art_buffer+1)+repmat(-art_buffer:art_buffer,length(bad_idx), 1);
            bad_idx = unique(bad_idx(:));
            bad_idx = bad_idx( (bad_idx>0) & (bad_idx<=nS) );
        end
        bt(bad_idx) = 1;
        BT(1,:,ep) = bt;
    end
end

% Remove too short periods with non artifacts
if ~isempty(P.minGoodTime)
    P.minGoodTime = P.minGoodTime * EEG.srate;
    if P.minGoodTime~=0 && P.minGoodTime<(size(EEG.data,2)-2)
        for e=1:nEp
            isbadsmpl = BT(1,:,e);
            tgood = ~isbadsmpl(1,:);
            goodi = [tgood(1) ( ( tgood(2:end) - tgood(1:end-1) ) == 1 ) ];
            goodf = [ ( ( tgood(1:end-1) - tgood(2:end) ) == 1 ) tgood(end) ];
            goodi = find( goodi );
            goodf = find( goodf );
            
            lim = [goodi' goodf'];
            dur = lim(:,2) - lim(:,1);
            idx = find(dur < P.minGoodTime);
            for i=1:numel(idx)
                isbadsmpl(1,goodi(idx(i)):goodf(idx(i))) = 1;
            end
            BT(1,:,e) = BT(1,:,e) | logical(isbadsmpl);
        end
    end
end

%% ------------------------------------------------------------------------
%% Display rejected data
d = 1*nS*nEp;
newd = sum(BT(:) & ~BTold(:));
alld = sum(BT(:));
fprintf('Total new bad times ______________ %3.2f %%\n', newd/d*100 )
fprintf('Total bad times __________________ %3.2f %%\n', alld/d*100 )

d = nEl*1*nEp;
newd = sum(BC(:) & ~BCold(:));
alld = sum(BC(:));
fprintf('Total new bad channels per epoch _ %3.2f %%\n', newd/d*100 )
fprintf('Total bad channels per epoch _____ %3.2f %%\n', alld/d*100 )

newd = sum(all(BC,3) & all(~BCold,3));
alld = sum(all(BC,3));
fprintf('Total new bad channels ___________ %d\n', newd )
fprintf('Total bad channels _______________ %d\n', alld )
% fprintf('Total manually rejected channels _ %d\n', length(EEG.artifacts.BCmanual) )

%% ------------------------------------------------------------------------
%% Update the rejection matrix
EEG.artifacts.BT = BT;
EEG.artifacts.BC = BC;
% EEG.reject.rejmanualE = permute(EEG.artifacts.BC,[1 3 2]);
fprintf('\n')
if exist('eega_summarypp','file')==2
    EEG = eega_summarypp(EEG);
end

%% ------------------------------------------------------------------------
%% Plot the rejection matrix
if P.plot
    eega_plot_artifacts(EEG);
end
end
