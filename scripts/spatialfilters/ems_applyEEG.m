function [EEGsf, T0, Tx] = ems_applyEEG(EEG, FactorCnd, CNDs, varargin)

fprintf('### EMS ###\n')

%% Default parameters
P.DataField     = 'data';
P.FactorField   = 'F';
P.BadData       = 'none';   % 'replacebymean' / 'none'
P.method        = 'ems_2conddiff_oneout'; % 'ems_2cond_diff' / 'ems_2cond_perm' / 'ems_2cond_oneout' / 'ems_2conddiff_perm' / 'ems_2conddiff_oneout' / 'ems_1cond_perm' / 'ems_1cond_oneout'
P.noiseCov      = 0;
P.T0            = [];
P.Tx            = [];
P.Np            = 100;
P.N2bias        = 0.7;
P.times         = [];
P.channels      = [];
P.trimmean      = [];

%% Optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('ems_applyEEG: Non recognized inputs')
end

%% Check the input
if any(strcmp(P.method,{'ems_2cond_perm' ; 'ems_2cond_oneout' ; 'ems_2conddiff_perm' ; 'ems_2conddiff_oneout'}))
    if length(CNDs)~=2
        error('ems_applyEEG: two conditions are requiered')
    end
elseif any(strcmp(P.method,{'ems_1cond_perm'; 'ems_1cond_oneout'}))
    if length(CNDs)~=1
        error('ems_applyEEG: one condition is requiered')
    end
end
if ~any(strcmp(P.method,{'ems_2cond_diff'; 'ems_2cond_perm' ; 'ems_2cond_oneout' ; 'ems_2conddiff_perm' ; 'ems_2conddiff_oneout';'ems_1cond_perm'; 'ems_1cond_oneout'}))
    error('ems_applyEEG: not recognized method')
end

if ischar(FactorCnd)
    FactorCnd = {FactorCnd};
end

if ischar(P.noiseCov)
    noiseCov = EEG.(P.noiseCov);
elseif P.noiseCov==1 || P.noiseCov==0
    noiseCov = P.noiseCov;
end

T0 = P.T0;
Tx = P.Tx;
sz = size(EEG.(P.DataField));

if isempty(P.times)
    P.times = [EEG.times(1) EEG.times(end)];
end
P.times = (EEG.times>=P.times(1)) & (EEG.times<=P.times(end));
if isempty(P.channels)
    P.channels = (1:sz(1));
end

%% Remove bad data
if ~isempty(EEG.(P.DataField))
   [~, Y ] = eega_rmvbaddata(EEG, 'BadData', P.BadData, 'DataField', P.DataField);
    % Y=EEG.(P.DataField);
    F = EEG.(P.FactorField);
else
    Y=[];
    F=[];
end

%% Spatial filter

if length(CNDs)==1
    theFilt = CNDs{1};
elseif length(CNDs)==2
    theFilt = sprintf('%s-%s',CNDs{1},CNDs{2});
end

if ~isempty(Y)
    
    
    %% Build the table defying the conditions
    TableCNDs = cnd_buildtablecond(FactorCnd, EEG.(P.FactorField));
    
    %% Get the conditions
    F = EEG.(P.FactorField);
    [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(TableCNDs, F);
    [~,CNDsIND] = ismember(CNDs, theCND);
    
    %% Spatial filter
    
    [~,CNDsIND] = ismember(CNDs, theCND);
    Favg{1}.name = 'diff';
    Favg{1}.val = {theFilt};
    Favg{1}.g = 1;

    % build a functioon handle for the spatial filter to apply
    fh = str2func(P.method);
    fprintf('Filter: %s \n', theFilt)
    
    % trials in each condition to compare
    COND = bsxfun(@eq, condition(:), CNDsIND(:)');
    COND = logical(COND);
    if isempty(T0), T0=any(COND,2); end
    if isempty(Tx), Tx=any(COND,2); end
    
%     if length(CNDs)==1
%         g1 = condition==CNDsIND(1);
%         if ~any(g1)
%             error('Conditions not found')
%         end
%         COND(:,1)=logical(g1);
%     elseif length(CNDs)==2
%         g1 = condition==CNDsIND(1);
%         g2 = condition==CNDsIND(2);
%         if ~any(g1) || ~any(g2)
%             error('Conditions not found')
%         end
%         COND(:,1)=logical(g1);
%         COND(:,2)=logical(g2);
%     end
%     if isempty(T0), T0=any(COND,2); end
%     if isempty(Tx), Tx=any(COND,2); end
    
    
%     
%     
%     % find the trials in the condtions
%     idxFcnd=[];
%     i=0;
%     while isempty(idxFcnd) && (i<=length(F))
%         i=i+1;
%         if strcmp(F{i}.name,FactorCnd)
%             idxFcnd=i;
%         end
%     end
%     if isempty(idxFcnd)
%         error('Condition factor not found')
%     end
%     
%     if length(CNDs)==1
%         g1 = find(strcmp(F{idxFcnd}.val,CNDs{1}));
%         if isempty(g1)
%             error('Conditions not found')
%         end
%         COND(:,1)=F{idxFcnd}.g==g1;
%         COND=logical(COND);
%     elseif length(CNDs)==2
%         g1 = find(strcmp(F{idxFcnd}.val,CNDs{1}));
%         g2 = find(strcmp(F{idxFcnd}.val,CNDs{2}));
%         if isempty(g1) || isempty(g2)
%             error('Conditions not found')
%         end
%         COND(:,1)=F{idxFcnd}.g==g1;
%         COND(:,2)=F{idxFcnd}.g==g2;
%         COND=logical(COND);
%     end
%     if isempty(T0), T0=any(COND,2); end
%     if isempty(Tx), Tx=any(COND,2); end
    
    % apply the function
    if all(sum(COND,1)>2)
        if any(strcmp(P.method,{'ems_2cond_perm' 'ems_2conddiff_perm' 'ems_1cond_perm'})) 
            [ x, w ] = fh(Y(P.channels,P.times,:), COND, T0, Tx, noiseCov, P.N2bias, P.Np);
        elseif any(strcmp(P.method,{'ems_2cond_oneout' 'ems_2conddiff_oneout' 'ems_1cond_oneout'}))
            [ x, w ] = fh(Y(P.channels,P.times,:), COND, T0, Tx, noiseCov);
        elseif any(strcmp(P.method,{'ems_2cond_diff'}))  % Sebastian function 
            [ x, w ] = ems_2cond_diff(Y(P.channels,P.times,:),COND,[],[],noiseCov);
            x = permute(x,[3 2 1]);
        end
        x1 = x(:,:,COND(Tx,1)); 
        x2 = x(:,:,COND(Tx,2)); 
        if isempty(P.trimmean) || (P.trimmean==0) 
            x1 = mean(x1,3);
            x2 = mean(x2,3);
        else
            x1 = trimmean(x1,P.trimmean,3);
            x2 = trimmean(x2,P.trimmean,3);
        end
        X = nan( 1, sz(2) );
        W = nan( sz(1), sz(2) );
        X(1, P.times, :) = x1-x2;
        W(P.channels, P.times) = w;
        for idf=1:length(F)
            F{idf}.g = F{idf}.g(Tx);
        end
    else
        warning('Not enougth good data')
        X = [];
        W = [];
        F = [];
    end
else
    X = [];
    W = [];
    Favg = [];
end


%% ------------------------------------------------------------------------
%% Fill the structure

if isfield(EEG, 'filename')
    [~,nameSubj,~] = fileparts(EEG.filename);
    EEGsf.subject = nameSubj;
end
EEGsf.Filter        = theFilt;
if isfield(EEG, 'FilterCNDs')
    EEGsf.FilterCNDs    = CNDs;
end
if isfield(EEG, 'chanlocs')
    EEGsf.chanlocs      = EEG.chanlocs;
end
if isfield(EEG, 'chaninfo')
    EEGsf.chaninfo      = EEG.chaninfo;
end
if isfield(EEG, 'srate')
    EEGsf.srate         = EEG.srate;
end
if isfield(EEG, 'times')
    EEGsf.times         = EEG.times;
end
if isfield(EEG, 'freq')
    EEGsf.freq         = EEG.freq;
end
EEGsf.nbchan        = size(X,1);
EEGsf.pnts          = size(X,2);
EEGsf.trials        = size(X,3);
EEGsf.W             = W;
EEGsf.X             = X;
EEGsf.F             = Favg;


end
