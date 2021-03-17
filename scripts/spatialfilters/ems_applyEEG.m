function [EEGsf, T0, Tx] = ems_applyEEG(EEG, FactorCnd, CNDs, varargin)

fprintf('### EMS ###\n\n')

st = dbstack;
fnnamestr = st.name;

%% Parameters

% Default parameters
P.DataField     = 'data';
P.FactorField   = 'F';
P.kPCA          = [];
P.method        = 'ems_2conddiff_perm_multitrials'; % 'ems_2cond_perm' | 'ems_2cond_oneout' | 'ems_2conddiff_perm' | 'ems_2conddiff_oneout' | 'ems_1cond_perm' | 'ems_1cond_oneout'
P.noiseCov      = 0;
P.timeCov0      = [];
P.T0            = [];
P.Tx            = [];
P.Np            = 100;
P.N2bias        = 0.75;
P.times         = [];
P.channels      = [];
P.ptrimmean     = 0;
P.keepdatafield = 0;

% Get the optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('%s: Non recognized inputs',fnnamestr)
end

% Check the input
if any(strcmp(P.method,{'ems_2cond_perm' ; 'ems_2cond_oneout' ; 'ems_2conddiff_perm' ; 'ems_2conddiff_oneout'}))
    if length(CNDs)~=2
        error('%s: two conditions are requiered',fnnamestr)
    end
elseif any(strcmp(P.method,{'ems_1cond_perm'; 'ems_1cond_oneout'}))
    if length(CNDs)~=1
        error('%s: one condition is requiered',fnnamestr)
    end
end
% if ~any(strcmp(P.method,{'ems_2cond_perm' ; 'ems_2cond_oneout' ; 'ems_2conddiff_perm' ; 'ems_2conddiff_oneout';'ems_1cond_perm'; 'ems_1cond_oneout'}))
%     error('%s: not recognized method',fnnamestr)
% end
if ischar(FactorCnd)
    FactorCnd = {FactorCnd};
end

if ischar(P.noiseCov)
    noiseCov = EEG.(P.noiseCov);
elseif P.noiseCov==1 || P.noiseCov==0
    noiseCov = P.noiseCov;
end
if isempty(P.timeCov0)
    timeCov0 = [];
elseif length(P.timeCov0)==2
    timeCov0 = (EEG.times>P.timeCov0(1)) & (EEG.times<P.timeCov0(2));
elseif (length(P.timeCov0)~=size(EEG.data,2))
    error('Bad time to compute the covariance of the noise') 
end

if isempty(P.times)
    P.times = [EEG.times(1) EEG.times(end)];
end
P.times = (EEG.times>=P.times(1)) & (EEG.times<=P.times(end));
if isempty(P.channels)
    P.channels = (1:size(EEG.(P.DataField),1));
end


%% Spatial filter
T0 = P.T0;
Tx = P.Tx;
Y = EEG.(P.DataField)(P.channels,P.times,:); 
[Ne,Ns,Nt] = size(Y);
if length(CNDs)==1
    theFilt = CNDs{1};
elseif length(CNDs)==2
    theFilt = sprintf('%s-%s',CNDs{1},CNDs{2});
end

if ~isempty(Y)  % if there is good data go on
    
    % PCA of the data
    if ~isempty(P.kPCA)
        fprintf('Reducing dimensions (PCA)\n')
        for i=1:Ne  % normalize the data
            d = Y(i,:,:);
            d = d(:)/norm(d(:));
            Y(i,:,:) = reshape(d,[1 Ns Nt]);
        end
        dd = reshape(Y,Ne,[])';
        dd = bsxfun(@minus,dd,mean(dd,1));
        c0 =(1/Nt).*(dd'*dd);
        c0 = c0/norm(c0);
        [V, S] = eig(c0);
        V = real(V);
        S = real(S);
        S = abs(S);
        [L, idx] = sort(diag(S), 'descend') ;
        V = V(:,idx);
        exppca1 = cumsum(L)/sum(L)*100;
        if P.kPCA<1
            ki = find(exppca1>=(P.kPCA*100),1);
            exppca1 = exppca1(ki);
        else
            ki = P.kPCA;
            if abs((exppca1(ki)-100))<1e-3
                ki = find(exppca1>=(100-1e-3),1);
            end
            exppca1 = exppca1(ki);
        end
        fprintf('PCA: the explained variance by the %i first components is %4.2f\n',ki,exppca1)
        Y = dd * V;
        Y = Y(:,1:ki);
%         % prooject back to snesors space
%         Y = Y * pinv(V(:,1:ki));
%         Y = Y';
        % reshape
%         Y = reshape(Y,[Ne, Ns, Nt ]);
        Y = reshape(Y',[ki, Ns, Nt ]);
        fprintf('\n');
    end

    % Get the conditions
    F = EEG.(P.FactorField);
    TableCNDs = cnd_buildtablecond(FactorCnd, EEG.(P.FactorField));
    [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(TableCNDs, F);
        
    [~,CNDsIND] = ismember(CNDs, theCND);
    
    % average factor
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
    
    % apply the function to the trials in that conditions
    [Ne, Ns, Nt] = size(Y);
    if all(sum(COND,1)>2)
        if any(strcmp(P.method,{'ems_2cond_oneout'; 'ems_2conddiff_oneout'}))
            [X, W] = fh(Y, COND, T0, Tx, noiseCov);
            X = mean(X(:,:,COND(:,1)),3) -  mean(X(:,:,COND(:,2)),3);
        elseif any(strcmp(P.method,{'ems_2cond_perm' ; 'ems_2conddiff_perm' }))
            [X, W] = fh(Y, COND, T0, Tx, noiseCov, P.N2bias, P.Np);  
            X = mean(X(:,:,COND(:,1)),3) -  mean(X(:,:,COND(:,2)),3);
        elseif any(strcmp(P.method,{'ems_2conddiff_perm_multitrials'}))
            [X, W] = fh(Y, COND, T0, Tx, noiseCov, timeCov0, P.N2bias, P.Np, P.ptrimmean);
        end
%         X = nan( 1, Ns, size(x,3) );
%         W = nan( Ne, Ns );
%         X(:,P.times,:) = x;
%         W(:,P.times) = w;
%         X = sum(x,1);
%         W = w;
    else
        warning('Not enougth good data')
        X = [];
        W = [];
        Favg = [];
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
if isfield(EEG,'InfoSBJ')
    EEGsf.InfoSBJ = EEG.InfoSBJ;
end
if isfield(EEG, 'chanlocs')
    EEGsf.chanlocs      = EEG.chanlocs(P.channels);
end
if isfield(EEG, 'srate')
    EEGsf.srate         = EEG.srate;
end
if isfield(EEG, 'times')
    EEGsf.times         = EEG.times(P.times);
end

EEGsf.nbchan        = size(X,1);
EEGsf.pnts          = size(X,2);
EEGsf.trials        = size(X,3);
EEGsf.W             = W;
EEGsf.X             = X;
EEGsf.F             = Favg;

if P.keepdatafield
    EEGsf.(P.DataField) = Y;
    EEGsf.(['F' P.DataField]) = EEG.(P.FactorField);
end


end
