function Dout = ems_stats_ttest(D,CNDs, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.ttesttail = 'both'; %'both' / 'right' / 'left'
P.latency = [D.times(1) D.times(end)];
P.DataField = 'X';
P.FactorField = 'F';
P.FactorCnd = 'diff';
P.FactorSBJ = 'SBJ';

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_sfiltcomparisonMarti: Non recognized inputs')
end

%% ------------------------------------------------------------------------

F=D.(P.FactorField);

% indexes to the subject factor
idxFcnd=[];
i=0;
while isempty(idxFcnd) && (i<=length(F))
    i=i+1;
    if strcmp(F{i}.name,P.FactorCnd)
        idxFcnd=i;
    end
end
if isempty(idxFcnd)
    error('Condition factor not found')
end
idxFsbj=[];
i=0;
while isempty(idxFsbj) && (i<=length(F))
    i=i+1;
    if strcmp(F{i}.name,P.FactorSBJ)
        idxFsbj=i;
    end
end
if isempty(idxFsbj)
    error('Subject factor not found')
end

% get the data

idxtimes = (D.times>=P.latency(1)) & (D.times<=P.latency(2));
D.(P.DataField)     = D.(P.DataField)(:,idxtimes,:);
if isfield(D,'W')
    D.W = D.W(:,idxtimes,:);
end
D.times = D.times(idxtimes);

gA = find(strcmp(F{idxFcnd}.val{1},CNDs{1}));
idxcnd_A = F{idxFcnd}.g==gA;
D.(P.DataField)     = D.(P.DataField)(:,:,idxcnd_A);
for i=1:length(F)
    F{i}.g = F{i}.g(idxcnd_A);
end
F{idxFcnd}.val = CNDs;
F{idxFcnd}.g(F{idxFcnd}.g==gA) = 1;
D.(P.FactorField)=F;

% calculate the grand average per condition
Dgavg = eega_avgdatabyfactors(D,{P.FactorCnd},'DataField',P.DataField,'FactorField',P.FactorField,'Stats',1);
if isfield(D,'W')
    Dgavg.W = mean(Dgavg.W,3);
    [Dgavg.W, Dgavg.W_stats.sd, Dgavg.W_stats.n, Dgavg.W_stats.error] = summarystats(Dgavg.W,3);
end

% t-test for the comparisons
Cmp.p = nan(size(D.(P.DataField),2),1);
Cmp.h = nan(size(D.(P.DataField),2),1);
d1=squeeze(D.(P.DataField)(1,:,F{idxFcnd}.g==1));
for i=1:size(D.(P.DataField),2)
    [Cmp.h(i,1),Cmp.p(i,1)] = ttest(d1(i,:)',0, 'Tail',P.ttesttail);
end
psort = sort(Cmp.p);
h1 = psort<=0.05;
m = length(h1);
m0 = sum(h1);
Cmp.pFDR = Cmp.p * (m/m0);


% Collect the data
Dout.Filter      = D.Filter;
if isfield('InfoSBJ',D)
    Dout.InfoSBJ     = D.InfoSBJ;
end
if isfield(D,'chanlocs')
    Dout.chanlocs    = D.chanlocs;
end
if isfield(D,'chaninfo')
    Dout.chaninfo    = D.chaninfo;
end
Dout.srate       = D.srate;
Dout.times       = D.times;
Dout.nbchan      = size(Dgavg.(P.DataField),1);
Dout.pnts        = size(Dgavg.(P.DataField),2);
Dout.trials      = size(Dgavg.(P.DataField),3);
if isfield(Dgavg,'W')
    Dout.W           = Dgavg.W;
end
Dout.chanlocs    = D.chanlocs;
% Dout.W_stats     = Dgavg.W_stats;
Dout.(P.DataField)           = Dgavg.(P.DataField);
Dout.([P.DataField '_stats'])     = Dgavg.([P.DataField '_stats']);
Dout.F           = Dgavg.F;
Dout.ttest       = Cmp;
end
