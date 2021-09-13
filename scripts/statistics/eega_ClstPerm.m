% This function perfoms the fieldtrip cluster permutation analysis

function [Avg, GAvg, stat, CNDsIDX] = eega_ClstPerm(EEGALL, varargin)

fprintf('### Cluster Based Permutation Analysis ### \n')

% check that Fieldtrip is in the path
if exist('ft_timelockstatistics','file')~=2
    error('The FieldTrip toolbox is not in the path. Add it to the path (https://www.fieldtriptoolbox.org/)')
end

%% ------------------------------------------------------------------------
%% Parameters

P.DataField         = 'data';
P.TimeField         = 'times';
P.FreqField         = 'freqs';
P.FField            = 'F';
P.method            = 'montecarlo';
P.latency           = 'all'; % time where the cluster permutation is run: interval or 'all'
P.frequency         = 'all';
P.idxelec           = [];
P.idxsbj            = [];
P.doabsolute        = 0;
P.conditions        = [];
P.sbjfactor         = 'SBJ';
P.cndfactor         = 'Condition';
P.design            = 'within';      % design type within or between
P.numrandomization  = 1000;
P.clusteralpha      = 0.05; % alpha level of the sample-specific test statistic that will be used for thresholding
P.clusterstatistic  = 'maxsum'; % test statistic that will be evaluated under the permutation distribution.
P.correctm          = 'cluster'; %string, apply multiple-comparison correction, 'no', 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
P.alpha             = 0.025; % alpha level of the permutation test
P.tail              = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
P.correcttail       = 'no'; %string, correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
P.clustertail       = 0;
P.avgoverchan       = 'no'; %'yes' or 'no' (default = 'no')
P.avgovertime       = 'no'; %'yes' or 'no' (default = 'no')
P.neighbourdist     = [];
P.neighbourmethod   = 'triangulation'; % 'distance', 'triangulation' or 'template'
P.chancmbneighbstructmat = [];
P.minnbchan         = 0; % minimum number of neighborhood channels that is required for a selected sample to be included in the clustering algorithm (default=0).

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_ClstPerm: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Check some usefull stuff
electrodes = size(EEGALL.(P.DataField),1);
F = EEGALL.(P.FField);
if isempty(EEGALL.chanlocs) && electrodes>1
    error('eega_ClstPerm: No channels localization')
end

%% ------------------------------------------------------------------------
%% Select the condtion, time and electrodes
isbjF = [];
j=1;
while isempty(isbjF)
    if strcmp(P.sbjfactor,F{j}.name)
        isbjF = j;
    end
    j=j+1;
end

icndF = [];
j=1;
while isempty(icndF)
    if strcmp(P.cndfactor,F{j}.name)
        icndF = j;
    end
    j=j+1;
end

if isempty(P.idxelec)
    P.idxelec = (1:electrodes);
end
if isempty(P.idxsbj)
    P.idxsbj = (1:length(F{isbjF}.val));
end
if isempty(P.conditions)
    P.conditions = F{icndF}.val;
end

nSBJ = length(P.idxsbj);
nCND = length(P.conditions);
nSBJCND = nan(nCND,1);

Avg = cell(1, nCND);
CNDsIDX = nan(1, nCND);
for cnd = 1:nCND
    thecnd = find( strcmp(P.conditions{cnd},F{icndF}.val) );
    idxcnd = F{icndF}.g == thecnd;
    CNDsIDX(cnd) = thecnd;
    sbjcounter = 0;
    for s=1:length(P.idxsbj)
        thesbj = s;
        idxsbj = F{isbjF}.g == thesbj;
        
        idxepoch = idxcnd(:) & idxsbj(:);
        if any(idxepoch)
            sbjcounter = sbjcounter+1;
            Avg{cnd}{sbjcounter} = eega_egglab2ft(EEGALL,'idxepoch',idxepoch,...
                'DataField',P.DataField,'TimeField',P.TimeField,'FreqField',P.FreqField);
        end
        
    end
    nSBJCND(cnd) = sbjcounter;
end

%% ------------------------------------------------------------------------
%% Neightbours
if electrodes<2
%     layout.label{1} = EEGALL.chanlocs(1).labels;
%     neighbours(1).label = EEGALL.chanlocs(1).labels;
%     neighbours(1).neighblabel = {''};
    layout.label{1} = 'El';
    neighbours(1).label = 'El';
    neighbours(1).neighblabel = {''};
else
    X = [EEGALL.chanlocs.X];
    Y = [EEGALL.chanlocs.Y];
    Z = [EEGALL.chanlocs.Z];
    ChLabels = cell(electrodes,1);
    Chanpos = cat(2,X(:),Y(:),Z(:));
    for el=1:electrodes
        ChLabels{el} = EEGALL.chanlocs(el).labels;
    end
    layout.label = ChLabels;
    layout.chanpos = Chanpos;
    neighbours = eega_neighborchannels(layout,'method',P.neighbourmethod,'neighbourdist',P.neighbourdist);
end

% tmpdata = Avg{1}{1};
% tmpcfg = [];
% tmpcfg.feedback = 'yes';
% tmpcfg.method = P.distmethod; 
% tmpcfg.neighbourdist = P.distneighbour;
% if electrodes<2
%     neighbours(1).label = ChLabels{1};
%     neighbours(1).neighblabel = '';
% else
%     neighbours = ft_prepare_neighbours(tmpcfg, tmpdata);
% end

%% ------------------------------------------------------------------------
%% Build the cgf structure

cfg = [];

% General
cfg.channel             = layout.label;
cfg.latency             = P.latency; 
if isfield(Avg{1}{1},'freq')
    cfg.frequency       = P.frequency; 
end
cfg.method              = P.method; % use the Monte Carlo Method to calculate the significance probability
cfg.correctm            = P.correctm;
cfg.parameter           = 'avg';
cfg.clusteralpha        = P.clusteralpha; 
cfg.clusterstatistic    = P.clusterstatistic; 
cfg.alpha               = P.alpha;
cfg.tail                = P.tail;
cfg.correcttail         = P.correcttail; 
cfg.clustertail         = P.clustertail;
cfg.numrandomization    = P.numrandomization; 
cfg.avgoverchan         = P.avgoverchan;
cfg.avgovertime         = P.avgovertime;

if nCND==2
    if strcmp(P.design,'between')
        
        design1 = [ones(1,nSBJCND(1)) 2*ones(1,nSBJCND(2))];
    
        cfg.statistic       = 'indepsamplesT'; % use the independent samples T-statistic as a measure to evaluate the effect at the sample level
        cfg.design          = design1;
        cfg.ivar            = 1; % the raw of the design matrix containing the independent variable
        
    elseif strcmp(P.design,'within')
        
        design1 = [ones(1,nSBJCND(1)) 2*ones(1,nSBJCND(2))];
        design2 = repmat((1:nSBJ), [1 nCND]);
        
        cfg.statistic       = 'depsamplesT';
        cfg.design          = [design2; design1];
        cfg.uvar            = 1; % raw of the design matrix containing the unit variable
        cfg.ivar            = 2; % raw of the design matrix containing the independent variable
        
    end
elseif nCND>2
    design1 = [];
    for i=1:nCND
        design1 = [design1 i*ones(1,nSBJ)];
    end
    if strcmp(P.design,'between')
        
        cfg.statistic       = 'indepsamplesF'; % use the independent samples F-statistic as a measure to evaluate the effect at the sample level
        cfg.design          = design1;
        cfg.ivar            = 1; % the raw of the design matrix containing the independent variable
        
    elseif strcmp(P.design,'within')
        
        design2 = repmat((1:nSBJ), [1 nCND]);
        
        cfg.statistic       = 'depsamplesFmultivariate';
        cfg.design          = [design2; design1];
        cfg.uvar            = 1; % raw of the design matrix containing the unit variable
        cfg.ivar            = 2; % raw of the design matrix containing the independent variable
        
    end
    
end
cfg.keepindividual      = 'yes';
cfg.normalizevar        = 'N-1';
cfg.subjects            = P.idxsbj;
cfg.neighbours          = neighbours;
cfg.minnbchan           = P.minnbchan;
cfg.chancmbneighbstructmat = P.chancmbneighbstructmat;

%% ------------------------------------------------------------------------
%% Permutation
data = [];
for i=1:nCND
    data=cat(1,data,Avg{i}(:));
end
if ~isfield(Avg{1}{1},'freq')
    [stat] = ft_timelockstatistics(cfg, data{:});
else
    [stat] = ft_freqstatistics(cfg, data{:});
end

% data1 = Avg{1};
% data2 = Avg{2};
% if ~isfield(Avg{1}{1},'freq')
%     [stat] = ft_timelockstatistics(cfg, data1{:}, data2{:});
% else
%     [stat] = ft_freqstatistics(cfg, data1{:}, data2{:});
% end

%% ------------------------------------------------------------------------
%% Grand Average
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
if ~isfield(Avg{1}{1},'freq')
    GAvg = cell(1,nCND);
    for i=1:nCND
        GAvg{i} = ft_timelockgrandaverage(cfg,Avg{i}{:});
        GAvg{i}.elec = Avg{1}{1}.elec;
    end
else
    GAvg = cell(1,nCND);
    for i=1:nCND
        GAvg{i} = ft_freqgrandaverage(cfg,Avg{i}{:});
        GAvg{i}.elec = Avg{1}{1}.elec;
    end
end

% %% ------------------------------------------------------------------------
% %% Grand Average Difference
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% GAvgDiff = ft_math(cfg,GAvg{1},GAvg{2});
% % GAvgDiff.dof = GAvg{1}.dof;

fprintf('\n')

end
