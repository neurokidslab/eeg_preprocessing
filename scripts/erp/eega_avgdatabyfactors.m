% -------------------------------------------------------------------------
% This function calculate the average across conditions, for the conditions
% identified in the eventField of the Epochs, when it takes the values
% eventValues. It also gives an average across all trials
%
% INPUTS
% EEG           EEG structure
% TableCNDs     it specifies what to average. 
%               
%               * It can be a table specifing how to average trials
%       - In each column it has the factors that will be usded to define
%       conditions (the name has to be same than in F{k}.name, and the
%       values have to be some of F{k}.val). The first coloumn is the
%       condition
%       - In each row each the possible combinations for each condition are
%       defined
%
%       For example
%
%            -----------------------------------------------
%           | Condition     |   Img1   |   Aud1   |   Aud2   |
%           |------------------------------------------------|
%           | Acongruent    |    a     |    x1    |    x1    |
%           | Acongruent    |    a     |    x2    |    x2    |
%           | Bcongruent    |    b     |    x1    |    x1    |
%           | Bcongruent    |    b     |    x2    |    x2    |
%           | Aincongruent  |    a     |    x1    |    x2    |
%           | Aincongruent  |    a     |    x2    |    x1    |
%           | Bincongruent  |    b     |    x1    |    x2    |
%           | Bincongruent  |    b     |    x2    |    x1    |
%            -----------------------------------------------
%
%               * It can also be a cell array indicating the factors to
%               take into acount to define conditions. In this case the
%               table will be built
%               * If it is empty all trials are average. 
%
% OPTIONAL INPUTS
% DataField     fields in EEG containig the data. It can be a cell with 
%               multiple fields. In that case the baseline correction is 
%               applied to all fileds. Default {'data'}
% FactorField   field in EEG containg the factors characterizing epochs. 
%               Default 'F'
% dim2avg       dimention across which data should be average. Default the
%               last dimension 
%
% OUTPUTS
% EEG           EEG structure with average data
%
% -------------------------------------------------------------------------

function EEGavg = eega_avgdatabyfactors(EEG, TableCNDs, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.DataField     = {'data'};
P.FactorField   = {'F'};
P.EEGlabFomat   = 1;
P.Stats         = 0;
P.Silent        = 0;
P.Weights       = [];
P.Plot          = 0;
P.dim2avg       = [];
P.trimmean      = [];
P.factors2keep  = [];

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_avgdata: Non recognized inputs')
end

if ~P.Silent; fprintf('### Average ERP ###\n'); end

if ~iscell(P.DataField)
    P.DataField={P.DataField};
end
if ~iscell(P.FactorField)
    P.FactorField={P.FactorField};
end
goodsbj=1;
for i=1:length( P.DataField)
    if isempty(EEG.(P.DataField{i}))
        goodsbj=0;
    end
end
if ~isempty(P.Weights)
    P.Weights = permute(P.Weights(:),[2 3 1]);
end
if isempty(P.factors2keep)
    P.factors2keep = (1:length(EEG.(P.FactorField{1})));
end

%% ------------------------------------------------------------------------
%% Build the table defying the conditions
if ~isempty(TableCNDs) && ~istable(TableCNDs)
    TableCNDs = cnd_buildtablecond(TableCNDs, EEG.(P.FactorField{1}));
end
if isempty(TableCNDs)
    goodsbj=0;
else
    F = EEG.(P.FactorField{1});
    [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(TableCNDs, F);
    ntheCND = length(theCND);
end   

%% ------------------------------------------------------------------------
%% Average
if ~isempty(P.dim2avg)
    dim2avg = P.dim2avg;
else
    dim2avg = ndims(EEG.(P.DataField{i}));
end
Dmean=cell(1,length(P.DataField));
Dsd=cell(1,length(P.DataField));
Dn=cell(1,length(P.DataField));
Derror=cell(1,length(P.DataField));
if goodsbj
    for i=1:length( P.DataField)
        s = size(EEG.(P.DataField{i}));
        ss = s; 
        if length(ss)==2; ss = [ss 1]; end
        ss(dim2avg) = [];
        Dmean{i} = nan([ss ntheCND]);
        Dsd{i} = nan([ss ntheCND]);
        Dn{i} = nan([ss ntheCND]);
        Derror{i} = nan([ss ntheCND]);
        if ~isempty(P.Weights)
            if size(P.Weights,P.dim2avg)~=size(EEG.(P.DataField{i}),P.dim2avg)
                error('wrong number of weights')
            end
            EEG.(P.DataField{i})=bsxfun(@mtimes,EEG.(P.DataField{i}),P.Weights*length(P.Weights));
        end
        for c=1:ntheCND
            
            if dim2avg==3
                dat = EEG.(P.DataField{i})(:,:,condition==c);
            elseif dim2avg==4
                dat = EEG.(P.DataField{i})(:,:,:,condition==c);
%             else
%                 dat = takelstdim(EEG.(P.DataField{i}), find(condition==c));
            end
            
            if isempty(P.trimmean) || (P.trimmean==0) || size(dat,3)<3
                [avg, sd, n, err] = meanstats(dat, dim2avg);
            else
                [avg, sd, n, err] = trimmeanstats(dat, P.trimmean, dim2avg);
            end
            if dim2avg==3
                Dmean{i}(:,:,c) = avg;
                Dsd{i}(:,:,c) = sd;
                Dn{i}(:,:,c) = n;
                Derror{i}(:,:,c) = err;
            elseif dim2avg==4
                Dmean{i}(:,:,:,c) = avg;
                Dsd{i}(:,:,:,c) = sd;
                Dn{i}(:,:,:,c) = n;
                Derror{i}(:,:,:,c) = err;
            end
        end
    end
else
    for i=1:length( P.DataField)
        Dmean{i} = [];
        Dsd{i} = [];
        Dn{i} = [];
        Derror{i} = [];
    end
end

%% ------------------------------------------------------------------------
%% Define the new factors
if goodsbj
    Favg = cnd_buildfactorsfromtable(TableCNDs,theCND);
else
    Favg = {};
end

% add factors that are consistent
favgnames = TableCNDs.Properties.VariableNames;
for fi=1:length(EEG.F)
    if ~any(strcmp(EEG.F{fi}.name, favgnames))
        vals = nan(size(TableCNDs,1),1);
        for icnd=1:size(TableCNDs,1)
            ok  = any(P.factors2keep == fi);
            idxcnd = condition==icnd;
            if any(idxcnd)
                valsi = EEG.F{fi}.g(idxcnd);
                if length(unique(valsi))==1
                    vals(icnd) = unique(valsi);
                else
                    ok = 0;
                end
            else
                vals(icnd) = 1;
            end
        end
        if ok
            fact = length(Favg)+1;
            Favg{fact}.name = EEG.F{fi}.name;
            Favg{fact}.val = EEG.F{fi}.val;
            Favg{fact}.g = vals;
        end
   end
end

%% ------------------------------------------------------------------------
%% Define the as empty is it was a bad subject
if goodsbj
    newname = cell(size(theCND));
    for i=1:length(theCND)
        newname{i} = theCND{i};
%         if length(theCND{i})>15
%             newname{i} = theCND{i}(end-15:end);
%         else
%             newname{i} = theCND{i};
%         end
    end
    newname = matlab.lang.makeValidName(newname);
    newname = matlab.lang.makeUniqueStrings(newname);
    TrialsxCND = array2table(trialsxCND,'VariableNames',newname');
else    
    TrialsxCND = [];
end

%% ------------------------------------------------------------------------
%% Fill the structure

% Fill the fields that are the same than EEG
fff=fieldnames(EEG);
fff(ismember(fff,P.FactorField{1})) = [];
for i=1:length( P.DataField)
    fff(ismember(fff,P.DataField{i})) = [];
end
for j=1:length(fff)
    EEGavg.(fff{j}) = EEG.(fff{j});
end

% Fill the fields that are different than EEG
s = size(Dmean{1});
EEGavg.nbchan = s(1);
EEGavg.pnts = s(2);
EEGavg.trials = s(end);
EEGavg.(P.FactorField{1}) = Favg;
for i=1:length( P.DataField)
    EEGavg.(P.DataField{i}) = Dmean{i};
    if P.Stats
        EEGavg.([P.DataField{i} '_stats']).sd = Dsd{i};
        EEGavg.([P.DataField{i} '_stats']).n = Dn{i};
        EEGavg.([P.DataField{i} '_stats']).error = Derror{i};
    end
end
EEGavg.TrialsxCND = TrialsxCND;
if strcmp('data',P.DataField)
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'summary')
        EEGavg.RejectedData = EEG.artifacts.summary;
    end
end

% Remove some fields
if P.EEGlabFomat
    f2rmv = {'epoch' 'reference' 'filter' 'artifacts'...
        'eventdescription' 'epochdescription' 'reject' 'stats' 'specdata' 'specicaact'...
        'splinefile' 'icasplinefile' 'dipfit' 'event' 'urevent'};
    for i=1:length(f2rmv)
        if isfield(EEGavg,f2rmv{i})
            EEGavg = rmfield(EEGavg,f2rmv{i});
        end
    end
else
    f2rmv = {'epoch' 'reference' 'filter' 'artifacts'...
        'eventdescription' 'epochdescription' 'reject' 'stats' 'specdata' 'specicaact'...
        'splinefile' 'icasplinefile' 'dipfit' 'event' 'urevent'...
        'icachansind' 'icaweights' 'icasphere' 'icawinv' 'icaact'};
    for i=1:length(f2rmv)
        if isfield(EEGavg,f2rmv{i})
            EEGavg = rmfield(EEGavg,f2rmv{i});
        end
    end
end

fprintf('\n')

%% ------------------------------------------------------------------------
%% Plot
if P.Plot
    plotavg(Dmean,Favg,times,trialsxCND)
end

end

%% ========================================================================

% function Txavg = buildtablecnds(Fxavg, F)
% 
% tb={};
% vnames={};
% k=0;
% for i=1:length(F)
%     if any(ismember(Fxavg,F{i}.name))
%         k=k+1;
%         vnames{k}=F{i}.name;
%         n=size(tb,1);
%         if n==0, n=1; end
%         val=F{i}.val;
%         tb=repmat(tb,[length(val) 1]);
%         val=repmat(val,[n 1]);
%         tb=cat(2,tb,val(:));
%     end
% end
% if ~isempty(tb)
%     cnds=tb(:,1)';
%     if size(tb,2)>1
%         for i=2:size(tb,2)
%             cnds=strcat(cnds,tb(:,i)');
%         end
%     end
%     vnames=cat(1,{'Condition'},vnames{:})';
%     Txavg = cell2table(cat(2,cnds(:),tb),'VariableNames',vnames);
% else
%     Txavg=[];
% end
% end
% 
% 
% function Favg = definefactorsfromtableCNDs(TableCNDs,theCND)
% ntheCND = length(theCND);
% Favg=cell(1,size(TableCNDs,2));
% for f=1:size(TableCNDs,2)
%     Favg{f}.name = TableCNDs.Properties.VariableNames{f};
%     Favg{f}.val = unique(TableCNDs{:,f},'sorted')';
% end
% Favg{f}.g=nan(1,ntheCND);
% for c=1:ntheCND
%     idx=find(ismember(TableCNDs{:,1},theCND{c}));
%     for f=1:length(Favg)
%         idxf=ismember(Favg{f}.val,TableCNDs{idx(1),f});
%         Favg{f}.g(c)=find(idxf);
%     end
% end
% end
% 
% function Favg=definefactorsfromtableFnew(TableFnew,theCND)
% ntheCND = length(theCND);
% Fnew = TableFnew.Properties.VariableNames;
% nFbew = length(Fnew);
% Favg = cell(1,length(Fnew));
% for f=1:nFbew
%     Favg{f}.name = Fnew{f};
%     Favg{f}.val = unique(TableFnew{:,f},'sorted');
%     Favg{f}.g = nan(1,ntheCND);
%     for c=1:ntheCND
%         Favg{f}.g(c) = find(strcmp(TableFnew{c,f},Favg{f}.val));
%     end
% end
% 
% 
% end

function plotavg(Davg,Favg,times,trialsxCND)

nCND = size(Davg,P.dim2avg);
ax_avg = [500 300];
m_avg = [50 50];
fig_ = [ ax_avg(1)+m_avg(1)*2 nCND*ax_avg(2)+m_avg(2)*(nCND+1)];
ax_avg = ax_avg./fig_;
m_avg = m_avg./fig_;

yl = max(abs(Davg(:)));
yl = [-yl yl];


figure('Position',[1 1 fig_],'Name',sprintf('SBJ %s',nameSubj));

for c=1:nCND
    axes('Position',[m_avg(1) (nCND-c)*ax_avg(2)+(nCND-c+1)*m_avg(2) ax_avg]);
    
    dc=Davg(:,:,c);
    if any(~isnan(dc(:)))
        plot(times,Davg(:,:,c)','k')
        ylim(yl)
        
        text(times(10),yl(1)+0.1*diff(yl),sprintf('n trials : %d',trialsxCND(c)))
        s='';
        for f=1:length(Favg)
            s = [s, sprintf('%s - ',Favg{f}.val{Favg{f}.g(c)})];
        end
        s = s(1:end-3);
        title(s)
    end
end

end


