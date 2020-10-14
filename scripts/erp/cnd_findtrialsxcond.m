function [condition, theCND, trialsxCND, Fcnd] = cnd_findtrialsxcond(TableCNDs, F)

CND = TableCNDs{:,1};
nCND = length(CND);
theCND = unique(CND,'sorted');
ntheCND = length(theCND);
theF = TableCNDs.Properties.VariableNames(2:end);
nF = length(theF);

nT=length(F{1}.g);

% Index for each relevant factor
idxF = nan(1,nF);
for f=1:nF
    idxf = 1;
    while ~strcmp(F{idxf}.name,theF{f})
        idxf = idxf+1;
    end
    idxF(f) = idxf;
end

% Define the condition of each trial
FxT=nan(nT,length(F),ntheCND);
condition=zeros(nT,1);
trialsxCND=zeros(1,ntheCND);
for i=1:nT
    vals = cell(1,nF);
    for f=1:nF
        vals{f}=F{idxF(f)}.val{F{idxF(f)}.g(i)};
    end
    test=find(all(strcmp(repmat(vals,[nCND 1]),TableCNDs{:,2:end}),2));
    if isempty(test)
        warning('trial %i does not belong to any condition',i)
    else
        cndi=TableCNDs{test(1),1};
        condition(i)=find( strcmp(cndi,theCND) );
        trialsxCND(condition(i))=trialsxCND(condition(i))+1;
        for f=1:length(F)
            if isempty(F{f}.g)
                disp('opd')
            else
                FxT(trialsxCND(condition(i)),f,condition(i))=F{f}.g(i);
            end
        end
    end
end

% Define the factors
Fcnd{1}.name='CND';
Fcnd{1}.val=theCND;
Fcnd{1}.g=condition(:)';

end