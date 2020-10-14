function Favg = cnd_buildfactorsfromtable(TableCNDs,theCND)
ntheCND = length(theCND);
Favg=cell(1,size(TableCNDs,2));
for f=1:size(TableCNDs,2)
    Favg{f}.name = TableCNDs.Properties.VariableNames{f};
    Favg{f}.val = unique(TableCNDs{:,f},'sorted')';
end
Favg{f}.g=nan(1,ntheCND);
for c=1:ntheCND
    idx=find(ismember(TableCNDs{:,1},theCND{c}));
    for f=1:length(Favg)
        idxf=ismember(Favg{f}.val,TableCNDs{idx(1),f});
        Favg{f}.g(c)=find(idxf);
    end
end
end