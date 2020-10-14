function T = cnd_buildtablecond(factors, F)

tb={};
vnames={};
k=0;
for i=1:length(F)
    if any(ismember(factors,F{i}.name))
        k=k+1;
        vnames{k}=F{i}.name;
        n=size(tb,1);
        if n==0, n=1; end
        val=F{i}.val;
        tb=repmat(tb,[length(val) 1]);
        val=repmat(val,[n 1]);
        tb=cat(2,tb,val(:));
    end
end
if ~isempty(tb)
    cnds=tb(:,1)';
    if size(tb,2)>1
        for i=2:size(tb,2)
            cnds = strcat(cnds,tb(:,i)');
        end
    end
    vnames = cat(1, {'Condition'}, vnames{:})';
    T = cell2table(cat(2,cnds(:), tb ),'VariableNames',vnames);
else
    T=[];
end
end
