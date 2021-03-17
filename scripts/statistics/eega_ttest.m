function [tstat, pval] = eega_ttest(EEG, opt)


if ~isfield(opt,'Fsbj')
    error('No subject factor was provided: Fsbj')
end
if ~isfield(opt,'Fcnd')
    error('No subject factor was provided: Fcnd')
end
if isa(opt.Fsbj,'double')
    opt.Fsbj = EEG.F{opt.Fsbj};
elseif ischar(opt.Fsbj)
    ok=0;
    i=0;
    while ~ok
        i=i+1;
        ok = strcmp(opt.Fsbj, EEG.F{i}.name);
    end
    if ok
        opt.Fsbj = EEG.F{i};
    else
        error('Fsbj not found')
    end    
end
if isa(opt.Fcnd,'double')
    opt.Fcnd = EEG.F{opt.Fcnd};
elseif ischar(opt.Fcnd)
    ok=0;
    i=0;
    while ~ok
        i=i+1;
        ok = strcmp(opt.Fcnd, EEG.F{i}.name);
    end
    if ok
        opt.Fcnd = EEG.F{i};
    else
        error('Fcnd not found')
    end      
end
if ~isfield(opt,'ch')
    opt.ch = (1:size(EEG.data,1));
end
if ~isfield(opt,'times')
    opt.idxtimes = (1:size(EEG.data,2));
else
    opt.idxtimes = find(EEG.times>=opt.times(1) & EEG.times<=opt.times(end));
end
if ~isfield(opt,'idxsbj')
    opt.idxsbj = (1:length(opt.Fsbj.val));
end
if ~isfield(opt,'idxcnd')
    opt.idxcnd = [1 2];
end
if ~isfield(opt,'test')
    opt.test = 'ttest';
end

% get the data
Nsmp = length(opt.idxtimes);
Nch = length(opt.ch);
Nsbj = length(opt.idxsbj);
D1 = nan(Nch,Nsmp,Nsbj);
D2 = nan(Nch,Nsmp,Nsbj);
for isbj=1:Nsbj
    idxsbji = opt.Fsbj.g==opt.idxsbj(isbj);
    idxcnd1 = opt.Fcnd.g==opt.idxcnd(1);
    idxcnd2 = opt.Fcnd.g==opt.idxcnd(2);
    D1(:,:,isbj) = EEG.data(opt.ch,opt.idxtimes,idxcnd1 & idxsbji);
    D2(:,:,isbj) = EEG.data(opt.ch,opt.idxtimes,idxcnd2 & idxsbji);
end

% t-test
fh = str2func(opt.test); 
tstat = nan(Nch,Nsmp);
pval = nan(Nch,Nsmp);
for ich=1:Nch
    for ismpl=1:Nsmp
        [h,p,ci,stats]=fh(squeeze(D1(ich,ismpl,:)),squeeze(D2(ich,ismpl,:)));
        tstat(ich,ismpl) = stats.tstat;
        pval(ich,ismpl) = p;
    end
end

end