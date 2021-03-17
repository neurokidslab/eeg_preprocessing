function eega_plot_roi(EEG, opt)

if ~isfield(opt,'ch')
    error('No channels were provided: ch')
end
if ~isfield(opt,'Fsbj')
    error('No subject factor was provided: Fsbj')
end
if ~isfield(opt,'Fcnd')
    error('No subject factor was provided: Fcnd')
end
if ~isfield(opt,'roiname')
    opt.roiname = '';
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
if ~isfield(opt,'idxsbj')
    opt.idxsbj = (1:length(opt.Fsbj.val));
end
if ~isfield(opt,'idxcnd')
    opt.idxcnd = (1:length(opt.Fcnd.val));
end
if ~isfield(opt,'xlim')
    opt.xlim = [min(EEG.times) max(EEG.times)];
end
if ~isfield(opt,'ylim')
    opt.ylim = [];
end
if ~isfield(opt,'colors')
    opt.colors = lines(length(opt.idxcnd));
end
if ~isfield(opt,'transparent')
    opt.transparent = 0.2;
end
if ~isfield(opt,'lwidth')
    opt.lwidth = 1.5;
end
if ~isfield(opt,'lstyle')
    opt.lstyle = repmat({'-'},[1 length(opt.idxcnd)]);
end
if ~isfield(opt,'ploterror')
    opt.ploterror = 1;
end
if ~isfield(opt,'chanlocs')
    opt.chanlocs = [];
end
if ~isfield(opt,'xevt')
    opt.xevt = [];
end
if ~isfield(opt,'pval')
    opt.pval = [];
end
if ~isfield(opt,'pvalFDR')
    opt.pvalFDR = [];
end
if ~isfield(opt,'pvaltwind')
    opt.pvaltwind = [];
end


Nsbj = length(opt.idxsbj);
Ncnd = length(opt.idxcnd);
Nsmp = size(EEG.data,2);

% get the data
D = nan(Nsmp,Ncnd,Nsbj);
for isbj=1:Nsbj
    idxsbji = opt.Fsbj.g==opt.idxsbj(isbj);
    for icnd=1:Ncnd
        idxcndi = opt.Fcnd.g==opt.idxcnd(icnd);
        idx = idxsbji & idxcndi;
        D(:,icnd,isbj) = mean(EEG.data(opt.ch,:,idx),1);
    end
end
opt.times = EEG.times;
opt.mean = mean(D,3);
opt.error = std(D,0,3)/sqrt(Nsbj);

if (strcmp(opt.pval,'do')) && (Ncnd==2)
    if isempty(opt.pvaltwind)
        opt.pvalidxtw = true(Nsmp,1);
    else
        opt.pvalidxtw = (opt.times>=opt.pvaltwind(1) & opt.times<=opt.pvaltwind(2));
        opt.pvalidxtw = opt.pvalidxtw(:);
    end
    opt.pval = nan(Nsmp,1);
    for ismp=1:Nsmp
        if opt.pvalidxtw(ismp)
            [~,opt.pval(ismp)] = ttest(squeeze(D(ismp,1,:)),squeeze(D(ismp,2,:)));
        end
    end 
end
if (strcmp(opt.pvalFDR,'do')) && ~isempty(opt.pval)
    opt.pvalFDR = nan(Nsmp,1);
    pval = opt.pval(~isnan(opt.pval));
    opt.pvalFDR(~isnan(opt.pval)) = correctFDR(pval);
end

% plot
figure('position',[0 0 800 500],'Color',[1 1 1]),
if ~isempty(opt.chanlocs)
    axerp = axes('Position',[0.30 0.10 0.65 0.80]);
else
    axerp = axes('Position',[0.10 0.10 0.80 0.80]);
end
hold on
Htc = {};
for icnd=1:Ncnd
    if opt.ploterror
        h = shadedErrorBar(opt.times, opt.mean(:,icnd), opt.error(:,icnd), ...
            {'color', opt.colors(icnd,:),'LineStyle',opt.lstyle{icnd},'LineWidth',opt.lwidth}, opt.transparent);
        Htc{icnd}=h.mainLine;
    else
        h = plot(opt.times, opt.mean(:,icnd),'color',opt.colors(icnd,:),'LineStyle',opt.lstyle{icnd},'LineWidth',opt.lwidth);
        Htc{icnd}=h;
    end
end
if isempty(opt.ylim)
    opt.ylim = ylim;
end
if ~isempty(opt.xevt)
    for il=1:length(opt.xevt)
        line(opt.xevt(il)*[1 1],opt.ylim,'LineStyle',':','LineWidth',1,'Color',[0 0 0])
    end
end
if ~isempty(opt.pval)
    idx = opt.pvalidxtw & (opt.pval<=0.05);
    v = (opt.ylim(1)+0.1*diff(opt.ylim));
    v = v*ones(sum(idx),1);
    c = 0.8*[1 1 1];
    plot(opt.times(idx),v,'s','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerSize',6)
end
if ~isempty(opt.pvalFDR)
    idx = opt.pvalidxtw & (opt.pvalFDR<=0.05);
    v = (opt.ylim(1)+0.1*diff(opt.ylim));
    v = v*ones(sum(idx),1);
    c = 0.6*[1 1 1];
    plot(opt.times(idx), v,'s','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerSize',6)
end

xlim(opt.xlim)
ylim(opt.ylim)
legend([Htc{:}]',opt.Fcnd.val(opt.idxcnd),'Location','southwest','FontSize',12)
set(gca,'Box','on')
title(opt.roiname)
if ~isempty(opt.chanlocs)
    axtopo = axes('Position',[0.025 0.60 0.25 0.25]);
    set(axtopo, 'visible', 'off')
    topodata = ones(size(EEG.data,1),1);
    topoplot(topodata, opt.chanlocs(1:128), 'maplimits',[0 1],'numcontour',0,...
        'style','both','headrad',0.5,'colormap',[1 1 1],...
        'nosedir', '+X', 'emarker2', {opt.ch, 'x','k', 2, 1});
end
set(gcf,'Color',[1 1 1])

end