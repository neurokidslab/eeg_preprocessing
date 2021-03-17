function eega_plot_erp(EEG, opt)

if nargin<2
    opt = [];
end
if ~isfield(opt,'cnd')
    opt.cnd = [];
end
if ~isfield(opt,'gfp')
    opt.gfp = 1;
end
if ~isfield(opt,'gfpcol')
    opt.gfpcol = [0 1 0];
end
if ~isfield(opt,'gfplwidth')
    opt.gfplwidth = 2;
end
if ~isfield(opt,'chmrklwidth')
    opt.gfplwidth = 2;
end
if ~isfield(opt,'xlim')
    opt.xlim = [min(EEG.times) max(EEG.times)];
end
if ~isfield(opt,'ylim')
    opt.ylim = [];
end
if ~isfield(opt,'xevt')
    opt.xevt = [];
end
if ~isfield(opt,'chmrk')
    opt.chmrk = [];
end
if ~isfield(opt,'chmrklwidth')
    opt.chmrklwidth = 2;
end
if ~isfield(opt,'figposition')
    opt.figposition = [0 0 700 500];
end

%% Build the table defying the conditions
if ~isempty(opt.cnd) && ~istable(opt.cnd)
    opt.cnd = cnd_buildtablecond(opt.cnd, EEG.(P.FactorField{1}));
end
if isempty(opt.cnd)
    condition = ones([1 size(EEG.data,3)]);
    theCND = {'all'};
    ntheCND = 1;
else
    F = EEG.(P.FactorField{1});
    [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(TableCNDs, F);
    ntheCND = length(theCND);
end   

%% PLot
for icnd=1:ntheCND
    d = nanmean(EEG.data(:,:,condition==icnd),3);
    gfp = std(d);
    figure('position',opt.figposition,'Color',[1 1 1])
    hold on
    plot(EEG.times, d','k')
    if ~isempty(opt.chmrk)
        for k=1:length(opt.chmrk)
            plot(EEG.times,d(opt.chmrk{k}{1},:)','color',opt.chmrk{k}{2},'linewidth',opt.chmrklwidth)
        end
    end
    if opt.gfp
        plot(EEG.times,gfp,'color',opt.gfpcol,'linewidth',opt.gfplwidth)
    end
    title(theCND{icnd})
    if isempty(opt.ylim)
        opt.ylim = ylim;
    end
    if ~isempty(opt.xevt)
        for il=1:length(opt.xevt)
            line(opt.xevt(il)*[1 1],opt.ylim,'LineStyle',':','LineWidth',1,'Color',[0 0 0])
        end
    end
    xlim(opt.xlim)
    ylim(opt.ylim)
    
    set(gca,'Box','on')

end

end
