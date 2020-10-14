function eega_plot_timeseries(EEG, idEl, idTime, idEp, Ylim,plotBCT,plotCCT,plotBT)

TIME=EEG.times;
DATA=EEG.data;

sz=size(DATA);

if nargin<2 || isempty(idEl)
    idEl = 1:sz(1);
end
if nargin<3 || isempty(idTime)
    idTime = 1:sz(2);
end
if nargin<4 || isempty(idEp)
    idEp = 1;
    DATA=reshape(DATA,sz(1),[]);
end
if nargin<5 
    Ylim = [];
end
if nargin<6 
    plotBCT = 1;
end
if nargin<7 
    plotCCT = 1;
end
if nargin<8
    plotBT = 1;
end

if isfield(EEG.artifacts,'BCT')
    BCT=EEG.artifacts.BCT;
    BCT=BCT(idEl,idTime,idEp);
    BCT=BCT(:);
else
    BCT=[];
end
if isfield(EEG.artifacts,'CCT')
    CCT=EEG.artifacts.CCT;
    CCT=CCT(idEl,idTime,idEp);
    CCT=CCT(:);
else
    CCT=[];
end
if isfield(EEG.artifacts,'BT')
    BT=EEG.artifacts.BT;
    BT=BT(1,idTime,idEp);
    BT=BT(:);
else
    BT=[];
end

TIME=TIME(idTime)-TIME(idTime(1));
TIME=TIME(:)/1000;
TIME=bsxfun(@plus,repmat(TIME,[1 length(idEp)]),[0:TIME(end):TIME(end)*(length(idEp)-1)]);
TIME=TIME(:);
DATA=DATA(idEl,idTime,idEp);
DATA=DATA(:);

hold on
H=cell(1,1);
L=cell(1,1);

h=plot(TIME,DATA,'k');
H{1}=h;
L{1}=sprintf('E%3d',idEl);
xlim([min(TIME) max(TIME)])
if ~isempty(Ylim)
    ylim(Ylim)
else
    Ylim=ylim;
end

if ~isempty(BCT) && plotBCT
    yfillBCT = 2*BCT-1;
    yfillBCT = [yfillBCT; -ones(length(idTime)*length(idEp),1)];
    xfillBCT = [TIME; flipud(TIME)];
    h=fill(xfillBCT,yfillBCT*Ylim,[1 0 0],'LineStyle','none','FaceAlpha',0.2);
    H{end+1}=h(1);
    L{end+1}='BCT';
end
if ~isempty(CCT) && plotCCT
    yfillCCT = 2*CCT-1;
    yfillCCT = [yfillCCT; -ones(length(idTime)*length(idEp),1)];
    xfillCCT = [TIME; flipud(TIME)];
    h=fill(xfillCCT,yfillCCT*Ylim,[0 1 0],'LineStyle','none','FaceAlpha',0.2);
    H{end+1}=h(1);
    L{end+1}='CCT';
end
if ~isempty(BT) && plotBT
    bt=BT;
    h=plot(TIME(bt),0*bt(bt),'s','MarkerSize',6,'Color',[1 0 0]);
    H{end+1}=h;
    L{end+1}='BT';
end
for i=1:length(idEp)
    line(TIME((i-1)*sz(2)+1)*[1 1],ylim,'color',[0 0 0],'LineStyle',':')
end
legend([H{:}]',L{:})
text(0.01,0.05,sprintf('samples [%d %d]',idTime(1),idTime(end)), 'units','normalized','fontsize',8)


end