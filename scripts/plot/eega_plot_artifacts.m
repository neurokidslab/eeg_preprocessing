function eega_plot_artifacts(EEG)

nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

cmap_M = [0 0 0; 1 1 0; 1 0.6445 0; 1 0 1; 1 0 0];
cmap_BCT = [0 0 0; 1 1 0];

M = zeros(size(EEG.artifacts.BCT));
M(EEG.artifacts.BCT) = 1;
M(repmat(EEG.artifacts.BT,[nEl 1 1])) = 2;
M(repmat(EEG.artifacts.BC,[1 nS 1]))  = 3;
M(repmat(EEG.artifacts.BE,[nEl nS 1])) = 4;
cl = [0 4];

figure('Position',[1 1 600 600],'Color',[1 1 1]),

% plot BCT
hax_BCT = axes('Position',[0.1 0.55 0.8 0.35]);
imagesc(reshape(EEG.artifacts.BCT,[nEl nS*nEp]));
colormap(hax_BCT,cmap_BCT)
set(gca,'FontSize',10)
ylabel('channel')
if nEp>1
    set(gca,'XTick',nS:nS*floor(nEp/10):nS*nEp)
    set(gca,'XTickLabel',cellfun(@num2str,num2cell(1:floor(nEp/10):nEp),'uniformoutput',0)  )
    xlabel('epoch','FontSize',12)
else
    set(gca,'XTick',1:floor(nS/10):nS)
    t = num2cell(EEG.times(1:floor(nS/10):nS)/1000);
    for i=1:length(t); t{i} = num2str(t{i},'%4.1f'); end
    set(gca,'XTickLabel',t)
    xlabel('time (s)','FontSize',12)
end
title('BCT')
caxis([-0.5 1.5])
colorbar('XTick',(0:1),'XTickLabel',{'good','bad'},'FontSize',10)

% plot BCT and the rejection channels x samples
hax_m = axes('Position',[0.1 0.1 0.8 0.35]);
imagesc(reshape(M,[nEl nS*nEp]))
colormap(hax_m,cmap_M)
set(gca,'FontSize',10)
ylabel('channel')
if nEp>1
    set(gca,'XTick',nS:nS*floor(nEp/10):nS*nEp)
    set(gca,'XTickLabel',cellfun(@num2str,num2cell(1:floor(nEp/10):nEp),'uniformoutput',0)  )
    xlabel('epoch','FontSize',12)
else
    set(gca,'XTick',1:floor(nS/10):nS)
    t = num2cell(EEG.times(1:floor(nS/10):nS)/1000);
    for i=1:length(t); t{i} = num2str(t{i},'%4.1f'); end
    set(gca,'XTickLabel',t)
    xlabel('time (s)','FontSize',12)
end
title('BCT+BT+BC')
caxis([-0.5 4.5])
colorbar('XTick',(0:4),'XTickLabel',{'good','bad','BT','BC','BE'},'FontSize',10)

end