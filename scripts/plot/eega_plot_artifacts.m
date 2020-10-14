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

figure('Position',[1 1 600 600]),

% plot BCT
hax_BCT = axes('Position',[0.1 0.55 0.8 0.35]);
title('BCT')
imagesc(reshape(EEG.artifacts.BCT,[nEl nS*nEp]));
colormap(hax_BCT,cmap_BCT)
ylabel('channel')
if nEp>1
    set(gca,'XTick',nS:nS*floor(nEp/10):nS*nEp)
    set(gca,'XTickLabel',cellfun(@num2str,num2cell(1:floor(nEp/10):nEp),'uniformoutput',0)  )
    xlabel('epoch')
else
    set(gca,'XTick',1:floor(nS/10):nS)
    set(gca,'XTickLabel',EEG.times(1:floor(nS/10):nS))
    xlabel('time')
end
caxis([0 1])
colorbar

% plot BCT and the rejection channels x samples
hax_m = axes('Position',[0.1 0.1 0.8 0.35]);
title('BCT+BT+BC')
imagesc(reshape(M,[nEl nS*nEp]))
colormap(hax_m,cmap_M)
ylabel('channel')
if nEp>1
    set(gca,'XTick',nS:nS*floor(nEp/10):nS*nEp)
    set(gca,'XTickLabel',cellfun(@num2str,num2cell(1:floor(nEp/10):nEp),'uniformoutput',0)  )
    xlabel('epoch')
else
    set(gca,'XTick',1:floor(nS/10):nS)
    set(gca,'XTickLabel',EEG.times(1:floor(nS/10):nS))
    xlabel('time')
end
caxis(cl)
colorbar

end