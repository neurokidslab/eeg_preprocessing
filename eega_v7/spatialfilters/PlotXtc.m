function PlotXtc(D, xl, yl,idxcnd)
if nargin<4
    idxcnd=[1 2];
end

dGavg = squeeze(D.Xgavg.mean(:,:,idxcnd));
eGavg = squeeze(D.Xgavg.error(:,:,idxcnd));
time = D.times;

figure('Position',[50 50 800, 400])

% Plot the gran average
axes('Position', [0.05 0.35 0.9 0.60])
cline = [0 0 1; 1 0 0];
hold on
lob = [];
for i=1:size(dGavg,2)
   h = shadedErrorBar(time,dGavg(:,i),eGavg(:,i),{'color',cline(i,:)}, 1 );
   lob(i) = h.mainLine;
end
ylim(yl)
xlim(xl)
line([0 0],ylim, 'color',[0 0 0])
legend(lob, D.Favg{1}.val)
title(D.Filter)


% t-test
axes('Position', [0.05 0.05 0.9 0.20])
hold on
plot(time, D.Cmp.p,'-','Color',[0 0 0])
plot(time, D.Cmp.pFDR,'-','Color',[1 0 0])
xlim([min(time),max(time)])
ylim([0 0.2])
xlim(xl)
set(gca,'Ydir','reverse')
line(xlim, [0.05 0.05],'linestyle',':','color',[0 0 0])
line(xlim, [0.01 0.01],'linestyle',':','color',[0 0 0])
