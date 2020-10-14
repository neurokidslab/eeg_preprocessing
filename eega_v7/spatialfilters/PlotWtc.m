function PlotWtc(D, TW)

zlim = [-0.04 0.04];
colormap = jet;
plotlayout = 0;
alpha = 90;
filelayout = '/home/an251296/Documents/MATLAB/toolbox/myfun_eega_v4/chanloc/GSN-HydroCel-129-bb.sfp';

ax_t = [150 150];
ax_m0 = [5 20];
ax_m = [50 0];


% Time windows for the topologies -----------------------------------------
time = D.times;
IDX = nan(size(TW,1),2);
nT = size(TW,1);
for i=1:nT
    idx = find( (time>TW(i,1)) & (time<TW(i,2)) );
    IDX(i,1) = idx(1);
    IDX(i,2) = idx(2);
end

% Topologies to plot ------------------------------------------------------
W = D.Wavg;
W = mean(W,3);
nE = size(W,1);
Widx = cell(1, nT);
for i=1:nT
    Widx{i}(:,:) = W(:,IDX(i,1):IDX(i,2));
end
for e=1:nE
    labE{e} = ['E' num2str(e)];
end

% Figure and axis sizes ---------------------------------------------------
screensize = get( groot, 'Screensize' );

fig(1) = ax_t(1)*nT+ax_m(1)*(nT)+ax_m0(1)*2;
fig(2) = ax_t(2)+ax_m(2)*2+ax_m0(2)*2;

fig1 = figure('Position', [screensize(1) screensize(2) fig(1) fig(2)]);
set(fig1,'color', [1 1 1])

norm = fig1.Position(3:4);
ax_t    = ax_t  ./ norm;
ax_m0   = ax_m0 ./ norm;
ax_m    = ax_m  ./ norm;

% Plot the topologies -----------------------------------------------------

% prepare the structure
cfg = [];
cfg.zlim = zlim;
cfg.xlim = 'maxmin';
cfg.projection = 'polar';
cfg.viewpoint = [];
cfg.rotate = alpha;  % rotation around the z-axis
cfg.overlap = 'shift';
cfg.colormap = colormap;
cfg.elecfile = filelayout;
cfg.parameter = 'avg';

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Prepare layout
tmpcfg     = removefields(cfg, 'inputfile');
cfg.layout = ft_prepare_layout(tmpcfg);
if plotlayout
    figure
    ft_plot_lay(cfg.layout);
end
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% plot
for i=1:nT

    ax0(1) = ax_m0(1) + (i-1)*(ax_m(1)+ax_t(1));
    ax0(2) = ax_m0(2);
    ax(i) = axes('Position', [ax0 ax_t], 'units','normalize');
    
    text(0.5,-0.05,sprintf('[%d - %d]',TW(i,1), TW(i,2)),...
        'Units', 'Normalize','HorizontalAlignment','center', 'FontSize',8)
    
%     d.avg = Widx{1};
%     d.time = time(IDX(i,1):IDX(i,2));
    d.avg = W;
    d.time = time;
    d.dof = 1;
    d.dimord = 'chan_time';
    d.fsample = 250;
    d.elec = [];
    d.label = labE;
    cfg.xlim=IDX(i,:);
    ft_topoplotER(cfg, d);
    
end
    
