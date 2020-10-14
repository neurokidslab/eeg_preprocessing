%% ========================================================================
%% PLOT THE TOPOLOGIES

function plot_topologiesdiff(GAvg, stat, timestep, srate, zlim)

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
diff = ft_math(cfg,GAvg{1},GAvg{2});

% cfg.layout = chanloc;
% [layout, ~] = ft_prepare_layout(cfg);
layout = GAvg{1}.elec;

% define parameters for plotting
scount = length(stat.time);
sstep = floor(timestep*srate); % duration of the time window to plot topological maps (samples)
timestep = sstep/srate;

j = 1000*[ stat.time(1)/1000 : timestep : stat.time(end)/1000 ];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
if j(end)~=stat.time(end)/1000
    j=[j stat.time(end)];
end
m = [1: sstep : scount];  % temporal endpoints in EEG samples
if m(end)~=scount
    m=[m scount];
end

% get relevant (significant) values
% pos_cluster_pvals = [stat.posclusters(:).prob];
% pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% pos = ismember(stat.posclusterslabelmat, pos_signif_clust);


pos_signif = stat.stat>tinv(0.975,GAvg{1}.dof(1));
pos = pos_signif;

% plot
figure;
cfg = [];
cfg.highlight = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = layout;
cfg.zlim = zlim;
cfg.colormap = jet;
nc = 4;
nr = (size(j,2)-1)/nc;
if nr~=floor(nr), nr = floor(nr)+1; end
for k = 1:size(j,2)-1;
    subplot(nr,nc,k);
    pos_int = all(pos(:, m(k):m(k+1)), 2);
    if any(pos_int)
        disp('')
    end
    cfg.xlim=[j(k) j(k+1)];
    cfg.highlightchannel = find(pos_int);
    ft_topoplotER(cfg, diff);
end


