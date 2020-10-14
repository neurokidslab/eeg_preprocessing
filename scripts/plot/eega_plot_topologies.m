%% ========================================================================
%% PLOT THE TOPOLOGIES

function hf = eega_plot_topologies( GAvg, varargin )


%% ========================================================================
%% Read the optional parameters

% Default parameters
P.time          =[];      
P.trials        =1;      
P.stat          =[];
P.zlim          =[-10 10];
P.timestep      =0.1;
P.colormap      =jet;
P.plotlayout    =0;
P.alphalayout   =90;
P.thefreq       = [];
P.nc            = 14;

% Optional parameters
if mod(length(varargin),2)==1
    warning('Optional parameters come by pairs')
    return
end
if ~isempty(varargin)
    for f=1:2:length(varargin)
        Pop.(varargin{f}) = varargin{f+1};
    end
    fff = fieldnames(P);
    for f=1:numel(fff)
        if isfield(Pop,fff{f})
            P.(fff{f}) = Pop.(fff{f});
        end
    end
end

%% ========================================================================
%% Arreange the data and plot
if isfield(GAvg,'time')
    datatype = 'timelocked';
elseif isfield(GAvg,'freq')
    datatype = 'powerspectrum';
end
srate=GAvg.fsample;

if isempty(P.time) && isfield(GAvg, 'time')
    P.time = [GAvg.time(1) GAvg.time(end)];
end
if ischar(P.colormap) && strcmp('redblue',P.colormap)
    figure; P.colormap = eval(P.colormap); close gcf
end

% Define parameters for plotting
switch datatype
    case 'timelocked'
        
        % Select the data 
        GAvg.avg = GAvg.avg(:,:,P.trials);

        % Temporal endpoints (in seconds) of the ERP average computed in each subplot
        sstep = floor(P.timestep*srate); % duration of the time window to plot topological maps (samples)
        timestep = sstep/srate;
        
        j = 1000*[ P.time(1)/1000 : timestep : P.time(end)/1000 ];
        if j(end)~=P.time(end)/1000
            j=[j P.time(end)];
        end
        
    case 'powerspectrum'
        
        % Select the data 
        GAvg.powspctrm = GAvg.powspctrm(:,:,P.trials);

        j = repmat(P.thefreq,[1 2]);
end


% Prepare the structure
cfg             = [];
cfg.highlight   = 'on';
cfg.comment     = 'xlim';
cfg.commentpos  = 'title';
cfg.projection  = 'polar';
cfg.viewpoint   = [];
cfg.rotate      = P.alphalayout;  % rotation around the z-axis
cfg.overlap     = 'shift';
cfg.zlim        = P.zlim;
cfg.colormap    = P.colormap;

% -------- Prepare layout --------
tmpcfg     = removefields(cfg, 'inputfile');
cfg.layout = ft_prepare_layout(tmpcfg, GAvg);
if P.plotlayout
    figure
    ft_plot_lay(cfg.layout);
end
% ---------------------------------

%% ========================================================================
%% plot
ax_m = [20 20];
ax_0 = [50 50];
ax_t = [150 150];
nc = min(P.nc,size(j,2)-1);
nr = (size(j,2)-1)/nc;
if nr~=floor(nr), nr = floor(nr)+1; end

pos(1) = 2*ax_0(1) + nc*ax_t(1) +(nc-1)*ax_m(1);
pos(2) = 2*ax_0(2) + nr*ax_t(2) +(nr-1)*ax_m(2);
hf = figure('Position',[0 0 pos]);

ax_m = ax_m ./ pos;
ax_0 = ax_0 ./ pos;
ax_t = ax_t ./ pos;

r=1;
for k = 1:size(j,2)-1
    if floor(k/nc)>r
        r=r+1;
    end
    c = k - floor((k-1)/nc);
    cfg.xlim=[j(k) j(k+1)];
    
    x = ax_0(1) + (c-1)*ax_t(1) + (c-1)*ax_m(1);
    y = ax_0(2) + (nr-r)*ax_t(2) + (nr-r)*ax_m(2);
    axes('Position',[x y ax_t])
    ft_topoplotER(cfg, GAvg);
end


%% AUXILIARY FUNCTIONS

% Red blue colormap
function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 
