% This function provides neighbor channels. It is an addaptation of the
% Fieldtrip function ft_prepare_neighbours
%
% Ana Flo November 2017

function neighbours = eega_neighborchannels(layout,varargin)

%% ------------------------------------------------------------------------
%% Parameters

P.neighbourdist = [];
P.method    = 'triangulation'; % 'distance', 'triangulation' or 'template'

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_neighborchannels: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Interpolate

if isfield(layout,'chanpos')
    chanpos = layout.chanpos;
else
    error('the layout has to have a field chanpos')
end
if isfield(layout,'label')
    label = layout.label;
else
    error('the layout has to have a field label')
end
if all(size(chanpos,2)~=[2 3])
    error('chanpos has size channels x 2 (X Y) or channels x 3 (X Y Z)')
end
if size(chanpos,1)~=length(label)
    error('the number of labels has to be the number of channels')
end

switch lower(P.method)
    case 'distance'
        % use a smart default for the distance
        if isempty(P.neighbourdist)
            d = pdist(chanpos);
            P.neighbourdist = max(d)/5;
            fprintf('using a distance threshold of %g\n', P.neighbourdist);
        end
        neighbours = compneighbstructfromgradelec(chanpos, label, P.neighbourdist);
    case {'triangulation', 'tri'} % the latter for reasons of simplicity
        if size(chanpos, 2)==2 || all(chanpos(:,3)==0)
            % the sensor positions are already on a 2D plane
            prj = chanpos(:,1:2);
        else
            % project sensor on a 2D plane
            prj = elproj(chanpos);
        end
        % make a 2d delaunay triangulation of the projected points
        tri = delaunay(prj(:,1), prj(:,2));
        tri_x = delaunay(prj(:,1)./2, prj(:,2));
        tri_y = delaunay(prj(:,1), prj(:,2)./2);
        tri = [tri; tri_x; tri_y];
        neighbours = compneighbstructfromtri(chanpos, label, tri);
    otherwise
        error('Method ''%s'' not known', P.method);
end

k = 0;
for i=1:length(neighbours)
    if isempty(neighbours(i).neighblabel)
        warning('no neighbours found for %s\n', neighbours(i).label);
    end
    k = k + length(neighbours(i).neighblabel);
end
if k==0
    error('No neighbours were found!');
end

fprintf('there are on average %.1f neighbours per channel\n', k/length(neighbours));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that compute the neighbourhood geometry from the
% gradiometer/electrode positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours] = compneighbstructfromgradelec(chanpos, label, neighbourdist)

nsensors = length(label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
    dist(i,:) = sqrt(sum((chanpos(1:nsensors,:) - repmat(chanpos(i,:), nsensors, 1)).^2,2))';
end;

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));

% construct a structured cell array with all neighbours
neighbours=struct;
for i=1:nsensors
    neighbours(i).label       = label{i};
    neighbours(i).neighblabel = label(find(channeighbstructmat(i,:)));
    neighbours(i).idx         = find(channeighbstructmat(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the neighbourhood geometry from the
% triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours] = compneighbstructfromtri(chanpos, label, tri)

nsensors = length(label);

channeighbstructmat = zeros(nsensors,nsensors);
% mark neighbours according to triangulation
for i=1:size(tri, 1)
    channeighbstructmat(tri(i, 1), tri(i, 2)) = 1;
    channeighbstructmat(tri(i, 1), tri(i, 3)) = 1;
    channeighbstructmat(tri(i, 2), tri(i, 1)) = 1;
    channeighbstructmat(tri(i, 3), tri(i, 1)) = 1;
    channeighbstructmat(tri(i, 2), tri(i, 3)) = 1;
    channeighbstructmat(tri(i, 3), tri(i, 2)) = 1;
end

% construct a structured cell array with all neighbours
neighbours = struct;
alldist = [];
for i=1:nsensors
    neighbours(i).label       = label{i};
    neighbidx                 = find(channeighbstructmat(i,:));
    neighbours(i).dist        = sqrt(sum((repmat(chanpos(i, :), numel(neighbidx), 1) - chanpos(neighbidx, :)).^2, 2));
    alldist                   = [alldist; neighbours(i).dist];
    neighbours(i).neighblabel = label(neighbidx);
    neighbours(i).idx         = neighbidx;
end

% remove neighbouring channels that are too far away (imporntant e.g. in
% case of missing sensors)
neighbdist = mean(alldist)+3*std(alldist);
for i=1:nsensors
    idx = neighbours(i).dist > neighbdist;
    neighbours(i).dist(idx)         = [];
    neighbours(i).neighblabel(idx)  = [];
    neighbours(i).idx(idx)          = [];
end
neighbours = rmfield(neighbours, 'dist');
