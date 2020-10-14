% -------------------------------------------------------------------------
% The function interpolates for each epoch the channels marked as bad
%
% INPUT
% DATA      : data matrix (channels x samples x epochs)
% GoodCh    : matirx indicating good channels per epoch (channels x 1 x epochs)
% chanlocs  : structure with the channels layout
% p_int     : maximun proportion of bad channels in a segment to interpolate
%
% OPTIONAL INPUTS
% distmethod    : methos used to compute the neighbours ('triangulation' or 
%                 'distance', default 'triangulation')
% pneigh        : maximun proportion of bad neighbor channels
% distneighbour : distance to define neighbours
%
% OUTPUT
% DATA      : data matrix with the interpolated data
% Interp    : logical matrix indicating the data that was interpolated
%
% USAGE
% [ DATA, Interp ] = eega_tInterpSpatial(DATA, GoodCh, chanlocs, 0.4, 'distmethod', 'triangulation')
%
% -------------------------------------------------------------------------

function [ DATA, Interp ] = eega_tInterpSpatial(DATA, GoodCh, chanlocs, p_int, varargin )

%% ------------------------------------------------------------------------
%% Parameters

P.distneighbour     = [];
P.distmethod        = 'triangulation'; % 'distance', 'triangulation' or 'template'
P.pneigh            = 1;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tInterpSpatial: Non recognized inputs')
end

% Check the input
if strcmp(P.distmethod,{'distance'  'triangulation'  'template' })
    error('eega_tInterpSpatialSegment: distmethod can have values ''distance'' / ''triangulation'' / ''template'' ')
end

%% ------------------------------------------------------------------------
%% Get the data from the EEG structure and other info
Interp  = zeros(size(DATA));

nEp  = size(DATA,3);
nEl  = size(DATA,1);
nS   = size(DATA,2);


%% ------------------------------------------------------------------------
%% INTERPOLATE

% -------------------------------------------------------------------------
% Calculate some useful position stuff
Xel = [chanlocs.X];
Yel = [chanlocs.Y];
Zel = [chanlocs.Z];
rad = sqrt(Xel.^2+Yel.^2+Zel.^2);
Xel =  Xel./rad;
Yel =  Yel./rad;
Zel =  Zel./rad;
Gel = computeg(Xel,Yel,Zel,Xel,Yel,Zel);

% -------------------------------------------------------------------------
% Define neightbours
if ~isempty(P.pneigh) && (P.pneigh<1)
    ChLabels = cell(nEl,1);
    Chanpos = cat(2,Xel(:),Yel(:),Zel(:));
    for el=1:nEl
        ChLabels{el} = chanlocs(el).labels;
    end
    layout.label = ChLabels;
    layout.chanpos = Chanpos;
    layout.label = ChLabels;
    layout.chanpos = Chanpos;
    neighbours = eega_neighborchannels(layout,'neighbourdist',P.distneighbour,'method',P.distmethod);
else
    neighbours = [];
end

% -------------------------------------------------------------------------
% Interpolate
fprintf('### Spatial interpolatiom of electrodes not working during the whole epoch ###\n')
for e = 1:nEp
    badch = logical(~GoodCh(:,:,e));
    fprintf('Epoch% 4.0d: ', e)
    if any(badch)
        fprintf('bad electrodes %s\n', mat2str(find(badch)'))
            
            % interpolate
            % ----------------------
            dataij = DATA(:,:,e);
            [badchandata, badint, interpdone] = eega_InterpSphericSpline(dataij, badch, badch, Xel, Yel, Zel, 'G', Gel, 'p_int', p_int, 'pneigh', P.pneigh, 'neighbours',neighbours);
            
            % update the information
            % ----------------------
            if interpdone
                DATA(badint,:,e)         = badchandata;
                Interp(badint,:,e)       = 1;
                Interp(badint,:,e)       = 1;
                
                if all(badint==badch)
                    fprintf('%s--> All bad channels were interpolated\n',repmat(' ',[11 1]))
                elseif any(badint)
                    fprintf('%s--> Some bad channels were interpolated\n',repmat(' ',[11 1]))
                end
            else               
                fprintf('%s--> None bad channels could be interpolated\n',repmat(' ',[11 1]))
            end
            
    else
        fprintf('none\n')
    end % if is a bad electrode in that epoch
        
end % epochs


fprintf('\n')
end

