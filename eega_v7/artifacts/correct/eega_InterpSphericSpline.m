% -------------------------------------------------------------------------
% This functions performs a spherical spline interpolation of bad channels
%
% INPUT
% data      : data (channels x samples)
% ch2int    : logical indexes indicating bad channels to interpolate
%            (channels x 1)
% badch     : logical indexes indicating bad channels to do NOT 
%             interpolate (channels x 1)
% X         : x position for each channel (channels x 1)
% Y         : y position for each channel (channels x 1)
% Z         : z position for each channel (channels x 1)
%
% OPTIONAL INPUTS
% G             : G matrix used for interpolation. If not provided it will
%                 be computed (slower)
% p_int         : maximun proportion of bad channels allowed to interpolate
% pneigh        : maximun proportion of bad neighbor channels
% neighbours    : structure defining nighbours. If not provided it will be 
%                 computed base on the electrodes positions
% distmethod    : methos used to compute the neighbours ('triangulation' or 
%                 'distance', default 'triangulation')
% distneighbour : distance to define neighbours
%
% OUTPUT
% interpdata    : interpolated data
% ChInt         : logical indexes indicating the interpoleted channels
% interpdone    : logical index indicating whethe the interpolation took
%                 place
%
% USAGE
% [badchandata, ChInt, interpdone] = eega_InterpSphericSpline(data, ch2int, badch, X, Y, Z, 'p_int', 0.4, 'distmethod', 'triangulation')
%
% -------------------------------------------------------------------------

function [interpdata, ChInt, interpdone] = eega_InterpSphericSpline(data, ch2int, badch, X, Y, Z, varargin)

%% ------------------------------------------------------------------------
%% Parameters

P.G             = [];
P.p_int         = 1;
P.pneigh        = 1;
P.neighbours    = [];
P.distneighbour = [];
P.distmethod    = 'triangulation'; % 'distance', 'triangulation' or 'template'

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_InterpSphericSpline: Non recognized inputs')
end

G               = P.G;
neighbours      = P.neighbours;

%% ------------------------------------------------------------------------
%% Interpolate

badch = logical(badch);
ch2int = logical(ch2int);
AllGoodCh = ~badch;
ChInt = ch2int;
electrodes = length(badch);

% define neightbours
% ------------------
if isempty(neighbours) && ~isempty(P.pneigh) && (P.pneigh<1)
    ChLabels = cell(electrodes,1);
    Chanpos = cat(2,X(:),Y(:),Z(:));
    for el=1:electrodes
        ChLabels{el} = sprintf('E%d',el);
    end
    layout.label = ChLabels;
    layout.chanpos = Chanpos;
    layout.label = ChLabels;
    layout.chanpos = Chanpos;
    neighbours = eega_neighborchannels(layout,'neighbourdist',P.distneighbour,'method',P.distmethod);
end

% check there are enought good channels
% -------------------------------------
if ~isempty(P.p_int)
    if (sum(badch)/length(badch)) <= P.p_int
        ok_p = 1;
    else
        ok_p = 0;
    end
else
    ok_p = 1;
end

% check there are enought good neightbor channels
% -----------------------------------------------
if ~isempty(neighbours) && (P.pneigh<1)
    allgoodch  = find(AllGoodCh);
    thebadch = find(ch2int);
    p_gn = ones(length(thebadch),1);
    for j=1:length(thebadch)
        goodneighbors = intersect(allgoodch,neighbours(thebadch(j)).idx);
        p_gn(j) = length(goodneighbors)/length(neighbours(thebadch(j)).idx);
    end

    if all(p_gn >= (1-P.pneigh))
        ok_pneighbors = 1;
    else
        chnoneighbors = thebadch(p_gn < (1-P.pneigh));
        ChInt(chnoneighbors) = false;
        if any(ChInt)
            ok_pneighbors = 1;
        else
            ok_pneighbors = 0;
        end
    end
else
    ok_pneighbors = 1;
end

ok = ok_p && ok_pneighbors;

% spheric spline interpolation
% -----------------------------
if ok

    values = data(logical(AllGoodCh),:);
    newchans = sum(ChInt);
    numpoints = size(values,2);

    if isempty(G)
        xelec = X(AllGoodCh);
        yelec = Y(AllGoodCh);
        zelec = Z(AllGoodCh);
        xbad = X(ChInt);
        ybad = Y(ChInt);
        zbad = Z(ChInt);
        Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
        Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);
    else
        Gelec = G(AllGoodCh,AllGoodCh);
        Gsph = G(ChInt,AllGoodCh);
    end

    % compute solution for parameters C
    meanvalues = mean(values);
    values = values - repmat(meanvalues, [size(values,1) 1]); % make mean zero

    values = [values;zeros(1,numpoints)];
    C = pinv([Gelec;ones(1,length(Gelec))]) * values;
    clear values;
    interpdata = zeros(newchans, numpoints);

    % apply results
    for j = 1:size(Gsph,1)
        interpdata(j,:) = sum(C .* repmat(Gsph(j,:)', [1 size(C,2)]));
    end
    interpdata = interpdata + repmat(meanvalues, [size(interpdata,1) 1]);
    interpdone = 1;
else
    interpdata = data(ChInt,:);
    interpdone = 0;
end

end
