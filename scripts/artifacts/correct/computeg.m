% -------------------------------------------------------------------------
% This functions computes G requiered for the spherical spline interpolation
%
% INPUT
% x     : vector of length n with the x positions of the points to interpolate
% y     : vector of length n with the y positions of the points to interpolate
% z     : vector of length n with the z positions of the points to interpolate
% xelec : vector of length m with the x positions of the good electrodes
% yelec : vector of length m with the x positions of the good electrodes
% zelec : vector of length m with the x positions of the good electrodes
%
% OUTPUT
% g     : matrix of size n x m
%
% -------------------------------------------------------------------------

function g = computeg(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +...
    (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
    (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

g = zeros(length(x(:)),length(xelec));

m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    if ismatlab
        L = legendre(n,EI);
    else % Octave legendre function cannot process 2-D matrices
        for icol = 1:size(EI,2)
            tmpL = legendre(n,EI(:,icol));
            if icol == 1, L = zeros([ size(tmpL) size(EI,2)]); end;
            L(:,:,icol) = tmpL;
        end;
    end;
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);

end
