function [Lo_D,Hi_D,Lo_R,Hi_R] = eega_orthfilt(W,P)
%ORTHFILT Orthogonal wavelet filter set.
%   [LO_D,HI_D,LO_R,HI_R] = ORTHFILT(W) computes the
%   four filters associated with the scaling filter W 
%   corresponding to a wavelet:
%   LO_D = decomposition low-pass filter
%   HI_D = decomposition high-pass filter
%   LO_R = reconstruction low-pass filter
%   HI_R = reconstruction high-pass filter.
%
%   See also BIORFILT, QMF, WFILTERS.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 13-May-2003.
%   Copyright 1995-2004 The MathWorks, Inc.

% Check arguments.
if nargin<2 , P = 0; end

% Normalize filter sum.
W = W/sum(W);

% Associated filters.
Lo_R = sqrt(2)*W;
Hi_R = qmf(Lo_R,P);
Hi_D = flip(Hi_R);
Lo_D = flip(Lo_R);

end

function y = qmf(x,p)
%QMF    Quadrature mirror filter.
%   Y = QMF(X,P) changes the signs of the even index entries
%   of the reversed vector filter coefficients X if P is even.
%   If P is odd the same holds for odd index entries.
%
%   Y = QMF(X) is equivalent to Y = QMF(X,0).

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 20-Dec-2010.
%   Copyright 1995-2010 The MathWorks, Inc.

% Check arguments.
if nargin == 1 , p = 0; end
if (p~=fix(p)) || (p<0)
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'))
end

% Compute quadrature mirror filter.
y = x(end:-1:1);
first = 2-rem(p,2);
y(first:2:end) = -y(first:2:end);
end
