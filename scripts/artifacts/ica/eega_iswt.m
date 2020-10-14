function varargout = eega_iswt(y)

% Get inputs.
validateattributes(y, {'numeric'}, {'2d', 'real'}, ...
    'eega_iswt', 'SWC', 1);
p = size(y,1);
n = p-1;
d = y(1:n,:);
a = y(p,:);

% Compute reconstruction filters.
F = eega_coifwavf('coif5');
[Lo_D,Hi_D,Lo_R,Hi_R] = eega_orthfilt(F);
   
a = a(size(a,1),:);
[n,lx] = size(d);
for k = n:-1:1
    step = 2^(k-1);
    last = step;
    for first = 1:last
      ind = first:step:lx;
      lon = length(ind);
      subind = ind(1:2:lon);
      x1 = idwtLOC(a(subind),d(k,subind),Lo_R,Hi_R,lon,0);
      subind = ind(2:2:lon);
      x2 = idwtLOC(a(subind),d(k,subind),Lo_R,Hi_R,lon,-1);
      a(ind) = 0.5*(x1+x2);
    end
end
varargout{1} = a;


%===============================================================%
% INTERNAL FUNCTIONS
%===============================================================%
function y = idwtLOC(a,d,lo_R,hi_R,lon,shift)

y = upconvLOC(a,lo_R,lon) + upconvLOC(d,hi_R,lon);
if shift==-1
    y = y([end,1:end-1]);
end
%---------------------------------------------------------------%
function y = upconvLOC(x,f,l)

lf = length(f);
y  = eega_dyadup(x,0,1);
y  = eega_wextend('1D','per',y,lf/2);
y  = eega_wconv1(y,f);
y  = eega_wkeep1(y,l,lf);
%===============================================================%
