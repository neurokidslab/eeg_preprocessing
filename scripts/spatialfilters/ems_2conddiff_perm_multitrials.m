% This function used the de Cheveigné NoiseTool
%
% INPUTS
%
% - data:   electrodes x time x trial
%
% - COND:   A matrix where each column is a logical index into DATA.
%           Each column selects a subset of trials belonging to one
%           experimental condition or category.
%
% - T0:     vector of 0 as 1 indicating the trials to bias.
%           (if it is empty all trials are used).
%
% - Tx:     vector of 0 as 1 indicating the trials to apply.
%           (if it is empty all trials are used).
%
% - docov:  A binary flag (0 or 1) to multiply by the inverse noise
%           covariance matrix. The matrix is computed based on the pooled
%           residuals around the condition means. With this option enabled,
%           the function implements Fisher's linear discriminant (FLD).
%           Caveat: the topographies of the spatial filter(s) are not
%           directly interpretable with this option enabled, although it
%           may improve sensitivity to very weak effects.
%
% - N2bias: proportion of the trials used to build the spatial filter.
%           (if it is empty half of the trials are used).
%
% - Np:     number of permutations.
%           (default 100).
%
% OUTPUTS
%
% - Y:      Selected trials projected into the filter
%           1 X nTimePoints X nTrials
%
% - W:      Spatial filters filter
%           nSensors X nTimePoints
%
%
% --------------------------------
% Ana Flo, Frebruary 2019
% --------------------------------

function [Y, W] = ems_2conddiff_perm_multitrials(DATA, COND, T0, Tx, docov, timecov0, N2bias, Np, ptrimmean )

if nargin<3
    T0 = [];
end
if nargin<4
    Tx = [];
end
if nargin<5
    docov = 0;
end
if nargin<6
    timecov0 = [];
end
if nargin<7
    N2bias = 0.5;
end
if nargin<8
    Np = 100;
end
if nargin<9
    ptrimmean = 0;
end


%% implement the EMS filter

fprintf('EMS filter \n')

Ne = size(DATA,1);
Ns = size(DATA,2);
Nt = size(DATA,3);

if isempty(T0), T0 = true(Nt,1); end
if isempty(Tx), Tx = true(Nt,1); end

Y = zeros(1, Ns);
W = zeros(Ne,Ns);
Npdone = 0;

tA = T0 & COND(:,1);
tB = T0 & COND(:,2);
n2bias = round(min(N2bias*sum(tA),N2bias*sum(tB)));

fprintf('Permutation:    ')
for i=1:Np
    fprintf('\b\b\b%03d',i)
    
    % select a sub-set of trials to build the filter
    [tA0, ~] = peakrand(find(tA),n2bias);
    [tB0, ~] = peakrand(find(tB),n2bias);
    
    t2applyA = COND(:,1) & Tx;
    t2applyA(tA0) = 0;
    t2applyA = find(t2applyA);
    t2applyB = COND(:,2) & Tx;
    t2applyB(tB0) = 0;
    t2applyB = find(t2applyB);
    
    % build the filter based on tA0 and tB0
    if ~isempty(t2applyA) && ~isempty(t2applyB) && length(tA0)>1 && length(tB0)>1
        [y, w] = emsfilterdata(DATA, tA0, tB0, t2applyA, t2applyB, docov, timecov0, ptrimmean);
        Y = Y + y;
        W = W + w;
        Npdone = Npdone+1;
    end
    
end
Y = Y/Npdone;
W = W/Npdone;

fprintf('\n')

end

% ---------------------------------------------------------------------
function [id1, id2] = peakrand(T,N)
n=length(T);
if n<2
    id1=[];
    id2=[];
else
    idx=randperm(n);
    id1=T(idx(1:N));
    id2=T(idx(N+1:end));
end
end

function [y, w] = emsfilterdata(data, t2biasA, t2biasB, t2applyA, t2applyB, docov, timecov0, ptrimmean)

if nargin<7 || isempty(timecov0)
    timecov0 = true(1,size(data,2));
end
if nargin<8
    ptrimmean = 0;
end

c0=[];
if numel(docov)>1
    if (size(docov,1)~=size(data,1)) || (size(docov,2)~=size(data,1))
        error('emsfilterdata: the covarinaca has to be squared matrix with size equal to the number of electrodes')
    else
        c0=docov;
        docov=1;
    end
end

% normalize
normalizedata = 0;
if normalizedata
    mu = mean(reshape(data,size(data,1),[]),2);
    var = std(reshape(data,size(data,1),[]),0,2);
    datan = (data - mu) ./ var;
else
    datan = data;
end

% average response (bias to phase locked signal)
if ptrimmean~=0
    DA = trimmean(datan(:,:,t2biasA),ptrimmean,3);
    DB = trimmean(datan(:,:,t2biasB),ptrimmean,3);
else
    DA = mean(datan(:,:,t2biasA),3);
    DB = mean(datan(:,:,t2biasB),3);
end
DDIFF = DA - DB;

if docov
    
    if isempty(c0)
        % covariance matrix for the noise
        ddA=permute(datan(:,timecov0,t2biasA),[2 1 3]);
        ddA=bsxfun(@minus,ddA,mean(ddA,3));
        ddB=permute(datan(:,timecov0,t2biasB),[2 1 3]);
        ddB=bsxfun(@minus,ddB,mean(ddB,3));
        c0=real( nt_cov(cat(3,ddA,ddB)) );
    end
    
    % compute the filter
    invc0=pinv(c0);
    w=invc0*DDIFF;
    
else
    w=DDIFF;
end

% normalize the filter
w = bsxfun(@rdivide,w,sqrt(sum(w.^2,1)));

% projection
if ptrimmean~=0
    dA = trimmean(data(:,:,t2applyA),ptrimmean,3);
    dB = trimmean(data(:,:,t2applyB),ptrimmean,3);
else
    dA = mean(data(:,:,t2applyA),3);
    dB = mean(data(:,:,t2applyB),3);
end
davg = dA - dB;
y = sum( davg .* w, 1);

end






















