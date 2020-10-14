% This function projects trials in two filters, one built based on one
% condition and the other on the oder condition. Afterwards calculates the
% difference between the projections.
% The filters are built multiple times based on a sub-set of trials
%
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
% - WA:     Spatial filter for the first condition
%           nSensors X nTimePoints
%
% - WB:     Spatial filter for the second condition
%           nSensors X nTimePoints
%
%
% --------------------------------
% Ana Flo, Frebruary 2018
% --------------------------------

function [Y, W, COND] = ems_2cond_perm(DATA, COND, T0, Tx, docov, N2bias, Np )

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
    N2bias = 0.5;
end
if nargin<7
    Np = 100;
end

%% implement the EMS filter

fprintf('EMS filter 2 conditions (based on sub-sets)\n')

Ne = size(DATA,1);
Ns = size(DATA,2);
Nt = size(DATA,3);

if isempty(T0), T0 = true(Nt,1); end
if isempty(Tx), Tx = true(Nt,1); end

Y = zeros(1, Ns,Nt);
WA = zeros(Ne,Ns,1);
WB = zeros(Ne,Ns,1);
Ntp = zeros(Nt,1);
Nw = 0;

tA = T0 & COND(:,1);
tB = T0 & COND(:,2);
n2bias = round(min(N2bias*sum(tA),N2bias*sum(tB)));

fprintf('Permutation:    ')
for i=1:Np
    fprintf('\b\b\b\b%04d',i)
    
    % select a sub-set of trials to build the filters
    [tA0, ~] = peakrand(find(tA),n2bias);
    [tB0, ~] = peakrand(find(tB),n2bias);
    
    t2apply = Tx;
    t2apply(tA0) = 0;
    t2apply(tB0) = 0;
    t2apply = find(t2apply);
    
    % build the filters based on tA0 and tB0
    if ~isempty(t2apply) && length(tA0)>1 && length(tB0)>1
        [y, wA, wB] = emsfilterdata(DATA, tA0, tB0, t2apply, docov);
        Y(:,:,t2apply) = Y(:,:,t2apply) + y;
        Ntp(t2apply) = Ntp(t2apply)+1;
        WA = WA + wA;
        WB = WB + wB;
        Nw = Nw+1;
    end
    
end
Ntp = Ntp(Tx);
Y   = Y(:,:,Tx);
Y   = bsxfun(@rdivide,Y,permute(Ntp,[2 3 1]));
WA  = WA/Nw;
WB  = WB/Nw;
COND = COND(Tx,:);
W = cat(3,WA,WB);

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

function [y, wA, wB] = emsfilterdata(data, t2biasA, t2biasB, t2apply, docov)

c0=[];
if numel(docov)>1
    if (size(docov,1)~=size(data,1)) || (size(docov,2)~=size(data,1))
        error('emsfilterdata: the covarinaca has to be squared matrix with size equal to the number of electrodes')
    else
        c0=docov;
        docov=1;
    end
end

y=nan(1,size(data,2),length(t2apply));

% average response (bias to phase locked signal)
DDA=real(mean(data(:,:,t2biasA),3));
DDB=real(mean(data(:,:,t2biasB),3));

if docov
    
    if isempty(c0)
        % covariance matrix for the noise
        ddA=permute(data(:,:,t2biasA),[2 1 3]);
        ddA=bsxfun(@minus,ddA,mean(ddA,3));
        ddB=permute(data(:,:,t2biasB),[2 1 3]);
        ddB=bsxfun(@minus,ddB,mean(ddB,3));
        c0A=real( nt_cov(ddA) );
        c0B=real( nt_cov(ddB) );
    else
        c0A=c0;
        c0B=c0;
    end
    
    % compute the filters
    wA=pinv(c0A)*DDA;
    wB=pinv(c0B)*DDB;
else
    wA=DDA;
    wB=DDB;
end

% normalize the filters
for ii=1:size(wA,2)
    wA(:,ii)=wA(:,ii)/norm(wA(:,ii));
end
for ii=1:size(wB,2)
    wB(:,ii)=wB(:,ii)/norm(wB(:,ii));
end

% projection
for k=1:length(t2apply)
    ykA=sum(data(:,:,t2apply(k)) .* wA, 1);
    ykB=sum(data(:,:,t2apply(k)) .* wB, 1);
    y(1,:,k)=permute(ykA-ykB,[2 1]);
end


end
