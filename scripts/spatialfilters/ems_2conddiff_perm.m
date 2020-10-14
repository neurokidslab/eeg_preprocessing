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
% Ana Flo, Frebruary 2018
% --------------------------------

function [Y, W, COND] = ems_2conddiff_perm(DATA, COND, T0, Tx, docov, N2bias, Np )

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

fprintf('EMS filter \n')

Ne = size(DATA,1);
Ns = size(DATA,2);
Nt = size(DATA,3);

if isempty(T0), T0 = true(Nt,1); end
if isempty(Tx), Tx = true(Nt,1); end

Y = zeros(1, Ns,Nt);
W = zeros(Ne,Ns,1);
Ntp = zeros(Nt,1);
Nw = 0;

tA = T0 & COND(:,1);
tB = T0 & COND(:,2);
n2bias = round(min(N2bias*sum(tA),N2bias*sum(tB)));

fprintf('Permutation:    ')
for i=1:Np
    fprintf('\b\b\b\b%04d',i)
    
    % select a sub-set of trials to build the filter
    [tA0, ~] = peakrand(find(tA),n2bias);
    [tB0, ~] = peakrand(find(tB),n2bias);
    
    t2apply = Tx;
    t2apply(tA0) = 0;
    t2apply(tB0) = 0;
    t2apply = find(t2apply);
    
    % build the filter based on tA0 and tB0
    if ~isempty(t2apply) && length(tA0)>1 && length(tB0)>1
        [y, w] = emsfilterdata(DATA, tA0, tB0, t2apply, docov);
        Y(:,:,t2apply) = Y(:,:,t2apply) + y;
        Ntp(t2apply) = Ntp(t2apply)+1;
        W = W + w;
        Nw = Nw+1;
    end
    
end
Ntp = Ntp(Tx);
Y   = Y(:,:,Tx);
Y   = bsxfun(@rdivide,Y,permute(Ntp,[2 3 1]));
W   = W/Nw;
COND = COND(Tx,:);


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

function [y, w] = emsfilterdata(data, t2biasA, t2biasB, t2apply, docov)

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
DD = mean(data(:,:,t2biasA),3)-mean(data(:,:,t2biasB),3);

if docov
    
    if isempty(c0)
        
        % covariance matrix for the noise
        ddA = data(:,:,t2biasA);
        ddA = bsxfun(@minus,ddA,mean(ddA,3));
        ddB = data(:,:,t2biasB);
        ddB = bsxfun(@minus,ddB,mean(ddB,3));
        dd = cat(3,ddA,ddB);
        c0 = zeros(size(dd,1));
        for i=1:size(dd,3)
            c0 = c0 + (dd(:,:,i)*dd(:,:,i)');
        end
        c0 = c0 *(1/size(dd,3));
        c0 = c0/norm(c0);
        c0 = real(c0);
        
%         % covariance matrix for the noise
%         ddA=permute(data(:,:,t2biasA),[2 1 3]);
%         ddA=bsxfun(@minus,ddA,mean(ddA,3));
%         ddB=permute(data(:,:,t2biasB),[2 1 3]);
%         ddB=bsxfun(@minus,ddB,mean(ddB,3));
%         c0=real( nt_cov(cat(3,ddA,ddB)) );
    end
    
    % compute the filter
    invc0=pinv(c0);
    w=invc0*DD;
    
else
    w=DD;
end

% normalize the filter
for ii=1:size(w,2)
    w(:,ii)=w(:,ii)/norm(w(:,ii));
end

% projection
for k=1:length(t2apply)
    x = data(:,:,t2apply(k));
    yk = sum(x.*w, 1);
    y(1,:,k)=permute(yk,[2 1]);
end


end
