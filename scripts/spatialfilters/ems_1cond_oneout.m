% This function projects trials to a filters built based on the average of
% a subset of trials.
% 
% The filter is built multiple times based on a sub-set of trials
%
% This function used the de Cheveigné NoiseTool
%
% INPUTS
%
% - DATA:   electrodes x time x trial
%
% - COND:   A column vector of logical indexes into DATA.
%           It selects a subset of trials belonging to one experimental
%           condition or category to biased the filter.
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
% - donorm: A binary flag (0 or 1) to normalize or not the filter. 
%           Default 1
%
% OUTPUTS
%
% - Y:      Selected trials projected into the filter
%           1 X nTimePoints X nTrials
%
% - W:      Spatial filter 
%           nSensors X nTimePoints
%
% --------------------------------
% Ana Flo, Frebruary 2019
% --------------------------------

function [Y, W, COND] = ems_1cond_oneout(DATA, COND, T0, Tx, docov, donorm )

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
    donorm = 1;
end
%% implement the EMS filter

fprintf('EMS filter (based on sub-sets)\n')

Ne = size(DATA,1);
Ns = size(DATA,2);
Nt = size(DATA,3);

if isempty(T0), T0 = true(Nt,1); end
if isempty(Tx), Tx = true(Nt,1); end

Y = zeros(1,Ns,Nt);
W = zeros(Ne,Ns,1);
Ntp = zeros(Nt,1);
Nw = 0;

tA = T0 & COND(:,1);

for i=1:length(Tx)
    if Tx(i)
        % select the trials to bias
        t0 = tA; t0(i) = 0;
        t2apply = tA;
        t2apply(t0) = 0;
        t2apply = find(t2apply);
        
        % build the filter based on t0
        if ~isempty(t2apply) && length(t0)>1
            [y, w] = emsfilterdata(DATA, t0, t2apply, docov, donorm);
            Y(:,:,t2apply) = Y(:,:,t2apply) + y;
            Ntp(t2apply) = Ntp(t2apply)+1;
            W = W + w;
            Nw = Nw+1;
        end
    end
end
Ntp = Ntp(Tx);
Y   = Y(:,:,Tx);
Y   = bsxfun(@rdivide,Y,permute(Ntp,[2 3 1]));
W   = W/Nw;
COND = COND(Tx,:);



end


function [y, w] = emsfilterdata(data, t2bias, t2apply, docov, donorm)

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
DD = mean(data(:,:,t2bias),3);
% DD = real(mean(data(:,:,t2bias),3));

if docov
    
    if isempty(c0)
        % covariance matrix for the noise
        dd=permute(data(:,:,t2bias),[2 1 3]);
        dd=bsxfun(@minus,dd,mean(dd,3));
        c0=real( nt_cov(dd) );
    end
    
    % compute the filter
    if size(c0,3)==1
        w = pinv(c0)*DD;
    else
        w = nan(size(DD,1),size(DD,2));
        for ii=1:size(c0,3)
            w(:,ii) = pinv(c0(:,:,ii))*DD(:,ii);
        end
    end
else
    w=DD;
end
% normalize the filter
if donorm
    for ii=1:size(w,2)
        w(:,ii)=w(:,ii)/norm(w(:,ii));
    end
end

% projection 
for k=1:length(t2apply)
    yk = sum(data(:,:,t2apply(k)) .* w, 1);
    y(1,:,k) = permute(yk,[2 1]);
end


end
