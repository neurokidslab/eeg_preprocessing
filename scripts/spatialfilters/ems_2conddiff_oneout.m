% This function used the de Cheveigné NoiseTool
%
% INPUTS
%
% - data:   electrodes x time x trial
%
% - COND:   A matrix where each column is a logical index into DATA.
%           Each column selects a subset of trials belonging
%           to one experimental condition or category.
%
% - docov:  A binary flag (0 or 1) to multiply by the inverse noise
%           covariance matrix. The matrix is computed based on the pooled
%           residuals around the condition means. With this option enabled,
%           the function implements Fisher's linear discriminant (FLD).
% 			If a square matrix is provided, it is used as covariance of the 
% 			noise.
%           Caveat: the topographies of the spatial filter(s) are not
%           directly interpretable with this option enabled, although it
%           may improve sensitivity to very weak effects.
%
% - T0:     vector of 0 as 1 indicating the trials to bias.
%           (if it is empty all trials are used).
%
% - Tx:     vector of 0 as 1 indicating the trials to apply.
%           (if it is empty all trials are used).
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

function [Y, W, COND] = ems_2conddiff_oneout(DATA, COND, T0, Tx, docov)

if nargin<3
    T0 = [];
end
if nargin<4
    Tx = [];
end
if nargin<5
    docov = 0;
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

fprintf('Trial:    ')
for i=1:length(Tx)
    if Tx(i)
        fprintf('\b\b\b\b%04d',i)
        
        % select the trials to bias
        tA0 = tA; tA0(i) = 0;
        tB0 = tB; tB0(i) = 0;
        t2apply = T0;
        t2apply(tA0) = 0;
        t2apply(tB0) = 0;
        t2apply = find(t2apply);
        
        % build the filter based on tA0 and tB0
        if length(tA0)>1 && length(tB0)>1
            [y, w] = emsfilterdata(DATA, tA0, tB0, t2apply, docov);
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

% t=(1/250:1/250:size(Y,2)/250)*1000-500;
% figure, plot(t, mean(Y(1,:,tA),3)-mean(Y(1,:,tB),3))


fprintf('\n')

end

% 
% function [y, w] = emsfilterdata(data, t2biasA, t2biasB, t2apply, Ncmp)
% 
% y = nan(1,size(data,2),length(t2apply));
% 
% % average response (bias to phase locked signal)
% DD = real(mean(data(:,:,t2biasA),3))-real(mean(data(:,:,t2biasB),3));
% 
% 
% % covariance matrix for the noise
% ddA = data(:,:,t2biasA);
% ddA = bsxfun(@minus,ddA,mean(ddA,3));
% ddB = data(:,:,t2biasB);
% ddB = bsxfun(@minus,ddB,mean(ddB,3));
% dd = cat(3,ddA,ddB);
% c0 = zeros(size(dd,1));
% for i=1:size(dd,3)
%     c0 = c0 + (dd(:,:,i)*dd(:,:,i)');
% end
% c0 = c0 *(1/size(dd,3));
% c0 = c0/norm(c0);
% c0 = real(c0);
% 
% % compute the covariance of the mean (bias)
% c1 = DD*DD';
% c1 = c1/norm(c1);
% c1 = real(c1);
% 
% % eigenvalues decomposition
% [V, L] = eig(c1,c0);
% [L, idx] = sort(diag(L), 'descend') ;
% V = V(:,idx);
% 
% % determine the spacial filter and normalize it
% w = V(:,1:Ncmp);
% for k=1:size(w,2)
%     w(:,k) = w(:,k) / norm(w(:,k)); % force the filter having unit norm
% end
% 
% % project each trial
% for k=1:length(t2apply)
%     dk = data(:,:,t2apply(k));
%     yk = dk' * w;
%     y(1:Ncmp,:,k) = yk';
% end
% 
% 
% end

function [y, w] = emsfilterdata(data, t2biasA, t2biasB, t2apply, docov)

c0=[];
if numel(docov)>1
    if (size(docov,1)~=size(data,1)) || (size(docov,2)~=size(data,1))
        error('emsfilterdata: the covariance has to be squared matrix with size equal to the number of electrodes')
    else
        c0=docov;
        docov=1;
    end
end

y = nan(1,size(data,2),length(t2apply));

% average response (bias to phase locked signal)
DD=real(mean(data(:,:,t2biasA),3))-real(mean(data(:,:,t2biasB),3));

if docov
    
    if isempty(c0)
        % covariance matrix for the noise
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
    end
    
    % compute the filter
%     w = c0\DD;
    w = pinv( c0 ) * DD;
    
else
    w = DD;
end


% % check if there are channels at each sample for which the filter is too 
% % high and set them as zero
% for i=1:size(w,2)
%     wi = w(:,i);
%     p = prctile(wi(:),[25 75]);
%     thrsh = 1.75*diff(p);
%     idx = (wi>(p(2)+thrsh)) | (wi<(p(1)-thrsh));
%     w(idx,i) = 0;
% end

% normalize the filter
w = bsxfun(@rdivide,w,sqrt(sum(w.^2,1)));

% projection
for k=1:length(t2apply)
    yk = sum(data(:, :, t2apply(k)) .* w, 1);
    y(1,:,k) = permute(yk, [2 1]);
end


end


% function [y, w] = emsfilterdata(data, t2biasA, t2biasB, t2apply, docov)
% 
% c0=[];
% if numel(docov)>1
%     if (size(docov,1)~=size(data,1)) || (size(docov,2)~=size(data,1))
%         error('emsfilterdata: the covariance has to be squared matrix with size equal to the number of electrodes')
%     else
%         c0=docov;
%         docov=1;
%     end
% end
% 
% y = nan(1,size(data,2),length(t2apply));
% 
% % average response (bias to phase locked signal)
% DD=real(mean(data(:,:,t2biasA),3))-real(mean(data(:,:,t2biasB),3));
% 
% if docov
%     
%     if isempty(c0)
%         % covariance matrix for the noise
%         ddA = permute(data(:,:,t2biasA),[2 1 3]);
%         ddA = bsxfun(@minus,ddA,mean(ddA,3));
%         ddB = permute(data(:,:,t2biasB),[2 1 3]);
%         ddB = bsxfun(@minus,ddB,mean(ddB,3));
%         c0 = real( nt_cov(cat(3,ddA,ddB)) );
%     end
%     
%     % compute the filter
%     w = pinv( (c0 - mean(c0(:))) / std(c0(:))) * ((DD - mean(DD)) / std(DD(:)));
%     
% else
%     w = DD;
% end
% 
% % normalize the filter
% for ii=1:size(w,2)
%     w(:,ii)=w(:,ii)/norm(w(:,ii));
% end
% 
% % projection
% for k=1:length(t2apply)
%     yk = sum(data(:, :, t2apply(k)) .* w, 1);
%     y(1,:,k) = permute(yk, [2 1]);
% end
% 
% 
% end
