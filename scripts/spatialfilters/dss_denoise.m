% Denosing Spatial filter based on de Cheveigné & Simon 2008
%
% INPUTS
% - data : electrodes x time x trial
% - K  : number of components to keep in the first PCA (or percentage of variance 0-1)
% - Nc : number of components to keep in the second PCA (or percentage of variance 0-1)
% - tbias : trials to compute the bias filter (trials x number of filters)
% - tbias : trials to apply the filters (trials x number of filters)
% - tw : indexes where use to calculate the filter (if ampty all the epoch)
%
% OUTPUTS
% - X : denoised data
% - Y : projection of the data in the Nc
% - W : weights
%
% --------------------------------
% Ana Flo, July 2020, created
% --------------------------------


function [X, Y, W, summ] = dss_denoise(data, K, Nc, tbias, tapply, tw )

if nargin<6 || isempty(tw)
    tw = true(1, size(data,2));
end
if isempty(tbias)
    tbias = true(1, size(data,3));
end
if isempty(K)
    K = 50;
end
if size(tbias,1)~=size(data,3)
    error('dft_ems: Different number of trials in the data and biased vector')
end
if isempty(tapply)
    tapply = tbias;
end
if size(tapply,1)~=size(data,3)
    error('dft_ems: Different number of trials in the data and apply vector')
end
if size(tapply,2)~=size(tbias,2)
    error('dft_ems: Different number of filters in the bias and apply vector')
end

[Ne, Ns, Nt] = size(data);

%% ------------------------------------------------------------------------
%% normalize the data
normel = nan(Ne,1);
for i=1:Ne
    d = data(i,:,:);
    normel(i) = norm(d(:));
    d = d(:)/normel(i);
    data(i,:,:) = reshape(d,[1 Ns Nt]);
end

%% ------------------------------------------------------------------------
%% spatial filter
Nfilt = size(tbias,2);
% X = nan(Ne, Ns, Nt);
X = data;
if Nc<=1
    W = nan(Ne, Ne, Nt);
    Y = nan(Ne, Ns, Nt);
else
    W = nan(Ne, Nc, Nt);
    Y = nan(Nc, Ns, Nt);
end
Nstw = sum(tw);
pca1n = nan(1,Nfilt);
pca1exp = nan(1,Nfilt);
pca2n = nan(1,Nfilt);
pca2exp = nan(1,Nfilt);
for iflt=1:Nfilt
    ntbias = sum(tbias(:,iflt));
    ntapply = sum(tapply(:,iflt));
    fprintf('DSS %d: bias on %d trials, apply to %d trials...\n',iflt,ntbias, ntapply)
    
    if (ntbias>2) && (ntapply>1)
        % perform fisrt PCA
        dd = reshape(data(:,tw,tbias(:,iflt)),Ne,[])';
        dd = bsxfun(@minus,dd,mean(dd,1));
        c0 =(1/ntbias).*(dd'*dd);
        c0 = c0/norm(c0);
        [V0, S0] = eig(c0);
        V0 = real(V0);
        S0 = real(S0);
        S0 = abs(S0);
        [L0, idx] = sort(diag(S0), 'descend') ;
        V0 = V0(:,idx);
        if K<=1
            exppca1 = cumsum(L0)/sum(L0)*100;
            ki = find(exppca1>=(K*100),1);
            exppca1 = exppca1(ki);
        else
            ki = K;
            exppca1 = sum(L0(1:ki))/sum(L0)*100;
        end
        pca1n(iflt) = ki;
        pca1exp(iflt) = exppca1;
        fprintf('- PCA1, the explained variance by the %i first components is %4.2f\n',ki,exppca1)
        Y0 = dd * V0;
        % ddred = Y0(:,1:K) * pinv(V0(:,1:K));
        Y0 = Y0(:,1:ki);
        
        % bias the data by averaging and perform a second PCA
        dd = Y0';
        dd = reshape(dd,[ki, Nstw, ntbias]);
        dd = mean(dd,3);
        dd = permute(dd,[2 1 3]);
        dd = bsxfun(@minus,dd,mean(dd,1));
        c1 =(1/ntbias).*(dd'*dd);
        c1 = c1/norm(c1);
        [V1, S1] = eig(c1);
        V1 = real(V1);
        S1 = real(S1);
        [L1, idx] = sort(diag(S1), 'descend') ;
        V1 = V1(:,idx);
        if Nc<=1
            exppca2 = cumsum(L1)/sum(L1)*100;
            nci = find(exppca2>=(Nc*100),1);
            exppca2 = exppca2(nci);
        else
            nci = Nc;
            exppca2 = sum(L1(1:nci))/sum(L1)*100;
        end
        pca2n(iflt) = nci;
        pca2exp(iflt) = exppca2;
        fprintf('- PCA2, the explained variance by the %i first components is %4.2f\n',nci,exppca2)
        
        % obtain the filter
        % Y1 = Y0 * V1;
        % Y = Y1(:,1:Nc);
        % X = Y1(:,1:Nc) * pinv(V1(:,1:Nc)) * pinv(V0(:,1:K));
        % Y = reshape(permute(Y,[2 1]),[Nc, Ns, Nt]);
        % X = reshape(permute(X,[2 1]),[Ne, Ns, Nt]);
        w = V0(:,1:ki) * V1;
        % N2 = diag(W'*c0*W);
        % W = W*diag(1./sqrt(N2)); % adjust so that components are normalized
        w = w(:,1:nci); % truncate to Nc components
        for i=1:nci
            w(:,i) = w(:,i) / norm(w(:,i)); % force the filter having unit norm
        end
            
        % apply to the data and project back to sensor space
        dd = reshape(data(:,:,tapply(:,iflt)),Ne,[])';
        y = dd * w;
        x = y * pinv(w);
        y = reshape(permute(y,[2 1]),[nci, Ns, ntapply]);
        x = reshape(permute(x,[2 1]),[Ne, Ns, ntapply]);
        X(:,:,tapply(:,iflt)) = x;
        Y(1:nci,:,tapply(:,iflt)) = y;
        W(:,1:nci,tapply(:,iflt)) = repmat(w,[1 1 ntapply]);
    end
end
if Nc<=1
    idcm = all(all(isnan(W),1),3);
    W(:,idcm,:) = [];
    Y(idcm,:,:) = [];
end

% multiply by the normalization factor
X = X .* repmat(normel,[1,Ns,Nt]);

% summary
summ.pca1n = pca1n;
summ.pca1exp = pca1exp;
summ.pca2n = pca2n;
summ.pca2exp = pca2exp;

end