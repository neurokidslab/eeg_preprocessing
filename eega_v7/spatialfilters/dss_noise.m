% This function used the de Cheveigné NoiseTool
%
% INPUTS
% - data : electrodes x time x trial
% - Nc : number of components to keep
% - triasl2bias : logical indexes indicating the trials to use to biased
% - K : number of components to keep in the first PCA
%
% OUTPUTS
% - X : denoised data
% - Y : projection of the data in the Nc
% - W : weights
%
% --------------------------------
% Ana Flo, October 2019, created
% --------------------------------

function [X, Y, W] = dss_noise(data, Nc, triasl2bias, K )


%% get optional inputs
if nargin<3 || isempty(triasl2bias)
    triasl2bias = true(1, size(data,3));
end
if nargin<4 || isempty(triasl2bias)
    K = 50;
end

%% implement the DSS filter

fprintf('DSS for denoising \n')

[Ne, Ns, Nt] = size(data);

% -----------------------------------
% calculate the covariance matrixes for the average response (bias to phase locked signal)
dd = mean(data(:,:,triasl2bias),3); 
dd = permute(dd,[2 1 3]);
dd = bsxfun(@minus,dd,mean(dd,1));
c1 =(1/sum(triasl2bias)).*(dd'*dd);
c1 = c1/norm(c1);

% -----------------------------------
% calculate the average of the covariance matrixes for the activity of all trials
dd = data(:,:,:);
dd = permute(dd,[2 1 3]);
dd = bsxfun(@minus,dd,mean(dd,1));
c0 = zeros(Ne);
for k=1:Nt
    c0 = c0 + (1/Nt).*(dd(:,:,k)'*dd(:,:,k));
end
c0 = c0/norm(c0);

% -----------------------------------
% compute dss filter
% 1st PCA
[V1, S1] = eig(c0);
V1 = real(V1);
S1 = real(S1);
S1 = abs(S1);
[L1, idx] = sort(diag(S1), 'descend') ;
V1 = V1(:,idx);
fprintf('PCA1: the explained variance by the %i first components is %4.2f\n',K,sum(L1(1:K))/sum(L1)*100)
N1 = diag(sqrt(1./(L1(1:K))));
c1 = N1'*V1(:,1:K)'*c1*V1(:,1:K)*N1;

% 2nd PCA
[V2, S2] = eig(c1);
V2 = real(V2);
S2 = real(S2);
[L2, idx] = sort(diag(S2), 'descend') ;
V2 = V2(:,idx);
fprintf('PCA2: the explained variance by the %i first components is %4.2f\n\n',Nc,sum(L2(1:Nc))/sum(L2)*100)

W = V1(:,1:K)*N1*V2;
N2 = diag(W'*c0*W);
W = W*diag(1./sqrt(N2)); % adjust so that components are normalized
W = W(:,1:Nc); % truncate to Nc components

% -----------------------------------
% projected the data to the components and project back to the sensor space
Y = nt_mmat(permute(data,[2 1 3]),W);
X = nt_mmat(Y,pinv(W));
Y = permute(Y,[2 1 3]);
X = permute(X,[2 1 3]);

end
