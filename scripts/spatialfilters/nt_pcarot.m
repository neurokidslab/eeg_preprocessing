function [topcs,eigenvalues]=nt_pcarot(cov,nkeep,threshold,minkeep,N)
% [topcs,eigenvalues]=pcarot(cov,nkeep,threshold,N) - PCA matrix from covariance
%
%  topcs: PCA rotation matrix
%  eigenvalues: PCA eigenvalues
%  
%  cov: covariance matrix
%  nkeep: number of component to keep
%  thresholds: discard components below this threshold
%  minkeep: I added an inferior limit for the number of components to keep 
%  N: eigs' K parameter (if absent: use eig)
%
% NoiseTools

if nargin<5; N=[]; end
if nargin<4; minkeep=1; end
if nargin<3; threshold=[]; end
if nargin<2; nkeep=[]; end

if ~isempty(N); 
    [V, S] = eigs(cov,N) ;  
else
    [V, S] = eig(cov) ;  
end

V=real(V);
S=real(S);
[eigenvalues, idx] = sort(diag(S)', 'descend') ;
topcs = V(:,idx);

% truncate
if ~isempty(threshold)
    ii=find(eigenvalues/eigenvalues(1)>threshold);
    if minkeep<length(ii)
        topcs=topcs(:,ii);
        eigenvalues=eigenvalues(ii);
    else
        topcs=topcs(:,1:minkeep);
        eigenvalues=eigenvalues(1:minkeep);
    end
end

if ~isempty(nkeep) & minkeep<nkeep
    nkeep=min(nkeep,size(topcs,2));
    topcs=topcs(:,1:nkeep);
    eigenvalues=eigenvalues(1:nkeep);
end
