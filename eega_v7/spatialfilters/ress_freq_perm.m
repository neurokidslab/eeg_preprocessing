% This function used the de Cheveigné NoiseTool
%
% INPUTS
%
% - DATA : electrodes x frequency x trial
%
% - T0:     vector of 0 as 1 indicating the trials to bias.
%           (if it is empty all trials are used).
%
% - Tx:     vector of 0 as 1 indicating the trials to apply.
%           (if it is empty all trials are used).
%
% - N2bias: proportion of the trials used to build the spatial filter.
%           (if it is empty half of the trials are used).
%
% - Np:     number of permutations.
%           (default 100).
%

% OUTPUTS
%
%
% --------------------------------
% Ana Flo, Frebruary 2018
% --------------------------------

function [Y, W] = ress_freq_perm(DATA, srate, ftarget, peakwidt, neighfreq, neighwidt, T0, Tx, N2bias, Np )

if nargin<7
    T0 = [];
end
if nargin<8
    Tx = [];
end
if nargin<9
    Np = 100;
end

%% implement the RESS filter

fprintf('RESS filter \n')

Ne = size(DATA,1);
Ns = size(DATA,2);
Nt = size(DATA,3);
Nf = length(ftarget);

if isempty(T0), T0 = true(Nt,1); end
if isempty(Tx), Tx = true(Nt,1); end

Y = zeros(Nf,Ns,Nt);
W = zeros(Ne,Nf);
Ntp = zeros(Nt,1);
Nw = 0;

n2bias = round(N2bias*sum(T0));

fprintf('Permutation:    ')
for i=1:Np
    fprintf('\b\b\b%03d',i)
    
    % select a sub-set of trials to build the filter
    [t0, ~] = peakrand(find(T0),n2bias);
    
    t2apply = Tx;
    t2apply(t0) = 0;
    t2apply = find(t2apply);
    
    % build the filters based on t0
    if ~isempty(t2apply) && length(t0)>1
        [y, w] = ressfilterdata(DATA, srate, t0, t2apply, ftarget, peakwidt, neighfreq, neighwidt );
        Y(:,:,t2apply) = Y(:,:,t2apply) + y;
        Ntp(t2apply) = Ntp(t2apply)+1;
        W = W + w;
        Nw = Nw+1;
    end   
       
end
Ntp=Ntp(Tx);
Y=Y(:,:,Tx);
Y = bsxfun(@rdivide,Y,permute(Ntp,[2 3 1]));
W=W/Nw;
fprintf('\n')

end

% ---------------------------------------------------------------------
function [idA, idB] = peakrand(T,N)
n=length(T);
if n<2
    idA=[];
    idB=[];
else
    idx=randperm(n);
    idA=T(idx(1:N));
    idB=T(idx(N+1:end));
end
end

% ---------------------------------------------------------------------
function [y, w] = ressfilterdata(data, srate, t2bias, t2apply, ftarget, peakwidt, neighfreq, neighwidt)

Ne=size(data,1);
y=nan(length(ftarget),size(data,2),length(t2apply));
w=nan(Ne,length(ftarget));

for f=1:length(ftarget)
    
    % compute covariance matrix at peak frequency
    fdatAt = filterFGx(data(:,:,t2bias),srate,ftarget(f),peakwidt);
    fdatAt = reshape(fdatAt, Ne,[]);
    fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
    covAt  = (fdatAt*fdatAt')/size(fdatAt,2);
    
    % compute covariance matrix for lower neighbor
    fdatLo = filterFGx(data(:,:,t2bias),srate,ftarget(f)-neighfreq,neighwidt);
    fdatLo = reshape(fdatLo, Ne,[]);
    fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
    covLo  = (fdatLo*fdatLo')/size(fdatLo,2);
    
    % compute covariance matrix for upper neighbor
    fdatHi = filterFGx(data(:,:,t2bias),srate,ftarget(f)+neighfreq,neighwidt);
    fdatHi = reshape(fdatHi, Ne,[]);
    fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
    covHi  = (fdatHi*fdatHi')/size(fdatHi,2);
    
    % perform generalized eigendecomposition. This is the meat & potatos of RESS
    [evecs,evals] = eig(covAt,(covHi+covLo)/2);
    [~,comp2plot] = max(diag(evals)); % find maximum component
    evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)
    w(:,f) = evecs(:,comp2plot);
        
    % reconstruct RESS component time series
    for ti=1:length(t2apply)
        y(f,:,ti) = w(:,f)'*squeeze(data(:,:,t2apply(ti)));
    end
        
end % loop over target frequencies

end
