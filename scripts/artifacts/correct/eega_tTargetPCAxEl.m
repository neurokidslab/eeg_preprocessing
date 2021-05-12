% -------------------------------------------------------------------------
% This function uses target PCA to correct brief artifacts
% It is applied by selecting the segments sequentially for each electrode
%
% INPUT
% DATA  : data (electrode x samples x epochs)
% BCT   : logical matrix indicating bad data (electrode x samples x epochs)
% IT    : logical matrix indicating interpolation (1 x samples x epochs)
% IC    : logical matrix indicating good channels (electrode x 1 x epochs)
% nSV   : number of principal components to remove from the data
% vSV   : proportion of variance to remove from the data
%
% OPTIONAL INPUTS
% maxTime       : maximun duration of a segment to apply terget PCA
% maskTIme      : mask aroun dthe artifacts (in seconds)
% splicemethod  : how the corrected segments are splice in the data
%       - 0 = none
%       - 1 = segments are aligened with the previous segment
%       - 2 = segments are aligned in the midle of the previous and the next
%       - 3 = segments are linearly fit between the two
%       - 4 = segments before and after are roboust detrended excluding bad
%           samples, then linearly fit between the rest of the data, and 
%           finally the bad segment is linearly fit inside
%
%
% OUTPUTS
% DATA      : data 
% Interp    : logical matrix indicating interpolated data
%
% USAGE
% [EEG] = eega_tTargetPCAxEl(DATA, BCT, IT, IC, nSV, vSV, 'maxTime',0.08)
%
% -------------------------------------------------------------------------


function [DATAgood, Interp] = eega_tTargetPCAxEl( DATA, BCT, IT, IC, nSV, vSV, varargin)

if nargin<3
    IT=[];
end
if nargin<4
    IC=[];
end
if nargin<5
    nSV=[];
end
if nargin<6
    vSV=[];
end

%% ------------------------------------------------------------------------
%% Parameters
P.maxTime       = 25;
P.maskTime      = 0;
P.order         = 3; % order of polynomial to fit the trend if P.splicemethod = 4
P.silent        = 0;
P.splicemethod  = 1; % 
P.wsize         = 1000; % number of samples to mask the bad segment to detrend if P.splicemethod = 4

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_tTargetPCAxEl: Non recognized inputs')
end
if isempty(IT)
    IT=true(1,size(DATA,2),size(DATA,3));
end
if isempty(IC)
    IC=true(size(DATA,1),1,size(DATA,3));
end
if isempty(nSV) && isempty(vSV)
    vSV=0.95;
    nSV=1;
end
if isempty(nSV)
    nSV=1;
end
if isempty(vSV)
    vSV=0.001;
end
if ~any(P.splicemethod==[0 1 2 3 4])
    error('eega_tTargetPCAxEl: wrong splicing method. Options 0 | 1 | 2| 3 | 4')
end

%% ------------------------------------------------------------------------
%% Apply target PCA

nEl     = size(DATA,1);
nS      = size(DATA,2);
nEp     = size(DATA,3);

DATA    = reshape(DATA,[nEl nS*nEp]);
BCT     = reshape(BCT,[nEl nS*nEp]);
IT      = reshape(IT,[1 nS*nEp]);
Interp  = zeros([nEl nS*nEp]);
if nEp>1
    IC=true(nEl,1);
end

% normalize so all electrodes have equal variance
% -----------------------------------------------
normEl=nan(nEl,1);
for el=1:nEl
    if sum(~BCT(el,:))>2
        normEl(el) = std(DATA(el,~BCT(el,:)));
    end
end
normEl(isnan(normEl) | (normEl<1e-3))=nanmean(normEl);
DATA = DATA./repmat(normEl,[1 nS*nEp]);

% correct
% -------
DATAgood = DATA;
for el=1:nEl
    
    if IC(el)
        
        elgood = sum(IC(1:el));
        
        if ~P.silent
            fprintf('Electrode % 4.0d. ', el)
        end
        
        % find the segments to interpolate
        % --------------------------------
        bad_if = badsegments(BCT(el,:), IT, nEp, nS, P.maskTime, P.maxTime);
        
        if ~isempty(bad_if)
            
            if ~P.silent
                fprintf('Total bad data % 8d (% 7.2f%%). Data to apply PCA % 8d (% 7.2f%%).\n',...
                    sum(sum(BCT(el,:,:))),sum(sum(BCT(el,:,:)))/(nS*nEp)*100, sum(bad_if(:,2)-bad_if(:,1)), sum(bad_if(:,2)-bad_if(:,1))/(nS*nEp)*100)
            end
            
            % target PCA
            % ----------
            [d,tC] = targetpca(DATA(IC,:), bad_if, nSV, vSV, elgood);
            d = d(elgood,:);
            Interp(el,tC) = 1;
            
            % splice the segments of data together
            % ------------------------------------
            if P.splicemethod==1 % align to the previous
                epoch_if = [(1:nS:nS*nEp)' (nS:nS:nS*nEp)'];
                d = eega_splicesgments1(d, bad_if, epoch_if, [], [], []);
            elseif P.splicemethod==2 % align in the midle of the previous and the next
                epoch_if = [(1:nS:nS*nEp)' (nS:nS:nS*nEp)'];
                d = eega_splicesgments2(d, bad_if, epoch_if, [], [], []);
            elseif P.splicemethod==3 % linearly fit between the two
                epoch_if = [(1:nS:nS*nEp)' (nS:nS:nS*nEp)'];
                d = eega_splicesgments3(d, bad_if, epoch_if, [], [], []);
            elseif P.splicemethod==4 % use detrending to detrend the data before and after the artifatc and the place the artifacted data linearly alignining
                epoch_if = [(1:nS:nS*nEp)' (nS:nS:nS*nEp)'];
                d = eega_splicesgments4(d, bad_if, epoch_if, [], [], [], P.wsize, P.order);
            end
            
            % store the data
            % --------------
            DATAgood(el,:) = d;

        else
            fprintf('No segments to apply target PCA were found.\n')
        end
        
        
    else
        fprintf('Electrode % 4.0d. Bad channel.\n', el)
    end
end
% DATA(IC,:) = DATAgood;

% rescale the data back
% ---------------------
DATAgood = DATAgood.*repmat(normEl, [1 nS*nEp]);
DATA = DATA.*repmat(normEl, [1 nS*nEp]);

%% ------------------------------------------------------------------------
%% Reshape the data
DATAgood = reshape(DATAgood, [nEl nS nEp]);
Interp = logical(reshape(Interp,[nEl nS nEp]));

end


% -------------------------------------------------------------------------
function bad_if = badsegments(baddata, intertime, nEp, nS, mask, maxtime)

bad_if = [];

for ep = 1:nEp
    
    tointer_el_ep = baddata(1,(ep-1)*nS+1:ep*nS) & intertime(1,(ep-1)*nS+1:ep*nS);
    
    % begining end of each segment
    badi = [ tointer_el_ep(1) diff(tointer_el_ep)==1];
    badf = [ diff(tointer_el_ep)==-1 tointer_el_ep(end) ];
    
    % mask
    badif_ep = [find(badi)'-mask find(badf)'+mask];
    badif_ep(badif_ep(:,1)<1,1) = 1;
    badif_ep(badif_ep(:,2)>nS,2) = nS;
    
    if ~isempty(badif_ep)
        
        % put together segments which overlap
        if size(badif_ep,1)>1
            id = find(badif_ep(1:end-1,2)>=badif_ep(2:end,1));
            for iii=1:length(id)
                badif_ep(id(iii),2) = badif_ep(id(iii)+1,2);
            end
            badif_ep(id+1,:) = [];
        end
        
        % remove too long segments
        baddur_ep = badif_ep(:,2) - badif_ep(:,1) + 1;
        baddur_ep = baddur_ep - 2*mask;
        idrmv = baddur_ep > maxtime;
        badif_ep(idrmv,:) = [];
        
        % go to total samples
        badif_ep = badif_ep + (ep-1)*nS;
        
        bad_if = [bad_if; badif_ep];
    end
    
end

end

% -------------------------------------------------------------------------
function [data,tC] = targetpca(data,badseg,nSV,vSV,el)

if isempty(el)
    el=true(size(data,1),1);
end

nlim=Inf; % target pca is apply on segments of data no longer that this
tC=false(1,size(data,2));

while ~isempty(badseg)
    
    tCi=false(1,size(data,2));
    
    d=cumsum(badseg(:,2)-badseg(:,1));
    idx=d<=nlim;
    if ~any(idx), idx(1)=1; end
    badsegi=badseg(idx,:);
    badseg(idx,:)=[];
    
    % get the data to correct (and remove the mean)
    y=data(:,badsegi(1,1):badsegi(1,2))';
    y=bsxfun(@minus,y,mean(y,1));
    tCi(badsegi(1,1):badsegi(1,2))=1;
    for iseg=2:size(badsegi,1)
        yi=data(:,badsegi(iseg,1):badsegi(iseg,2))';
        yi=bsxfun(@minus,yi,mean(yi,1));
        y=cat(1,y,yi);
        tCi(badsegi(iseg,1):badsegi(iseg,2))=1;
    end
    y=bsxfun(@minus,y,mean(y,1));
    
    % pca
    covMatrix = cov(y);
    [V,S] = eig(covMatrix);
    [eigenvalues, idx] = sort(diag(S), 'descend') ;
    eigenvectors = V(:,idx);
    score = y*eigenvectors;
    expvar=var(score)';
    expvar=expvar/sum(expvar);
    
    evexp = (cumsum(expvar)>=vSV);
    vSV=find(evexp);
    vSV=vSV(1);
    
    nSV=max(nSV,vSV);
    
    ev = zeros(size(eigenvectors,1),1);
    ev(1:nSV) = 1;
    ev = diag(ev);
    
    yc = y - y*eigenvectors*ev*eigenvectors';
    
    % store the corrected data
    data(el,tCi)=permute(yc(:,el),[2 1]);
    tC(tCi)=1;
    
end

end

