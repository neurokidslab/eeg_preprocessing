%function [wIC] = wICA(IC, mult=1, L=5, threshtype='global', sorh='s', applyonlevel=(1:L), aproxpart=0)
function [wIC] = eega_wt_thresh_xlevel(IC, varargin)
%--------------- function [wIC,A,W] = wICA(data,varargin) -----------------
%
% Wavelet thresholding to remove low-amplitude activity from the computed ICs.
% This is useful for extracting artifact-only ICs in EEG (for example), and
% then subtracting the artifact-reconstruction from the original data. 
%
%               >>> INPUTS >>>
% Required: 
%   IC = data matrix in row format
% Optional:
%   mult = threshold multiplier...multiplies the computed threshold from
%       "ddencmp" by this number. Higher thresh multipliers = less
%       "background" (or low amp. signal) is kept in the wICs.
%   L = level set for stationary wavelet transform. Higher levels give
%       better frequency resolution, but less temporal resolution. 
%       Default = 5
%
%               <<< OUTPUTS <<<
%   ICw = non-wavelet ICs (optional)
%   
%       * you can reconstruct the artifact-only signals as:
%               artifacts = A*wIC;
%       - upon reconstruction, you can then subtract the artifacts from your
%       original data set to remove artifacts, for instance.
%
% The IC are transfomed into the wavelet space using the stationary
% wavelet transform. A threshold is estimated within each level and the
% coeficients smaller than the threshold are reduce usign soft thresholding
% The obtained components should be mainly dominated by artifacts. If they 
% are transformed back to data space the artifacts are obtained and can be
% remove from the data before applying ICA another time
% It uses "coif5" as wavelet family 
%
% By Ana Flo June 23th 2020
%---------------------------------------------------------------------------------------

% check inputs
if nargin>1 && ~isempty(varargin{1})
mult=varargin{1};else mult=1;end
if nargin>2 && ~isempty(varargin{2})
L=varargin{2}; else L=5;end
if nargin>3 && ~isempty(varargin{3})
threshtype=varargin{3}; else threshtype='global';end  % 'global': universal global threshold | 'xlevel': applied per level
if nargin>4 && ~isempty(varargin{4})
sorh=varargin{4}; else sorh='s';end  %% soft thresholding
if nargin>5 && ~isempty(varargin{5})
applyonlevel=varargin{5}; else applyonlevel=(1:L);end  %% soft thresholding
if nargin>6 && ~isempty(varargin{6})
aproxpart=varargin{6}; else aproxpart=0;end  %% remove the aproximation part

% padding data for proper wavelet transform...data must be divisible by
% 2^L, where L = level set for the stationary wavelet transform
modulus = mod(size(IC,2),2^L); %2^level (level for wavelet)
if modulus ~=0
    extra = zeros(1,(2^L)-modulus);
else
    extra = [];
end
      
% loop through ICs and perform wavelet thresholding
for s = 1:size(IC,1)
    if ~isempty(extra)
        sig = [IC(s,:),extra]; % pad with zeros
    else
        sig = IC(s,:);
    end
    
    % use stationary wavelet transform (SWT) to wavelet transform the ICs
    swc = eega_swt(sig,L);
    
    % compute the thresholds
    thresh = estimatethresh(swc, threshtype);
    if length(thresh)==L
        thresh = thresh(applyonlevel);
    end
    
    % threshold the wavelet to remove small values
    thresh = thresh*mult;   % multiply threshold by scalar
    Y = nan(size(swc));
    Y(applyonlevel,:) = wthreshold(swc(applyonlevel,:),sorh,thresh);
    if aproxpart==0
        Y(end,:) = zeros(1,size(Y,2));  % the aproximation part is completely removed
    else
        Y(end,:) = swc(end,:);
    end
    
    % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
    wIC(s,:) = eega_iswt(Y);

    clear Y sig thresh swc 
end

% remove extra padding
if ~isempty(extra)
    wIC = wIC(:,1:end-numel(extra));
end

% % plot the ICs vs. wICs
% if plotting>0
%     Fs = 1;
%     disp('Plotting');
%     figure,
%     subplot(3,1,1);
%         multisignalplot(IC,Fs,'r');
%         title('ICs');
%     subplot(3,1,2);
%         multisignalplot(wIC,Fs,'r');
%         title('wICs')
%     subplot(3,1,3);
%         multisignalplot(IC-wIC,Fs,'r');
%         title('Difference (IC - wIC)');
% end

end

% -------------------------------------------------------------------------
% AUXILIARY FUNCTIONS

function y = wthreshold(x,sorh,t)
%WTHRESH Perform soft or hard thresholding. 
%   Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's')
%   or hard (if SORH = 'h') T-thresholding  of the input 
%   vector or matrix X. T is the threshold value.
%
%   Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft 
%   thresholding is shrinkage.
%
%   Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard
%   thresholding is cruder.
%
%   See also WDEN, WDENCMP, WPDENCMP.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 24-Jul-2007.
%   Copyright 1995-2018 The MathWorks, Inc.

if isStringScalar(sorh)
    sorh = convertStringsToChars(sorh);
end
if numel(t)==1
    t=repmat(t,[1 size(x,1)]);
end

switch sorh
    case 's'
        y = nan(size(x));
        for i=1:size(x,1)
            tmp = (abs(x(i,:))-t(i));
            tmp = (tmp+abs(tmp))/2;
            y(i,:) = sign(x(i,:)).*tmp;
        end
        
    case 'h'
        y = nan(size(x));
        for i=1:size(x,1)
            y(i,:) = x(i,:).*(abs(x(i,:))>t(i));
        end
    otherwise
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'))
end
end

function thresh = estimatethresh(x, type)

switch type
    
    case 'global'
        % using Donoho,Universal threshold
        normaliz = median(abs(x(1,:)))/0.6745;
        n = size(x,2);
        thr = sqrt(2*log(n));
        thresh = thr*normaliz;
        
    case 'xlevel'
        % using Jhonstone and Silverman 1997 (threshold for each level)
        L = size(x,1);
        n = size(x,2);
        normaliz = nan(L,1);
        for j=1:L
            normaliz(j) = median(abs(x(j,:)))/0.6745;
        end
        thr = sqrt(2*log(n));
        thresh = thr*normaliz;
end
end