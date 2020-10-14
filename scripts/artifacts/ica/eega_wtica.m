% In order to have a good ICA decomposition the length of the recording
% should be > 30 * (number of channles)^2
% This limit is used in order to determine the number of channles that are 
% considered for each ICA decomposition. All channles will be analized in
% different runs.

% See Gheetha & Geethalakshmi 2011
%   1. wavelet Coif5, level 6
%   2. identify spikes at each level
%   3. identify occular and muscular artifacted zones using coeficient of variation
%   4. apply denosing to fix the threshold value and function for the artifactes zone
%   5. inverse stationary wavelet tranform

function EEG = eega_wtica(EEG, stdlatbelsch, subspacech, varargin)

%% Optional parameters
if nargin>3 && ~isempty(varargin{1})
filthighpass=varargin{1}; else filthighpass=1;end
if nargin>4 && ~isempty(varargin{2})
filtlowpass=varargin{2}; else filtlowpass=[];end
if nargin>5 && ~isempty(varargin{3})
level=varargin{3}; else level=6;end
if nargin>6 && ~isempty(varargin{4})
mult=varargin{4}; else mult=1;end
if nargin>7 && ~isempty(varargin{5})
threshtype=varargin{5}; else threshtype='global';end  % 'global': universal global threshold | 'xlevel': applied per level
if nargin>8 && ~isempty(varargin{6})
usealphaband=varargin{6}; else usealphaband=1;end  % alpha band for MARA


%% Make copy of the data
eeg = EEG;
eeg = rmfield(eeg,'artifacts');
bc = EEG.artifacts.BC;
bt = EEG.artifacts.BT;

%% Low pass filter
if ~isempty(filtlowpass)
    eeg = eega_filter(eeg, [], filtlowpass);
end

%% High pass filter
if ~isempty(filthighpass)
    eeg = eega_filter(eeg, filthighpass, []);
end

%% Load the channels sub-space
FID = fopen(subspacech);
CH = {};
tline = fgetl(FID);
l = 1;
while ischar(tline)     
    CH = [CH; strsplit(tline,'\t','CollapseDelimiters',0)];
    tline = fgetl(FID);
    l=l+1;
end
fclose(FID);

%% New label based on the 10-20 system
FID = fopen(stdlatbelsch);
data = textscan(FID,'%s\t%s');
fclose(FID);
newLab = [data{1} data{2}];
L = {eeg.chanlocs(:).labels};
labchange = false(1,size(newLab,1));
for i=1:size(newLab,1)
    id = strcmpi(newLab{i,1},L);
    if any(id)
        eeg.chanlocs(id).labels = newLab{i,2};
        labchange(id) = 1;
    end
end
for i=1:numel(CH)
    id = strcmp(newLab(:,1),CH(i));
   if any(id)
       CH{i} = newLab{id,2};
   end
end

%% Common channels in the sub-spaces
if size(CH,2)>1
    fixChLab = CH(:,1);
    for i=1:size(CH,2)
        fixChLab = intersect(fixChLab,CH(:,i));
    end
    fixCh = false(length(eeg.chanlocs),1);
    L = {eeg.chanlocs(:).labels};
    for i=1:length(fixChLab)
        id = strcmp(fixChLab{i},L);
        if any(id)
            fixCh(id) = 1;
        end
    end
else
    fixCh = false(length(eeg.chanlocs),1);
end

%% Interpolate missing chanells from the common channels
if any(fixCh & bc)
    d = eega_tInterpSpatial( eeg.data, ~bc, eeg.chanlocs, 0.35);
    eeg.data(fixCh & bc,:,:) = d(fixCh & bc,:,:);
    bc(fixCh & bc,:,:) = 0;
end

%% Remove bad samples
eeg.data = eeg.data(:,~bt,:);
eeg.pnts = size(eeg.data,2);
eeg.times = eeg.times(~bt);
eeg.xmin = eeg.times(1);
eeg.xmax = eeg.times(end);

%% Remove bad channels
goodch = find(~bc);
eeg = pop_select( eeg,'channel', goodch);
eeg = eeg_checkset( eeg );

%% See if the number of channels to analyse per time is not too high
nS = size(eeg.data,2)*size(eeg.data,3);
nICAlim = floor(sqrt(nS / 30));
if size(CH,1)>nICAlim
    warning('The number of components should be smaller than %d',nICAlim); 
end

%% Re-arrange the channels 
L = {eeg.chanlocs(:).labels};
chxround = cell(1,size(CH,2));
for i=1:size(CH,2)
    [ch, iL, ~]  = intersect(upper(L),upper(CH(:,i)));
    chxround{i} = iL;
end

%% Compute ICA for each channels sub-space and remove the artifacts
nrounds = length(chxround);
ICA = struct();
for i=1:nrounds
    
    %% Take the subspace channels 
    eegi = pop_select( eeg,'channel', chxround{i});
    eegi = eeg_checkset( eegi );
    
    %% Channels sub-space in terms of the original EEG structure
    chi = false(size(EEG.data,1),1);
    chi(goodch(chxround{i})) = 1;
    
    %% First ICA
    [OUTEEG, ~] = pop_runica(eegi, 'extended',1,'interupt','on'); %runica for parametric, default extended for finding subgaussian distributions
    W = OUTEEG.icaweights*OUTEEG.icasphere;
    A = inv(W);
    IC = reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    clear OUTEEG
    
    %% Wavelet-thresholding on the IC to remove big artifacts
    try
        [wIC] = eega_wt_thresh_xlevel(IC, mult, level, threshtype, 's');
    catch wica_err
        rethrow(wica_err)
    end
        
    %reconstruct artifact signal as channelsxsamples format from the wavelet coefficients
    artifacts = A*wIC;
    artifacts = reshape(artifacts,size(eegi.data));

    %subtract out wavelet artifact signal from EEG signal
    eegi.data = eegi.data-artifacts;   
    
    %% Second ICA on the clean data
    eegi = pop_runica(eegi, 'extended',1,'interupt','on');
    
    %% MARA to identify components corrisponding to artifacts
    eegi = eega_icarunmara(eegi, 0, usealphaband);
%     [~,eegi,~] = processMARA_babies( eegi, eegi, eegi, usealphaband,[0, 0, 1, 1, 1] );
    cmp2rmv = find(eegi.reject.gcompreject == 1);
    
    [~, varrmv] = compvar(eegi.data, eegi.icaact, eegi.icawinv, cmp2rmv);
    perccmprmv = length(cmp2rmv)/size(eegi.data,1)*100;
    fprintf('IC removed by MARA ......... %4.2f%% (%d out of %d) \n', perccmprmv,length(cmp2rmv),size(eegi.data,1))
    fprintf('Variance removed by MARA ... %4.2f%% \n',varrmv)
        
    %% store the ICa decomposition matrixes 
    ICA(i).filthighpass = filthighpass;
    ICA(i).filtlowpass = filtlowpass;
    ICA(i).ch = chi;
    ICA(i).smpls = ~bt;
    ICA(i).icaweights = eegi.icaweights;
    ICA(i).icasphere = eegi.icasphere;
    ICA(i).cmp2rmv = cmp2rmv;
    ICA(i).perccmprmv = perccmprmv;
    ICA(i).varrmv = varrmv;
    
end

%% Store the ICA decomposition and rejection data
EEG.ica = ICA;

end
