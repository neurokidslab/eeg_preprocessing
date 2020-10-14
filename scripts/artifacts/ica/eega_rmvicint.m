function EEG = eega_rmvicint(EEG, filthighpass, filtlowpass)
if nargin<2 
    filthighpass=[];
end
if nargin>3 
    filtlowpass=[];
end

% compute the grand sum-squared data
var_org  = sum(sum(EEG.data(~EEG.artifacts.BC,~EEG.artifacts.BT).^2)); 

% do a copy
eeg = EEG;

% Low pass filter
if ~isempty(filtlowpass)
    eeg = eega_filter(eeg, [], filtlowpass);
end

% High pass filter
if ~isempty(filthighpass)
    eeg = eega_filter(eeg, filthighpass, []);
end

% interpolate the bad channels
% databad = eeg.data(eeg.artifacts.BC,:,:);
eeg.data = eega_tInterpSpatial( eeg.data, ~eeg.artifacts.BC, eeg.chanlocs, 0.35);
    
% loop over ica decompositions
nrounds = length(eeg.ica);
dataclean = zeros(size(eeg.data));
dataart = zeros(size(eeg.data));
N = zeros(size(eeg.data));
for i=1:nrounds

    %data restricted to the channels subset
    chi = eeg.ica(i).ch;
    eegi = pop_select( eeg,'channel', find(chi));
    eegi = eeg_checkset( eegi );
    
    %data restricted to the samples
    smplsi = eeg.ica(i).smpls;
    eegi.data = eegi.data(:,smplsi,:);
    eegi.pnts = size(eegi.data,2);
    eegi.times = eegi.times(smplsi);
    eegi.xmin = eegi.times(1);
    eegi.xmax = eegi.times(end);

    %ICA information
    eegi.icaweights     = eeg.ica(i).icaweights;
    eegi.icasphere      = eeg.ica(i).icasphere;
    eegi.icawinv        = pinv( eegi.icaweights*eegi.icasphere );
    eegi.icaact         = eegi.icaweights *eegi.icasphere * eegi.data;
    eegi.icachansind    = (1:size(eegi.data,1));
    cmp2rmv             = eeg.ica(i).cmp2rmv;
    
    %percentage of components removed
    icrmv = 100*length(cmp2rmv)/size(eegi.data,1);
    fprintf('IC removed ......... %4.2f%% (%d out of %d) \n', icrmv,length(cmp2rmv),size(eegi.data,1))
    
    %removed variance
    [~, varrmv] = compvar(eegi.data, eegi.icaact, eegi.icawinv, cmp2rmv); 
    fprintf('Variance removed ... %4.2f%% \n',varrmv)
    
    %remove the components identied as artifacts
    dorg = eegi.data;
    if length(cmp2rmv)==size(eegi.data,1)
        eegi.data = zeros(size(eegi.data));
    elseif ~isempty(cmp2rmv)
        eegi = pop_subcomp( eegi, cmp2rmv, 0);
    end
    dart = dorg - eegi.data; % the artifatcs removed
    
%     % interpolate the artifacts for the rest of the channels 
%     D = zeros(size(eeg.data,1),sum(smplsi));
%     D(chi,:) = dart;
%     D = eega_tInterpSpatial( D, chi, eeg.chanlocs, 1);
%     
%     %store the artifacts  
%     dataart(:,smplsi) = dataart(:,smplsi) + D;
%     N(:,smplsi) = N(:,smplsi)+1;
    
    % interpolate the rest of the channels 
    d = eeg.data;
    d(chi,smplsi) = eegi.data;
    d = eega_tInterpSpatial( d, chi, eeg.chanlocs, 1);
    
    %store the clean data  
    dataclean(:,smplsi) = dataclean(:,smplsi) + d(:,smplsi);
    N(:,smplsi) = N(:,smplsi)+1;
      
end

% % divide the artifacted data for the number of estimations
% dataart(EEG.artifacts.BC,EEG.artifacts.BT) = 0;
% N(N==0) = 1;
% dataart = dataart ./ N;
% 
% % remove the artifacts from the original data
% EEG.data = EEG.data - dataart;

% divide the clean data for the number of estimations
dataclean(dataclean==0) = nan;
dataclean(EEG.artifacts.BC,:) = nan;
dataclean = dataclean ./ N;

% remove the artifacts from the original data
idx = ~isnan(dataclean);
EEG.data(idx) = EEG.data(idx) - ( eeg.data(idx) - dataclean(idx) );
% EEG.data(~isnan(dataclean)) = dataclean(~isnan(dataclean));
% EEG.data(EEG.artifacts.BC,:) = databad;

% compute the grand sum-squared data for the clean data
var_clean  = sum(sum(EEG.data(~EEG.artifacts.BC,~EEG.artifacts.BT).^2)); 
varrmv = 100*(1- var_clean/var_org);
EEG.icainfo.varrmv = varrmv;
fprintf('--> Total variance removed: %4.2f%% \n\n',varrmv)

end