% this function uses the info in ICA to remove the artifacts
% It was not written to work on epoched data...
%
% Ana Flo July 2020

function EEG = eega_rmvic(EEG, ICA, varargin)

%% Optional parameters

% Default parameters
P.applymara     = 1;
P.usealphaband  = 1;
P.stdlatbelsch  = [];
P.rmvart        = 1;  % 0= do not remove; 1= remove IC; 2= remove IC and wavelet thresholding

% get the optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_rmvic: Non recognized inputs')
end

%% Compute the grand sum-squared data
var_org  = sum(sum(EEG.data(~EEG.artifacts.BC,~EEG.artifacts.BT).^2)); 
    
%% Loop over ica decompositions
nrounds = length(ICA);
dataart = zeros(size(EEG.data));
N = zeros(size(EEG.data));
varrmv = nan(1,nrounds);
icrmv = nan(1,nrounds);
for i=1:nrounds

    %channels subset
    chi = ICA(i).ch;
    
    %samples subset
    smplsi = find(~EEG.artifacts.BT);
    useart = 1;
    if length(intersect(smplsi, ICA(i).smplsi))~=length(smplsi)
        warning('Artifacts obtained from wavelet-thresholding will not be taken into account!')
        useart = 0;
    end
    
    %take the data to remove the components 
    tmpdata = eeg.data(chi,smplsi,:);
    tmpdata = reshape( tmpdata, length(chi), []);
    
    %interpolate bad channels
    bad = EEG.artifacts.BC(chi);
    tmpdata = eega_tInterpSpatial( tmpdata, ~bad, EEG.chanlocs(chi), 1);
    
    %zero mean
    tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); 
    
    %remove the wavelet-thresholding artifacts
    if useart
        tmpdata = tmpdata - ICA(i).artifacts;
    end
    
    %create the eeglab structure
    eegi.setname        = EEG.setname;
    eegi.filename       = EEG.filename;
    eegi.filepath       = EEG.filepath;
    eegi.data           = tmpdata;
    eegi.nbchan         = size(eegi.data,1);
    eegi.trials         = size(eegi.data,3);
    eegi.pnts           = size(eegi.data,2);
    eegi.srate          = EEG.srate;
    eegi.times          = (1:eegi.pnts)/eegi.srate;
    eegi.xmin           = eegi.times(1);
    eegi.xmax           = eegi.times(end);
    eegi.chanlocs       = EEG.chanlocs(chi);
    eegi.icaweights     = ICA(i).icaweights;
    eegi.icasphere      = ICA(i).icasphere;
    eegi.icawinv        = pinv( eegi.icaweights*eegi.icasphere );
    eegi.icaact         = eegi.icaweights *eegi.icasphere * eegi.data;
    eegi.icachansind    = (1:size(eegi.data,1));
    cmp2rmv             = ICA(i).cmp2rmv;
    
    %apply MARA
    if applymara
        %change the labels of the channels to agree with the MARA labels
        eegi.chanlocs = eega_changechlable(eegi.chanlocs,stdlatbelsch);
        
        %run MARA
        eegi = eega_icarunmara(eegi, 0, usealphaband);
        cmp2rmv = find(eegi.reject.gcompreject == 1);  % [~,eegi,~] = processMARA_babies( eegi, eegi, eegi, usealphaband,[0, 0, 1, 1, 1] );
        ICA(i).cmp2rmv = cmp2rmv;
    end
    
    %percentage of components removed
    ndat = size(eegi.icaweights,1);
    icrmv(i) = 100*length(cmp2rmv)/ndat;
    fprintf('IC removed ......... %4.2f%% (%d out of %d) \n', icrmv(i),length(cmp2rmv),ndat)
    
    %removed variance
    [~, varrmv(i)] = compvar(eegi.data, eegi.icaact, eegi.icawinv, cmp2rmv); 
    fprintf('Variance removed ... %4.2f%% \n',varrmv(i))
    
    %remove the components identied as artifacts
    if length(cmp2rmv)==size(eegi.data,1)
        eegirmv.data = zeros(size(eegi.data));
    elseif ~isempty(cmp2rmv)
        eegirmv = pop_subcomp( eegi, cmp2rmv, 0);
    end
    
    %store the artifacts 
    if useart && (P.rmvart==2)
        dataart(chi,smplsi) = dataart(chi,smplsi) + (eegi.data - eegirmv.data) + ICA(i).artifacts;
        N(chi,smplsi) = N(chi,smplsi)+1;
    else
        dataart(chi,smplsi) = dataart(chi,smplsi) + (eegi.data - eegirmv.data);
        N(chi,smplsi) = N(chi,smplsi)+1;
    end
    
end

%% divide the artifacts for the number of estimations
dataart(EEG.artifacts.BC,:) = 0;
N(N==0) = 1;
dataart = dataart ./ N;

%% remove the artifacts from the original data
EEG.data = EEG.data - dataart;

%% compute the grand sum-squared data for the clean data
var_clean  = sum(sum(EEG.data(~EEG.artifacts.BC,~EEG.artifacts.BT).^2)); 
varrmvtot = 100*(1- var_clean/var_org);
EEG.icainfo.varrmvtot = varrmvtot;
EEG.icainfo.varrmv = varrmv;
EEG.icainfo.icrmv = icrmv;
fprintf('--> Total variance removed: %4.2f%% \n\n',varrmvtot)

%% Remove the field ICA
EEG = rmfield(EEG,'ica');

end