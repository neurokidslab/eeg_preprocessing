% In order to have a good ICA decomposition the length of the recording
% should be > 30 * (number of channles)^2
% This limit is used in order to determine the number of channles that are
% considered for each ICA decomposition. All channles will be analized in
% different runs.

% See Geetha & Geethalakshmi 2011
%   1. wavelet Coif5, level 6
%   2. identify spikes at each level
%   3. identify occular and muscular artifacted zones using coeficient of variation
%   4. apply denosing to fix the threshold value and function for the artifactes zone
%   5. inverse stationary wavelet tranform
%
% EEG = eega_pcawtica(EEG, 'opt1', val1,...)
function EEG = eega_pcawtica(EEG, varargin)

%% Optional parameters

% Default parameters
P.subspacech    = [];
P.filthighpass  = 2;        % High pass filter applied before ICA. Default 2
P.filtlowpass   = 40;       % Low pass filter applied before ICA. Default 40
P.npc           = 0;        % Number of components to keep in the PCA. If 0 it is not applied. Default 0
P.level         = 6;
P.mult          = 1;
P.threshtype    = 'xlevel'; % 'xlevel' | 'global'
P.classifyIC    = 1;        % Apply an authomatic classification of IC components. Default (1)
P.classifyICfun = 'iMARA';  % Function use to classify IC. Default 'iMARA'
P.changelabelch = 0;        % Change the labels of the channels to be conistent with the classification algorithm. Default False
P.labelch       = [];       % Cell of size n x 2 with the labels in chanlocs and the new names. It can also be the name of a text file with the old and new names. If empty, the default is the convertion from a EGI 129 layout
P.rmvart        = 1;        % 0= do not remove; 1= remove IC; 2= remove IC and wavelet thresholding
P.rmvartfilt    = 1;        % Estimate the artifacts based on filtered data (1) or on original data (0). Default 1 
P.saveica       = 1;        % Save a file with the ICA weights. Default True
P.icaname       = 'ica_';   % Name added to save the ICA decomposition. Default 'ica_'
P.icapath       = [];       % Path were to to save the ICA decomposition. By default in EEG.filepath
P.artthreshold = 0.9;

% get the optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_pcawtica: Non recognized inputs')
end

% If an authomatic ICA rejection has to be applied check that it is install
if P.classifyIC
    if exist(P.classifyICfun,'file')~=2
		switch P.classifyICfun
			case 'iMARA'
				error('The %s toolbox is not in the path. Install the plugin (https://github.com/Ira-marriott/iMARA/)',P.funmara)
			case 'MARA'
				error('The %s toolbox is not in the path. Install the plugin (https://github.com/irenne/MARA/)',P.funmara)
			otherwise
				error('The %s toolbox is not in the path.',P.funmara)
		end
    end
end

%% Compute the grand sum-squared data
var_org  = sum(sum(EEG.data(~EEG.artifacts.BC,~EEG.artifacts.BT).^2));

%% Make copy of the data
eeg = EEG;
eeg = rmfield(eeg,'artifacts');
bc = EEG.artifacts.BC;
bt = EEG.artifacts.BT;
bct = EEG.artifacts.BCT;
goodch = find(~bc);

%% Check that there is any good data
gooddata = 1;
if all(bc) || all(bt)
    gooddata = 0;
end

%% If there is good data run ICA
if gooddata
    
    %% 1. Low pass filter
    if ~isempty(P.filtlowpass)
        minphase = 0;
        eeg = pop_eegfiltnew(eeg, [], P.filtlowpass,  [], 0, [], [], minphase);
    end
    
    %% 2. High pass filter
    if ~isempty(P.filthighpass)
        minphase = 0;
        eeg = pop_eegfiltnew(eeg, P.filthighpass, [],  [], 0, [], [], minphase);
    end
    
    %% 3. Steps relevant only if ICA has to be applied on a sub-group of channels
	
	%  3.1 Load the channels sub-space 
    if ~isempty(P.subspacech)
        FID = fopen(P.subspacech);
        CH = {};
        tline = fgetl(FID);
        l = 1;
        while ischar(tline)
            CH = [CH; strsplit(tline,'\t','CollapseDelimiters',0)];
            tline = fgetl(FID);
            l=l+1;
        end
        fclose(FID);
    else
        CH = {EEG.chanlocs(:).labels}';
    end
    
    %  3.2 Search for the common channels in the sub-spaces
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
    
    %  3.3 Interpolate missing channels from the common channels
    if any(fixCh & bc)
        d = eega_tInterpSpatial( eeg.data, ~bc, eeg.chanlocs, 1);
        eeg.data(fixCh & bc,:,:) = d(fixCh & bc,:,:);
        bc(fixCh & bc,:,:) = 0;
    end
    
	%% 4. Remove bad data before appliyng ICA to obtain a better decomposition
	
    %  4.1 Remove bad samples
    eeg.data = eeg.data(:,~bt,:);
    eeg.pnts = size(eeg.data,2);
    eeg.times = eeg.times(~bt);
    eeg.xmin = eeg.times(1);
    eeg.xmax = eeg.times(end);
    
    %  4.2 Remove bad channels
    if sum(~bc)<EEG.nbchan
        eeg = pop_select( eeg,'channel', find(~bc));
        eeg = eeg_checkset( eeg );
    end
    
    %  4.3 Set as zero the remaining bad data
    eeg.data(bct(~bc,~bt)) = 0;
    
    %% 5. Check if the number of channels to analyse per time is not too high relative to the leght of the recording
    nS = size(eeg.data,2)*size(eeg.data,3);
    nICAlim = floor(sqrt(nS / 30));
    war=0;
    if (P.npc==0 && size(CH,1)>nICAlim); war=1;
    elseif (P.npc~=0 && P.npc>nICAlim); war=1;
    end
    if war
        warning('The number of components should be smaller than %d',nICAlim);
    end
    
    %% 6. Re-arrange the channels labels relative to the channel sub-space
    L = {eeg.chanlocs(:).labels};
    chxround = false(length(L),size(CH,2));
    for i=1:size(CH,2)
        [ch, indL, indCHi]  = intersect(upper(L),upper(CH(:,i)));
        chxround(indL,i) = 1;
    end
    
    %% 7. Compute ICA for each channels sub-space and remove the artifacts
    nrounds = size(CH,2);
    ICA = struct();
    varrmv = nan(1,nrounds);
    icrmv = nan(1,nrounds);
    dataart = zeros(size(EEG.data));
    N = zeros(size(EEG.data));
    for i=1:nrounds
        
		% 7.1 Prepare the data 
		
        %channels sub-space in terms of the original EEG structure
        chi = goodch(chxround(:,i));
        
        %samples subspace
        smplsi = find(~bt);
        
        %take and reshape the data (channels x samples)
        tmpdata = eeg.data(chxround(:,i),:,:);
        tmpdata = reshape( tmpdata, sum(chxround(:,i)), eeg.pnts*eeg.trials);
		
		%zero mean
        tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); 
        
        %  7.2 First PCA+ICA
        if P.npc~=0
            [eigenvec, eigenval, score] = dopca(tmpdata');
            explained = eigenval/sum(eigenval)*100;
            fprintf('The first %d principal components carying  %4.2f%% of the varianceare are kept... \n', P.npc, sum(explained(1:P.npc)))
            fprintf('Performing first ICA decomposition... \n')
            [icaweights,icasphere] = runica( score(1:P.npc,:), 'lrate', 0.001, 'extended',1 );
            W = icaweights*icasphere;
            A = pinv(W);
            IC = W * score(1:P.npc,:);
        else
            fprintf('Performing first ICA decomposition... \n')
            [icaweights,icasphere] = runica( tmpdata, 'lrate', 0.001, 'extended',1 );
            W = icaweights*icasphere;
            A = pinv(W);
            IC = W * tmpdata;
        end
        
        %  7.3 Wavelet-thresholding on the IC to remove big artifacts
        fprintf('Performing wavelet-thresholding... \n')
        
        %wavelet decomposition and thresholding
        [wIC] = eega_wt_thresh_xlevel(IC, P.mult, P.level, P.threshtype, 's');
        
        %reconstruct artifact signal from the wavelet coefficients
        if P.npc~=0
            artifacts = eigenvec(:,1:P.npc)*A*wIC;
        else
            artifacts = A*wIC;
        end
        
        clear eigenvec eigenval score explained icaweights icasphere W A IC
        
        %  7.4 Take the data again and substract the artifact signal
        tmpdata = eeg.data(chxround(:,i),:,:);
        tmpdata = reshape( tmpdata, sum(chxround(:,i)), eeg.pnts*eeg.trials);
        tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); %zero mean
        tmpdata = tmpdata - artifacts;
        
        %  7.5 Second PCA+ICA on the clean data
        if P.npc~=0
            [eigenvec, eigenval, score] = dopca(tmpdata');
            explained = eigenval/sum(eigenval)*100;
            %         tmpdata = eigenvec(:,1:P.npc) * score(1:P.npc,:);
            fprintf('The first %d principal components carying  %4.2f%% of the varianceare are kept... \n', P.npc, sum(explained(1:P.npc)))
            fprintf('Performing second ICA decomposition on denoised data... \n')
            [icaweights,icasphere] = runica( score(1:P.npc,:), 'lrate', 0.001, 'extended',1 );
            W = icaweights*icasphere*eigenvec(:,1:P.npc)';
            S = eye(size(tmpdata,1));
        else
            explained = [];
            fprintf('Performing second ICA decomposition on denoised data... \n')
            [icaweights,icasphere] = runica( tmpdata, 'lrate', 0.001, 'extended',1 );
            W = icaweights;
            S = icasphere;
        end
        
        %  7.6 authomatically identify components corrisponding to artifacts
        if P.classifyIC
            fprintf('Applying %s to identify artifacts IC... \n',P.classifyICfun)
            
            %create the eeglab structure
            eegi.setname        = eeg.setname;
            eegi.filename       = eeg.filename;
            eegi.filepath       = eeg.filepath;
            eegi.data           = tmpdata;
            eegi.nbchan         = size(eegi.data,1);
            eegi.trials         = size(eegi.data,3);
            eegi.pnts           = size(eegi.data,2);
            eegi.srate          = eeg.srate;
            eegi.times          = eeg.times;
            eegi.xmin           = eegi.times(1);
            eegi.xmax           = eegi.times(end);
            eegi.chanlocs       = eeg.chanlocs(chxround(:,i));
            eegi.icaweights     = W;
            eegi.icasphere      = S;
            eegi.icawinv        = pinv( eegi.icaweights*eegi.icasphere );
            eegi.icaact         = eegi.icaweights *eegi.icasphere * eegi.data;
            eegi.icachansind    = (1:size(eegi.data,1));
            eegi.etc            = eeg.etc;
            
            %change the labels of the channels to agree with the MARA labels.
            if P.changelabelch
                eegi.chanlocs = eega_changechlable(eegi.chanlocs,P.labelch);
            end
            
            %run MARA/iMARA
            switch P.classifyICfun
                case 'MARA'
                    [artcomps, info] = MARA(eegi);
                case 'iMARA'
                    [artcomps,info] = iMARA(eegi,P.artthreshold);
                case 'MARA_babies'
                    usealphaband = 0;
                    [artcomps, info] = MARA_babies(eegi,usealphaband);
            end
            eegi.reject.gcompreject = false(1,eegi.nbchan);
            eegi.reject.gcompreject(artcomps) = 1;
            cmp2rmv = find(eegi.reject.gcompreject == 1);
            
            %percentage of components removed
            ndat = size(W,1);
            icrmv(i) = 100*length(cmp2rmv)/ndat;
            fprintf('IC removed ......... %4.2f%% (%d out of %d) \n', icrmv(i),length(cmp2rmv),ndat)
            
            %removed variance
            [~, varrmv(i)] = compvar(eegi.data, eegi.icaact, eegi.icawinv, cmp2rmv);
            fprintf('Variance removed ... %4.2f%% \n',varrmv(i))
            
        else
            cmp2rmv = [];
            icrmv(i) = 0;
            varrmv(i) = 0;
        end
        
        %  7.7 Estimate the artifacts
        if P.rmvart~=0
            
            if P.rmvartfilt
                %remove the components identied as artifacts
                if P.npc~=0 && length(cmp2rmv)==P.npc~=0
                    eegirmv.data = zeros(size(eegi.data));
                elseif P.npc==0 && length(cmp2rmv)==size(eegi.data,1)
                    eegirmv.data = zeros(size(eegi.data));
                elseif ~isempty(cmp2rmv)
                    eegirmv = pop_subcomp( eegi, cmp2rmv, 0);
                else
                    eegirmv = eegi;
                end
                
                %store the artifacts
                if P.rmvart==1
                    dataart(chi,smplsi) = dataart(chi,smplsi) + (eegi.data - eegirmv.data);
                    N(chi,smplsi) = N(chi,smplsi)+1;
                elseif P.rmvart==2
                    dataart(chi,smplsi) = dataart(chi,smplsi) + (eegi.data - eegirmv.data) + artifacts;
                    N(chi,smplsi) = N(chi,smplsi)+1;
                end
                
                clear eegi eegirmv
        
            else
                clear eegi
                
                EEGi = pop_select( EEG,'channel', find(~bc));
                EEGi = eeg_checkset(EEGi);
                EEGi.icaweights     = W;
                EEGi.icasphere      = S;
                EEGi.icawinv        = pinv( EEGi.icaweights*EEGi.icasphere );
                EEGi.icaact         = EEGi.icaweights *EEGi.icasphere * EEGi.data;
                EEGi.icachansind    = (1:size(EEGi.data,1));
                EEGi.reject.gcompreject = false(1,EEGi.nbchan);
                EEGi.reject.gcompreject(artcomps) = 1;
                EEGi = eeg_checkset(EEGi);
                
                %remove the components identied as artifacts
                if P.npc~=0 && length(cmp2rmv)==P.npc~=0
                    EEGirmv.data = zeros(size(EEGirmv.data));
                elseif P.npc==0 && length(cmp2rmv)==size(EEGirmv.data,1)
                    EEGirmv.data = zeros(size(eegirmv.data));
                elseif ~isempty(cmp2rmv)
                    EEGirmv = pop_subcomp( EEGi, cmp2rmv, 0);
                else
                    EEGirmv = EEGi;
                end
                
                %store the artifacts
                if P.rmvart==1
                    dataart(chi,smplsi) = dataart(chi,smplsi) + (EEGi.data(:,smplsi) - EEGirmv.data(:,smplsi));
                    N(chi,smplsi) = N(chi,smplsi)+1;
                elseif P.rmvart==2
                    dataart(chi,smplsi) = dataart(chi,smplsi) + (EEGi.data(:,smplsi) - EEGirmv.data(:,smplsi)) + artifacts;
                    N(chi,smplsi) = N(chi,smplsi)+1;
                end
                
                clear EEGi EEGirmv
        
            end
                    
        end
        
        %  7.8 store the ICa decomposition matrixes
        ICA(i).filthighpass = P.filthighpass;
        ICA(i).filtlowpass = P.filtlowpass;
        ICA(i).ch = chi;
        ICA(i).smpls = smplsi;
        ICA(i).artifacts = artifacts;
        ICA(i).pcanpc = P.npc;
        ICA(i).pcaexplained = sum(explained(1:P.npc));
        ICA(i).icaweights = W;
        ICA(i).icasphere = S;
        ICA(i).cmp2rmv = cmp2rmv;
        ICA(i).icrmv = icrmv(i);
        ICA(i).varrmv = varrmv(i);
        
        
    end
    
    %% 8. Remove the artifacts
    if P.rmvart~=0
        
        %divide the artifacts for the number of estimations
        dataart(EEG.artifacts.BC,:) = 0;
        N(N==0) = 1;
        dataart = dataart ./ N;
        
        %remove the artifacts extracted by ICA
        EEG.data = EEG.data - dataart;
        
        %compute the grand sum-squared data for the clean data
        var_clean  = sum(sum(EEG.data(~EEG.artifacts.BC,~EEG.artifacts.BT).^2));
        varrmvtot = 100*(1- var_clean/var_org);
        EEG.icainfo.varrmvtot = varrmvtot;
        EEG.icainfo.varrmv = varrmv;
        EEG.icainfo.icrmv = icrmv;
        fprintf('--> Total variance removed: %4.2f%% \n\n',varrmvtot)
        
    end
    
    %% 9. Save the ICA decomposition
    if P.saveica
        [~,name,~] = fileparts(EEG.filename);
        if isempty(P.icapath)
            filename = fullfile(EEG.filepath, [P.icaname name '.mat']);
        else
            filename = fullfile(P.icapath, [P.icaname name '.mat']);
        end
        save(filename, 'ICA')
    end
else
    warning('No good data to apply ICA. Nothing done\n')
end

end

function [eigenvec, eigenval, score] = dopca(d)
[m,~] = size(d);
sigma = (1/m) * (d' * d);
[eigenvec,S] = eig(sigma);
[eigenval,index] = sort(diag(S));
index = rot90(rot90(index));
eigenval = rot90(rot90(eigenval))';
eigenvec = eigenvec(:,index);
score = eigenvec'*d';
end

function tmprank2 = getrank(tmpdata);

tmprank = rank(tmpdata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here: alternate computation of the rank by Sven Hoffman
%tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
covarianceMatrix = cov(tmpdata', 1);
[E, D] = eig (covarianceMatrix);
rankTolerance = 1e-7;
tmprank2=sum (diag (D) > rankTolerance);
if tmprank ~= tmprank2
    fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
    tmprank2 = max(tmprank, tmprank2);
end
end