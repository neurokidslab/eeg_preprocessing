function EEG = eega_spatialsmoothing(EEG, varargin)

fprintf('### Spatial smoothing of the data ###\n')

%% ------------------------------------------------------------------------
%% Parameters

P.nbvoisins     = 5; % Number of neighbouring electrodes to consider
P.alpha         = 0.5; % Weight given to the original channel and  1-alpha to the neighbours
P.ref           = 1; % SI the data average referenced?
P.chmapfile     = 'SplineInterpolationBBKT129chan.mat';

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_normalization: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%% Load the channels map and arreange it
load(P.chmapfile)
if ~exist('neighbourchan')
    for e=1:size(G,1)
        [a,b]=sort(G(e,:),'descend');
        neighbourchan(e,:)=b(2:end);
        weighthchan(e,:)=a(2:end);
    end
    weighthchan=weighthchan/max(max(weighthchan));
end
nCh = size(EEG.data,1);
if nCh < size(H,1) %not reference avg in the data
    H = H(1:nCh,1:nCh);
    G = G(1:nCh,1:nCh);
end
if ~P.ref  %%reconstruct data from 129 for when it is used as neighbourgh
    rmv129=1;
    neighbours=neighbourchan(129,1:P.nbvoisins)';
    EEG.data(129,:,:)=nanmean(EEG.data(neighbours,:,:),1);
end
if P.ref && nCh==128
    rmv129=1;
    EEG.data(129,:,:)=zeros(1, size(EEG.data,2),size(EEG.data,3));
end


%% ------------------------------------------------------------------------
%% Smooth
for ch=1:nCh
    neighbours=neighbourchan(ch,1:P.nbvoisins)';
    voisins=nanmean(EEG.data(neighbours,:,:),1);  
    EEG.data(ch,:,:) =  EEG.data(ch,:,:)*P.alpha + voisins*(1-P.alpha);  
end
if rmv129
    EEG.data(129,:,:)=[];
end

fprintf('\n')

end