% Run the MARA lagorithm to reject components fomr the ICA

function EEG = eega_icarunmara(EEG, keeppre, usealphaband)

if nargin<2
    keeppre = 0;
end
if nargin<3
    usealphaband = 1;
end
% if nargin<3
%     alphaband = [8 13];
% end

%% ========================================================================
%% Copy the data
eeg = EEG;

%% ========================================================================
%% Run MARA
[artcomps, info] = MARA_babies(eeg, usealphaband);
% [artcomps, info] = MARA_babies(eeg, alphaband);

if keeppre
    EEG.reject.gcompreject(artcomps) = 1;
else
    EEG.reject.gcompreject = false(1,EEG.nbchan);
    EEG.reject.gcompreject(artcomps) = 1;
end


end