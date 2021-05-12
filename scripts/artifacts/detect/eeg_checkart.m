function EEG = eeg_checkart(EEG)

[nEl, nS, nEp] = size(EEG.data);
if ~isfield(EEG, 'artifacts')
    EEG.artifacts.algorithm.parameters = [];
    EEG.artifacts.algorithm.stepname = [];
    EEG.artifacts.algorithm.rejxstep = [];
    EEG.artifacts.BCT = false(nEl,nS,nEp);
    EEG.artifacts.BC = false(nEl,1,nEp);
    EEG.artifacts.BCmanual = [];
    EEG.artifacts.BT = false(1,nS,nEp);
    EEG.artifacts.BE = false(1,1,nEp);
    EEG.artifacts.BS = false(1);
end
if ~isfield(EEG.artifacts, 'BCT')
    EEG.artifacts.BCT = false(nEl,nS,nEp);
end
end