% This function updates the strucutr with the summary information of the
% preprocessing 


function EEG = eega_summarypp(EEG)

[Nch,Nsmpl,Nep] = size(EEG.data);

% for non epoched data
if isfield(EEG,'epoch') && isempty(EEG.epoch)
    EEG.summary.cont.nsmpl = Nsmpl;
    EEG.summary.cont.nch = Nch;
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
        EEG.summary.cont.BCTn = sum(EEG.artifacts.BCT(:));
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
        EEG.summary.cont.CCTn = sum(EEG.artifacts.CCT(:));
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
        EEG.summary.cont.BTn = sum(EEG.artifacts.BT);
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
        EEG.summary.cont.BCn = sum(EEG.artifacts.BC);
        EEG.summary.cont.BC = find(EEG.artifacts.BC);
    end
end

% for epoched data
if isfield(EEG,'epoch') && ~isempty(EEG.epoch)
    EEG.summary.epoch.nsmpl = Nsmpl*Nep;
    EEG.summary.epoch.nep = Nep;
    EEG.summary.epoch.nch = Nch;
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
        EEG.summary.epoch.BCTn = sum(EEG.artifacts.BCT(:));
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
        EEG.summary.epoch.CCTn = sum(EEG.artifacts.CCT(:));
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
        EEG.summary.epoch.BTn = sum(EEG.artifacts.BT(:));
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
        EEG.summary.epoch.BCn = sum(EEG.artifacts.BC(:));
    end
end

end