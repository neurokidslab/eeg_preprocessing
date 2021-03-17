function EEG = eega_rmvelectrodes(EEG,Ch2Rmv)

if ~isempty(Ch2Rmv)
    ChAll = (1:size(EEG.data,1));
    Ch2Keep = ChAll;
    Ch2Keep(Ch2Rmv) = [];
    
    % remove the channels
    EEG = pop_select(EEG,'channel',Ch2Keep);
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BCT')
        EEG.artifacts.BCT = EEG.artifacts.BCT(Ch2Keep,:,:);
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BC')
        EEG.artifacts.BC = EEG.artifacts.BC(Ch2Keep,:,:);
    end
    if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'CCT')
        EEG.artifacts.CCT = EEG.artifacts.CCT(Ch2Keep,:,:);
    end
end

end