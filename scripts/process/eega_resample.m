function EEG = eega_resample(EEG, newfrq)

orgfrq = EEG.srate;
szorg = size(EEG.data);

EEG = pop_resample( EEG, newfrq);
sznew = size(EEG.data);

if isfield(EEG,'artifacts')
    if isfield(EEG.artifacts,'BCT')
        x = EEG.artifacts.BCT;
        x = reshape(x, szorg(1),[])';
        x = resample(double(x),newfrq,orgfrq);
        x = logical(round(x));
        EEG.artifacts.BCT = reshape(x', szorg(1), []);
    end
    if isfield(EEG.artifacts,'CCT')
        x = EEG.artifacts.CCT;
        x = reshape(x, szorg(1),[])';
        x = resample(double(x),newfrq,orgfrq);
        x = logical(round(x));
        EEG.artifacts.CCT = reshape(x', szorg(1), []);
    end
    if isfield(EEG.artifacts,'BT')
        x = EEG.artifacts.BT;
        x = reshape(x, 1,[])';
        x = resample(double(x),newfrq,orgfrq);
        x = logical(round(x));
        EEG.artifacts.BT = reshape(x', 1, []);
    end
end

% remove the urevent field
if isfield(EEG.event,'urevent')
    EEG.event = rmfield(EEG.event,'urevent');
end

end