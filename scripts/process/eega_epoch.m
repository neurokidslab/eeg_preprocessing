function [EEG] = eega_epoch(EEG, events, timelim)

fprintf('### Epoching ###\n')

nEli    = size(EEG.data,1);
nSi     = size(EEG.data,2);
nEpi    = size(EEG.data,3);

Nev = numel(EEG.event);
evtype = cell(Nev,1);
evlatency = nan(Nev,1);
evidx = zeros(Nev,1);
for e=1:Nev
    if ~ischar(EEG.event(e).type)
        evtype{e} = strtrim(num2str(EEG.event(e).type));
    else
        evtype{e} = strtrim(EEG.event(e).type);
    end
    evlatency(e) = EEG.event(e).latency;
    evidx(e) = any(strcmp(evtype{e}, events));
end

timelimsmpl = round(timelim * EEG.srate);

%% Epoch the data
[EEG, indices] = pop_epoch( EEG, events, timelim );

%% Epoch structures related with the artifacts
evidx = find(evidx);
evidx = evidx(indices);
evtime = evlatency(evidx);

if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BCT')
    [ EEG.artifacts.BCT ] = epoch( EEG.artifacts.BCT, evtime, timelimsmpl );
    EEG.artifacts.BCT = logical(EEG.artifacts.BCT);
end

if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BC')
    [ EEG.artifacts.BC ] = epoch( repmat(EEG.artifacts.BC,[1 nSi 1]), evtime, timelimsmpl );
    EEG.artifacts.BC = logical(EEG.artifacts.BC(:,1,:));
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BT')
    [ EEG.artifacts.BT ] = epoch( EEG.artifacts.BT, evtime, timelimsmpl );
    EEG.artifacts.BT = logical(EEG.artifacts.BT);
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'CCT')
    [ EEG.artifacts.CCT ] = epoch( EEG.artifacts.CCT, evtime, timelimsmpl );
    EEG.artifacts.CCT = logical(EEG.artifacts.CCT);
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BE')
    [ EEG.artifacts.BE ] = epoch( repmat(EEG.artifacts.BE,[1 nSi 1]), evtime, timelimsmpl );
    EEG.artifacts.BE = logical(EEG.artifacts.BE(:,1,:));
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts,'BEm')
    [ EEG.artifacts.BEm ] = epoch( repmat(EEG.artifacts.BEm,[1 nSi 1]), evtime, timelimsmpl );
    EEG.artifacts.BEm = logical(EEG.artifacts.BEm(:,1,:));
end

%% Epoch other structures
if isfield(EEG,'reference') && (size(EEG.reference,2)==nSi) && (size(EEG.reference,3)==nEpi)
    [ EEG.ref ] = epoch( EEG.ref, evtime, timelimsmpl );
end

%% Add as an event the orignal epoch number
if isfield(EEG, 'epoch')
    for e=1:length(EEG.epoch)
        EEG.epoch(e).eventEpoch = repmat({sprintf('%04d',e)},[1 length(EEG.epoch(e).event)]);
    end
end

%% Update the summary structure
if exist('eega_summarypp','file')==2
    EEG = eega_summarypp(EEG);
end

fprintf('\n')
end

