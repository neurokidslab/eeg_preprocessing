% This function adds events after a reference event
%
%   INPUTS
%   EEG         =
%   newevent    = cellarray with the new event's name
%   refevent    = event from where the new events are added
%   eventspos   = position of the events in seconds (n x 1)
%
% -------------------------------------------------------------------------
% Ana Flo May 2017
% -------------------------------------------------------------------------

function EEGout = eega_addmultipleevents( EEG, newevent, refevent, eventspos )

EEGout = EEG;

fprintf('### Ading events %s after the event %s ...\n',newevent,refevent)


% Find the reference event
evtype = cell(1,size(EEG.event,2));
for j=1:size(EEG.event,2)
    evtype{j} = strtrim(EEG.event(j).type);
end
REfev = find(strcmp(evtype,refevent));
REflat = nan(1,numel(REfev));
for j=1:numel(REfev)
    REflat(j) = EEG.event(REfev(j)).latency;
end

% Add after each reference event
fprintf('Adding event %s, number   ', newevent)
for k=1:numel(REfev)
    for i=1:length(eventspos)
        fprintf('\b\b')
        fprintf('%02.0d', i)
        
        
        % Get the latencies of the existying events
        evlat = zeros(1,size(EEGout.event,2));
        for j=1:size(EEGout.event,2)
            evlat(j) = EEGout.event(j).latency;
        end
        
        % Get the position for the new event
        evpos = sum( (evlat / EEGout.srate ) <= (REflat(k)/EEGout.srate + eventspos(i)) ) + 1;
        if evpos>0
            
            % Create the event
            fff = fieldnames(EEGout.event);
            for f = 1:numel(fff)
                newev.(fff{f}) = [];
            end
            newev.type = newevent;
            newev.latency = round(REflat(k) + eventspos(i)*EEG.srate); % latency in samples
            newev.latency = newev.latency/EEG.srate; % latency in seconds
            newev.urevent = numel(EEGout.urevent) + 1;
            if isfield(newev,'Type')
                newev.Type = 'Stimulus Event';
            end
            if isfield(newev,'Track')
                newev.Track = '';
            end
            if isfield(newev,'Code')
                newev.Code = newevent;
            end
            if isfield(EEG.event(REfev(k)),'epoch')
                newev.epoch = EEG.event(REfev(k)).epoch;
            end
            
            fvals = cell(1,numel(fff));
            for j=1:numel(fff)
                fvals{j} =  newev.(fff{j});
            end
            EEGout = pop_editeventvals( EEGout, 'add', [evpos fvals]);
                 
        end
    end
end
fprintf('\n')

fprintf('\n')

end