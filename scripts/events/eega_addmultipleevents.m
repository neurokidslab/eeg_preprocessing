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
        fprintf('\b\b\b\b')
        fprintf('%04.0d', i)
        
        
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
%             newev.latency = newev.latency/EEG.srate; % latency in seconds
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
            if isfield(EEGout.event(REfev(k)),'epoch')
                newev.epoch = EEGout.event(REfev(k)).epoch;
            end
            
            fvals = cell(1,numel(fff));
            for j=1:numel(fff)
                fvals{j} =  newev.(fff{j});
            end
            
            % add the event - added by Ana Flo 06/12/2021
            fnames_i = fff;
            fvals_i = fvals;
            evn = length(EEGout.event)+1;
            for ff=1:length(fnames_i)
                EEGout.event(evn).(fnames_i{ff}) = fvals_i{ff};
            end
            
            % add the urevent - added by Ana Flo 06/12/2021
            fnames_i_ur = {'type' 'latency'};
            fvals_i_ur = cat(2,fvals_i(strcmp(fnames_i,'type')), fvals_i(strcmp(fnames_i,'latency')));
            evn = length(EEGout.urevent)+1;
            for ff=1:length(fnames_i_ur)
                EEGout.urevent(evn).(fnames_i_ur{ff}) = fvals_i_ur{ff};
            end
            
%             EEGout = pop_editeventvals( EEGout, 'add', [evpos fvals]);  - commented by Ana Flo 06/12/2021
%             EEGout = pop_editeventvals( EEGout, 'add', [1 fvals]);
%             EEGout = pop_editeventvals( EEGout, 'add', [length(EEGout.event)+1 fvals]);
                 
        end
    end
end
% Check the structure - added by Ana Flo 06/12/2021
EEGout = eeg_checkset(EEGout, 'eventconsistency');
EEGout = eeg_checkset(EEGout, 'checkur');

fprintf('\n')

fprintf('\n')

end