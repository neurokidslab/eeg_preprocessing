% This function chaages the targets events with certain field
%
%   INPUTS          -------------------------------------------------------
%   EEG             :
%   targevent       : event to change
%   field           : field to change
%   newfieldval     : new value
%   field2cond      : fields to condition the change
%   vals2cond       : values to condition the change
%
% -------------------------------------------------------------------------
% Ana Flo September 2017
% -------------------------------------------------------------------------

function EEGout = eega_changeeventcond( EEG, targevent, field, newfieldval, field2cond, vals2cond)

EEGout = EEG;


% Events type
evtype = cell(1,size(EEG.event,2));
for j=1:size(EEG.event,2)
    evtype{j} = strtrim(EEG.event(j).type);
end

% Find the target events
Targev = find(strcmp(evtype,targevent));

% change the evnts
nEvchange = 0;
fprintf('### Changing events of type %s... \n', targevent)
for e=1:length(Targev)
    theev = Targev(e);
    
    % check if ti respect all the conditions
    cond = false(1, length(field2cond));
    for i=1:length(field2cond)
        cond(i) = strcmp(EEG.event(theev).(field2cond{i}), vals2cond{i});
    end
    
    % Change
    if all(cond)
%         fprintf('Event number: %05.0d\n', theev)
        
        % create the field if it does not exist
        if ~isfield(EEGout.event(theev), field)
            EEGout.event(theev).(field) = [];
            if isfield(EEGout,'epoch')
                for iiep=1:length(EEGout.epoch)
                    EEGout.epoch(iiep).(['event' field]) = repmat({[]},[1 length(EEGout.epoch(iiep).event)]);
                end
            end
        end
        
        % modify the field
        EEGout.event(theev).(field) = newfieldval;
        if isfield(EEGout,'epoch') && ~isempty(EEGout.epoch)
            for iiep=1:length(EEGout.epoch)
                idxEpEv = EEGout.epoch(iiep).event==idxevent(i);
                if any(idxEpEv)
                    EEGout.epoch(iiep).(['event' field]){idxEpEv} = newval_i;
                end
            end
        end
        
        nEvchange = nEvchange+1; 
    end
end
fprintf('\n')
fprintf('The field %s of %d events were modified\n\n', field, nEvchange)

EEGout = eeg_checkset(EEGout, 'eventconsistency');
EEGout = eeg_checkset(EEGout, 'checkur');


end
