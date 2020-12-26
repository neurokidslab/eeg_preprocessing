% This function chaages the targets events indicated by the indexes after
% the reference event
%
%   INPUTS          -------------------------------------------------------
%   EEG             :
%   targevent       : event to change
%   field           : field to change
%   newfieldval     : new value
%   refevent        : event from where it is counted
%   idxrefevent     : indexes of the reference event
%   idxtargevent    : indexes of the events from all the possible ones
%   twtargevent     : time window where the target avents can be, measured  
%                     from the referent. In seconds.
%   field           : name of the field to chnage
%   newfieldval     : new field value
%                     If it is a sigle value the provided value is used
%                     If it is a cell array of size nx2 when the old value
%                     is equal to (i,1) it is replaced by (i,2)
%
% -------------------------------------------------------------------------
% Ana Flo September 2017
% -------------------------------------------------------------------------

function EEGout = eega_changeevents( EEG, targevent, field, newfieldval,...
    refevent, idxrefevent, idxtargevent, twtargevent)

if nargin<5 
    refevent = [];
end
if nargin<6 || isempty(idxrefevent)
    idxrefevent = 'all';
end
if nargin<7 || isempty(idxtargevent)
    idxtargevent = 'all';
end
if nargin<8 || isempty(twtargevent)
    twtargevent = [0 Inf];
end


EEGout = EEG;
if ischar(newfieldval)
    newfieldval = {newfieldval};
end

% Events type
evtype = cell(1,size(EEG.event,2));
for j=1:size(EEG.event,2)
    evtype{j} = strtrim(EEG.event(j).type);
end

% Find the reference events
if isempty(refevent)
    Refev = 1;
else
    Refev = find(strcmp(evtype,refevent));
    if ~strcmp(idxrefevent,'all')
        if idxrefevent<=length(Refev)
            Refev = Refev(idxrefevent);
        else
            warning('The reference event %s number %d was not found', refevent, idxrefevent)
            return
        end
    end
end
RefevLat = nan(size(Refev));
for i =1:length(Refev)
    RefevLat(i) = EEG.event(Refev(i)).latency;
end

% Find the target events
Targev = find(strcmp(evtype,targevent));
TargevLat = nan(size(Targev));
for i =1:length(Targev)
    TargevLat(i) = EEG.event(Targev(i)).latency;
end

% Change after each reference event
nEvchange =0;
fprintf('### Changing events of type %s... \n', targevent)
for k=1:numel(Refev)
    
    % Find the target events in the time window
    tw = RefevLat(k) + twtargevent*EEG.srate;
    Targev_ktime = Targev( (TargevLat <= tw(2)) & (TargevLat >= tw(1)) );
    Targev_k = Targev(Targev>=Refev);
    posevent = intersect(Targev_k,Targev_ktime);
        
    % Find the target events acording to the provided indexes
    if strcmp('all', idxtargevent)
        idxevent = posevent;
    elseif strcmp('first', idxtargevent)
        idxevent = posevent(1);
    else
        idxtargevent(idxtargevent>length(posevent)) = [];
        idxevent = posevent(idxtargevent);
    end
    
    % Change the events
    for i=1:length(idxevent)
        fprintf('Event number: %05.0d\n',idxevent(i))
        
        % Determine the new value
        if numel(newfieldval)==1
            newval_i = newfieldval{1};
        else
            oldval_i = EEG.event(idxevent(i)).(field);
            idx = strcmp(newfieldval(:,1), oldval_i);
            if any(idx)
                newval_i = newfieldval{idx,2};
            else
                newval_i = [];
            end
        end
        
        if ~isempty(newval_i)
            
            % Edit the event
            if ~isfield(EEGout.event,field)
                for iiev=1:length(EEGout.event)
                    EEGout.event(iiev).(field) = [];
                end
                if isfield(EEGout,'epoch')
                    for iiep=1:length(EEGout.epoch)
                        EEGout.epoch(iiep).(['event' field]) = repmat({[]},[1 length(EEGout.epoch(iiep).event)]);
                    end
                end
            end
            EEGout.event(idxevent(i)).(field) = newval_i;
            if isfield(EEGout,'epoch') && ~isempty(EEGout.epoch)
                for iiep=1:length(EEGout.epoch)
                    idxEpEv = EEGout.epoch(iiep).event==idxevent(i);
                    if any(idxEpEv)
                        EEGout.epoch(iiep).(['event' field]){idxEpEv} = newval_i;
                    end
                end
            end
            %             EEGout=pop_editeventvals(EEGout,'changefield',{idxevent(i) field newval_i});
            
            nEvchange = nEvchange+1;
            
        end
    end
end
EEGout = eeg_checkset(EEGout);
EEGout = eeg_checkset(EEGout, 'eventconsistency');
% EEGout = eeg_checkset(EEGout, 'checkur');

fprintf('\n')

fprintf('The field %s of %d events were modified\n\n', field, nEvchange)

end
