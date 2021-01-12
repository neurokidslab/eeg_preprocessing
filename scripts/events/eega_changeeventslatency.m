% This function adds to the latency of the target events indicated by the indexes after
% the reference event. It works fine on non epoched data
%
%   INPUTS          -------------------------------------------------------
%   EEG             :
%   targevent       : event to change
%   time2add        : time to add in seconds
%   refevent        : event from where it is counted
%   idxrefevent     : indexes of the reference event
%   idxtargevent    : indexes of the events from all the possible ones
%   twtargevent     : time window where the target avents can be, measured  
%                     from the referent. In seconds.
%
% -------------------------------------------------------------------------
% Ana Flo September 2017
% -------------------------------------------------------------------------

function EEGout = eega_changeeventslatency( EEG, targevent, time2add,...
    refevent, idxrefevent, idxtargevent, twtargevent)

if nargin<4 
    refevent = [];
end
if nargin<5 || isempty(idxrefevent)
    idxrefevent = 'all';
end
if nargin<6 || isempty(idxtargevent)
    idxtargevent = 'all';
end
if nargin<7 || isempty(twtargevent)
    twtargevent = [0 Inf];
end


EEGout = EEG;

% time to add to samples
smpl2add = round(time2add*EEG.srate);

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
    
    % Change the latency
    for i=1:length(idxevent)
        fprintf('Event number: %05.0d\n',idxevent(i))
        
        % Determine the new time
        newtime =  EEG.event(idxevent(i)).latency + smpl2add;
        
        % Edit the event
        EEGout = correctlatency(EEGout, idxevent(i), newtime);
        
        % Update the counter
        nEvchange = nEvchange+1;
        
        
    end
end
EEGout = eeg_checkset(EEGout);
EEGout = eeg_checkset(EEGout, 'eventconsistency');
% EEGout = eeg_checkset(EEGout, 'checkur');

fprintf('\n')

fprintf('The latency of %d events were modified\n\n', nEvchange)

end

function EEG = correctlatency(EEG, valnum, editval)

EEG.event(valnum).latency = editval;

if isfield(EEG, 'urevent') & isfield(EEG.event, 'urevent')
    urvalnum  = EEG.event(valnum).urevent;
    
    % latency case
    % ------------
    if ~isempty(editval)
        if isfield(EEG.urevent, 'epoch')
            urepoch = EEG.urevent(urvalnum).epoch;
            
            
            % find closest event latency
            % --------------------------
            if valnum<length(EEG.event)
                if EEG.event(valnum+1).epoch == urepoch
                    urlatency = EEG.urevent(EEG.event(valnum+1).urevent).latency;
                    latency   = EEG.event(valnum+1).latency;
                end;
            end;
            if valnum>1
                if EEG.event(valnum-1).epoch == urepoch
                    urlatency = EEG.urevent(EEG.event(valnum-1).urevent).latency;
                    latency   = EEG.event(valnum-1).latency;
                end;
            end;
            
            % update event
            % ------------
            if exist('urlatency') ~=1
                disp('Urevent not updated: could not find other event in the epoch');
            else
                editval = urlatency - ( latency - editval ); % new latency value
            end;
        else
            editval = eeg_urlatency(EEG.event, EEG.event(valnum).latency);
        end;
    else % empty editval
        EEG.event(valnum).latency = NaN;
    end;
    
    EEG.urevent(urvalnum).latency = editval;
end;


end

