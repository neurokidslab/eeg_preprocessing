% This function corrects the latency of the events indicated in "events"
% using another event type, "eventT".
%
% Usually "eventsT" are DINs. %
% INPUTS:
%
%   EEG         = EEGLab structure with the EEG data
%   eventT      = event indicating the correct latency
%   events      = cell array of strings with the events to correct
%   correctEvT  = number of requiered eventT after the event. if it's not
%               provided this is not checked
%
% OUTPUT:
%
%   EEG
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------


function EEG = eega_latencyevent(EEG, eventT, events, correctEvT)

if nargin < 4
    correctEvT = [];
end

fprintf('### Correct latencies ###\n')

nEv = size(EEG.event,2);
etype = cell(nEv,1); % Events type
elatency = nan(nEv,1); % Events latency
for i=1:nEv
    etype{i} = strtrim(EEG.event(i).type);
    elatency(i) = EEG.event(i).latency;
end

noEvT = find(~strcmp(etype,eventT));
iii=0;
correction = [];
for E = 1:numel(events)
    isE = find(strcmp(events{E},etype));
    if ~isempty(isE)
        for j=1:numel(isE)
            
            dochange = 1;
            
            next = find(noEvT>isE(j));
            if isempty(next)
                dinafter = ( elatency > elatency(isE(j)) );
                ndinafter = sum(dinafter);
                theEvT = find(dinafter);
            else
                next = noEvT(next(1));
                dinafter = ( elatency > elatency(isE(j)) ) & ( (elatency < elatency(next) ) );
                ndinafter = sum(dinafter);
                theEvT = find(dinafter);
            end
            if ~isempty(theEvT)
                theEvT = theEvT(1);
            else
                dochange = 0;
            end
            
            % Check the number of DINs
            if ~isempty(correctEvT) && (ndinafter~=correctEvT(E))
                warning('Incorrect number of %s events found !!!',eventT)
                fprintf(' Number of %s events after %s: %d\n',eventT, events{E}, sum(dinafter))
            end
            
            % Change
            if dochange
                iii = iii + 1;
                fprintf('Correcting the latency of event %s by %s - event %d - change %d...\n',events{E}, eventT, isE(j), iii)
                correction(iii) = EEG.event(theEvT).latency - EEG.event(isE(j)).latency;
                
                % the event
                EEG = correctlatency(EEG, isE(j), EEG.event(theEvT).latency);
%                 newlatency = (EEG.event(theEvT).latency - 1) / EEG.srate + EEG.xmin;  % must be in seconds
%                 EEG = pop_editeventvals(EEG, 'changefield', {isE(j) 'latency' newlatency});

            end 
        end       
    end
end
EEG = eeg_checkset(EEG);
EEG = eeg_checkset(EEG, 'eventconsistency');
% EEG = eeg_checkset(EEG, 'checkur');

correction = correction/EEG.srate*1000;
fprintf('The latency of %d events was corrected\n',iii)
fprintf('- Mean delay %f ms\n',mean(correction))
fprintf('- SD of the delay %f ms\n',std(correction))
fprintf('- min of the delay %f ms\n',min(correction))
fprintf('- max of the delay %f ms\n',max(correction))

fprintf('\n')

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
