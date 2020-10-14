% This function shows the number of DINs after each non DIN event
%
% INPUTS:
%
%   EEG         = EEGLab structure with the EEG data
%   Event       = string with the event which number will be checked
%   refEvent    = cell array with other events. The number of times that
%               Event appears after each of them will be checked
%   NcorrectEv  = correct number of Event after each refEvent. If it is
%               provided a warning message will be shown if there is not
%               agreement
%
% OUTPUT:
%
%   EEG
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------

function eega_shownumevents( EEG, Event, refEvent, NcorrectEv )

if nargin<3
    refEvent = {''};
end
if nargin<4
    NcorrectEv = [];
end
if ischar(refEvent)
    refEvent = {refEvent};
end
if ischar(Event)
    Event = {Event};
end

fprintf('### Number of events ###\n')

nEv = numel(EEG.event);
evtype = cell(nEv,1); % Events type
evlatency = nan(nEv,1); % Events latency
for i=1:nEv
    evtype{i} = strtrim(EEG.event(i).type);
    evlatency(i) = EEG.event(i).latency;
end
isEv = find(strcmp(evtype,Event));
noEv = find(~strcmp(evtype,Event));

fprintf('Total %s events : %d\n',Event{1},length(isEv))
for E = 1:numel(refEvent)
    isE = find(strcmp(refEvent{E},evtype));
    if ~isempty(isE)
        for j=1:numel(isE)
            next = find(noEv>isE(j));
            if ~isempty(next)
                next = noEv(next(1));
                dinafter = ( evlatency > evlatency(isE(j)) ) & ( (evlatency < evlatency(next) ) );
            else
                dinafter = ( evlatency > evlatency(isE(j)) );
            end
            ndinafter = sum(dinafter);
            
            if ~isempty(NcorrectEv) && ndinafter~=NcorrectEv(E)
                warning('Incorrect number of %s events found !!!',Event{1})
                fprintf('%s\n', repmat('-',[20 1]))
                fprintf('Event %s\n', refEvent{E})
                fprintf('   - Event after %s\n', evtype{next})
                fprintf('   - Number of %s events after : %d\n',Event{1}, sum(dinafter))
            else
                fprintf('%d %s events after the event %s\n',sum(dinafter), Event{1}, refEvent{E})
            end           
        end
    end
end
fprintf('\n')








