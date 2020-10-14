% This function removes DINs that are too close.
% It is useful to remove the second DIN mark that many times is send
% The first DIN has to be correct
%
% INPUTS:
%
%   EEG         = EEGLab structure with the EEG data
%   dinevent    = cell array with the DIN event
%   mindindist  = minimun possble distansse between DINs. If a closer than
%               mindindist is found it is remove. In seconds
%   dinnumber   = number of DINs that should be present
%
% OUTPUT:
%
%   EEG
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------

function EEG = eega_rmvcloseeventsA( EEG, Event, mindindist )

if ischar(Event)
    Event = {Event};
end


fprintf('### Removing very close events ###\n')

mindindist = EEG.srate * mindindist; % latency in samples

nEv = numel(EEG.event);
Ev2rmv = zeros(nEv,1); % Events to remove
evtype = cell(nEv,1); % Events type
evlatency = nan(nEv,1); % Events latency
for i=1:nEv
    evtype{i} = EEG.event(i).type;
    evlatency(i) = EEG.event(i).latency;
end
isdin = find(strcmp(evtype,Event));
nEvDINs = numel(isdin);

% Show the info for the first DIN
fprintf('For the first DIN\n')
if isdin(1)>1
    fprintf('	- Previous event type     : %s\n', evtype{isdin(1)-1})
    fprintf('	- Previous event latency  : %f\n', (evlatency(isdin(1)-1) - evlatency(isdin(1)))/EEG.srate)
else
    fprintf('	- Previous event type     : none\n')
end
fprintf('	- Next event type         : %s\n', evtype{isdin(1)+1})
fprintf('	- Next event latency      : %f\n', (evlatency(isdin(1)+1) - evlatency(isdin(1)))/EEG.srate)

runit = 'x';
while ~any(runit=='yn');
    runit = input('Is this the first DIN of the experiment? (y/n): ','s');
end

% Check DINs from the first DIN and remove
if strcmp(runit,'y')
    j = 1;
    while j<nEvDINs
        latnext = evlatency(isdin(j+1)) - evlatency(isdin(j));
        if latnext>mindindist
            j = j + 1;
        else
            Ev2rmv(isdin(j)+1) = 1;
            j = j + 2;
        end
    end
    
    % Remove the DINs   
    if any(Ev2rmv)
        Ev2rmv = find(Ev2rmv);
        EEG = pop_editeventvals(EEG, 'delete', Ev2rmv);
        fprintf('%d %s events where deleted:\n', numel(Ev2rmv), Event{1})
    end
else 
    fprintf('Nothing will be done\n')
end

fprintf('\n')





