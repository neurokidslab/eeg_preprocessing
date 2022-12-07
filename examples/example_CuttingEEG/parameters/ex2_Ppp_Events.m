%% -----------------------------------------------------------------------
% Correcting events
%
% This script sets the parameters for 
% - importing extra events
% - correcting events latencies
% - removing uneccessary events
%
% Ana Fló, December 2022
%
% ------------------------------------------------------------------------

function P = ex2_Ppp_Events

% Import extra information for the events
% -------------------------------------------------------------------------
%
% Add extra information in an even file (.evt) 
% This step should be avoided unless extra information in the event 
% file is not imported from the raw data to the EEGLAB structure. 
% The function adding events has been created to add events sent by 
% Psychtoolbox to EGI. The addition of new events will not work with other
% acquisition systems. Customized functions might be required

% (1) add | (0) do not add
P.import.apply = 1;  

% Event to use to align differences in time offsets between the event file 
% and the EEG structure
% Only necessary if events are added (Ppp.importevntapply = 1)
P.import.ev0 = 'STRT';   


% Correct events latencies using another event
% -------------------------------------------------------------------------
%
% Correct the latency of a given event using Digital Input events (DINs)
% Leave it empty to do not correct event's latencies
% In this example the latency of the events 'Icue' 'Iout' 'Ieye' are 
%  corrected by the latency of 'DIN6'
P.corrlatencies.ev = [];
P.corrlatencies.ev(1).event1 = {'Icue'};  % event which latencies will be corrected
P.corrlatencies.ev(1).event2 = 'DIN6';    % event to use to correct the latency
P.corrlatencies.ev(2).event1 = {'Iout'};  % event which latencies will be corrected
P.corrlatencies.ev(2).event2 = 'DIN6';    % event to use to correct the latency
P.corrlatencies.ev(3).event1 = {'Ieye'};  % event whdich latencies will be corrected
P.corrlatencies.ev(3).event2 = 'DIN6';    % event to use to correct the latency


% Delete unuseful events. 
% -------------------------------------------------------------------------
%
% Remove the indicated events
% Leave it empy to do not remove events
P.rmv.ev = {'DIN1' 'DIN6'};

end
