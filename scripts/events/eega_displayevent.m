% This function displays all the events indicates in events 
%
% INPUTS:
%
%   EEG         = EEGLab structure with the EEG data
%   events      = cell array of strings with specific event to display
%
% OUTPUT:
%
%   EEG
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------


function tbl = eega_displayevent(EEG, events)

if nargin < 2
    events = [];
end

ne = size(EEG.event,2);
fields = fieldnames(EEG.event);
itype = logical(strcmp('type',fields));

etypes = struct2cell(EEG.event);
etypes = squeeze(etypes(itype,:,:));
if isempty(events)
    events=unique(etypes);
end

isev = ones(ne,numel(events));
if ~isempty(events)
    for i=1:numel(events)
        isev(:,i) = strcmp(strtrim(etypes),strtrim(events{i}));
    end
    isev = any(isev,2);
end
nisev = sum(isev);

type = cell(nisev,1);
latency = nan(nisev,1);
urevent = nan(nisev,1);
nievent = nan(nisev,1);
i = 1;
for e=1:ne
    
    if isev(e)
        type{i} = EEG.event(e).type;
        latency(i) = EEG.event(e).latency/EEG.srate;
        urevent(i) = EEG.event(e).urevent;
        nievent(i) = e;
        i = i+1;
    end
    
end

latdiff = diff(latency);
latdiff = [0 ; latdiff];
rawnames = cellfun(@num2str,num2cell([1:nisev]'),'UniformOutpu',0);
tbl = table(nievent, urevent, type, latency, latdiff,... 
    'VariableNames', {'Index' 'UR' 'Type' 'Latency' 'LatencyDiff'},...
    'RowNames',rawnames);

fprintf('The list of events is: \n')
disp(tbl)




