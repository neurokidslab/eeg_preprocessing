% This function function shows the events of some types and asks which
% events you want to delete
%
% INPUTS:
%
%   EEG         = EEGLab structure with the EEG data
%   showevents  = events to show
%   exclevents  = events to avoid
%
% OUTPUT:
%
%   EEG
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------


function EEGout = eega_removeevent(EEG, targetEv, exclEv, tormv)


if nargin < 2
    targetEv = [];
end
if nargin < 3
    exclEv = [];
end
if nargin < 4
    tormv = [];
end

if ischar(targetEv)
    targetEv = {targetEv};
end
if ischar(exclEv)
    exclEv = {exclEv};
end

fprintf('### Remove events ###\n')


EEGout = EEG;

ne = size(EEG.event,2);
fields = fieldnames(EEG.event);
itype = logical(strcmp('type',fields));

etypes = struct2cell(EEG.event);
etypes = squeeze(etypes(itype,:,:));

isnoev = zeros(ne,numel(exclEv));
if ~isempty(exclEv)
    for i=1:numel(exclEv)
        isnoev(:,i) = strcmp(strtrim(etypes),exclEv{i});
    end
end
isnoev = any(isnoev,2);

iskey = ones(ne,1);
if ~isempty(targetEv)
    iskey = ones(ne,numel(targetEv));
    for i=1:numel(targetEv)
        iskey(:,i) = strcmp(strtrim(etypes),targetEv{i});
    end
    iskey = any(iskey,2);
end

isev = ~isnoev & iskey;
nisev = sum(isev);

if nisev~=0
    type = cell(nisev,1);
    latency = nan(nisev,1);
    urevent = nan(nisev,1);
    nievent = nan(nisev,1);
    i = 1;
    for e=1:ne
        
        if isev(e)
            type{i} = EEG.event(e).type;
            latency(i) = EEG.event(e).latency/EEG.srate;
            if isfield(EEG.event(e),'urevent')
                urevent(i) = EEG.event(e).urevent;
            end
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
    
%     fprintf('The list of events to remove is: \n')
%     disp(tbl)
    if ischar(tormv) && strcmp(tormv,'all')
        ev2delet = 1:nisev;
    elseif isempty(tormv)
        ev2delet = input('Which events do you want to delet? (numbers in column 1): ');
    elseif all(tormv<=nisev) && all(tormv>0)
        ev2delet = tormv;
    end
    
    % Remove the irrelevant events
    if ~isempty(ev2delet)
        ev2delet = nievent(ev2delet);
        EEGout = pop_editeventvals(EEG, 'delete', ev2delet);
        fprintf('%d events where deleted:\n',numel(ev2delet))
    end
    
end
fprintf('\n')


end




