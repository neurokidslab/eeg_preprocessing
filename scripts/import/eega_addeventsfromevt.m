% This function reads the event file (.txt), extract other relavant
% information and integrates it in the EEG structure
%
%   INPUTS
%   EEG         = EEGLab structure with the EEG data
%   pathIn      = path where the event file is
%   latlim      = tolerance in the latency (in samples) between the latnecy
%                   in the event-file and in the EEG structure
%
%   OUTPUT
%   EEG
% -------------------------------------------------------------------------
% Ana Flo June 2020
% -------------------------------------------------------------------------


function EEG = eega_addeventsfromevt( EEG, pathIn, latlim )

fprintf('### Importing information of the events ###\n')

%% Oprtional inputs
if nargin<3; latlim = 3; end

%% Read the event-file

%find an event file with the correct name
[~,name,~] = fileparts(EEG.filename);
listing = dir([pathIn filesep '*.evt']);
for i=1:numel(listing)
    [~,namei,~] = fileparts(listing(i).name);
    k = strfind(name,namei);
    if any(k)
        thename = name(k:length(name));
    end
end
fileventname = [thename '.evt'];
filevent = fullfile(pathIn, fileventname);

%read the event file
L = cell(3000,50);
formatSpec = '%s';
fid = fopen(filevent);
if fid==-1
    error('eega_infoevents: The events file could not be read')
else
    fprintf('Reading information from %s\n',filevent)
end
tline = fgets(fid);
ie = 0;
while ischar(tline)
    ie = ie+1;
    li = textscan(tline,formatSpec,100,'Delimiter','\t'); li = li{1};
    L(ie,1:numel(li)) = li;
    tline = fgets(fid);
end
fclose(fid);
idx_emp = cellfun(@isempty,L);
idx = all(idx_emp,2);
L(idx,:) = [];
idx = all(idx_emp,1);
L(:,idx) = [];

%read the line with the colum names
idx = find(strcmp('Code',L(:,1)));
Ename = L(idx,:);
EnameOK = ~cellfun(@isempty,Ename);

%take the data with the events
Ed = L(idx+1:end,:);

%position in the event file
id_onset = logical(strcmp('Onset',strtrim(Ename)));  % column with the onsets
id_dur = logical(strcmp('Duration',strtrim(Ename)));  % column with the durations
id_code = logical(strcmp('Code',strtrim(Ename)));  % column with the type of stimulus

%get the indexes of the extra fields in the event file
id_extra = find(~EnameOK);
id_extraname = id_extra(1:2:end);
id_extravals = id_extra(2:2:end);

%% Convert all time to samples
for ie=1:size(Ed,1)
    %onsets in samples
    t_txtx = Ed{ie,id_onset};
    tpos = strfind(t_txtx,':');
    t_sec = datenum( t_txtx(tpos(1)-2:end),  'HH:MM:SS.FFF' ) .* (24*60*60) - ...
        datenum( '00:00:00.000', 'HH:MM:SS.FFF' ) .* (24*60*60);
    Ed{ie,id_onset} = round(t_sec * EEG.srate);
    %durations in samples
    t_txtx = Ed{ie,id_dur};
    tpos = strfind(t_txtx,':');
    t_sec = datenum( t_txtx(tpos(1)-2:end),  'HH:MM:SS.FFF' ) .* (24*60*60) - ...
        datenum( '00:00:00.000', 'HH:MM:SS.FFF' ) .* (24*60*60);
    Ed{ie,id_dur} = round(t_sec * EEG.srate);
end

%% Add avents that do not exist in the EEGLAB structure


% % If there are not events add the first one in the event file
% if isempty(EEG.event)
%     EEG.event(1).type = Ed{1,id_code};
%     EEG.event(1).duration = Ed{1,id_dur};
%     EEG.event(1).latency =  Ed{1,id_onset};  % in samples
%     EEG = eeg_checkset(EEG);
% end
  
%get the type and latencies from the EEGLAB structure
type = cell(length(EEG.event),1);
lat = nan(length(EEG.event),1);
for i=1:length(EEG.event)
    type{i} = EEG.event(i).type;
    lat(i) = EEG.event(i).latency / EEG.srate;  % in seconds
end

%loop through the events in the event file
ncount = 0;
for i=1:size(Ed,1)
    
    %check if the event exists in the EELAB structure
    if ~isempty(type)
        idx_type = strcmp(type, Ed(i,id_code));
        idx_time = abs(lat - Ed{i,id_onset})<=latlim;
        exists = any(idx_type==1 & idx_time==1);
    else
        exists = 0;
    end
    if ~exists
        
        fprintf('Adding new event (n=%d)...\n',ncount+1);

        %event number
        evn = length(EEG.event)+1;
        urevn = length(EEG.urevent)+1;
        EEG.urevent(urevn).type = Ed{i,id_code};
        EEG.urevent(urevn).latency = Ed{i,id_onset};
        EEG.event(evn).type = Ed{i,id_code};
        EEG.event(evn).duration = Ed{i,id_dur};
        EEG.event(evn).latency =  Ed{i,id_onset};  % in samples
        EEG.event(evn).urevent = urevn;  
        EEG = eeg_checkset(EEG, 'eventconsistency');
        EEG = eeg_checkset(EEG, 'checkur');
        
%         %add an empty event a the end
%         if isempty(EEG.event)
%             fnames = fieldnames(EEG.event);
%             nfff = length(fnames);
%             fvals = cell(1,nfff);
%         else
%         end
%         EEG = pop_editeventvals( EEG, 'insert', [length(EEG.event) fvals]);
%         
%         
%         %get the extra values
%         ev_codes = cell(length(id_extraname),2);
%         for k=1:length(id_extraname)
%             ev_codes{k,1} = Ed{i,id_extraname(k)};
%             ev_codes{k,2} = Ed{i,id_extravals(k)};
%         end
%         
%         %store the values of the fields
%         ev_type = Ed{i,id_code};
%         ev_duration = Ed{i,id_dur} / EEG.srate;     % duration in seconds 
%         ev_latency = Ed{i,id_onset} / EEG.srate;    % latency in seconds 
%         
%         %modify the different fields (last latency to reorder at the end)
%         EEG.event(evn).codes = ev_codes;
%         EEG = pop_editeventvals( EEG, 'changefield', {evn 'type' ev_type});
%         EEG = pop_editeventvals( EEG, 'changefield', {evn 'duration' ev_duration});
%         EEG = pop_editeventvals( EEG, 'changefield', {evn 'latency' ev_latency});
%            
        % update the counter
        ncount = ncount+1;
    end
end

fprintf('--> %d new events were added\n',ncount);

% %% Check the structure
% EEG = eeg_checkset(EEG, 'eventconsistency');
% EEG = eeg_checkset(EEG, 'checkur');
% EEG = eeg_checkset(EEG);


end
