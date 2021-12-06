% This function reads the event file (.txt), extract other relavant
% information adn integrates it in the EEG structure
%
%   INPUTS
%   EEG         : EEGLab structure with the EEG data
%   filevent    : name of the txt file with the events
%   event0      : event use to align the timings if they are not the same
%               If not provided the first common event is used
%   eventtype   : event type for which extra information is collected
%               If nor provided, all events are updated
%   latlim      : tolerance in the latency (in samples) between the latency
%               in the event-file and in the EEG structure. Default value 3
%
%   OUTPUT
%   EEG
% -------------------------------------------------------------------------
% Ana Flo September 2018
% modify in September 2020.
% -------------------------------------------------------------------------


function EEG = eega_importinfoevents( EEG, pathIn, event0, eventtype, latlim )

fprintf('### Importing information of the events ###\n')

%% Oprtional inputs
if nargin<3; event0 = []; end
if nargin<4; eventtype = []; end
if nargin<5; latlim = 3; end

if isempty(eventtype)
    ev = cell(length(EEG.event),1);
    for e=1:length(EEG.event)
        ev{e} = strtrim(EEG.event(e).type);
    end
    eventtype = unique(ev);
else
    if ischar(eventtype)
        eventtype = {eventtype};
    end
end
if ~isempty(event0)
    if isnumeric(event0)
        event0 = {EEG.event(event0).type};
    end
    if ischar(event0)
        event0 = {event0};
    end
end


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

%read the line with the colum names
idx = find(strcmp('Code',L(:,1)));
Ename = L(idx,:);
EnameOK = ~cellfun(@isempty,Ename);

%take the data with the events
Ed = L(idx+1:end,:);

%% Convert all time to samples for the event file
id_time = logical(strcmp('Onset',strtrim(Ename)));  % column with the times
for ie=1:size(Ed,1)
    t_txtx = Ed{ie,id_time};
    tpos = strfind(t_txtx,':');
    t_sec = datenum( t_txtx(tpos(1)-2:end),  'HH:MM:SS.FFF' ) .* (24*60*60) - ...
        datenum( '00:00:00.000', 'HH:MM:SS.FFF' ) .* (24*60*60);
    Ed{ie,id_time} = round(t_sec * EEG.srate);
end
TIMEevt = cell2mat(Ed(:,id_time));
TYPEevt = Ed(:,logical(strcmp('Code', strtrim(Ename))));

%% Get all the timing and names of the events for the EEG structure
TIMEeeg = nan(length(EEG.event),1);
TYPEeeg = cell(length(EEG.event),1);
for ie=1:length(EEG.event)
    TIMEeeg(ie) = EEG.urevent(EEG.event(ie).urevent).latency;
    TYPEeeg{ie} = EEG.urevent(EEG.event(ie).urevent).type;
end

%% Determine the events to use for aligment
if isempty(event0)
    event0 = intersect(unique(TYPEeeg), unique(TYPEevt)); 
    if isempty(event0)
        error('Commond event were not found')
    end
end

%% Align the times to the same t0

% find the time zero events
t0_idx_evt = [];
t0_idx_eeg = [];
for i=1:length(event0)
    for j=1:length(EEG.urevent)
        if strcmp(strtrim(EEG.urevent(j).type), strtrim(event0{i}))
            t0_idx_eeg = [t0_idx_eeg; j];
        end
    end
    for j=1:length(TYPEevt)
        if strcmp(strtrim(TYPEevt{j}), strtrim(event0{i}))
            t0_idx_evt = [t0_idx_evt; j];
        end
    end
end
t0_idx_eeg = sort(t0_idx_eeg);
t0_idx_evt = sort(t0_idx_evt);

% check that all the events in the text file exist in the eeg structure
if length(t0_idx_evt)~=length(t0_idx_eeg)
    error('Missmatch in the number of t0 events')
end

% check that the difference in timing is constant
tt = [0; TIMEevt(t0_idx_evt) - TIMEeeg(t0_idx_eeg)];
idxjumps = diff(tt)>latlim;
if sum(idxjumps)==0
    fprintf('Timings are aligned, nothing will be done. \n');
elseif sum(idxjumps)==1
    for i=1:length(event0)
        fprintf('Timing will be aligned based on the first event %s \n', event0{i});
    end
elseif sum(idxjumps)>1
    warning('Dissaliging observed at %d points!! The data will try to be realigned', sum(idxjumps))
end

% bring to zero between jumps in the timing
idxjumps = find(idxjumps);
if ~isempty(idxjumps)
    
    evt2align = t0_idx_evt(idxjumps);
    eeg2align = t0_idx_eeg(idxjumps);
    
    % first aligment, apply it for all times before and after
    fprintf('Aligment 1\n')
    fprintf(' - Aligmant based on event %s, event %d in the event file, latency %d \n', TYPEevt{evt2align(1)}, evt2align(1), TIMEevt(evt2align(1)))
    fprintf(' - Aligmant based on event %s, event %d in the EEG structure, latency %d \n', TYPEeeg{evt2align(1)}, eeg2align(1), TIMEeeg(eeg2align(1)))
    t0 = TIMEevt(evt2align(1)) - TIMEeeg(eeg2align(1));
    if t0~=0
        fprintf('--> A time difference %d samples was found between the event.\n', t0)
        TIMEevt = TIMEevt - t0;
        for ie=1:size(Ed,1)
            Ed{ie,id_time} = Ed{ie,id_time} - t0;
        end
    end
    
    % further aligment if necessary
    if length(idxjumps)>1
        for ij=2:length(idxjumps)
            fprintf('Aligment %d\n',ij)
            fprintf(' - Aligmant based on event %s, event %d in the event file, latency %d \n', TYPEevt{evt2align(ij)}, evt2align(ij), TIMEevt(evt2align(ij)))
            fprintf(' - Aligmant based on event %s, event %d in the EEG structure, latency %d \n', TYPEeeg{evt2align(ij)}, eeg2align(ij), TIMEeeg(eeg2align(ij)))
            t0 = TIMEevt(evt2align(ij)) - TIMEeeg(eeg2align(ij));
            if t0~=0
                fprintf('--> A time difference %d samples was found between the event.\n', t0)
                TIMEevt(evt2align(ij):end) = TIMEevt(evt2align(ij):end) - t0;
                for ie=evt2align(ij):size(Ed,1)
                    Ed{ie,id_time} = Ed{ie,id_time} - t0;
                end
            end
        end
    end
end

%% Identify the events to which information will be added
idxEv = zeros(length(EEG.event),1);
for i=1:length(EEG.event)
    idxEv(i) = any(strcmp(eventtype, strtrim(EEG.event(i).type)));
end
idxEv = find(idxEv);

idxEd = zeros(length(TYPEevt), length(eventtype));
for i=1:length(TYPEevt)
    idxEd(i,:) = strcmp(eventtype, strtrim(TYPEevt{i}));
end
Edev = Ed(logical(any(idxEd,2)),:);


%% Search for extra information for each event in EEG
count = 0;
typeev = cellfun(@strtrim, Edev(:,1), 'UniformOutput', 0);
latev = cell2mat(Edev(:,id_time));
v_diff = [];
for i=1:length(idxEv)
    ieeg = idxEv(i);
    latieeg = EEG.urevent(EEG.event(ieeg).urevent).latency;
    typeeeg = strtrim(EEG.event(ieeg).type);
    
    % find the event of the correct type
    typeOK = find( strcmp( typeeeg, typeev) );
    if ~isempty(typeOK)
        [v, idx] = min( abs(latev(typeOK) - latieeg) );
        ie = typeOK(idx);
        v_diff = [v_diff; v];
        
        if v<=latlim
            % Edit the event
            j = 1;
            while j <= size(Edev,2)
                if all(~strcmp(Ename{j},{'Onset' 'Duration'}))
                    if ~EnameOK(j) && isempty(Edev{ie,j})
                        j = j+2;
                    else
                        if EnameOK(j)
                            jj = j;
                            thefield = Ename{j};
                            j = j+1;
                        elseif ~EnameOK(j) && ~isempty(Edev{ie,j})
                            jj = j+1;
                            thefield =strtrim(Edev{ie,j});
                            j = j+2;
                        end
                        
                        if ~isfield(EEG.event,thefield)
                            for iiev=1:length(EEG.event)
                                EEG.event(iiev).(thefield) = [];
                            end
                            if isfield(EEG,'epoch')
                                for iiep=1:length(EEG.epoch)
                                    EEG.epoch(iiep).(['event' thefield]) = repmat({''},[1 length(EEG.epoch(iiep).event)]);
                                end
                            end
                        end
                        
                        EEG.event(ieeg).(thefield) = strtrim(Edev{ie,jj});
                        if isfield(EEG,'epoch')
                            for iiep=1:length(EEG.epoch)
                                idxEpEv = EEG.epoch(iiep).event==ieeg;
                                if any(idxEpEv)
                                    EEG.epoch(iiep).(['event' thefield]){idxEpEv} = strtrim(Edev{ie,jj});
                                end
                            end
                        end
                        %                     EEGout=pop_editeventvals(EEGout,'changefield',{ieeg thefield strtrim(Edev{ie,jj})});
                    end
                else
                    j = j+1;
                end
                
            end
            
            count = count+1;
            % fprintf('-   New information was added for the event %d %s\n', ieeg, EEG.event(ieeg).type);
        else
            fprintf('WARNING!!! %d samples difference for the event %s - %d - urevent %d \n no information was added \n',...
                v, EEG.event(ieeg).type, ieeg, EEG.event(ieeg).urevent);
        end
    else
        fprintf('WARNING!!! No information was found for the event %s - %d - urevent %d \n',...
            EEG.event(ieeg).type, ieeg, EEG.event(ieeg).urevent);
    end
end
fprintf('\nNew information added to %d events \n\n', count)
% figure, histogram(v_diff)

%% Search for events in the event file no present in the EEG structure
% If there are not events add the first one in the event file
if isempty(EEG.event)
    EEG.event(1).type = Ed{1,id_code};
    EEG.event(1).duration = Ed{1,id_dur} * EEG.srate;
    EEG.event(1).latency =  Ed{1,id_onset} * EEG.srate;
end
  
%get the type and latencies from the EEGLAB structure
TIMEeeg = nan(length(EEG.event),1);
TYPEeeg = cell(length(EEG.event),1);
for ie=1:length(EEG.event)
    TIMEeeg(ie) = EEG.urevent(EEG.event(ie).urevent).latency;
    TYPEeeg{ie} = EEG.urevent(EEG.event(ie).urevent).type;
end

%get the indexes of the extra fields in the event file
id_extra = find(~EnameOK);
id_extraname = id_extra(1:2:end);
id_extravals = id_extra(2:2:end);

% 
% type = cell(length(EEG.event),1);
% lat = nan(length(EEG.event),1);
% for i=1:length(EEG.event)
%     type{i} = EEG.event(i).type;
%     lat(i) = EEG.event(i).latency;
% end

%loop through the events in the event file
ncount = 0;
for i=1:size(Ed,1)
    
    %check if the event exists in the EELAB structure
    idx_type = strcmp(TYPEeeg, TYPEevt{i});
    idx_time = abs(TIMEeeg - TIMEevt(i))<=latlim;
    exists = any(idx_type==1 & idx_time==1);
    
    if ~exists
         
        fnamesextra_i = Ed(i,id_extraname);
        fnamesextra_i(cellfun(@isempty,fnamesextra_i)) = [];
        fvalsextra_i = Ed(i,id_extravals);
        fvalsextra_i(cellfun(@isempty,fvalsextra_i)) = [];
        
        % add the event
        fnamesbase_i = {'type' 'latency' 'urevent'};
        fvalsbase_i = cat(2,TYPEevt(i), TIMEevt(i), length(EEG.urevent)+1);
        
        fnames_i = cat(2, fnamesbase_i, fnamesextra_i);
        fvals_i = cat(2, fvalsbase_i, fvalsextra_i);
        
        evn = length(EEG.event)+1;
        for ff=1:length(fnames_i)
            EEG.event(evn).(fnames_i{ff}) = fvals_i{ff};
        end
        
        % add the urevent - added by Ana Flo 06/12/2021
        fnames_i_ur = {'type' 'latency'};
        fvals_i_ur = cat(2,TYPEevt(i), TIMEevt(i));
        
        evn = length(EEG.urevent)+1;
        for ff=1:length(fnames_i_ur)
            EEG.urevent(evn).(fnames_i_ur{ff}) = fvals_i_ur{ff};
        end
        
        % update the counter
        ncount = ncount+1;
    end
end
fprintf('%d new events were added\n',ncount);

%% Check the structure
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'checkur');


end
