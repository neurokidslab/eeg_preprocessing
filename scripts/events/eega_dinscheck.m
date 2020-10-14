% This function check if some DINs were missing based on the distance
% between DINs inside an epoch and the numbe rof DINs that there should be
% in an epoch. It is asumed that each DIN generates two DIN events, one
% when the DIN starts and one when it finishes.
%
% INPUTS:
%
%   EEG         = EEGLab structure with the EEG data
%   dinevent    = cell array of strings with specific event to look for
%   dinduration = duration of the DIN in ms
%   dinnumber   = number of DINs that should be present
%
% OUTPUT:
%
%   EEG
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------

function EEGout = eega_dinscheck( EEG, dinevent, dinduration, dinnumber, dindist )

if nargin<5
    dindist = [];
end

tolerance = 5*(1/EEG.srate*1000); % tolerance for checking the time (ms)

nep = numel(EEG.epoch);

% Get the fields of the events
fff = fieldnames(EEG.event);
for f = 1:numel(fff);
    dinevport.(fff{f}) = [];
end
dinevport.type = dinevent{1};
dinevport.Code = 'DIN1';
dinevport.Label = 'AUDIO';
dinevport.Type = 'Stimulus Event';
dinevport.Track = 'DIN 1';


ievaddtot = 0; % counter for the total number of events added
EEGout = EEG;
if nep==0
    nblocks = 1;
else
    nblocks = nep;
end
for ep = 1:nblocks
    
    fprintf('%s\n',repmat('.',50,1))
    fprintf('Epoch %d \n',ep);
    fprintf('%s\n',repmat('.',50,1))
    
    % --------------------------------------
    % Identify the DINs and find its latency
    if nep==0,
        etypes = cell(size(EEG.event,2),1);
        elatency = nan(size(EEG.event,2),1);
        for ie=1:size(EEG.event,2)
            etypes{ie} = EEG.event(ie).type;
            elatency(ie) = EEG.event(ie).latency/EEG.srate*1000;
        end
    else
        etypes      = EEG.epoch(ep).eventtype;
        elatency    = cell2mat(EEG.epoch(ep).eventlatency);
    end
    isdin       = strcmp(strtrim(etypes),dinevent);
    etypes      = etypes(isdin);
    elatency    = elatency(isdin);
    
    
    % ----------------------------------------
    % Calculate the DINs duration and distance
    ndinsev = numel(etypes);                    % total number of DIN events
    ndins = numel(etypes)/2;                    % total number of DINs
    idin_lat = elatency(1:2:ndins*2);    % begining of DINs
    fdin_lat = elatency(2:2:ndins*2);    % end of DINs
    dins_d = (fdin_lat(1:floor(ndins)) - idin_lat(1:floor(ndins)));             % duration of DINs
    dins_diff = diff(idin_lat);                 % distance between consecutive DINs
    diffbig = (dins_d > (dinduration+tolerance));
    diffsmall = (dins_d < (dinduration-tolerance));
    
    % ---------------------------
    % Check if averything is fine
    if any(diffbig) || any(diffsmall) || ndins ~= dinnumber
        fprintf('WARNING!!! There is a problem with the DINs \n');
        problem = 1;
    else problem = 0;
    end
    
    % -------------------------------------
    % Try to identify and solve the problem
    newdintimes = nan(dinnumber,1);
    newev_latency = nan(ndinsev,1);
    if problem
        fprintf('Trying to identify and solve the problem... \n');
        
        idin        = 1;    % counter of DINs in the epoch
        iev         = 1;    % counter of events in the epoch
        ievaddep    = 0;    % counter of the events added in the epoch
        
        if ndinsev==1
            
            % Show the information and ask what to do
            fprintf('## Only 1 DIN event was found!!\n');
            fprintf('   Event time %d\n', elatency(iev))
            
            add = 'x';
            while ~any(add=='yn');
                add = input('   Do you want to add a DIN ? (y/n) ','s');
            end
            
            if strcmp(add,'y')
                evis = 0;
                while ~any(evis==[1 2]);
                    evis = input('   Is the event the onset (1) or the offset (2) of the DIN ? ');
                end
                if evis==1
                    newev_latency = elatency(iev) + dinduration;
                    fprintf('   An event will be added %d ms after \n', dinduration);
                elseif evis==2
                    newev_latency = elatency(iev) - dinduration;
                    fprintf('   An event will be added %f ms before \n', dinduration);
                end
                
            else
                fprintf('   Nothing will be done \n');
            end
            
        else
            
            % Show all the information
            fprintf('- %f were found \n',ndins);
            fprintf('- Maximun DIN duration: %f \n',max(dins_d));
            fprintf('- Minimun DIN duration: %f \n',min(dins_d));
            fprintf('- Maximun time between consecutive DINs: %f \n',max(dins_diff));
            fprintf('- Minimun time between consecutive DINs: %f \n',min(dins_diff));
            
            % Check that there are not missing DINs
            while iev < ndinsev
                
                diffi = elatency(iev+1) - elatency(iev);
                
                if (diffi > (dinduration+tolerance))
                    
                    % Calculete the distance to the closer events and DINs
                    if iev > 1, tev_pre = elatency(iev-1) - elatency(iev);
                    else tev_pre = nan;
                    end
                    if iev < ndinsev, tev_next = elatency(iev+1) - elatency(iev);
                    else tev_next = nan;
                    end
                    if idin > 1, tdin_pre = newdintimes(idin-1) - elatency(iev);
                    else tdin_pre = nan;
                    end
                    
                    % Show the information and ask what to do
                    fprintf('## Problem with event %d DIN %d !!\n', iev, idin );
                    fprintf('   Closer events:     %d ms and %d ms \n', tev_pre, tev_next);
                    fprintf('   Previous good DIN: %d ms \n', tdin_pre);
                    
                    add = 'x';
                    while ~any(add=='yn');
                        add = input('   Do you want to add a DIN ? (y/n) ','s');
                    end
                    if strcmp(add,'y')
                        ievaddtot    = ievaddtot + 1;
                        ievaddep  = ievaddep + 1;
                        
                        evis = 0;
                        while ~any(evis==[1 2]);
                            evis = input('   Is the event the onset (1) or the offset (2) of the DIN ? ');
                        end
                        if evis==1
                            newev_latency(ievaddep) = elatency(iev) + dinduration;
                            newdintimes(idin) = elatency(iev);
                            fprintf('   An event will be added %d ms after \n', dinduration);
                        elseif evis==2
                            newev_latency(ievaddep) = elatency(iev) - dinduration;
                            newdintimes(idin) = elatency(iev) - dinduration;
                            fprintf('   An event will be added %f ms before \n', dinduration);
                        end
                        
                        idin = idin + 1;
                        iev = iev + 1;
                    else
                        fprintf('   Nothing will be done \n');
                        
                        newdintimes(idin) = elatency(iev);
                        idin = idin + 1;
                        iev = iev + 1;
                    end
                    
                else
                    newdintimes(idin) = elatency(iev);
                    idin = idin + 1;
                    iev = iev + 2;
                end
                
            end
            
            newev_latency(isnan(newev_latency)) = [];
            
            % Add the events
            if ~isempty(newev_latency)             
                eegnew = EEG;
                if nep==0
                    lat0 = 0;
                else
                    i_ev_1 = EEG.epoch(ep).event(1,1);
                    lat0 = EEG.event(i_ev_1).latency;
                end
                
                for i=1:numel(newev_latency)
                    newlatency = newev_latency(i)/1000*EEG.srate + lat0;
                    
                    i_new = size(eegnew.event,2)+1;
                    eegnew.event(i_new) = dinevport;
                    eegnew.event(i_new).latency = newlatency;
                    if nep~=0
                        eegnew.event(i_new).epoch = ep; 
                    end
                end
                eegnew = eeg_checkset(eegnew, 'eventconsistency');           
            end
            
            % Recalculate everything
            if ~isempty(newev_latency)
                if nep==0,
                    etypesnew = cell(size(eegnew.event,2),1);
                    elatencynew = nan(size(eegnew.event,2),1);
                    for ie=1:size(eegnew.event,2)
                        etypesnew{ie} = eegnew.event(ie).type;
                        elatencynew(ie) = eegnew.event(ie).latency/eegnew.srate*1000;
                    end
                else
                    etypesnew      = eegnew.epoch(ep).eventtype;
                    elatencynew    = cell2mat(eegnew.epoch(ep).eventlatency);
                end
                isdinnew       = strcmp(strtrim(etypesnew),dinevent);
                etypesnew      = etypesnew(isdinnew);
                elatencynew    = elatencynew(isdinnew);
                ndinsevnew     = numel(etypesnew);                      % total number of DIN events
                ndinsnew       = numel(etypesnew)/2;                    % total number of DINs
                idin_latnew    = elatencynew(1:2:ndinsnew*2);    % begining of DINs
                fdin_latnew    = elatencynew(2:2:ndinsnew*2);    % end of DINs
                dins_dnew      = (fdin_latnew(1:floor(ndinsnew)) - idin_latnew(1:floor(ndinsnew)));           % duration of DINs
                dins_diffnew   = diff(idin_latnew);                     % distance between consecutive DINs
                diffbignew     = (dins_dnew > (dinduration+tolerance));
                diffsmallnew   = (dins_dnew < (dinduration-tolerance));
                
                %             % If the number of DINs is bad and the distance between DINs is
                %             % provided, it tries to solve the problem by adding DINs
                %             if (ndinsnew ~= dinnumber) && (~isempty(dindist))
                %                 fprintf('   It seems there are missing DINs.\n DINs will be added according to the DINs distance\n')
                %                 idx_1 = dins_diffnew > 2*(dindist-tolerance);
                %                 idx_1 = find(idx_1);
                %                 for i=1:numel(idx_1)
                %                     iev = idx_1(i) - 1;
                %                     newev_latency(ievaddep) = elatency(iev) + dinduration;
                %                 end
                
                
                % Check if averything is fine
                if any(diffbignew) || any(diffsmallnew) || ndinsnew ~= dinnumber
                    fprintf('   ## The problem could not be solved. Nothing was modified.\n');
                else
                    EEGout = eegnew;
                    fprintf('   ## The problem was solved!! \n');
                    elatency = elatencynew;
                    ndins = ndinsnew;
                    dins_d = dins_dnew;
                    dins_diff = dins_diffnew;
                    problem = 0;
                end
            end
            
        end
        
    else
        fprintf('DINs seem to be fine!! \n');
    end
    
    if ndinsev <= 2
        fprintf('- %d DIN was found \n', ndins);
        fprintf('- Event time %d\n', elatency(1))
    else
        fprintf('- %d DINs were found \n',ndins);
        fprintf('- Maximun DIN duration: %f \n',max(dins_d));
        fprintf('- Minimun DIN duration: %f \n',min(dins_d));
        fprintf('- Maximun time between consecutive DINs: %f \n',max(dins_diff));
        fprintf('- Minimun time between consecutive DINs: %f \n',min(dins_diff));
    end
    
    if problem
        disp('   WARNING!!! The problem was not solved')
        
    end
end




