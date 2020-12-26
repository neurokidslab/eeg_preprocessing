% -------------------------------------------------------------------------
% This function identifies and defines the factors for each epoch
% The defined factors can later be used to determine conditions and average
% across trials
%
% INPUTS
% EEG           EEG structure
% EvFieldFactor cell array with the names of the events used to define 
%               factors. Its length defines the total number of factors, nF
%
% OPTIONAL INPUTS
% EvNumber      number of the event to look in each epoch in order to
%               determine the value for each factor. It can be empty, it 
%               can be a number, or 'end'. If it is empty (default )the 
%               event at latancy = 0 is used. If it is a vector of length 1, 
%               the same event is used for all factors, if it is a vector 
%               of length nF, different events are used. If it is 'end',
%               the last event in the epoch is used
% EvName        string with an event field to identify the relevant 
%               event in the epoch to define the factors. If it is provided
%               EvNumber is ignored. 
%                   It come together with EvValue
% EvValue       value of the EvName to identify the relevant event
% FactorsName   name of the field where the factors are stored in EEF.
%               Default 'F'
% Factors       table with the factors
%
% OUTPUTS
% EEG           EEG structure
%
% -------------------------------------------------------------------------

function EEG = eega_definefactors(EEG, EvFieldFactor, varargin)

%% ------------------------------------------------------------------------
%% Parameters
P.Factors = [];
P.FactorsName = 'F';
P.EvNumber = [];
P.EvName = 'eventlatency';
P.EvValue = 0;
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_definefactors: Non recognized inputs')
end
if ~isempty(P.EvNumber) && ~iscell(P.EvNumber)
    P.EvNumber = {P.EvNumber};
end
if ~isempty(P.EvName) && ~iscell(P.EvName)
    P.EvName = {P.EvName};
end
if ~isempty(P.EvValue) && ~iscell(P.EvValue)
    P.EvValue = {P.EvValue};
end
if iscell(P.EvNumber) && ~all(strcmp('end',P.EvNumber))
    sprintf('The EvNumber has to be a number or ''end''')
end
if ~isempty(P.EvNumber) && ~isempty(P.EvName) && ~isempty(P.EvValue)
    sprintf('The EvNumber input will be ignored')
    P.EvName = [];
    P.EvValue = [];
end
if (~isempty(P.EvName) && isempty(P.EvValue)) || (isempty(P.EvName) && ~isempty(P.EvValue))
    sprintf('Both EvName and EvValue have to be provided')
end
if (length(P.EvName)~=length(P.EvValue))
    sprintf('EvName and EvValue have to have the same length')
end

%% ------------------------------------------------------------------------
%% Define the subjects based on the info of each epoch
fprintf('### Defining factors ###\n')

nEp = length(EEG.epoch);
nF = numel(EvFieldFactor);
if ~isempty(P.EvNumber) && (length(P.EvNumber)==1)
    P.EvNumber = repmat(P.EvNumber(:),[nF 1]);
end
if ~isempty(P.EvName) && (length(P.EvName)==1)
    P.EvName = repmat(P.EvName(:),[nF 1]);
end
if ~isempty(P.EvValue) && (length(P.EvValue)==1)
    P.EvValue = repmat(P.EvValue(:),[nF 1]);
end

F = cell(nF+1,1);
for k=1:nF
    F{k}.name = EvFieldFactor{k};
    G = cell(1,nEp);
    for e = 1:nEp
        if ~iscell(EEG.epoch(e).eventtype)
            v_ek = EEG.epoch(e).(EvFieldFactor{k});
        else
            % find the relevents event
            if isempty(P.EvNumber)
                if length(EEG.epoch(e).eventlatency)==1
                    eventnumber = 1;
                else
                    idxev = false(size(EEG.epoch(e).(P.EvName{k})));
                    for i=1:length(EEG.epoch(e).(P.EvName{k}))
                        if ischar(P.EvValue{k})
                            idxev(i) = strcmp(strtrim(EEG.epoch(e).(P.EvName{k}){i}), strtrim(P.EvValue{k}));
                        else
                            idxev(i) = EEG.epoch(e).(P.EvName{k}){i}==P.EvValue{k};
                        end
                        
                    end
                    eventnumber = find(idxev);
                end
            elseif iscell(P.EvNumber) && strcmp(P.EvNumber{k},'end')
                eventnumber = length(EEG.epoch(e).event);
            else
                eventnumber = P.EvNumber(k);
            end
            
            
%             
%             if isempty(P.EvNumber)
%                 if length(EEG.epoch(e).eventlatency)==1
%                     eventnumber = 1;
%                 else
%                     eventnumber = find( cellfun(@eq, EEG.epoch(e).eventlatency, repmat({0},size(EEG.epoch(e).eventlatency)) ) );
%                 end
%             elseif iscell(P.EvNumber) && strcmp(P.EvNumber{k},'end')
%                 eventnumber = length(EEG.epoch(e).event);
%             else
%                 eventnumber = P.EvNumber(k);
%             end
            if length(EEG.epoch(e).(EvFieldFactor{k}))>=eventnumber
                v_ek = EEG.epoch(e).(EvFieldFactor{k})(eventnumber);
            elseif ~isempty(EEG.epoch(e).(EvFieldFactor{k}))
                v_ek = EEG.epoch(e).(EvFieldFactor{k})(end);
            else
                warning('no events in epoch %d',e)
                v_ek = [];
            end
        end
        if iscell(v_ek)
            v_ek =v_ek{1};
        end
        if isnumeric(v_ek)
            v_ek = num2str(v_ek);
        end
        v_ek = matlab.lang.makeValidName(strtrim(v_ek));
        G{e} = v_ek;
        
        
    end
    F{k}.val = unique(G);
    F{k}.g = nan(1,nEp);
    for e = 1:nEp
        F{k}.g(e) = find(strcmp(G{e},F{k}.val));
    end
    
end
F{nF+1}.name = 'Epoch';
F{nF+1}.val = cell(1,nEp);
F{nF+1}.g = (1:nEp);
for e = 1:nEp
    if isfield(EEG.epoch(e), 'eventEpoch')
        if ~isempty(EEG.epoch(e).eventEpoch)
            F{nF+1}.val{e} = genvarname(EEG.epoch(e).eventEpoch{1});
        end
    else
        F{nF+1}.val{e} = genvarname(num2str(e));
    end
end

% Defining new factors if provided
if ~isempty(P.Factors)
    nfnames = P.Factors.Properties.VariableNames;
    idEvent = ismember(nfnames(:),EvFieldFactor);
    tEvent = P.Factors(:,idEvent);
    tNF = P.Factors(:,~idEvent);
    Fnew=cell(size(tNF,2),1);
    for f=1:length(Fnew)
        Fnew{f}.name = tNF.Properties.VariableNames(f);
        Fnew{f}.val = unique(tNF{:,f})';
        Fnew{f}.g = zeros(1,nEp);
    end
    
    % find the events and the new factors
    idFE = zeros(1,size(tEvent,2));
    for j =1:size(tEvent,2)
        i=1;
        while ~idFE(j) && i<=length(F)
            if strcmp(F{i}.name,tEvent.Properties.VariableNames(j))
                idFE(j) = i;
            else
                i=i+1;
            end         
        end 
        if ~idFE(j)
            error('event %s not found',tEvent.Properties.VariableNames(j))
        end
    end
    
    % determine for each epoch
    for e = 1:nEp
        events = cell(1,size(tEvent,2));
        for j =1:size(tEvent,2)
            events(j) = F{idFE(j)}.val(F{idFE(j)}.g(e));
        end
        id = logical(all(strcmp(repmat(events,[size(tEvent,1) 1]),tEvent{:,:}),2));
        if any(id)
            for f=1:length(Fnew)
                gi = find(strcmp(tNF{id,f},Fnew{f}.val));
                Fnew{f}.g(e) = gi;
            end
        else
            warning('Trial %d could not be categorized',e)
            for f=1:length(Fnew)
                Fnew{f}.g(e) = nan;
                
            end
        end
    end
    
    F = cat(1,F(:),Fnew(:));
end

EEG.(P.FactorsName) = F;

fprintf('\n')

end


