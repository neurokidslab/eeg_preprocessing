function EEG = eega_eventextra2fields(EEG, extras, readformat)

if nargin<2 || isempty(extras)
    extras = 'codes';
end
if nargin<3
    readformat = [];
end
    

%% Make the codes in EEG.event(i).codes a field
for i=1:length(EEG.event)
    fnames = fieldnames(EEG.event(1));
    if any(strcmp(fnames,extras))
        for k=1:size(EEG.event(i).(extras),1)
            newfff = EEG.event(i).(extras){k,1};
            newval = EEG.event(i).(extras){k,2};
            if ~isempty(newval) && ischar(newval) && ~isempty(readformat)
                newval = textscan(newval,readformat);
                newval = newval{1};
            end
            if ~isempty(newfff)
                EEG.event(i).(newfff) = newval;
            end
        end
        EEG.event(i).(extras) = [];
    end
end

%% Check the structure
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'checkur');

end