function EEG = eega_addchinfo(EEG, filechanloc, varargin)

try
    chanlocs = readlocs( filechanloc, varargin{:} );
    check = 1:size(chanlocs,2);
    rm = zeros(1, size(chanlocs,2));
    for k=1:numel(check)
        if strcmp(chanlocs(check(k)).type,'FID')
            rm(check(k)) = 1;
        end
    end
    chanlocs(logical(rm)) = []; % remove the localization
catch
    warning('eega_addchinfo: The channels layout could not be loaded')
    chanlocs = [];
end
if size(EEG.data,1)==(size(chanlocs,2)-1)
    chanlocs(end) = [];
end % remove the reference
EEG.chaninfo.filename = filechanloc;
EEG.chanlocs = chanlocs;

fprintf('\n')

end