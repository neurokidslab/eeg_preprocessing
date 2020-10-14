% This function imports data
%
%   INPUTS
%   pathIn      = path where the input files are
%
% -------------------------------------------------------------------------
% Ana Flo June 2020
% -------------------------------------------------------------------------

function EEG = eega_importdata( filename )

%file parts
[~,nameSubj,ext] = fileparts(filename);

%read the data
fprintf('Importing file %s ...\n', filename)
fprintf('File type %s \n', ext)
switch ext
    case '.raw'
        EEG = pop_readegi(filename);
    case '.mff'
        EEG = pop_readegimff(filename);
end

%add the subject name to the EEGLAB structure
EEG.filename = nameSubj;

%check the set
EEG = eeg_checkset(EEG);

end