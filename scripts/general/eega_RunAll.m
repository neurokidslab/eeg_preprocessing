% -------------------------------------------------------------------------
% This function runs the series of functions specified as part of the input
% for all the files found in the input folder and saves the output in the 
% output folder
%
% INPUTS
% filenames     string with the files
% pathIn        path where the input files are
% pathOut       path where the outpu files are saved. If the folder does
%               not exist, it creastes it
%
% OPTIONAL INPUTS
% prename       string with something to add before the filename
% runto         If 1, the functions are run to all the files found. 
%               If 2, the functions are run to the files still not analyzed
% saveformat    'set' | 'mat'
% function_i    (varargin{i}) name of the function to run 
% inputs_i      (varargin{i+1}) inputs for function i 
%
% USAGE
% eega_RunAll(filenames, pathIn, pathOut,...
%       'eega_refavg', {},...
%       'eega_normalization', {'epochs','each','electrodes','all'},...
%       'prename','refn_')
% the functions to apply have to be of the form
%
% EEG = function(EEG, other inputs)
%
% -------------------------------------------------------------------------

function eega_RunAll( filenames, pathIn, pathOut, varargin )

%% =========================================================================
%% Parameters

% Defalut values
prename = '';
saveformat = 'set';
doeegall = 0;
runto = 0;

% Take the optional parameters
if mod(length(varargin),2)==1
    error('eega_RunAll: Optional parameters come by pairs')
else
    idop = zeros(size(varargin));
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'runto'
                runto = varargin{i+1};
                idop(i:i+1) = 1; 
            case 'prename'
                prename = varargin{i+1};
                idop(i:i+1) = 1; 
            case 'saveformat'
                saveformat = varargin{i+1};
                idop(i:i+1) = 1;
            case 'saveeegall'
                doeegall = varargin{i+1};
                idop(i:i+1) = 1;
        end
    end
end

% Take the functions
thefunctions = varargin(~idop);

%% =========================================================================
%% Ask which on which files de action wants to be perfomed
fprintf( 'What would you like to do?\n')
fprintf( '	1. Run it for all the files found and overwrite previous output files (1)\n')
fprintf( '	2. Run it only for the new files (2)\n')
fprintf( '	3. Ask in each case if a previous analysis is found (3)\n')
while ~any(runto==[1 2 3]);
    runto = input(sprintf('Chose an option 1, 2 or 3: '));
end

%% ========================================================================
%% List of files
addpath(pathIn);
ssList = dir(fullfile(pathIn,filenames));
rem = false(length(ssList),1);
for tf = 1:length(ssList)
    if ssList(tf).name(1)=='.';
        rem(tf) = 1;
    end
end
ssList(rem) = [];
subjects    = 1:length(ssList);

if isempty(ssList)
    warning('eega_RunAll: No files were found')
    return
end

%% =========================================================================
%% Create the output folder if it does not exis
if ~isempty(pathOut)
    if exist(pathOut,'dir')==0
        mkdir(pathOut)
    end
    addpath(pathOut);
end

%% ========================================================================
%% Do it for all the files

for s = 1:length(subjects)  % loop over subjects
    
    t0 = tic;
    clear EEG
    do = 1;
    
    % ----------------------------------------------------------------------
    % Subject and file names
    [~,nameSubj,ext] = fileparts(ssList(subjects(s)).name);
    if strcmp(saveformat,'set')
        nameOut = sprintf('%s%s.set',prename,nameSubj);
    elseif strcmp(saveformat,'mat')
        nameOut = sprintf('%s%s.mat',prename,nameSubj);
    else
        nameOut = [];
    end
    
    fprintf('\n%s\nSubject: %s\n%s\n',repmat('-',[1 50]),nameSubj,repmat('-',[1 50]))
    
    if runto==3 || runto==2
        % ----------------------------------------------------------------------
        % Check if a file with the same name already exists
        fullnameOut = fullfile(pathOut, nameOut);
        if exist(fullnameOut,'file') && runto==3
            fprintf( 'ATTENTION!! a file was found for this subject \n' )
            resp = input(sprintf('  Do you want to over-write it? YES(y)/NO(n):'), 's');
            if strcmp('n',resp); do = 0; end
        elseif exist(fullnameOut,'file') && runto==2
            do = 0;
        end
    end
    if isempty(nameOut)
        do = 1;
    end
    
    if do
        % ----------------------------------------------------------------------
        % Load the data
        if strcmp(thefunctions{1},'eega_importdata') ||...
                strcmp(thefunctions{1},'pop_readegi') ||...
                strcmp(thefunctions{1},'pop_readegimff')
            
            thefun = str2func(thefunctions{1});
            inputs = thefunctions{2};           
            EEG = thefun( fullfile(pathIn,ssList(subjects(s)).name), inputs{:} );
            fun_i = 3;
            dosave = 1;
            
        else
            if strcmp(ssList(subjects(s)).name(end-2:end),'set')
                EEG = pop_loadset( ssList(subjects(s)).name, pathIn );
                fun_i = 1;
            elseif strcmp(ssList(subjects(s)).name(end-2:end),'mat')
                EEG = load( fullfile(pathIn, ssList(subjects(s)).name) );
                EEG = EEG.EEG;
                fun_i = 1;
            end
            dosave = 0;
        end
        fprintf('\n')

        % ----------------------------------------------------------------------
        % Run all the functions
        while fun_i<length(thefunctions)
            
            thefun = str2func(thefunctions{fun_i});
            inputs = thefunctions{fun_i+1};
            
            if ~isempty(EEG.data)
                if nargout(thefun)==0
                    thefun(EEG, inputs{:});
                else
                    EEG = thefun(EEG, inputs{:});
                    dosave = 1;
                end
            else
                warning('Empty data!')
            end
            
            fun_i = fun_i+2;            
        end
        
        t = toc(t0);
        hours   = floor(t / 3600);
        t = t - hours * 3600;
        mins    = floor(t / 60);
        t = t - mins * 60;
        secs    = floor(t);
        fprintf('\nTotal time fo subject %s : %dh : %dm : %ds\n',nameSubj, hours,mins,secs)
        fprintf('%s\n',repmat('-',[1 50]))
        
        % ----------------------------------------------------------------------
        % Save
        if dosave && ~isempty(pathOut)
            if strcmp(saveformat,'set')
                pop_saveset( EEG, nameOut, pathOut);
            elseif strcmp(saveformat,'mat')
                save(fullfile(pathOut,nameOut),'EEG')
            end
        end
        
    end
    
end % loop over subjects
if doeegall
    save(fullfile(pathOut,sprintf('%sEEGALL',prename),EEGALL))
end
end
