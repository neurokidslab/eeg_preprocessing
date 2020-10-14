% This function imports to Brainstorm
%
%   INPUTS
%   filenames       = string with the files
%   pathIn          = path where the input files are
%   pathOut         = path where the output file is saved. If the folder
%                   does not exist, it creastes it
%   CNDs            = name of the conditions to import
%   filenum         = verctor indicating the file numbers to import
%
% -------------------------------------------------------------------------
% Ana Flo March 2017
% -------------------------------------------------------------------------

function eega_Import2Brainstorm( filenames, pathIn, factor, varargin )

%% ========================================================================
%% Parameters

% Default parameters
P.DataField = 'data';
P.GenericSubjectName = 'S';
P.IncSbj = [];

% Optional parameters
if mod(length(varargin),2)==1
    warning('Optional parameters come by pairs')
    return
end
if ~isempty(varargin)
    for f=1:2:length(varargin)
        Pop.(varargin{f}) = varargin{f+1};
    end
    fff = fieldnames(P);
    for f=1:numel(fff)
        if isfield(Pop,fff{f})
            P.(fff{f}) = Pop.(fff{f});
        end
    end
end
DataField = P.DataField;
GenericSubjectName = P.GenericSubjectName;
IncSbj = P.IncSbj;

%% ========================================================================
%% List of files
ssList = dir(fullfile(pathIn,filenames));
rem = false(length(ssList),1);
for tf = 1:length(ssList)
    if ssList(tf).name(1)=='.';
        rem(tf) = 1;
    end
end
ssList(rem) = [];
subjects    = 1:length(ssList);
if isempty(IncSbj)
    for s = 1:length(subjects)
        [~,nameSubj,Ext] = fileparts(ssList(subjects(s)).name);
        IncSbj{s} = nameSubj;
    end
end

%% ========================================================================
%% Do it for all the files
S = 0;
for s = 1:length(subjects)
    
    % ---------------------------------------------------------------------
    % Subject and file names
    [~,nameSubj,Ext] = fileparts(ssList(subjects(s)).name);
    
    Si=0;
    si=0;
    while Si==0 && si<length(IncSbj)
        si=si+1;
        k = strfind(nameSubj,IncSbj{si});
        if ~isempty(k)
            Si=si;
        end
    end
    
    if Si
        S = S+1;
        clear EEG
        fprintf('\n%s\nSubject %s\n --> added\n%s\n',repmat('-',[1 50]),nameSubj,repmat('-',[1 50]))
        
        % ----------------------------------------------------------------------
        % Load the data
        if strcmp(Ext,'.set')
            EEG = pop_loadset( ssList(s).name, pathIn );
        elseif strcmp(Ext,'.mat')
            EEG = load(fullfile(pathIn,  ssList(s).name));
            fff=fieldnames(EEG);
            EEG = EEG.(fff{1});
        end
        
        % ----------------------------------------------------------------------
        % Remove the reference
        if size(EEG.data,1)==129
            EEG.data = EEG.data(1:128,:,:);
            EEG.chanlocs = EEG.chanlocs(1:128);
            EEG.nbchan = 128;
        end
        
        % ---------------------------------------------------------------------
        % Get the conditions
        idxF=1;
        while ~strcmp(EEG.F{idxF}.name,factor)
            idxF=idxF+1;
        end
        CNDs = EEG.F{idxF}.val;
        
        % ---------------------------------------------------------------------
        % Add the subject
        SbjNewName = [GenericSubjectName num2str(S,'%02.0f')];
        db_add_subject(SbjNewName);
        
        Comment = nameSubj;
        
        nT = size(EEG.(DataField),3);
        for t = 1:nT
            
            matimport = squeeze(EEG.(DataField)(:,:,t));
            matimport(isnan(matimport)) = 0;
            matimport = double(matimport);
            dimchan = find(size(matimport)== EEG.nbchan);
            if dimchan==2
                matimport = matimport';
            end
            
            DataMat.F           = matimport;
            DataMat.Channel     = EEG.chanlocs;
            DataMat.ChannelFlag = ones(EEG.nbchan, 1);
            DataMat.Time        = EEG.times/1000;
            DataMat.Comment     = Comment;
            DataMat.DataType    = 'recordings';
            DataMat.Device      = 'Unknown';
            DataMat.nAvg        = 1;
            
            iStudies = db_add_condition(SbjNewName, CNDs{t});
            db_add_data(iStudies, CNDs{t});
            db_add(iStudies,DataMat);
            
        end
        
    end
end



