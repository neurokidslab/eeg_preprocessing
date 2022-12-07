% This functions gives a summary of the preprocessing for a dataset
%
% INPUTS:
%
% filenames     It can be a cell array with the complete file names of all 
%               the files to include in the report. In this case, pathIn 
%               should be empty.
%               Otherwise, if pathIn is provided, it should be a string 
%               that serves to retrive all the files in pathIn.
%
% path2data     Path to the files. Only necessary if a string pattern is
%               provided in filenames. It is also the path were the summary
%               is saved
%
% patternname   Pattern to find in the names to create the table 
%               (subjetcs identifier). If empty, the whole files name are 
%               used
% patternn      Number of caracters to use from the begiging of the patter 
%               to create the table. If empty, the name is extracted from 
%               the begign of the patter till the end of the name

function [T] = eega_printsummary(filenames, path2data, patternname, patternn, silent)
if nargin<5
    silent=1;
end
if nargin<4
    patternn=[];
end
if nargin<3
    patternname=[];
end
if nargin<2
    path2data=[];
end
if ischar(filenames) && isempty(path2data)
    error('eega_printsummary: not path to the files was provided')
end

%% List of files
if ischar(filenames)
    ssList = dir(fullfile(path2data,filenames));
    rem = false(length(ssList),1);
    for tf = 1:length(ssList)
        if ssList(tf).name(1)=='.';
            rem(tf) = 1;
        end
    end
    ssList(rem) = [];
    nsbj = length(ssList);
    filenames = cell(nsbj,1);
    for i=1:length(ssList)
        filenames{i} = fullfile(ssList(i).folder, ssList(i).name);
    end
else
    nsbj = length(filenames);
end

%% Do it for all the files
if nsbj>0
    
    SBJ = cell(nsbj,1);
    
    summaryIDX = nan(3,1);
    what = cell(3,1);
    
    Nch = nan(nsbj,length(summaryIDX));
    Nspl = nan(nsbj,length(summaryIDX));
    Nep = nan(nsbj,length(summaryIDX));
    BCT = nan(nsbj,length(summaryIDX));
    CCT = nan(nsbj,length(summaryIDX));
    BT = nan(nsbj,length(summaryIDX));
    BC = nan(nsbj,length(summaryIDX));
    
    for s = 1:nsbj
        
        % -----------------------------------------------------------------
        % Load the data
        [folderSubj,nameSubj,Ext] = fileparts(filenames{s});
        fprintf('\n%s\nSubject %s\n',repmat('-',[1 50]),nameSubj)
        clear EEG
        if strcmp(Ext,'.set')
            EEG = pop_loadset( [nameSubj Ext], folderSubj );
        elseif strcmp(Ext,'.mat')
            EEG = load(fullfile(folderSubj,  [nameSubj Ext]));
            feeg = fieldnames(EEG);
            EEG = EEG.(feeg{1});
        end
        
        % -----------------------------------------------------------------
        % Take the name of the subject
        
        if isempty(patternname)
            SBJ{s} = EEG.filename;
        else
            k = strfind(EEG.filename,patternname);
            if ~isempty(patternn)
                SBJ{s} = EEG.filename(k:k+patternn-1);
            else
                SBJ{s} = EEG.filename(k:end);
            end
        end
        
        % -----------------------------------------------------------------
        % Find the info for the different stages during pre-processing
        
        if isfield(EEG,'summaryhistory')
            
            sumhist = 1;
            
            % Last stage before epoching
            idx = find(EEG.summaryhistory.nep==1);
            if ~isempty(idx)
                summaryIDX(1) = idx(end);
                what{1} = 'continuos data';
            end
            
            % If the data is epoched
            if any(EEG.summaryhistory.nep>1)
                
                % Before rejecting bad epochs
                idx = find(diff(EEG.summaryhistory.nep)<0);
                if ~isempty(idx)
                    summaryIDX(2) = idx(end);
                    what{2} = 'epoched data';
                end
                
                % Finale stage
                %         idx = find(diff(EEG.summaryhistory.nep)<0);
                idx  = length(EEG.summaryhistory.nep);
                if ~isempty(idx)
                    summaryIDX(3) = idx(1);
                    what{3} = 'epoched data without bad epochs';
                end
            end
            
        else
            sumhist = 0;
        end
       
        % -----------------------------------------------------------------
        % Take the important info
        
        if ~sumhist && isfield(EEG,'summary')
            Nch(s,1) = EEG.summary.nch;
            Nspl(s,1) = EEG.summary.nsmpl;
            Nep(s,1) = EEG.summary.nep;
            BCT(s,1) = 100* EEG.summary.BCTn / (Nspl(s)*Nch(s)*Nep(s));
            CCT(s,1) = 100* EEG.summary.CCTn / (Nspl(s)*Nch(s)*Nep(s));
            BT(s,1) = 100* EEG.summary.BTn / (Nspl(s)*Nep(s));
            BC(s,1) = 100* EEG.summary.BCn / (Nch(s)*Nep(s));
        elseif sumhist
            for j=1:length(summaryIDX)
                if ~isnan(summaryIDX(j))
                    Nch(s,j) = EEG.summaryhistory.nch(summaryIDX(j));
                    Nspl(s,j) = EEG.summaryhistory.nsmpl(summaryIDX(j));
                    Nep(s,j) = EEG.summaryhistory.nep(summaryIDX(j));
                    BCT(s,j) = 100* EEG.summaryhistory.BCTn(summaryIDX(j)) / (Nspl(s,j)*Nch(s)*Nep(s,j));
                    CCT(s,j) = 100* EEG.summaryhistory.CCTn(summaryIDX(j)) / (Nspl(s,j)*Nch(s)*Nep(s,j));
                    BT(s,j) = 100* EEG.summaryhistory.BTn(summaryIDX(j)) / (Nspl(s,j)*Nep(s,j));
                    BC(s,j) = 100* EEG.summaryhistory.BCn(summaryIDX(j)) / (Nch(s,j)*Nep(s,j));
                    
                end
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % Remove empty tables
    idxrmv = all(isnan(Nch),1);
    summaryIDX(idxrmv) = [];
    what(idxrmv) = [];
    Nch(:,idxrmv) = [];
    Nspl(:,idxrmv) = [];
    Nep(:,idxrmv) = [];
    BCT(:,idxrmv) = [];
    CCT(:,idxrmv) = [];
    BT(:,idxrmv) = [];
    BC(:,idxrmv) = [];
    
    % ---------------------------------------------------------------------
    % Print the table
    VarNames = {'Ch' 'Smpl' 'Ep' 'BCTx100' 'CCTx100' 'BTx100' 'BCx100'};
    T = cell(1,length(summaryIDX));
    for j=1:length(summaryIDX)
        T{j} = table(Nch(:,j),Nspl(:,j),Nep(:,j),BCT(:,j),CCT(:,j),BT(:,j),BC(:,j), 'VariableNames',VarNames,'RowNames',SBJ);
        writetable(T{j},fullfile(path2data,sprintf('%s.txt',what{j})),'Delimiter','\t','WriteRowNames',1);
        if ~silent
            fprintf('\n%s\nSummary for all subject: %s \n\n',repmat('-',[1 50]),what{j})
            disp(T{j})
        end
    end
else
    warning('No data files were found!')
end

end