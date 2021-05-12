% -------------------------------------------------------------------------
% This functions merges different EEG structures into one by concatenating
% the data in the third dimension
%
% INPUTS
% filenames     string with the files names
% pathIn        path where the input files are
%
% OPTIONAL INPUTS
% DataField     fields in EEG containig the data. It can be a cell with 
%               multiple fields. In that case the baseline correction is 
%               applied to all fileds. Default {'data'}
% FactorField   field in EEG containg the factors characterizing epochs. 
%               Default 'F'
% FactorFieldnames cell with the names of the factors to keep in the merged
%               structure.
% dim2cat       dimension along which data is concatenated. Default last
% IncSbj        defines the subjects to concatenate. It can be a string
%               'all' then all files compatible with filenames will be
%               concatenated. It can be a cell array with the names (or
%               part of the names) of the subjects to include
%
% OUTPUTS
% EEGALL        EEG structure with all the data
%
% -------------------------------------------------------------------------

function EEGALL = eega_MergeData( filenames, pathIn, varargin )

%% ========================================================================
%% Parameters

P.DataField     = {'data'};
P.FactorField   = {'F'};
P.FactorFieldnames   = [];
P.IncSbj        = 'all';
P.dim2cat       = [];
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_MergeData: Non recognized inputs')
end

if ischar(P.DataField)
    P.DataField = {P.DataField};
end

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
if ischar(P.IncSbj) && strcmp(P.IncSbj,'all')
    allsbj=1;
else
    allsbj=0;
    Sbjnames = P.IncSbj;
end

%% ========================================================================
%% Do it for all the files
EEGALL  = [];
S = 0;
SummaryHistory = [];
if ~isempty(subjects)
    for s = 1:length(subjects)
        
        % ---------------------------------------------------------------------
        % Subject and file names
        [~,nameSubj,Ext] = fileparts(ssList(subjects(s)).name);
        fprintf('\n%s\nSubject %s\n',repmat('-',[1 50]),nameSubj)
        
        Si=0;
        si=0;
        while Si==0 && si<length(P.IncSbj)
            si=si+1;
            if allsbj
                Si=si;
                Sbjnames{si} = nameSubj;
            else
                k = strfind(nameSubj,P.IncSbj{si});
                if ~isempty(k)
                    Si=si;
                end
            end
        end
        
        if Si
            
            % ----------------------------------------------------------------------
            % Load the data
            clear EEG
            if strcmp(Ext,'.set')
                EEG = pop_loadset( ssList(s).name, pathIn );
            elseif strcmp(Ext,'.mat')
                EEG = load(fullfile(pathIn,  ssList(s).name));
                feeg = fieldnames(EEG);
                EEG = EEG.(feeg{1});
            end
            
            oktoadd = 1;
            for k=1:length(P.DataField)
                oktoadd = oktoadd && ~isempty(EEG.(P.DataField{k}));
            end
            
            if Si && oktoadd
                S = S+1;
                
                % ----------------------------------------------------------------------
                % Add the data
                
                for k=1:length(P.DataField)
                    thefield = P.DataField{k};
                    if S==1
                        if numel(EEG.(thefield))>0
                            EEGALL.(thefield) = EEG.(thefield);
                        end
                    else
                        if numel(EEG.(thefield))>0
                            if isempty(P.dim2cat) || length(P.dim2cat)<k
                                dim2cat = ndims(EEG.(thefield));
                            else
                                dim2cat = P.dim2cat(k);
                            end
                            EEGALL.(thefield) = cat(dim2cat,EEGALL.(thefield),EEG.(thefield));
                        end
                    end
                end
                
                % ----------------------------------------------------------------------
                % Add the factors
                for k=1:length(P.FactorField)
                    thefield = P.FactorField{k};
                    if isempty(P.dim2cat) || length(P.dim2cat)<k
                        dim2cat = ndims(EEG.(P.DataField{1}));
                    else
                        dim2cat = P.dim2cat(k);
                    end
                    
                    if S==1
                        if isempty(P.FactorFieldnames)
                            EEGALL.(thefield) = EEG.(thefield);
                        else
                            aa = cell(1,length(EEG.(thefield)));
                            for jj=1:length(EEG.(thefield))
                                aa{jj} = EEG.(thefield){jj}.name;
                            end
                            for ii=1:length(P.FactorFieldnames)
                                EEGALL.(thefield){ii}.name = P.FactorFieldnames{ii};
                                if any(strcmp(aa,P.FactorFieldnames{ii}))
                                    jj = find(strcmp(aa,P.FactorFieldnames{ii}));
                                    EEGALL.(thefield){ii}.val = EEG.(thefield){jj}.val;
                                    EEGALL.(thefield){ii}.g = EEG.(thefield){jj}.g;
                                else
                                    EEGALL.(thefield){ii}.val = {};
                                    EEGALL.(thefield){ii}.g = nan(1,length(EEG.(thefield){1}.g));
                                end
                            end
                            
                        end
                        isbjf = length(EEGALL.(thefield))+1;
                        EEGALL.(thefield){isbjf}.name = 'SBJ';
                        EEGALL.(thefield){isbjf}.val{1,S} = Sbjnames{Si};
                        EEGALL.(thefield){isbjf}.g = S*ones(1,size(EEG.(P.DataField{1}),dim2cat));
                        
                    else
                        aa = cell(1,length(EEG.(thefield)));
                        for jj=1:length(EEG.(thefield))
                            aa{jj} = EEG.(thefield){jj}.name;
                        end
                        for ii=1:isbjf-1
                            jj = find(strcmp(aa,EEGALL.(thefield){ii}.name));
                            if ~isempty(jj)
                                newval = unique([EEG.(thefield){jj}.val, EEGALL.(thefield){ii}.val], 'sorted');
                                g_eegall = nan(size(EEGALL.(thefield){ii}.g));
                                g_eeg = nan(size(EEG.(thefield){jj}.g));
                                for c=1:length(newval)
                                    idx = strcmp(newval{c}, EEGALL.(thefield){ii}.val(EEGALL.(thefield){ii}.g));
                                    g_eegall(idx) = c;
                                    idx = strcmp(newval{c}, EEG.(thefield){jj}.val(EEG.(thefield){jj}.g));
                                    g_eeg(idx) = c;
                                end
                                EEGALL.(thefield){ii}.val = newval;
                                EEGALL.(thefield){ii}.g = cat(2, g_eegall(:)', g_eeg(:)');
                            else
                                EEGALL.(thefield){ii}.g = cat(2, EEGALL.(thefield){ii}.g, nan(1,length(EEG.(thefield){1}.g)));
                            end
                                
                        end
                        EEGALL.(thefield){isbjf}.val{1,S} = Sbjnames{Si};
                        EEGALL.(thefield){isbjf}.g = cat(2,EEGALL.(thefield){ii+1}.g,S*ones(1,size(EEG.(P.DataField{1}),dim2cat)));
                    end
                end
                
                
                % ----------------------------------------------------------------------
                % Add the summary history
                if isfield(EEG,'summaryhistory') && ~isempty(EEG.summaryhistory)
                    if S==1
                        SummaryHistory = EEG.summaryhistory ;
                    else
                        fff = fieldnames(SummaryHistory);
                        for isum=1:length(fff)
                            SummaryHistory.(fff{isum}) = cat(1,SummaryHistory.(fff{isum}),EEG.summaryhistory.(fff{isum}));
                        end
                    end
                end
                
                % ----------------------------------------------------------------------
                % Add the number of trials
                if isfield(EEG,'TrialsxCND') && ~isempty(EEG.TrialsxCND)
                    if S==1
                        TrialsxCND = EEG.TrialsxCND ;
                    else
                        Vp = TrialsxCND.Properties.VariableNames;
                        Vn = EEG.TrialsxCND.Properties.VariableNames;
                        V = unique([Vp Vn]);
                        m = cellfun(@strcmp,repmat(V,[length(Vp) 1]),repmat(Vp',[1 length(V)]));
                        m = ~any(m,1);  % new conditions
                        if any(m)  % add new condition
                            newcnd = find(m);
                            oldcnd = (1:length(V));
                            oldcnd(newcnd) = [];
                            vals = zeros(size(TrialsxCND{:,:},1),length(V));
                            vals(:,oldcnd) = TrialsxCND{:,:};
                            TrialsxCND = array2table(vals,'VariableNames',V);
                        end
                        vals = [TrialsxCND{:,:}; zeros(1, length(V))];
                        for v=1:length(V) % add the new values
                            in = strcmp(V{v},Vn);
                            if any(in)
                                vals(end,v) = EEG.TrialsxCND{1,in};
                            end
                        end
                        TrialsxCND = array2table(vals,'VariableNames',V) ;
                        
%                         trials = [];
%                         for v=1:length(V) % add the new values
%                             ip = strcmp(V{v},Vp);
%                             in = strcmp(V{v},Vn);
%                             if any(in)
%                                 trials(:,v) = [TrialsxCND{:,ip};EEG.TrialsxCND{:,in}];
%                             else
%                                 trials(:,v) = [TrialsxCND{:,ip}; zeros(1,size(TrialsxCND,2))];
%                             end
%                         end
%                         TrialsxCND = array2table(trials,'VariableNames',V) ;
                    end
                end
                if isfield(EEG,'RejectedData')
                    
                    if S==1
                        RejectedData{1}=EEG.RejectedData;
                    else
                        RejectedData{length(RejectedData)+1} = EEG.RejectedData;
                    end
                end
                
                % ----------------------------------------------------------------------
                % Show the output
                fprintf('--> added\n%s\n',repmat('-',[1 50]))
            else
                % ----------------------------------------------------------------------
                % Show the output
                fprintf('not included\n%s\n',repmat('-',[1 50]))
            end
            
        else
            % ----------------------------------------------------------------------
            % Show the output
            fprintf('--> not included\n%s\n',repmat('-',[1 50]))
        end
    end
    if exist('TrialsxCND','var')
%         Fthefield = [];
%         i=1;
%         while isempty(Fthefield)
%             if isfield(EEGALL,'F')
%                 Fthefield = 'F';
%             elseif isfield(EEGALL,['F' P.DataField{i}])
%                 Fthefield = ['F' P.DataField{i}];
%             else
%                 i=i+1;
%             end
%         end
        for k=1:length(P.FactorField)
            Fthefield = P.FactorField{k};
            EEGALL.InfoSBJ.TrialsxCND{k} = array2table(TrialsxCND{:,:},...
                'VariableNames',TrialsxCND.Properties.VariableNames,'RowNames', EEGALL.(Fthefield){end}.val);
        end
    end
    if exist('RejectedData','var')
        EEGALL.InfoSBJ.RejectedData = RejectedData;
    end
    
    EEGALL.nbchan = size(EEGALL.(P.DataField{1}),1);
    EEGALL.trials = size(EEGALL.(P.DataField{1}),3);
    EEGALL.filename = 'all_avg';
    
    if isfield(EEG,'Filter')
        EEGALL.Filter = EEG.Filter;
    end
    if isfield(EEG,'srate')
        EEGALL.srate = EEG.srate;
    end
    if isfield(EEG,'times')
        EEGALL.times = EEG.times;
    end
    if isfield(EEG,'freqs')
        EEGALL.freqs = EEG.freqs;
    end
    if isfield(EEG,'freq')
        EEGALL.freq = EEG.freq;
    end
    if isfield(EEG,'freqs')
        EEGALL.freqs = EEG.freqs;
    end
    if isfield(EEG,'timestf')
        EEGALL.timestf = EEG.timestf;
    end
    if isfield(EEG,'freqstf')
        EEGALL.freqstf = EEG.freqstf;
    end
    if isfield(EEG,'chanlocs')
        EEGALL.chanlocs = EEG.chanlocs;
    end
    if isfield(EEG,'urchanlocs')
        EEGALL.urchanlocs = EEG.urchanlocs;
    end
    if isfield(EEG,'chaninfo')
        EEGALL.chaninfo = EEG.chaninfo;
    end
    if ~isempty(SummaryHistory)
        EEGALL.summaryhistory = SummaryHistory;
    end
    
else
    fprintf('Not file were found\n')
    return
end

fprintf('\n')
end
