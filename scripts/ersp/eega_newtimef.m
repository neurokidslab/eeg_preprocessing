% This function calculates the Event Related Spectral Perturbation using the EEGLAB function pop_timef
%
% INPUTS:
%
% OUTPUTS:
%
% -------------------------------------------------------------------------
% Ana Flo November 2018
% -------------------------------------------------------------------------

function EEGtf = eega_newtimef(EEG, TableCNDs, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.Silent        = 0;

P.DataField     = {'data'};
P.FactorField   = {'F'};
P.BadData       = 'none';   % 'replacebynan' / 'replacebymean' / 'none'
P.EEGlabFomat   = 1;

P.rmverp        = 0;

%P.typeproc      = 1;
P.channels      = (1:size(EEG.data,1));
P.tlimits       = [EEG.times(1) EEG.times(end)];
P.cycles        = [3 0.5];

P.plotersp      = 'off';
P.plotitc       = 'off';
P.plotphase     = 'off';

% Get the optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
% if ~OK
%     for i=1:2:length(extrainput)
%         fprintf('--> Extra otional inputs %s = %s', extrainput{i}, extrainput{i});
%     end
% end

if ~P.Silent; fprintf('### Average ERP ###\n'); end

%% ------------------------------------------------------------------------
%% Build the table defying the conditions 
goodsbj = 1;
if ~isempty(TableCNDs) && ~istable(TableCNDs)
    TableCNDs = cnd_buildtablecond(TableCNDs, EEG.(P.FactorField{1}));
end
if isempty(TableCNDs)
    goodsbj=0;
end


%% ------------------------------------------------------------------------
%% Remove bad data and get the conditions
if goodsbj
    
    F = EEG.(P.FactorField{1});
    [ EEG ] = eega_rmvbaddata(EEG, 'BadData', P.BadData, 'DataField', P.DataField, 'Silent',P.Silent);
    [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(TableCNDs, F);
    idxcndrmv = trialsxCND==0;
    theCND(idxcndrmv) = [];
    trialsxCND(idxcndrmv) = [];
    ntheCND = length(theCND);
    
    for i=1:length(theCND)
        theCND{i}= matlab.lang.makeValidName(theCND{i});
    end
    TrialsxCND = array2table(trialsxCND,'VariableNames',theCND');
    
    if goodsbj && any(trialsxCND==0)
        goodsbj = 0;
    end
    
end

if goodsbj
    
    %% --------------------------------------------------------------------
    %% Average in the time domain and substract the evoke response
    if P.rmverp
        for i=1:length( P.DataField)
            for c=1:ntheCND
                d = EEG.(P.DataField{i})(:,:,condition==c);
                Dmean = nanmean(d,3);
                EEG.(P.DataField{i})(:,:,condition==c) = d - Dmean;
            end   
        end
    end
    
    %% --------------------------------------------------------------------
    %% Define the new factors
    Favg = cnd_buildfactorsfromtable(TableCNDs,theCND);
    
    %% --------------------------------------------------------------------
    %% Time-Frequency analysis
    ERSP = [];
    ITC = [];
    POWBASE = [];
    ERSPBOOT = [];
    ITCBOOT = [];
    times =[];
    freqs =[];
    idx_time = (P.tlimits(1)<=EEG.times) & (P.tlimits(2)>=EEG.times);
    for cnd=1:length(theCND)
        idx_cnd = (condition == cnd);
        ersp_cnd = [];
        itc_cnd = [];
        powbase_cnd = [];
        erspboot_cnd = [];
        itcboot_cnd = [];
        for ch = 1:length(P.channels)
            thech = P.channels(ch);
            dat = EEG.(P.DataField{1})(thech, idx_time, idx_cnd);
            frames = size(dat, 2);
            dat = dat(:)';
            if ~isempty(dat)
                [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(dat,...
                    frames, P.tlimits, EEG.srate, P.cycles,...
                    'plotersp', P.plotersp, 'plotitc', P.plotitc, 'plotphase', P.plotphase,...
                    extrainput{:} );
                
                ersp_cnd = cat(1, ersp_cnd, permute(ersp, [3 2 1]));
                itc_cnd = cat(1, itc_cnd, permute(itc, [3 2 1]));
                powbase_cnd = cat(1, powbase_cnd, permute(powbase, [3 2 1]));
                erspboot_cnd = cat(1, erspboot_cnd, permute(erspboot, [3 2 1]));
                itcboot_cnd = cat(1, itcboot_cnd, permute(itcboot, [3 2 1]));
            else
                goodsbj = 0;
            end
        end
        ERSP = cat(4, ERSP, ersp_cnd);
        ITC = cat(4, ITC, itc_cnd);
        POWBASE = cat(4, POWBASE, powbase_cnd);
        ERSPBOOT = cat(4, ERSPBOOT, erspboot_cnd);
        ITCBOOT = cat(4, ITCBOOT, itcboot_cnd);
    end
    
    %% --------------------------------------------------------------------
    %% Alocate the data
    EEGtf.originaldata.filename = EEG.filename;
    EEGtf.originaldata.filepath = EEG.filepath;
    EEGtf.nbchan = size(ERSP,1);
    EEGtf.pnts = length(times);
    EEGtf.trials = size(ERSP,4);
    EEGtf.ersp = ERSP;
    EEGtf.itc = ITC;
    EEGtf.powbase = POWBASE;
    EEGtf.erspboot = ERSPBOOT;
    EEGtf.itcboot = ITCBOOT;
    EEGtf.times = times;
    EEGtf.freqs = freqs;
    EEGtf.srate = 1/(times(2)-times(1))*1000;
    EEGtf.chanlocs = EEG.chanlocs;
    EEGtf.chaninfo = EEG.chaninfo;
    EEGtf.(P.FactorField{1}) = Favg;
    EEGtf.TrialsxCND = TrialsxCND;
    
else
    warning('No good data!!')
    EEGtf = [];
end


end