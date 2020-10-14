function EEGtf = eega_timefreq(EEG, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.Silent        = 0;
P.DataField     = 'data';
P.FactorField   = 'F';

P.TableCNDs     = [];
P.rmverp        = 0;
P.channels      = (1:size(EEG.data,1));
P.tlimits       = [EEG.times(1) EEG.times(end)];
P.cycles        = [3 0.5];

% Get the optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
% if ~OK
%     for i=1:2:length(extrainput)
%         fprintf('--> Extra otional inputs %s = %s', extrainput{i}, extrainput{i});
%     end
% end

if ~P.Silent; fprintf('### Average ERP ###\n'); end

%% --------------------------------------------------------------------
%% Substract the evoke response
if P.rmverp
    
    % Identify the conditions
    if ~isempty(P.TableCNDs)  % conditions based on factors   
        if ~istable(P.TableCNDs)
            P.TableCNDs = cnd_buildtablecond(P.TableCNDs, EEG.(P.FactorField));
        end
        [condition, theCND, trialsxCND, Favg] = cnd_findtrialsxcond(P.TableCNDs, EEG.(P.FactorField));
        idxcndrmv = trialsxCND==0;
        theCND(idxcndrmv) = [];
        ntheCND = length(theCND);
               
    else  % average across all trials
        ntheCND = 1;
        condition = ones(1,size(EEG.(P.DataField),3));
    end
    
    % remove the evoke response in each condition
    for c=1:ntheCND
        d = EEG.(P.DataField)(:,:,condition==c);
        Dmean = nanmean(d,3);
        EEG.(P.DataField)(:,:,condition==c) = d - Dmean;
    end
     
end

%% --------------------------------------------------------------------
%% Time-Frequency analysis
TF = [];
times =[];
freqs =[];
str = repmat('-',[1 20]);
for ch = 1:length(P.channels)
    thech = P.channels(ch);
    fprintf('%s\nChannel %d\n%s\n',str,thech,str)
    dat = squeeze(EEG.(P.DataField)(thech, :, :));
    if size(dat,1)==1; dat = dat'; end
    [tf,freqs,times] = timefreq(dat, EEG.srate,...
        'cycles',P.cycles,'tlimits',P.tlimits,...
        extrainput{:} );
    tf = permute(tf,[4 1 2 3]);
    TF = cat(1, TF, tf);
end

%% --------------------------------------------------------------------
%% Alocate the data
EEGtf.setname = EEG.setname;
EEGtf.filename = EEG.filename;
EEGtf.filepath = EEG.filepath;
EEGtf.nbchan = size(TF,1);
EEGtf.trials = size(TF,4);
EEGtf.pnts = length(times);
EEGtf.srate = 1/(times(2)-times(1))*1000;
EEGtf.xmax = times(end);
EEGtf.xmin = times(1);
EEGtf.times = times;
EEGtf.freqs = freqs;
EEGtf.data = TF;
EEGtf.chanlocs = EEG.chanlocs;
EEGtf.chaninfo = EEG.chaninfo;
if isfield(EEG,P.FactorField)
    EEGtf.(P.FactorField) = EEG.(P.FactorField);
end
if isfield(EEG,'artifacts')
    EEGtf.artifacts = tf_resampleart(EEG.artifacts, EEG.times, EEGtf.times);
end

end

% resample the artifacts structure
function art_ = tf_resampleart(art, timesorg, timesnew)

    art_ = art;
    Tx = nan(1,length(timesnew));
    for i=1:length(timesnew)
        [v,idx] = min(abs(timesnew(i)-timesorg));
        Tx(i) = idx;
    end
    
    if isfield(art,'BCT')
        art_.BCT = art.BCT(:,Tx,:);
    end
    if isfield(art,'BT')
        art_.BT = art.BT(:,Tx,:);
    end
end

% % resample the artifacts structure
% function art_ = tf_resampleart(art, timesorg, timesnew)
% 
%     art_ = art;
%     
%     timesorg_ = repmat(timesorg',[1 length(timesnew)]);
%     timesnew_ = repmat(timesnew,[length(timesorg)-1 1]);
%     linktimes = timesnew_>=timesorg_(1:end-1,:) & timesnew_<timesorg_(2:end,:);
%     linktimes_ = linktimes .* repmat((1:length(timesorg)-1)',[1 length(timesnew)]);
%     Tx = sum(linktimes_,1);
%     
%     if isfield(art,'BCT')
%         art_.BCT = art.BCT(:,Tx,:);
%     end
%     if isfield(art,'BT')
%         art_.BT = art.BT(:,Tx,:);
%     end
% end
