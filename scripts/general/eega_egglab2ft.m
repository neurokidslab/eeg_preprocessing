% This functions transform an structure generated using EEGlab to a
% fieldtrip structure for time lock data

function FT = eega_egglab2ft(EEG,varargin)

%% ------------------------------------------------------------------------
%% Read the optional parameters

% Default parameters

P.DataField='data';
P.TimeField='times';
P.FreqField='freqs';
P.idxtime=[];      % design type within or between
P.idxfreq=[];
P.idxelec=[];
P.idxepoch=[];
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_egglab2ft: Non recognized inputs')
end

%% ------------------------------------------------------------------------
%%
sz = size(EEG.(P.DataField));
nEl = sz(1);
nEp = sz(end);
if (length(sz)==3) && ~isempty(EEG.(P.TimeField))
    nS = sz(2);
    dS = 2;
    nF = 0;
    dF = 0;
elseif (length(sz)==3) && ~isempty(EEG.(P.FreqField))
    nF = sz(2);
    dF = 2;
    nS = 0;
    dS = 0;
elseif (length(sz)==4) && ~isempty(EEG.(P.TimeField)) && ~isempty(EEG.(P.FreqField))
    if (sz(2)==length(EEG.(P.TimeField))) && (sz(3)==length(EEG.(P.FreqField)))
        nS = sz(2);
        dS = 2;
        nF = sz(3);
        dF = 3;
    elseif sz(2)==length(EEG.(P.FreqField)) && (sz(3)==length(EEG.(P.TimeField)))
        nS = sz(3);
        dS = 3;
        nF = sz(2);
        dF = 2;
    else
        error('wrong dimentions')
    end
end

% indexes 
if isempty(P.idxelec)
    P.idxelec = (1:nEl);
end
if isempty(P.idxtime)
    P.idxtime = (1:nS);
end
if isempty(P.idxfreq)
    P.idxfreq = (1:nF);
end
if isempty(P.idxepoch)
    P.idxepoch = (1:nEp);
end

% sampling rate
FT.fsample  = EEG.srate; 

% the timepoints in seconds and the frequency in hertz
if ~isempty(P.idxtime)
    FT.time = EEG.(P.TimeField)(P.idxtime);
end
if ~isempty(P.idxfreq)
    FT.freq = EEG.(P.FreqField)(P.idxfreq);
end

% the numeric data
if nS~=0 && nF==0
    FT.avg = EEG.(P.DataField)(P.idxelec, P.idxtime, P.idxepoch); 
    FT.dimord   = 'chan_time'; % defines how the numeric data should be interpreted
elseif nF~=0 && nS==0
    FT.avg = EEG.(P.DataField)(P.idxelec,P.idxfreq,P.idxepoch); % the numeric data
    FT.dimord   = 'chan_freq'; % defines how the numeric data should be interpreted
elseif nS~=0 && nF~=0
    if nF==3 && nS==2  
        FT.avg = EEG.(P.DataField)(P.idxelec,P.idxtime,P.idxfreq,P.idxepoch); % the numeric data
        FT.avg = permute(FT.avg, [1 3 2 4]); % permute time and freqeuncy
    else
        FT.avg = EEG.(P.DataField)(P.idxelec,P.idxfreq,P.idxtime,P.idxepoch); % the numeric data
    end
    FT.dimord   = 'chan_freq_time'; % defines how the numeric data should be interpreted
end


% Electrodes layout
ChLabels = cell(nEl,1);
Elecpos = nan(nEl,3);
for el=1:nEl
    ChLabels{el} = EEG.chanlocs(el).labels;
    Elecpos(el,1) = EEG.chanlocs(el).X;
    Elecpos(el,2) = EEG.chanlocs(el).Y;
    Elecpos(el,3) = EEG.chanlocs(el).Z;
end
Chanpos = Elecpos;
layout.label = ChLabels;
layout.elecpos = Elecpos;
layout.chanpos = Chanpos;
FT.label    = ChLabels(P.idxelec); % the channel labels
FT.elec     = layout; % information about the sensor array

end