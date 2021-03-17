function Tsummary = eega_summaryartifacts(EEG)

checkfields = {'BCT' 'BC' 'BT' 'BE' 'CCT' 'BS'};
fieldslabels = {'BCT (channels x time)'...
                'BC  (channels)'...
                'BT  (times)'...
                'BE  (epochs)'...
                'CCT (channels x time)',...
                'BS  (bad subject)'};

nEl = size(EEG.data,1);
nS = size(EEG.data,2);
nEp = size(EEG.data,3);

Tot = nan(length(checkfields),1);
TotPos = nan(length(checkfields),1);
PercPos = nan(length(checkfields),1);

% BCT
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
    Tot(1) = nEl*nS*nEp;
    TotPos(1) = sum(EEG.artifacts.BCT(:));
    PercPos(1) = TotPos(1)/Tot(1)*100;
end

% BC
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
    Tot(2) = nEl*nEp;
    TotPos(2) = sum(EEG.artifacts.BC(:));
    PercPos(2) = TotPos(2)/Tot(2)*100;
end

% BT
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
    Tot(3) = nS*nEp;
    TotPos(3) = sum(EEG.artifacts.BT(:));
    PercPos(3) = TotPos(3)/Tot(3)*100;
end


% BE
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BE')
    Tot(4) = nEp;
    TotPos(4) = sum(EEG.artifacts.BE(:));
    PercPos(4) = TotPos(4)/Tot(4)*100;
end


% CCT
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
    Tot(5) = nEl*nS*nEp;
    TotPos(5) = sum(EEG.artifacts.CCT(:));
    PercPos(5) = TotPos(5)/Tot(5)*100;
end

if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BS')
    Tot(6) = 1;
    TotPos(6) = EEG.artifacts.BS;
    PercPos(6) = TotPos(6)/Tot(6)*100;
end

% Summary table
Tsummary = table(TotPos, Tot, PercPos,...
    'VariableNames',{'Positive' 'Total' 'Percentage'},...
    'RowNames',fieldslabels);

end


