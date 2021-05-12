% This function updates the strucutr with the summary information of the
% preprocessing 


function EEG = eega_summarypp(EEG)

[Nch,Nsmpl,Nep] = size(EEG.data);

summary.nsmpl = Nsmpl;
summary.nch = Nch;
summary.nep = Nep;
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
    summary.BCTn = sum(EEG.artifacts.BCT(:));
else
    summary.BCTn = 0;
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
    summary.CCTn = sum(EEG.artifacts.CCT(:));
else
    summary.CCTn = 0;
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
    summary.BTn = sum(EEG.artifacts.BT(:));
else
    summary.BTn = 0;
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
    summary.BCn = sum(EEG.artifacts.BC(:));
else
    summary.BCn = 0;
end
if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BE')
    summary.BEn = sum(EEG.artifacts.BE(:));
else
    summary.BEn = 0;
end

% save the summary
EEG.summary = summary;
% save a summary history
if isfield(EEG,'summaryhistory')
    fff = fieldnames(summary);
    for j=1:length(fff)
        EEG.summaryhistory.(fff{j}) = cat(2,EEG.summaryhistory.(fff{j}), summary.(fff{j}));
    end
else
    EEG.summaryhistory = summary;
end

% 
% 
% % for non epoched data
% if isfield(EEG,'epoch') && isempty(EEG.epoch) %&& ~isfield(EEG.summary,'epoch')
%     EEG.summary.cont.nsmpl = Nsmpl;
%     EEG.summary.cont.nch = Nch;
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
%         EEG.summary.cont.BCTn = sum(EEG.artifacts.BCT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
%         EEG.summary.cont.CCTn = sum(EEG.artifacts.CCT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
%         EEG.summary.cont.BTn = sum(EEG.artifacts.BT);
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
%         EEG.summary.cont.BCn = sum(EEG.artifacts.BC);
%         EEG.summary.cont.BC = find(EEG.artifacts.BC);
%     end
% end
% 
% % for epoched data
% if isfield(EEG,'epoch') && ~isempty(EEG.epoch) 
%     EEG.summary.epoch.nsmpl = Nsmpl*Nep;
%     EEG.summary.epoch.nep = Nep;
%     EEG.summary.epoch.nch = Nch;
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
%         EEG.summary.epoch.BCTn = sum(EEG.artifacts.BCT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
%         EEG.summary.epoch.CCTn = sum(EEG.artifacts.CCT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
%         EEG.summary.epoch.BTn = sum(EEG.artifacts.BT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
%         EEG.summary.epoch.BCn = sum(EEG.artifacts.BC(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BE')
%         EEG.summary.epoch.BEn = sum(EEG.artifacts.BE(:));
%     end
% end
% if isfield(EEG,'epoch') && isempty(EEG.epoch) && isfield(EEG.summary,'epoch')
%     EEG.summary.epoch.nsmpl = Nsmpl*Nep;
%     EEG.summary.epoch.nep = Nep;
%     EEG.summary.epoch.nch = Nch;
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BCT')
%         EEG.summary.epoch.BCTn = sum(EEG.artifacts.BCT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'CCT')
%         EEG.summary.epoch.CCTn = sum(EEG.artifacts.CCT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BT')
%         EEG.summary.epoch.BTn = sum(EEG.artifacts.BT(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BC')
%         EEG.summary.epoch.BCn = sum(EEG.artifacts.BC(:));
%     end
%     if isfield(EEG,'artifacts') && isfield(EEG.artifacts, 'BE')
%         EEG.summary.epoch.BEn = sum(EEG.artifacts.BE(:));
%     end
% end

end