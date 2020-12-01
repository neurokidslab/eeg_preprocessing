function EEG = eega_demean(EEG)
fprintf('### Zero mean ###\n')
% [~, data_good] = eega_rmvbaddata(EEG, 'BadData', 'replacebynan');
[nEl, nS, nEp ] = size(EEG.data);
mu = nanmean(reshape(EEG.data,[nEl nS*nEp]),2);
mu(isnan(mu)) = 0;
EEG.data = bsxfun(@minus,EEG.data,mu);
fprintf('\n')
end