function [EEG, mu, sd] = eega_ZscoreForArt(EEG)

[nEl, nS, nEp] = size(EEG.data);

% take the d
if ~isfield(EEG, 'artifacts') && ~isfield(EEG.artifacts, 'BCT')
    d = EEG.data;
    d(EEG.artifacts.BCT) = NaN;
    d(repmat(EEG.artifacts.BT,[nEl 1 1])) = NaN;
else
    d = EEG.data;
end

% compute the mean and standard desviation
d = reshape(d,[nEl nS*nEp]);
mu = nanmean(d,2);
sd = nanstd(d,[],2);
mu(isnan(mu)) = 0;
sd(isnan(sd)) = 1;

% zscore
EEG.data = (EEG.data - repmat( mu, [1 nS nEp])) ./ repmat( sd, [1 nS nEp]);

end