function [EEG, Y, W] = dss_denoise_EEG(EEG, keepnPCA1, nComps, fbias, fapply, varargin)

fprintf('### Denoising using DSS filter ###\n')

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.DataField = 'data';
P.FactorField = 'F';
P.timewind = [];

% Optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dss_denoise_EEG: Non recognized inputs')
end

if isempty(fapply)
    fapply = fbias;
end
if size(fapply,1)~=size(fbias,1)
    error('dft_ems_EEG: Different number of filters in the bias and fapply')
end

%% ------------------------------------------------------------------------
%% Build vectors defining how to apply the filter
[Ne, Ns, Nt] = size(EEG.(P.DataField));
Nfilt = size(fbias,1);
tbias = false(Nt,Nfilt);
tapply = false(Nt,Nfilt);
F = EEG.(P.FactorField);
for ifilt=1:Nfilt
    
    % define the trials to bias
    TableCNDs = cnd_buildtablecond(fbias{ifilt,1}, F);
    [condition, theCND, trialsxCND, ~] = cnd_findtrialsxcond(TableCNDs, F);
    [thecnd, idxcnd] = intersect(theCND, fbias{ifilt,2});
    idxbias = any(repmat(idxcnd(:)',[Nt 1])==condition(:), 2);
    tbias(idxbias,ifilt) = 1;
    
    % define the trials to apply
    if ~isempty(fapply{ifilt,1}) && ~isempty(fapply{ifilt,2})
        TableCNDs = cnd_buildtablecond(fapply{ifilt,1}, F);
        [condition, theCND, trialsxCND, ~] = cnd_findtrialsxcond(TableCNDs, F);
        [thecnd, idxcnd] = intersect(theCND, fapply{ifilt,2});
        idxbias = any(repmat(idxcnd(:)',[Nt 1])==condition(:), 2);
        tapply(idxbias,ifilt) = 1;
    else
        tapply(:,ifilt) = 1;
    end
     
end

%% ------------------------------------------------------------------------
%% apply the DSS for denoising
d = EEG.(P.DataField);
if isempty(P.timewind)
    twidx = [];
else
    twidx = EEG.times>=P.timewind(1) & EEG.times<=P.timewind(2);
end
[X, Y, W] = dss_denoise(d, keepnPCA1, nComps, tbias, tapply, twidx );
EEG.(P.DataField) = X;
fprintf('\n')


end