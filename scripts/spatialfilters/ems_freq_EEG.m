function [EEG, freq, Xdft, W, E] = ems_freq_EEG(EEG, nComps, fbias, fnoise, fapply, varargin)

fprintf('### Effect Match Spatial Filter ###\n')

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.DataField = 'data';
P.FactorField = 'F';
P.emsoneout = 0;
P.fband = [0 100];
P.emseig = 1;

% Optional parameters
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dss_denoise_EEG: Non recognized inputs')
end

if isempty(fnoise)
    fnoise = repmat({{} {}},[size(fbias,1) 1]);
end
if isempty(fapply)
    fapply = fbias;
end

if size(fnoise,1)~=size(fbias,1)
    error('dft_ems_EEG: Different number of filters in the bias and noise')
end
if size(fapply,1)~=size(fbias,1)
    error('dft_ems_EEG: Different number of filters in the bias and fapply')
end


srate = EEG.srate;
d = EEG.(P.DataField);

%% ------------------------------------------------------------------------
%% Build vectors defining how to apply the filter
[Ne, Ns, Nt] = size(EEG.(P.DataField));
Nfilt = size(fbias,1);
tbias = false(Nt,Nfilt);
tnoise = false(Nt,Nfilt);
tapply = false(Nt,Nfilt);
F = EEG.(P.FactorField);
for ifilt=1:Nfilt
    
    % define the trials to bias
    TableCNDs = cnd_buildtablecond(fbias{ifilt,1}, F);
    [condition, theCND, trialsxCND, ~] = cnd_findtrialsxcond(TableCNDs, F);
    [thecnd, idxcnd] = intersect(theCND, fbias{ifilt,2});
    idxbias = any(repmat(idxcnd(:)',[Nt 1])==condition(:), 2);
    tbias(idxbias,ifilt) = 1;
    
    % define the trials for noise
    if ~isempty(fnoise{ifilt,1}) && ~isempty(fnoise{ifilt,2})
        TableCNDs = cnd_buildtablecond(fnoise{ifilt,1}, F);
        [condition, theCND, trialsxCND, ~] = cnd_findtrialsxcond(TableCNDs, F);
        [thecnd, idxcnd] = intersect(theCND, fnoise{ifilt,2});
        idxbias = any(repmat(idxcnd(:)',[Nt 1])==condition(:), 2);
        tnoise(idxbias,ifilt) = 1;
    else
        tnoise(:,ifilt) = 1;
    end
    
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

% remove trial with nans
trlnan = any(any(isnan(d),1),2);
trlnan = trlnan(:);
if any(trlnan)
    warning('%d trials contain NaNs',sum(trlnan))
end
tbias(trlnan,:) = 0;
tnoise(trlnan,:) = 0;

%% ------------------------------------------------------------------------
%% apply the EMS
if P.emsoneout
    if P.emseig
        [freq, Xdft, W, E] = ems_freq_oneout(d, srate, nComps, tbias, tnoise, tapply,'fmin',P.fband(1),'fmax',P.fband(2));
    else
        [freq, Xdft, W] = ems1_freq_oneout(d, srate, tbias, tnoise, tapply,'fmin',P.fband(1),'fmax',P.fband(2));
        E = [];
    end
else
    if P.emseig
        [freq, Xdft, W, E] = ems_freq(d, srate, nComps, tbias, tnoise, tapply, 'fmin',P.fband(1),'fmax',P.fband(2));
    else
        [freq, Xdft, W] = ems1_freq(d, srate, tbias, tnoise, tapply, 'fmin',P.fband(1),'fmax',P.fband(2));
        E = [];
    end
end


fprintf('\n')
end