% -------------------------------------------------------------------------
% This function function normalize the data by z-score it.
% The average and standar desviation are calculated during the time and the
% epochs specified in the input
%
% INPUTS
% EEG           EEG structure
%
% OPTIONAL INPUTS
% BadData       how to treat bad data 
%               Possible values: 'replacebynan' | 'replacebymean' | 'none'
%               Default 'none'
% epochs        compute the mean and standard deviation for normalization 
%               across each epoch ('single') or all epochs ('all')
%               Default 'all'
% electrodes    compute the mean and standard deviation for normalization 
%               across each electrode ('single') or all  electrodes ('all')
%               Default 'all'
% latency       time window where the mean and standar desviaiton are
%               computed. Can be [t1 t2] or 'all'. Default 'all'
% mean          mean to use. If it is empty is computed from the data. If
%               it is a value, this value is used. If it is description,
%               then it is obtained from EEG.description.mean
% sd            standar desviaiton to use. If it is empty is computed from
%               the data. If it is a value, this value is used. If it is 
%               description, then it is obtained from EEG.description.sd
%
% OUTPUTS
% EEG           EEG structure
%
% -------------------------------------------------------------------------

function [EEG] = eega_normalization( EEG, varargin )

fprintf('### Z-score data normalization ###\n')

%% ------------------------------------------------------------------------
%% Parameters

P.BadData       = 'none';   % 'replacebynan' / 'replacebymean' / 'none'
P.epochs        = 'all'; % 'all' / 'single' / 'timewindow'
P.electrodes    = 'all'; % 'all' / 'single'
P.latency       = 'all';  % 'all' / interval (if epochs is set to sinlge then the normalization is done based on the interval relative to the epoch time, if epoch is set to timewindow the interval is centerd at each sample)
P.mean          = [];   % 'description' / [] / value
P.sd            = [];   % 'description' / [] / value
P.rescale       = 1;   % rescale the data to have a distribtion with this standar desviation

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_normalization: Non recognized inputs')
end

% check the inputs
if ~any(strcmp(P.epochs,{'single' 'all' 'timewindow'}))
    error('eega_normalization: The option ''epochs'' can have values ''single'' / ''all''')
end
if ~any(strcmp(P.BadData,{'replacebynan' 'replacebymean' 'none'}))
    error('eega_normalization: The option ''BadData'' can have values ''replacebynan'' / ''replacebymean'' / ''none''')
end
if ~strcmp(P.latency,'all') && ~(isnumeric(P.latency) && length(P.latency)==2) && ~isempty(P.latency)
    error('eega_normalization: Invalid input ''latency''')
end
if ~strcmp(P.mean,'description') && ~isnumeric(P.mean) && ~isempty(P.mean)
    error('eega_normalization: Invalid input ''mean''')
end
if ~strcmp(P.sd,'description') && ~isnumeric(P.sd) && ~isempty(P.sd)
    error('eega_normalization: Invalid input ''sd''')
end


nEl = size(EEG.data,1);
nS  = size(EEG.data,2);
nEp = size(EEG.data,3);

%% ------------------------------------------------------------------------
%% Normalization

% Obtaine the mean and stamdar desviation
if isempty(P.mean) || isempty(P.sd)    
    if isempty(P.latency) || strcmp(P.latency,'all')
        P.latency = [EEG.times(1) EEG.times(end)];
    end
    switch P.epochs
        
        case 'all'
            fprintf(' - Normalization applied over %s electrodes \n',P.electrodes)
            fprintf(' - Normalization applied over %s epochs \n',P.epochs)
            fprintf(' - Normalization time [%d, %d] \n',P.latency(1), P.latency(2))
            idxt = (EEG.times >= P.latency(1)) & (EEG.times <= P.latency(2));
            [~,D] = eega_rmvbaddata(EEG, 'BadData', P.BadData);
            D = D(:,idxt,:);
            D = reshape(D,[nEl size(D,2)*nEp]);
            if strcmp(P.electrodes,'all')
                D = reshape(D,[1 nEl*size(D,2) size(D,3)]);
            end
            if isempty(P.mean); theM = nanmean(D,2); 
            else; theM = P.mean; end
            if isempty(P.sd); theSD = nanstd(D,[],2); 
            else; theSD = P.sd; end
            
        case 'single'
            fprintf(' - Normalization applied over %s electrodes \n', P.electrodes)
            fprintf(' - Normalization applied over %s epochs \n', P.epochs)
            fprintf(' - Normalization time [%d, %d] \n',P.latency(1), P.latency(2))
            idxt = (EEG.times >= P.latency(1)) & (EEG.times <= P.latency(2));
            [~,D] = eega_rmvbaddata(EEG, 'BadData', P.BadData);
            D = D(:,idxt,:);
            if strcmp(P.electrodes,'all')
                D = reshape(D,[1 nEl*size(D,2) size(D,3)]);
            end
            if isempty(P.mean); theM = nanmean(D,2);
            else; theM = P.mean; end
            if isempty(P.sd); theSD = nanstd(D,[],2); 
            else; theSD = P.sd; end
            
        case 'timewindow'
            fprintf(' - Normalization applied over %s electrodes \n',P.electrodes)
            fprintf(' - Normalization applied over sliding time windows \n')
            fprintf(' - Normalization time [%d, %d] \n',P.latency(1), P.latency(2))
            if strcmp(P.electrodes,'all')
                if isempty(P.mean); theM = nan(1,nS,nEp);
                else; theM = P.mean; end
                if isempty(P.sd); theSD = nan(1,nS,nEp);
                else; theSD = P.sd; end
            else
                if isempty(P.mean); theM = nan(nEl,nS,nEp);
                else; theM = P.mean; end
                if isempty(P.sd); theSD = nan(nEl,nS,nEp);
                else; theSD = P.sd; end
            end
            [~,D] = eega_rmvbaddata(EEG, 'BadData', P.BadData);
            for ismpl=1:nS
                idxt = ( EEG.times >= (EEG.times(ismpl) + P.latency(1)) ) & ( EEG.times <= (EEG.times(ismpl) + P.latency(2)) );
                d = D(:,idxt,:);
                if strcmp(P.electrodes,'all')
                    d = reshape(d,[1 nEl*size(d,2) size(d,3)]);
                end
                if isempty(P.mean); theM(:,ismpl,:) = nanmean(d,2); end
                if isempty(P.sd); theSD(:,ismpl,:) = nanstd(d,[],2); end
            end
    end
end               
if strcmp(P.mean,'description')
    try
        fprintf(' - Mean obtained from the description \n')
        P.mean = EEG.description.mean;
        sigma = EEG.description.sd;
    catch
        error('eega_normalization: The provided mean was not found')
    end
end
if strcmp(P.sd,'description')
    try
        fprintf(' - Standar desviation obtained from the description \n')
        P.sd = EEG.description.sd;
    catch
        error('eega_normalization: The provided standar desviation was not found')
    end
end

EEG.data = bsxfun(@minus, EEG.data, theM);
EEG.data = bsxfun(@rdivide, EEG.data, theSD);
EEG.data = P.rescale * EEG.data;

fprintf('\n')

end
