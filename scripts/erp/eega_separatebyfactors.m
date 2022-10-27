% -------------------------------------------------------------------------
% This function separates trials by conditions, for the conditions
% identified in the eventField of the Epochs, when it takes the values
% eventValues. 
%
% INPUTS
% EEG           EEG structure
% TableCNDs     it specifies what to average. 
%               
%               * It can be a table specifing how to average trials
%       - In each column it has the factors that will be usded to define
%       conditions (the name has to be same than in F{k}.name, and the
%       values have to be some of F{k}.val). The first coloumn is the
%       condition
%       - In each row each the possible combinations for each condition are
%       defined
%
%       For example
%
%            -----------------------------------------------
%           | Condition     |   Img1   |   Aud1   |   Aud2   |
%           |------------------------------------------------|
%           | Acongruent    |    a     |    x1    |    x1    |
%           | Acongruent    |    a     |    x2    |    x2    |
%           | Bcongruent    |    b     |    x1    |    x1    |
%           | Bcongruent    |    b     |    x2    |    x2    |
%           | Aincongruent  |    a     |    x1    |    x2    |
%           | Aincongruent  |    a     |    x2    |    x1    |
%           | Bincongruent  |    b     |    x1    |    x2    |
%           | Bincongruent  |    b     |    x2    |    x1    |
%            -----------------------------------------------
%
%               * It can also be a cell array indicating the factors to
%               take into acount to define conditions. In this case the
%               table will be built
%               * If it is empty all trials are average. 
%
% OPTIONAL INPUTS
% DataField     fields in EEG containig the data. It can be a cell with 
%               multiple fields. In that case the baseline correction is 
%               applied to all fileds. Default {'data'}
% FactorField   field in EEG containg the factors characterizing epochs. 
%               Default 'F'
% dim2avg       dimention across which data should be average. Default the
%               last dimension 
%
% OUTPUTS
% EEG           EEG structure with average data
%
% -------------------------------------------------------------------------

function EEGcnd = eega_separatebyfactors(EEG, TableCNDs, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.FactorField   = {'F'};
P.EEGlabFomat   = 1;
P.Silent        = 0;
P.dimtrl        = [];

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_separatebyfactors: Non recognized inputs')
end

if ~P.Silent; fprintf('### Separating epochs by conditions ###\n'); end

if ~iscell(P.FactorField)
    P.FactorField={P.FactorField};
end
goodsbj=1;
if isempty(EEG.data)
    goodsbj=0;
end

%% ------------------------------------------------------------------------
%% Build the table defying the conditions
if ~isempty(TableCNDs) && ~istable(TableCNDs)
    TableCNDs = cnd_buildtablecond(TableCNDs, EEG.(P.FactorField{1}));
end
if isempty(TableCNDs)
    goodsbj=0;
else
    F = EEG.(P.FactorField{1});
    [condition, theCND, trialsxCND, ~] = cnd_findtrialsxcond(TableCNDs, F);
    ntheCND = length(theCND);
end   

%% ------------------------------------------------------------------------
%% Separate
EEGcnd = [];

for c=1:ntheCND
    cndname = matlab.lang.makeValidName(theCND{c});
    cndname = matlab.lang.makeUniqueStrings(cndname);
    
    if goodsbj

        idxtrl = find(condition==c);
        EEGcnd.(cndname) = pop_select(EEG, 'trial', idxtrl);
        
        % Factor
        f = EEG.F;
        for j=1:length(EEG.F)
            f{j}.g = EEG.F{j}.g(idxtrl);
        end
        EEGcnd.(cndname).F = f;

    else
        for i=1:length( P.DataField)
            EEGcnd.(theCND) = [];
        end
    end

end

fprintf('\n')


end






