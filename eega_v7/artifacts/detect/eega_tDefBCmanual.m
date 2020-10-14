% -------------------------------------------------------------------------
% The functions asks for the bad channels. The input has to be provided as
% a vector with the bad channels 
% For example
%       EEG = eega_DefBCmanual(EEG)
%       Bad channels: [25 86 114]
%
% INPUTS
%   EEG : EEGLab structure
%
% OUTPUTS
%   EEG
% 
% -------------------------------------------------------------------------


function EEG = eega_tDefBCmanual(EEG, BCmanual)

fprintf('### Manually define bad channels ###\n')

if nargin<2
    % Ask for the bad channels
    BCmanual = input('Bad channels: ');
end

EEG.artifacts.BCmanual = BCmanual;

end

