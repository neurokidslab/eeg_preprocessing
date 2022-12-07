function P = ex2_Perp_DSS

% components to keep in the first PCA
P.k          = 50; 
% components to keep in the second PCA
P.n          = 15; 
% Define the trials to use to bias the filter. 
% Different filters can be used for different trials sets
% If empty, one filter is created using all trials
P.fbias      = [];
% Define the trials to which the DSS is applied 
% Different filters can be applied for different trials sets
% If empty, the single filter is applied to all trials
P.fapply     = [];


end