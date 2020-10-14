% This function makes some basic statistics with the data
%
% INPUTS
% DATA : matrxi with the data
% DIMENTION : dimention along which the statictics are computed
%
% OUTPUTS
% MEAN  : mean
% SD    : standar desviation
% N     : number of observations
% ERROR : standar error

function [mean, sd, n, error] = meanstats(data, dimention)

if nargin < 2
    dimention = 1;
end


mean = nanmean(data,dimention);
sd = nanstd(data,0,dimention);
n = sum(~isnan(data),dimention);
error = sd./sqrt(n);

end