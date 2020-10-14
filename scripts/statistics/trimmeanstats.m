% This function makes some basic statistics with the trim data
%
% INPUTS
% DATA : matrxi with the data
% PRCT : percent of the data to remove from ecach side
% DIMENTION : dimention along which the statictics are computed
%
% OUTPUTS
% MEAN  : mean
% SD    : standar desviation
% N     : number of observations
% ERROR : standar error

function [avg, sd, n, error] = trimmeanstats(data, prct, dimention)

if nargin < 3
    dimention = 1;
end

% reshape the data in order to be a two dimensional matrix with the
% dimension to compute the mean in the first place
dmorg = (1:ndims(data));
dmrsh = dmorg;
dmrsh(1) = dmorg(dimention);
dmrsh(dimention) = dmorg(1);
data = permute(data,dmrsh);
szrsh = size(data);
data = reshape(data,[szrsh(1) prod(szrsh(2:end))]);

% remove the outliers
n2rmv = floor(size(data,1)*prct/100);
data = sort(data);
data = data(n2rmv+1:end-n2rmv,:);

% reshape the data back
data = reshape(data, [size(data,1) szrsh(2:end)]);
data = permute(data,dmrsh);

% calculate the mean standar, desviation, standar error 
avg = nanmean(data,dimention);
sd = nanstd(data,0,dimention);
n = sum(~isnan(data),dimention);
error = sd./sqrt(n);

end