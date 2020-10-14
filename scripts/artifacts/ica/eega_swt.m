function varargout = eega_swt(x,n)

% Get inputs.
validateattributes(n, {'numeric'}, ...
    {'integer', 'positive', 'scalar'}, 'eega_swt', 'N', 2);

validateattributes(x, {'numeric'}, {'vector', 'real'}, 'eega_swt', 'X', 1);

% Use row vector.
x = x(:)';

s = length(x);
pow = 2^n;
if rem(s,pow)>0
    sOK = ceil(s/pow)*pow;
    error(message('Wavelet:moreMSGRF:SWT_length_MSG',n,s,sOK))
end

% Compute decomposition filters.
F = eega_coifwavf('coif5');
[Lo_D,Hi_D,Lo_R,Hi_R] = eega_orthfilt(F);

% Compute stationary wavelet coefficients.
evenoddVal = 0;
evenLEN    = 1;
swd = zeros(n,s);
swa = zeros(n,s);
for k = 1:n

    % Extension.
    lf = length(Lo_D);
    x  = eega_wextend('1D','per',x,lf/2);

    % Decomposition.
    swd(k,:) = eega_wkeep1(eega_wconv1(x,Hi_D),s,lf+1);
    swa(k,:) = eega_wkeep1(eega_wconv1(x,Lo_D),s,lf+1);

    % upsample filters.
    Lo_D = eega_dyadup(Lo_D,evenoddVal,evenLEN);
    Hi_D = eega_dyadup(Hi_D,evenoddVal,evenLEN);

    % New value of x.
    x = swa(k,:);

end

if nargout==2
    varargout = {swa,swd};
else
    varargout{1} = [swd ; swa(n,:)];
end

