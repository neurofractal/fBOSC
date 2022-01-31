function value = calc_prctile(vector, percentile)
% Origianl author: Martin Cousineau, 2020

% Try to use toolbox function
try
    value = prctile(vector, percentile);
    return;
catch
end

if ~isvector(vector)
    error('Only vectors supported.');
end
