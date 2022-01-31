function peak_params = fit_peak_guess(guess, freqs, flat_spec, peak_type, guess_weight, std_limits, hOT)
%     Fits a group of peak guesses with a fit function.
%
%     Parameters
%     ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       peakType : {'gaussian', 'cauchy', 'best'}
%           Which types of peaks are being fitted.
%       guess_weight : 'none', 'weak', 'strong'
%           Parameter to weigh initial estimates during optimization.
%       std_limits: 1x2 array
%           Minimum and maximum standard deviations for distribution.
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       peak_params : mx3, where m =  No. of peaks.
%           Peak parameters post-optimization.

    
    if hOT % Use OptimToolbox for fmincon
        options = optimset('Display', 'off', 'TolX', 1e-3, 'TolFun', 1e-5, ...
        'MaxFunEvals', 3000, 'MaxIter', 3000); % Tuned options
        lb = [guess(:,1)-guess(:,3)*2,zeros(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(1)];
        ub = [guess(:,1)+guess(:,3)*2,inf(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(2)];
        peak_params = fmincon(@error_model_constr,guess,[],[],[],[], ...
            lb,ub,[],options,freqs,flat_spec, peak_type);
    else % Use basic simplex approach, fminsearch, with guess_weight
        options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-5, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
        peak_params = fminsearch(@error_model,...
            guess, options, freqs, flat_spec, peak_type, guess, guess_weight);
    end
end
