function aperiodic_params = simple_ap_fit(freqs, power_spectrum, aperiodic_mode)
%       Fit the aperiodic component of the power spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

%       Set guess params for lorentzian aperiodic fit, guess params set at init
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-2, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
    
    exp_guess = abs(power_spectrum(end) - power_spectrum(1) / ...
                            freqs(end)) - log10(freqs(1));
    
    switch (aperiodic_mode)
        case 'fixed'  % no knee
            guess_vec = [power_spectrum(1), exp_guess];
            aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs, power_spectrum);
        case 'knee'
            guess_vec = [power_spectrum(1),1, exp_guess];
            aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs, power_spectrum);
    end

end
