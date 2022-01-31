function ap_vals = gen_aperiodic(freqs,aperiodic_params,aperiodic_mode)
%       Generate aperiodic values, from parameter definition.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%       	Frequency vector to create aperiodic component for.
%       aperiodic_params : 1x3 array
%           Parameters that define the aperiodic component.
%       aperiodic_mode : {'fixed', 'knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       ap_vals : 1d array
%           Generated aperiodic values.

    switch aperiodic_mode
        case 'fixed'  % no knee
            ap_vals = expo_nk_function(freqs,aperiodic_params);
        case 'knee'
            ap_vals = expo_function(freqs,aperiodic_params);
    end
end