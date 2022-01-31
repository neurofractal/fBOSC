function spectrum_flat = flatten_spectrum(freqs, power_spectrum, robust_aperiodic_params, aperiodic_mode)
%       Flatten the power spectrum by removing the aperiodic component.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       robust_aperiodic_params : 1x2 or 1x3 array (see aperiodic_mode)
%           Parameter estimates for aperiodic fit.
%       aperiodic_mode : 1 or 2
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       spectrum_flat : 1xn array
%           Flattened (aperiodic removed) power spectrum.


spectrum_flat = power_spectrum - gen_aperiodic(freqs,robust_aperiodic_params,aperiodic_mode);

end

function [model_params,peak_function] = fit_peaks(freqs, flat_iter, max_n_peaks, peak_threshold, min_peak_height, gauss_std_limits, proxThresh, peakType, guess_weight, hOT)
%       Iteratively fit peaks to flattened spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       max_n_peaks : double
%           Maximum number of gaussians to fit within the spectrum.
%       peak_threshold : double
%           Threshold (in standard deviations of noise floor) to detect a peak.
%       min_peak_height : double
%           Minimum height of a peak (in log10).
%       gauss_std_limits : 1x2 double
%           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
%       proxThresh : double
%           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%       peakType : {'gaussian', 'cauchy'}
%           Which types of peaks are being fitted
%       guess_weight : {'none', 'weak', 'strong'}
%           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       gaussian_params : mx3 array, where m = No. of peaks.
%           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].
    switch peakType 
        case 'gaussian' % gaussian only
            peak_function = @gaussian; % Identify peaks as gaussian
            % Initialize matrix of guess parameters for gaussian fitting.
            guess_params = zeros(max_n_peaks, 3);
            % Save intact flat_spectrum
            flat_spec = flat_iter;
            % Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian.
            % Stopping procedure based on either the limit on # of peaks,
            % or the relative or absolute height thresholds.
            for guess = 1:max_n_peaks
                % Find candidate peak - the maximum point of the flattened spectrum.
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);

                % Stop searching for peaks once max_height drops below height threshold.
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end

                % Set the guess parameters for gaussian fitting - mean and height.
                guess_freq = freqs(max_ind);
                guess_height = max_height;

                % Halt fitting process if candidate peak drops below minimum height.
                if guess_height <= min_peak_height
                    break
                end

                % Data-driven first guess at standard deviation
                % Find half height index on each side of the center frequency.
                half_height = 0.5 * max_height;

                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height);

                % Keep bandwidth estimation from the shortest side.
                % We grab shortest to avoid estimating very large std from overalapping peaks.
                % Grab the shortest side, ignoring a side if the half max was not found.
                % Note: will fail if both le & ri ind's end up as None (probably shouldn't happen).
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate std from FWHM. Calculate FWHM, converting to Hz, get guess std from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_std = fwhm / (2 * sqrt(2 * log(2)));

                % Check that guess std isn't outside preset std limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_std < gauss_std_limits(1)
                    guess_std = gauss_std_limits(1);
                end
                if guess_std > gauss_std_limits(2)
                    guess_std = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq, guess_height, guess_std];

                % Subtract best-guess gaussian.
                peak_gauss = gaussian(freqs, guess_freq, guess_height, guess_std);
                flat_iter = flat_iter - peak_gauss;

            end
            % Remove unused guesses
            guess_params(guess_params(:,1) == 0,:) = [];

            % Check peaks based on edges, and on overlap
            % Drop any that violate requirements.
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 1, guess_weight, gauss_std_limits,hOT);
            else
                model_params = zeros(1, 3);
            end
            
        case 'cauchy' % cauchy only
            peak_function = @cauchy; % Identify peaks as cauchy
            guess_params = zeros(max_n_peaks, 3);
            flat_spec = flat_iter;
            for guess = 1:max_n_peaks
                max_ind = find(flat_iter == max(flat_iter));
                max_height = flat_iter(max_ind);
                if max_height <= peak_threshold * std(flat_iter)
                    break
                end
                guess_freq = freqs(max_ind);
                guess_height = max_height;
                if guess_height <= min_peak_height
                    break
                end
                half_height = 0.5 * max_height;
                le_ind = sum(flat_iter(1:max_ind) <= half_height);
                ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height);
                short_side = min(abs([le_ind,ri_ind]-max_ind));

                % Estimate gamma from FWHM. Calculate FWHM, converting to Hz, get guess gamma from FWHM
                fwhm = short_side * 2 * (freqs(2)-freqs(1));
                guess_gamma = fwhm/2;
                % Check that guess gamma isn't outside preset limits; restrict if so.
                % Note: without this, curve_fitting fails if given guess > or < bounds.
                if guess_gamma < gauss_std_limits(1)
                    guess_gamma = gauss_std_limits(1);
                end
                if guess_gamma > gauss_std_limits(2)
                    guess_gamma = gauss_std_limits(2);
                end

                % Collect guess parameters.
                guess_params(guess,:) = [guess_freq(1), guess_height, guess_gamma];

                % Subtract best-guess cauchy.
                peak_cauchy = cauchy(freqs, guess_freq(1), guess_height, guess_gamma);
                flat_iter = flat_iter - peak_cauchy;

            end
            guess_params(guess_params(:,1) == 0,:) = [];
            guess_params = drop_peak_cf(guess_params, proxThresh, [min(freqs) max(freqs)]);
            guess_params = drop_peak_overlap(guess_params, proxThresh);

            % If there are peak guesses, fit the peaks, and sort results.
            if ~isempty(guess_params)
                model_params = fit_peak_guess(guess_params, freqs, flat_spec, 2, guess_weight, gauss_std_limits,hOT);
            else
                model_params = zeros(1, 3);
            end
    end
            
end