function guess = drop_peak_cf(guess, bw_std_edge, freq_range)
%       Check whether to drop peaks based on center's proximity to the edge of the spectrum.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.

    cf_params = guess(:,1)';
    bw_params = guess(:,3)' * bw_std_edge;

    % Check if peaks within drop threshold from the edge of the frequency range.

    keep_peak = abs(cf_params-freq_range(1)) > bw_params & ...
        abs(cf_params-freq_range(2)) > bw_params;

    % Drop peaks that fail the center frequency edge criterion
    guess = guess(keep_peak,:);

end