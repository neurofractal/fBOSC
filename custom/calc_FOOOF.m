function fooof_results = calc_FOOOF(freqs,mean_pow)

use_matlab_fit_function = 1;
cfg.fBOSC.fooof.verbose = 0;


diff_freq = diff(freqs);
if ~all(diff_freq == diff_freq(1))
    warning(['The frequency resolution of cfg.fBOSC.F is ',...
        'not constant. Results could be weird. Proceed with caution']);
end


% specparam opts
opt                     = [];
opt.freq_range          = freqs([1 end]);
opt.peak_width_limits   = [2 6];
opt.max_peaks           = 3;
opt.min_peak_height     = 0.1;
opt.aperiodic_mode      = 'knee'
opt.peak_threshold      = 2.0;   % 2 std dev: parameter for interface simplification
% Matlab-only options
opt.peak_type           = 'gaussian'; % alternative: cauchy
opt.proximity_threshold = 2;
opt.guess_weight        = 'none';
opt.thresh_after        = true;
if license('test','optimization_toolbox') % check for optimization toolbox
    opt.hOT = 1;
    if cfg.fBOSC.fooof.verbose
        disp('FOOOF: Using constrained optimization, Guess Weight ignored.');
    end
else
    opt.hOT = 0;
    if cfg.fBOSC.fooof.verbose
        disp('FOOOF: Using unconstrained optimization, with Guess Weights.');
    end
end
opt.rmoutliers          = 'yes';
opt.maxfreq             = 2.5;
opt.maxtime             = 6;
opt.minnear             = 3;

% log10 the power
spec = log10(mean_pow);

% Fit aperiodic using robust_ap_fit
aperiodic_pars = robust_ap_fit(freqs, spec, opt.aperiodic_mode);

if cfg.fBOSC.fooof.verbose
    figure;
    ap_fit1 = gen_aperiodic(freqs, aperiodic_pars, opt.aperiodic_mode);
    plot(log10(freqs),spec); hold on;
    plot(log10(freqs),ap_fit1); hold on;
    title('Initial 1/f fit');
    set(gca,'FontSize',18);
    xlabel('Log Freqs (Hz)','FontSize',20);
    ylabel('Log Power a.u.','FontSize',20);
end

% Remove aperiodic component
flat_spec = flatten_spectrum(freqs, spec, aperiodic_pars, opt.aperiodic_mode);

% figure;
% plot(log10(freqs), flat_spec, '-ro');

% % Fit peaks
% % Where frequency resolution is linear use default peak fitting
if all(diff_freq == diff_freq(1))
    [peak_pars, peak_function] = fit_peaks(freqs, flat_spec, ...
        opt.max_peaks, opt.peak_threshold, opt.min_peak_height, ...
        opt.peak_width_limits/2, opt.proximity_threshold, ...
        opt.peak_type, opt.guess_weight,opt.hOT,...
        use_matlab_fit_function,cfg.fBOSC.fooof.verbose);
else
    % Where frequency resolution is not linear, increase it via
    % spline interpolation.
    %
    % EXPERIMENTAL - proceed with caution
    %
    if cfg.fBOSC.fooof.verbose
        disp(['Increasing frequency resolution to '...
            num2str(min(diff_freq)) ' Hz']);
    end
    [peak_pars, peak_function] = fit_peaks_interp(freqs,...
        flat_spec, opt.max_peaks, opt.peak_threshold, ...
        opt.min_peak_height, opt.peak_width_limits/2, ...
        opt.proximity_threshold, opt.peak_type, ...
        opt.guess_weight,opt.hOT,...
        use_matlab_fit_function,cfg.fBOSC.fooof.verbose);
end


if opt.thresh_after && ~opt.hOT  % Check thresholding requirements are met for unbounded optimization
    peak_pars(peak_pars(:,2) < opt.min_peak_height,:)     = []; % remove peaks shorter than limit
    peak_pars(peak_pars(:,3) < opt.peak_width_limits(1)/2,:)  = []; % remove peaks narrower than limit
    peak_pars(peak_pars(:,3) > opt.peak_width_limits(2)/2,:)  = []; % remove peaks broader than limit
    peak_pars = drop_peak_cf(peak_pars, opt.proximity_threshold, opt.freq_range); % remove peaks outside frequency limits
    peak_pars(peak_pars(:,1) < 0,:) = []; % remove peaks with a centre frequency less than zero (bypass drop_peak_cf)
    peak_pars = drop_peak_overlap(peak_pars, opt.proximity_threshold); % remove smallest of two peaks fit too closely
end

% Remove peaks and Refit aperiodic
aperiodic = spec;
for peak = 1:size(peak_pars,1)
    aperiodic = aperiodic - peak_function(freqs,peak_pars(peak,1), peak_pars(peak,2), peak_pars(peak,3));
end

if cfg.fBOSC.fooof.verbose
    figure;
    plot(log10(freqs), aperiodic, 'red');
    title('Flattened Spectrum with No Peaks');
    set(gca,'FontSize',18);
    xlabel('Log Freqs (Hz)','FontSize',20);
    ylabel('Log Power a.u.','FontSize',20);
end

aperiodic_pars = simple_ap_fit(freqs, aperiodic, opt.aperiodic_mode);
% Generate model fit
ap_fit = gen_aperiodic(freqs, aperiodic_pars, opt.aperiodic_mode);
model_fit = ap_fit;
for peak = 1:size(peak_pars,1)
    model_fit = model_fit + peak_function(freqs,peak_pars(peak,1),...
        peak_pars(peak,2),peak_pars(peak,3));
end

% Calculate model error
MSE = sum((spec - model_fit).^2)/length(model_fit);
rsq_tmp = corrcoef(spec,model_fit).^2;

if cfg.fBOSC.fooof.verbose
    figure;
    plot(log10(freqs),log10(mean_pow),'k','LineWidth',3); hold on;
    plot(log10(freqs), ap_fit, '--r','LineWidth',2);
    plot(log10(freqs), model_fit, 'green');
    title(['Overall fit. Error: ', num2str( rsq_tmp(2))]);
    legend({'','1/f fit','model fit'});
    set(gca,'FontSize',18);
    xlabel('Log Freqs (Hz)','FontSize',20);
    ylabel('Log Power a.u.','FontSize',20);
end

% Return FOOOF results
fooof_results                   = [];
fooof_results.r_squared         = rsq_tmp(2);
fooof_results.error             = MSE;
fooof_results.peak_params       = peak_pars;
fooof_results.ap_fit            = ap_fit;
fooof_results.aperiodic_params  = aperiodic_pars;
fooof_results.freqs             = freqs;
fooof_results.power_spectrum    = spec;
fooof_results.fooofed_spectrum  = model_fit;

