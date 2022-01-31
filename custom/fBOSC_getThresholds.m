%%
% fBOSC_getThresholds - Calculate 1/f fit using FOOOF and then 
%                       get the thresholds for oscillation detection
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________
%    
%    This file is part of the fBOSC library.
%    License: The GNU General Public License v3.0
%
%    Built on work by Kosciessa and colleagues
%    https://github.com/jkosciessa/eBOSC
%__________________________________________________________________________


function [fBOSC, pt, dt] = fBOSC_getThresholds(cfg, TFR, fBOSC)
% Inputs: 
%           cfg | config structure with cfg.fBOSC field
%           TFR | time-frequency matrix
%           fBOSC | main eBOSC output structure; will be updated
%
% Outputs: 
%           fBOSC   | updated w.r.t. background info (see below)
%                   | bg_pow: overall power spectrum
%                   | bg_log10_pow: overall power spectrum (log10)
%                   | pv: intercept and slope of fit
%                   | mp: linear background power
%                   | pt: power threshold
%           pt | empirical power threshold
%           dt | duration threshold
    
    

    % Concatenate TFRs into a 3D array
    % disp('Calculating mean of all TFRs');
    TFR_pad = [];
    for tr = 1:length(TFR.trial)
        TFR_pad{tr} = TFR.trial{tr}(:,cfg.fBOSC.pad.background_sample+1:end-cfg.fBOSC.pad.background_sample);
    end
    
    BG = [TFR_pad{:}];
    
    % Get Freqs and  Power
    freqs = cfg.fBOSC.F;
    % For consistency with the rest of the toolbox, the mean is
    % calculated on logged version of the data.. which is then unlogged to
    % pass to fooof
    mean_pow = 10.^(mean(log10(BG(:,:)),2))';
    
    %figure; plot(log10(freqs),log10(mean_pow));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE TO SELF: Is this the right order (log then mean)?
    % BOSC and eBOSC both seem to do this, but then the threshold
    % might not match with the original data?
    %
    % Potentially this is a bug?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch cfg.fBOSC.fooof.version
        %% Run FOOOF Using MATLAB wrapper for Python
        case 'python'
            % FOOOF settings
            if strcmp(cfg.fBOSC.fooof.aperiodic_mode,'old')
                settings = struct();
                setting.verbose = 0; % Use defaults
            else
                settings = cfg.fBOSC.fooof;
            end
            f_range = freqs([1 end]);
            
            fooof_results = fooof(freqs, mean_pow, f_range, settings,true);
            
        case 'matlab'
            %% Run FOOOF using MATLAB Code
            
            % NOTE: This will only work well for situations in which
            % frequency resolution is constant (i.e. 1,2,3,4...Hz). Future
            % work will be needed to improve this. For now warn the user
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
            opt.aperiodic_mode      = cfg.fBOSC.fooof.aperiodic_mode;
            opt.peak_threshold      = 2.0;   % 2 std dev: parameter for interface simplification
            % Matlab-only options
            opt.peak_type           = 'gaussian'; % alternative: cauchy
            opt.proximity_threshold = 2;
            opt.guess_weight        = 'none';
            opt.thresh_after        = true;
            if license('test','optimization_toolbox') % check for optimization toolbox
                opt.hOT = 1;
                disp('Using constrained optimization, Guess Weight ignored.')
            else
                opt.hOT = 0;
                disp('Using unconstrained optimization, with Guess Weights.')
            end
            opt.rmoutliers          = 'yes';
            opt.maxfreq             = 2.5;
            opt.maxtime             = 6;
            opt.minnear             = 3;
            
            % log10 the power
            spec = log10(mean_pow);
            
            %figure; plot(log10(freqs),spec); hold on;
            
            % Fit aperiodic using robust_ap_fit
            aperiodic_pars = robust_ap_fit(freqs, spec, opt.aperiodic_mode);
            
%             figure;
%             ap_fit1 = gen_aperiodic(freqs, aperiodic_pars, opt.aperiodic_mode);
%             plot(log10(freqs),spec); hold on;
%             plot(log10(freqs),ap_fit1); hold on;
            
            % Remove aperiodic component
            flat_spec = flatten_spectrum(freqs, spec, aperiodic_pars, opt.aperiodic_mode);
            
%             figure;
%             plot(log10(freqs), flat_spec, '-ro');
            
            % Fit peaks
            [peak_pars, peak_function] = fit_peaks(freqs, flat_spec, opt.max_peaks, opt.peak_threshold, opt.min_peak_height, ...
                opt.peak_width_limits/2, opt.proximity_threshold, opt.peak_type, opt.guess_weight,opt.hOT);
            
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
            
%             figure;
%             plot(log10(freqs),model_fit); hold on;
%             plot(log10(freqs),fooof_results2.fooofed_spectrum); hold on;
%             plot(log10(freqs), aperiodic, 'red');
            
            aperiodic_pars = simple_ap_fit(freqs, aperiodic, opt.aperiodic_mode);
            % Generate model fit
            ap_fit = gen_aperiodic(freqs, aperiodic_pars, opt.aperiodic_mode);
            model_fit = ap_fit;
            for peak = 1:size(peak_pars,1)
                model_fit = model_fit + peak_function(freqs,peak_pars(peak,1),...
                    peak_pars(peak,2),peak_pars(peak,3));
            end
%             
%             figure;
%             plot(log10(freqs),log10(mean_pow)); hold on;
%             plot(log10(freqs), ap_fit, 'red');
%             plot(log10(freqs), fooof_results.ap_fit, 'green');
            
            % Calculate model error
            MSE = sum((spec - model_fit).^2)/length(model_fit);
            rsq_tmp = corrcoef(spec,model_fit).^2;
            
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
            
    end
    
    %% Process Results
    mp = 10.^fooof_results.ap_fit;
    
    % Instead of FOOOF, use original BOSC ordinary least squares regression
    % to find the
    [pv_BOSC,~]=BOSC_bgfit(cfg.fBOSC.F,BG);
    
%     b = robustfit(log10(freqs),mean(log10(BG(:,:)),2)'); clear fitInput;
%     pv(1) = b(2); pv(2) = b(1);
    
    mp_old = 10.^(polyval(pv_BOSC,log10(cfg.fBOSC.F))); 
    
    % Force use of original BOSC mp value rather than FOOOF
    if strcmp(cfg.fBOSC.fooof.aperiodic_mode,'old')
        mp = mp_old;
        warning('NOT RECOMMENDED - not using FOOOF');
    end
    
    % compute fBOSC power (pt) and duration (dt) thresholds: 
    % power threshold is based on a chi-square distribution with df=2 and mean as estimated above
    pt=chi2inv(cfg.fBOSC.threshold.percentile,2)*mp/2; % chi2inv.m is part of the statistics toolbox of Matlab and Octave
    pt_old=chi2inv(cfg.fBOSC.threshold.percentile,2)*mp_old/2; % chi2inv.m is part of the statistics toolbox of Matlab and Octave
    % duration threshold is the specified number of cycles, so it scales with frequency
    dt=(cfg.fBOSC.threshold.duration*cfg.fBOSC.fsample./cfg.fBOSC.F)';

    % save multiple time-invariant estimates that could be of interest:
    % overall wavelet power spectrum (NOT only background)
    fBOSC.static.bg_pow(cfg.tmp.channel(1),:)               = ...
        mean(BG(:,cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample),2);
    % log10-transformed wavelet power spectrum (NOT only background)
    fBOSC.static.bg_log10_pow(cfg.tmp.channel(1),:)         = ...
        mean(log10(BG(:,cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample)),2);
    % aperiodic parameters from fooof
    fBOSC.static.aperiodic_params(cfg.tmp.channel(1),:)     = ...
        fooof_results.aperiodic_params;
    % aperiodic fir from fooof
    fBOSC.static.ap_fit(cfg.tmp.channel(1),:)     = ...
        fooof_results.ap_fit;
    % r squared value from fooof
    fBOSC.static.r_squared(cfg.tmp.channel(1),:)     = ...
        fooof_results.r_squared;
    % error term from fooof
    fBOSC.static.error(cfg.tmp.channel(1),:)     = ...
        fooof_results.error;
    % Fooofed spectrum
    fBOSC.static.fooofed_spectrum(cfg.tmp.channel(1),:)     = ...
        fooof_results.fooofed_spectrum;
    % linear background power at each estimated frequency
    fBOSC.static.mp(cfg.tmp.channel(1),:)            = mp;
    fBOSC.static.mp_old(cfg.tmp.channel(1),:)        = mp_old;
    % statistical power threshold
    fBOSC.static.pt(cfg.tmp.channel(1),:)            = pt;
    fBOSC.static.pt_old(cfg.tmp.channel(1),:)        = pt_old;
end
