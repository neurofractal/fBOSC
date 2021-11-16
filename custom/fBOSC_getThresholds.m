%    This file is part of the extended Better OSCillation detection (eBOSC) library.
%
%    The eBOSC library is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The eBOSC library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2020 Julian Q. Kosciessa, Thomas H. Grandy, Douglas D. Garrett & Markus Werkle-Bergner

function [fBOSC, pt, dt] = fBOSC_getThresholds(cfg, TFR, fBOSC)
% This function estimates the static duration and power thresholds and
% saves information regarding the overall spectrum and background.
%
% Inputs: 
%           cfg | config structure with cfg.eBOSC field
%           TFR | time-frequency matrix
%           eBOSC | main eBOSC output structure; will be updated
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
    disp('Calculating mean of all TFRs');
    TFR_pad = [];
    for tr = 1:length(TFR.trial)
        TFR_pad{tr} = TFR.trial{tr}(:,cfg.fBOSC.pad.background_sample+1:end-cfg.fBOSC.pad.background_sample);
    end
    
    BG = [TFR_pad{:}];
    
    % Get Freqs and  Power
    freqs = cfg.fBOSC.F;
    % For consistency with the rest of the eBOSC toolbox, the mean is
    % calculated on logged version of the data.. which is then unlogged to
    % pass to fooof
    mean_pow = 10.^(mean(log10(BG(:,:)),2))';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE TO SELF: Is this the right order (log then mean)?
    % BOSC and eBOSC both seem to do this, but then the threshold
    % might not match with the original data?
    %
    % Potentially this is a bug?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Run fooof
    % FOOOF settings
    if strcmp(cfg.fBOSC.fooof.aperiodic_mode,'old')
        settings = struct();
        setting.verbose = 0; % Use defaults
    else
        settings = cfg.fBOSC.fooof;  
    end
    f_range = [freqs(1), freqs(end)];

    % Run FOOOF
    tic;
    fooof_results = fooof(freqs, mean_pow, f_range, settings,true);
    elapsed = toc;
    %fprintf('TIC TOC: %g\n', elapsed);
    
    mp = 10.^fooof_results.ap_fit;
    
    % perform the robust linear fit, only including putatively aperiodic components (i.e., peak exclusion)
    b = robustfit(log10(freqs),mean(log10(BG(:,:)),2)'); clear fitInput;
    pv(1) = b(2); pv(2) = b(1);
    
    mp_old = 10.^(polyval(pv,log10(cfg.fBOSC.F))); 
    
    % Force use of old mp value
    if strcmp(cfg.fBOSC.fooof.aperiodic_mode,'old')
        mp = mp_old;
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
