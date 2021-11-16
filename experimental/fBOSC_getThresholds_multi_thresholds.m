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
%    Copyright 2021 Robert Seymour

function [fBOSC, pt_all, dt,mean_pow,mp] = fBOSC_getThresholds_multi_thresholds(cfg, TFR, fBOSC)
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
        
    pt_all = [];
    
    for tr = 1:length(cfg.fBOSC.threshold.percentile)
        % compute fBOSC power (pt) and duration (dt) thresholds:
        % power threshold is based on a chi-square distribution with df=2 and mean as estimated above
        pt_all(tr,:) = (chi2inv(cfg.fBOSC.threshold.percentile(tr),2)*mp/2); % chi2inv.m is part of the statistics toolbox of Matlab and Octave
    end
    
    % duration threshold is the specified number of cycles, so it scales with frequency
    dt=(cfg.fBOSC.threshold.duration*cfg.fBOSC.fsample./cfg.fBOSC.F)';

    
end
