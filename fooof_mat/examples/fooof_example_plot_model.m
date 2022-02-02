%% FOOOF Matlab Wrapper Example - Plot a FOOOF model
%
% This example computes an example power spectrum model fit, 
% and then plots the result.
%

%% Run Example

% Load data
load('data/ch_dat_one.mat');

% Calculate a power spectrum with Welch's method
[psd, freqs] = pwelch(ch_dat_one, 500, [], [], s_rate);

% Transpose, to make inputs row vectors
freqs = freqs';
psd = psd';

% FOOOF settings
settings                = struct();  % Use defaults
settings.aperiodic_mode = 'knee';
settings.verbose        = 0;

f_range = [1, 30];

% Run FOOOF, also returning the model
fooof_results = fooof(freqs, psd, f_range, settings, true);

% Plot the resulting model
log_freqs = 1;
fooof_plot(fooof_results,log_freqs)


