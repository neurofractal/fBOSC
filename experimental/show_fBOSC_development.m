%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show an example of where BOSC/eBOSC fails and where fBOSC works!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('/Volumes/Robert T5/RS_HCP_VE');

cfg                     = [];
cfg.fBOSC.F             = 2.^[1:.125:5.4];    % frequency sampling
cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample       = 250;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

% fooof
cfg.fBOSC.fooof.aperiodic_mode    = 'knee';

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.fBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.fBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.fBOSC.postproc.edgeOnly = 'no';        % Deconvolution only at on- and offsets of fBOSC.episodes? (default = 'yes')
cfg.fBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.fBOSC.channel = [44]; % select posterior channels (default: all)
cfg.fBOSC.trial = []; % select trials (default: all)
cfg.fBOSC.trial_background = []; % select trials for background (default: all)

%% run fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, VE_for_HMM);

%%
% Plot for sanity
%fBOSC.label = {'Left Hippocampus'};
cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

cfg.log_freqs = 1;
cfg.plot_old = 0;
fBOSC_fooof_plot(cfg,fBOSC)
ylim([-0.5 2]);
title('Fitted with FOOOF');
print('example_FOOOF','-dpng','-r300');
cfg.plot_old = 1;
fBOSC_fooof_plot(cfg,fBOSC)
title('Fitted with BOSC');
ylim([-0.5 2]);
print('example_BOSC','-dpng','-r300');

%% Plot Mean Power Spectrum
figure;
set(gcf,'Position',[1 1 1000 600]);
plot(log10(cfg.fBOSC.F),fBOSC.static.bg_log10_pow,'k','LineWidth',5); hold on;
set(gca,'FontSize',20);
xlabel('log10 Frequency (Hz)','FontSize',30);
ylabel('log10 Power (a.u.)','FontSize',30);
ylim([2.5 4.5]);

print([condition{c} 'PSD_with_' num2str(number_of_osc(n)) 'osc'],'-dpng','-r300');
