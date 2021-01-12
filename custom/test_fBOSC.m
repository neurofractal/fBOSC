%%

data_path     = '/Volumes/Robert T5/Robert_HCP/';
save_path     = '/Volumes/Robert T5/RS_HCP_VE';

cd(data_path);

subject = {'102816'};
cd(fullfile(save_path,subject{1}));

load('VE_for_HMM.mat');

%%
start_fBOSC

%% fBOSC parameters

% general setup
cfg.fBOSC.F             = [3:1:50];    % frequency sampling
cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample       = 500;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 1;       % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = .5;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 1;       % padding of segments for BG (only avoiding edge artifacts)

% fooof
cfg.fBOSC.fooof.aperiodic_mode    = 'knee';

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(2, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.fBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.fBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.fBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.fBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.fBOSC.channel = [58:60]; % select posterior channels (default: all)
cfg.fBOSC.trial = []; % select trials (default: all)
cfg.fBOSC.trial_background = []; % select trials for background (default: all)

%% run fBOSC
clear fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, data);

%%
cfg.log_freqs = 1;
cfg.plot_old = 0;
fBOSC_fooof_plot(cfg,fBOSC)

%%
indChan = 1; indTrial = 1; % Here we select the first trial and first channel we encoded (see cfg.eBOSC.channel).

disp(['Results are for trial ', num2str(cfg.fBOSC.trial(indTrial)), ' at channel ', VE_for_HMM.label{cfg.fBOSC.channel(indChan)}])

%% get original time series for plotting
origData = VE_for_HMM.trial{indTrial}(cfg.fBOSC.channel(indChan),...
    cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);

origData_time = VE_for_HMM.time{indTrial}(...
    cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);

%%
h = figure('units','normalized','position',[.1 .1 .6 .3]);
hold on; 
plot(origData_time,squeeze(origData), 'k');
tmpDetected = single(squeeze(nanmean(fBOSC.detected(indChan, ...
    indTrial,cfg.fBOSC.F > 8 & cfg.fBOSC.F < 13,:),3))>0); ...
    tmpDetected(tmpDetected==0) = NaN;
plot(origData_time,squeeze(origData).*tmpDetected', 'r');
xlim([12 15])
xlabel('Time (s)'); ylabel('Power (a.u)');
[~, hobj, ~, ~] = legend({'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'northoutside'); legend('boxoff')
hl = findobj(hobj,'type','line');
ht = findobj(hobj,'type','text')
set(ht,'FontSize',20);
set(hl,'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',26)








