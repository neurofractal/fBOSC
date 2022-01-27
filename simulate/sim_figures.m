%% Make Figures to Demonstrate Simulations

load('aperiodic3');
load('aperiodic4');

%% Sim some data
cfg                         = [];
cfg.freq                    = 10; % Simulated a 10Hz oscillation
cfg.amplitude               = 24; % SNR of the simulated oscillation
cfg.cycles                  = 32; % How many cycles?
cfg.time                    = 10; % 20s of simulated data
cfg.trial                   = ntrials; % How many 'trials'?
[aperiodic_out_alpha, ...
    osc_alpha, ...
    data_alpha, ...
    alpha_osc_array]      = sim_fBOSC(cfg,250,aperiodic4);

cfg                         = [];
cfg.freq                    = 4; % Simulated a 10Hz oscillation
cfg.amplitude               = 24; % SNR of the simulated oscillation
cfg.cycles                  = 13; % How many cycles?
cfg.time                    = 10; % 20s of simulated data
cfg.trial                   = ntrials; % How many 'trials'?
[aperiodic_out_theta, ...
    osc_theta, ...
    data_theta,...
    theta_osc_array]        = sim_fBOSC(cfg,250,aperiodic4);

% Plot the theta for sanity
figure; subplot(2,1,1); plot(osc_alpha(1,:));
subplot(2,1,2);plot(osc_theta(1,:));


% Combine all trials into one array
data_comb_osc       = zeros(cfg.trial,(2*cfg.time*Fs));
for i = 1:cfg.trial
    data_comb_osc(i,:) = horzcat(data_alpha.trial{i},data_theta.trial{i});
end

% Plot Figure to Show Example Trial
hFig=figure; 
set(gcf,'Position',[1 1 1200 600]);
plot(0:(1/250):((size(data_comb_osc,2)/250)-1/250),...
    ft_preproc_highpassfilter(data_comb_osc(120,:),250,2),'k');
set(gca,'YTickLabel',[]);
color = get(hFig,'Color');
set(gca,'YColor',color,'TickDir','out')
set(gca,'FontSize',25)
xlabel('Time (s)','FontSize',35);
print('example_sim','-dpng','-r300');

print('example_osc','-dpng','-r300');

%% Sim some data with no oscillations
cfg                         = [];
cfg.freq                    = 10; % Simulated a 10Hz oscillation
cfg.amplitude               = 0; % SNR of the simulated oscillation
cfg.cycles                  = 32; % How many cycles?
cfg.time                    = 10; % 20s of simulated data
cfg.trial                   = 100; % How many 'trials'?
[~, ~, data_no_osc_knee,~]  = sim_fBOSC(cfg,250,aperiodic3);
[~, ~, data_no_osc,~]       = sim_fBOSC(cfg,250,aperiodic4);


cfg                     = [];
% general setup
cfg.fBOSC.F             = [1:1:80];    % frequency sampling
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
cfg.fBOSC.postproc.use      = 'yes';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.fBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.fBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of fBOSC.episodes? (default = 'yes')
cfg.fBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.fBOSC.channel = [1]; % select posterior channels (default: all)
cfg.fBOSC.trial = []; % select trials (default: all)
cfg.fBOSC.trial_background = []; % select trials for background (default: all)

%% run fBOSC & Plot for knee
clear fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, data_no_osc_knee);


figure();
set(gcf,'Position',[100 100 800 600])
hold on

plt_freqs = log10(cfg.fBOSC.F);

% Plot the original data
data = plot(plt_freqs, fBOSC.static.bg_log10_pow(1,:),...
    'black','LineWidth',5);

set(gca,'Fontsize',25)
xlabel({'log10(Frequency)';'(Hz)'},'FontSize',40);
ylabel({'log10(Power)';'(a.u.)'},'FontSize',40);
print('sims_knee_no_osc','-dpng','-r300');

%% run fBOSC and plot for no knee
clear fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, data_no_osc);

figure();
set(gcf,'Position',[100 100 800 600])
hold on

plt_freqs = log10(cfg.fBOSC.F);

% Plot the original data
data = plot(plt_freqs, fBOSC.static.bg_log10_pow(1,:),...
    'black','LineWidth',3);

set(gca,'Fontsize',25)
xlabel({'log10(Frequency)';'(Hz)'},'FontSize',40);
ylabel({'log10(Power)';'(a.u.)'},'FontSize',40);
print('sims_no_osc','-dpng','-r300');


%% Sim some theta/alpha data on

cfg                         = [];
cfg.freq                    = 4; % Simulated a 10Hz oscillation
cfg.amplitude               = 14; % SNR of the simulated oscillation
cfg.cycles                  = 13; % How many cycles?
cfg.time                    = 10; % 20s of simulated data
cfg.trial                   = 100; % How many 'trials'?
[~, ~, data_theta_knee,~]   = sim_fBOSC(cfg,250,aperiodic3);
[~, ~, data_theta,~]        = sim_fBOSC(cfg,250,aperiodic4);
cfg.freq                    = 10; % Simulated a 10Hz oscillation
cfg.cycles                  = 32; % How many cycles?
[~, ~, data_alpha_knee,~]   = sim_fBOSC(cfg,250,aperiodic3);
[~, ~, data_alpha,~]        = sim_fBOSC(cfg,250,aperiodic4);


cfg                     = [];
% general setup
cfg.fBOSC.F             = [1:1:80];    % frequency sampling
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
cfg.fBOSC.postproc.use      = 'yes';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.fBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.fBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of fBOSC.episodes? (default = 'yes')
cfg.fBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.fBOSC.channel = [1]; % select posterior channels (default: all)
cfg.fBOSC.trial = []; % select trials (default: all)
cfg.fBOSC.trial_background = []; % select trials for background (default: all)

%% run fBOSC and plot for no knee
clear fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, data_theta_knee);

figure();
set(gcf,'Position',[100 100 800 600])
hold on

plt_freqs = log10(cfg.fBOSC.F);

% Plot the original data
data = plot(plt_freqs, fBOSC.static.bg_log10_pow(1,:),...
    'black','LineWidth',3);

set(gca,'Fontsize',25)
xlabel({'log10(Frequency)';'(Hz)'},'FontSize',40);
ylabel({'log10(Power)';'(a.u.)'},'FontSize',40);
print('data_theta_knee','-dpng','-r300');

ft_databrowser([],data_theta_knee);
