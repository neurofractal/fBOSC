cd(fullfile(cd,'simulate'));

%% Make simulated data with oscillations embedded within a 
%  non-linear 1/f aperiodic signal

% Load the previously computed non-linear 1/f aperiodic signal
load('aperiodic3.mat');

% Settings
ntrials                     = 10;  % Number of trials 
Fs                          = 250; % Sampling Rate
SNR                         = 10;   % SNR of oscillation
len_of_trial                = 2;   % Length of each trial in s

% Trim aperiodic data down to size of ntrials*len_of_trial
aperiodic                   = aperiodic3(1:ntrials,[1:(len_of_trial*Fs)]);

cfg                         = [];
cfg.freq                    = 10; % Simulated a 10Hz oscillation
cfg.amplitude               = SNR; % SNR of the simulated oscillation
cfg.cycles                  = 10; % How many cycles?
cfg.time                    = 2; % 20s of simulated data
cfg.trial                   = ntrials; % How many 'trials'?
[aperiodic_out_alpha,osc_alpha, data_alpha,~]...
                            = sim_fBOSC(cfg,Fs,aperiodic);

%% Plot the combined data
ft_databrowser([],data_alpha)

%% Plot the aperiodic and periodic data separately
figure; subplot(3,1,1);plot(data_alpha.time{1},osc_alpha(1,:));
subplot(3,1,2);plot(data_alpha.time{1},aperiodic_out_alpha(1,:));
subplot(3,1,3);plot(data_alpha.time{1},data_alpha.trial{1});


%% fBOSC parameters

% general setup
cfg.fBOSC.F             = [3:1:80];    % frequency sampling
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

%% run fBOSC
clear fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, data_alpha);

%%
cfg.log_freqs = 1;
cfg.plot_old = 0;
fBOSC_fooof_plot(cfg,fBOSC)

%%
indChan = 1; indTrial = 1; % Here we select the first trial and first channel we encoded (see cfg.eBOSC.channel).

disp(['Results are for trial ', num2str(cfg.fBOSC.trial(indTrial)), ' at channel ', VE_for_HMM.label{cfg.fBOSC.channel(indChan)}])

%% get original time series for plotting

figure;

for indTrial = 1:10
    
    origData = data_alpha.trial{indTrial}(cfg.fBOSC.channel(indChan),...
        cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
    
    origData_time = data_alpha.time{indTrial}(...
        cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
    
    %%
    subplot(10,1,indTrial); hold on;
    plot(origData_time,squeeze(origData), 'k');
    tmpDetected = single(squeeze(nanmean(fBOSC.detected_ep(indChan, ...
        indTrial,cfg.fBOSC.F > 8 & cfg.fBOSC.F < 13,:),3))>0); ...
        tmpDetected(tmpDetected==0) = NaN;
    plot(origData_time,squeeze(origData).*tmpDetected', 'r');
    %xlim([12 15])
    xlabel('Time (s)'); ylabel('Power (a.u)');
%     [~, hobj, ~, ~] = legend({'Original signal'; 'Rhythmic signal'}, ...
%         'orientation', 'horizontal', 'location', 'northoutside'); legend('boxoff')
%     hl = findobj(hobj,'type','line');
%     ht = findobj(hobj,'type','text')
%     set(ht,'FontSize',20);
%    set(hl,'LineWidth',3);
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    
end


%%
figure;
for indTrial = 1:10
    TFR_i = zscore(TFR.trial{indTrial}(:,...
        cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample));
    
    tmpDetected = logical(squeeze(fBOSC.detected_ep(1,indTrial,:,:)));
    
    subplot(10,1,indTrial);
    imagesc(TFR_i,'AlphaData', tmpDetected); hold on;
    imagesc(TFR_i, 'AlphaData', 0.2);
    set(gca,'YDir','normal')
    
    set(gca,'YLim',[4 16]);% something like this
    %set(gca,'XTicks',[4 12]);% something like this
   
    if indTrial==indTrial(end)
        xlabel('Time (s)');
    end
    
end





