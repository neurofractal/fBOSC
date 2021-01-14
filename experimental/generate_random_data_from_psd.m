%%
start_fBOSC
cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

%% 
load('aperiodic1.mat');

% alpha amplitudes
amplitude = [0 2 4 6 8 12 16 24];

% alpha cycles
cycles = [2 4 8 16 32 64 128 200]; % total of 14 seconds

%%
cfg                         = [];
cfg.freq                    = 10; % Simulated a 10Hz oscillation
cfg.amplitude               = 24; % SNR of the simulated oscillation
cfg.cycles                  = 64; % How many cycles?
cfg.time                    = 20; % 20s of simulated data
cfg.trial                   = 500; % How many 'trials'?
[aperiodic_out_alpha, ...
osc_alpha, data_alpha]      = sim_fBOSC(cfg,500,aperiodic1);

cfg                         = [];
cfg.freq                    = 30; % Simulated a 10Hz oscillation
cfg.amplitude               = 24; % SNR of the simulated oscillation
cfg.cycles                  = 192; % How many cycles?
cfg.time                    = 20; % 20s of simulated data
cfg.trial                   = 500; % How many 'trials'?
[aperiodic_out_theta, ...
osc_theta, data_theta]      = sim_fBOSC(cfg,500,aperiodic1);


figure; plot(osc_theta(1,:));

% Combine all trials into one array
data_comb_osc       = zeros(cfg.trial,(2*cfg.time*500));
for i = 1:cfg.trial
    data_comb_osc(i,:) = horzcat(data_alpha.trial{i},data_theta.trial{i});
end

% Combine all aperiodic data into one array
data_comb_aperiodic = horzcat(aperiodic_out_alpha(:,:),...
    aperiodic_out_theta(:,:));

% Set up BOSC parameters
cfg                     = [];
cfg.fBOSC.F             = 2.^[1:.125:5.25];    % frequency sampling
cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample       = 500;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 1;       % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = .5;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 1;       % padding of segments for BG (only avoiding edge artifacts)

% fooof
cfg.fBOSC.fooof.aperiodic_mode    = 'knee';

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

%%
cfg.fBOSC.trial_background = 1:size(data_comb_osc,1);
% calculate the sample points for paddding
cfg.fBOSC.pad.tfr_sample = cfg.fBOSC.pad.tfr_s.*cfg.fBOSC.fsample;                          % automatic sample point calculation
cfg.fBOSC.pad.detection_sample = cfg.fBOSC.pad.detection_s.*cfg.fBOSC.fsample;              % automatic sample point calculation
cfg.fBOSC.pad.total_s = cfg.fBOSC.pad.tfr_s + cfg.fBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
cfg.fBOSC.pad.total_sample = cfg.fBOSC.pad.tfr_sample + cfg.fBOSC.pad.detection_sample;
cfg.fBOSC.pad.background_sample = cfg.fBOSC.pad.tfr_sample;
cfg.tmp.channel = 1;

% Compute wavelet transform for the combined data
TFR = [];
for indTrial = 1:size(data_comb_osc,1)
    TFR.trial{indTrial} = BOSC_tf(data_comb_osc(indTrial,:),...
        cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
end; clear indTrial

% Get thresholds + fit line using FOOOF
[fBOSC, pt, dt] = fBOSC_getThresholds(cfg, TFR, []);

% Get thresholds using eBOSC with 3-12Hz excluded
cfg.eBOSC = cfg.fBOSC;
cfg.eBOSC.F_fit = cfg.fBOSC.F(~(cfg.fBOSC.F >= 8 & cfg.fBOSC.F <= 12))
% Remove beta freqs
cfg.eBOSC.F_fit = cfg.eBOSC.F_fit(~(cfg.eBOSC.F_fit >= 25 & cfg.eBOSC.F_fit <= 35));

[eBOSC, pt, dt] = eBOSC_getThresholds_multipeak(cfg, TFR, []);

fBOSC.label = {'chan1'};
cfg.log_freqs = 1;
cfg.plot_old = 1;
fBOSC_fooof_plot(cfg,fBOSC)

% Compute wavelet transform for just the aperiodic data
TFR_aperiodic = [];
for indTrial = 1:1:size(data_comb_osc,1)
    TFR_aperiodic.trial{indTrial} = BOSC_tf(data_comb_aperiodic(indTrial,:),...
        cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
end; clear indTrial

% Calculate RMSE
RMSE_all    = zeros(500,3);

for i = 1:500
    % First BOSC
    trial_aperiodic_mean = mean(log10(TFR_aperiodic.trial{i}(:,:)),2);
    fit_to_use = log10(fBOSC.static.mp_old)';
    RMSE_trial = sqrt(mean((trial_aperiodic_mean-fit_to_use)).^2); % Root mean square error 
    RMSE_all(i,1) = RMSE_trial;
    
    % Now eBOSC with frequencies excluded
    fit_to_use = log10(eBOSC.static.mp)';
    RMSE_trial = sqrt(mean((fit_to_use-trial_aperiodic_mean)).^2); % Root mean square error 
    RMSE_all(i,2) = RMSE_trial;    
    
    % Now FOOOF
    fit_to_use = log10(fBOSC.static.mp)';
    RMSE_trial = sqrt(mean((fit_to_use-trial_aperiodic_mean)).^2); % Root mean square error 
    RMSE_all(i,3) = RMSE_trial;    
end

figure;boxplot(RMSE_all,{'BOSC','eBOSC','fBOSC'});

RMSE                        = reshape(RMSE_all,[1500,1]);
condition                   = [];
[condition{1,1:500}]        = deal('BOSC');
[condition{1,501:1000}]     = deal('eBOSC');
[condition{1,1001:1500}]    = deal('fBOSC');
condition                   = condition';

t = table(RMSE,condition);
writetable(t,'alpha_beta_knee.csv')



















B(k,:,:) = BOSC_tf(signal,cfg.eBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
BGonly_FT(k,:,:) = BOSC_tf(aperiodic1(k,:),cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
freq = 19;
FT_SNR_local(a,c,k) = mean(B(k,freq,AlphaPlace))./mean(BGonly_FT(k,freq,AlphaPlace));
FT_SNR_total(a,c,k) = mean(B(k,19,AlphaPlace))./mean(BGonly_FT(k,19,:));


figure; plot(signal);
    
data = [];
data.trial{1} = signal;
data.label = {'chan1'};
data.time{1} = [0.002:0.002:20];
    
    
%% fBOSC parameters

% general setup
cfg.fBOSC.F             = 2.^[1:.125:5.25];    % frequency sampling
cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample       = 500;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 1;       % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = .5;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 1;       % padding of segments for BG (only avoiding edge artifacts)

% fooof
cfg.fBOSC.fooof.aperiodic_mode    = 'fixed';

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.fBOSC.postproc.use      = 'yes';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.fBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.fBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.fBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.fBOSC.channel = [1]; % select posterior channels (default: all)
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
indChan = 1; indTrial = 1; % Here we select the first trial and first channel we encoded (see cfg.fBOSC.channel).

disp(['Results are for trial ', num2str(cfg.fBOSC.trial(indTrial)), ' at channel ', data.label{cfg.fBOSC.channel(indChan)}])

%% get original time series for plotting
origData = data.trial{indTrial}(cfg.fBOSC.channel(indChan),...
    cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);

origData_time = data.time{indTrial}(...
    cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);

%%
h = figure('units','normalized','position',[.1 .1 .6 .3]);
hold on; 
plot(origData_time,squeeze(origData), 'k');
tmpDetected = single(squeeze(nanmean(fBOSC.detected(indChan, ...
    indTrial,cfg.fBOSC.F > 8 & cfg.fBOSC.F < 12,:),3))>0); ...
    tmpDetected(tmpDetected==0) = NaN;
plot(origData_time,squeeze(origData).*tmpDetected', 'r');
xlim([0 20])
xlabel('Time (s)'); ylabel('Power (a.u)');
[~, hobj, ~, ~] = legend({'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'northoutside'); legend('boxoff')
hl = findobj(hobj,'type','line');
ht = findobj(hobj,'type','text')
set(ht,'FontSize',20);
set(hl,'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',26)

%%
figure; plot(alpha);

episodes = fBOSC.episodes;

%% Calculate hit rate and abundance
if ~isempty(episodes) % if episodes
    alpha_eps = find([episodes{:,3}]>= 8 & [episodes{:,3}]<= 12)';
    if ~isempty(alpha_eps) % if rhythmic episode
        % create unique detected vector
        alpha_locs = cat(1,episodes{alpha_eps,1});
        alpha_locs = unique(alpha_locs(:,2))'; % detected rhythm SPs
        % extract abundance
        abundance_ep(a,c,k,1) = numel(alpha_locs)/numberOfEffectivePoints;
        % check for overlap with simulated rhythm
        if a > 1
            % identify data points in final episode space that were
            % simulated as rhythmic. Note that this gives an
            % accurate amount of total signals as well.
            alpha_sim_locs = zeros(1,amountTimePoints);
            alpha_sim_locs(rhythmIdxVector) = 1;
            alpha_sim_locs = find(alpha_sim_locs(cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad) == 1);
        else
            alpha_sim_locs = [];
        end
        
        % encode amount of rhythm/non-rhythm
        Amount.Alpha(c) = numel(alpha_sim_locs);
        Amount.NoAlpha(c) = numberOfEffectivePoints-numel(alpha_sim_locs);
        
        Hits = intersect(alpha_locs, alpha_sim_locs); % detected & simulated
        




