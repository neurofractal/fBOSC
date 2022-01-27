
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now let's apply these thresholds to single trials using both fBOSC
% and eBOSC approaches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For hit-rate vs false alarm rates we need to use multiple thresholds
cfg.fBOSC.threshold.percentile  = horzcat([0.05:0.05:0.95],[0.91:0.01:0.99]);
[fBOSC, pt_all, dt]             = fBOSC_getThresholds_multi_thresholds(cfg, TFR, []);
[fBOSC, pt_all_eBOSC, dt]       = eBOSC_getThresholds_multi_thresholds(cfg, TFR, []);

%% Get logical array for known times of alpha and theta oscillations

alpha_osc_array;
theta_osc_array = theta_osc_array + size(osc_alpha,2);

osc_combined = [osc_alpha(1,:) osc_theta(1,:)]

% Create logical array for theta
theta_logical = zeros(1,size(osc_combined,2));
theta_logical(theta_osc_array) = 1;
figure; plot(osc_combined); hold on;
plot(theta_logical);
theta_logical = theta_logical(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample)
theta_logical = theta_logical(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
theta_logical = logical(theta_logical);

% Create logical array for alpha
alpha_logical = zeros(1,size(osc_combined,2));
alpha_logical(alpha_osc_array) = 1;
figure; plot(osc_combined); hold on;
plot(alpha_logical);
alpha_logical = alpha_logical(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample)
alpha_logical = alpha_logical(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
alpha_logical = logical(alpha_logical);

%%
hit_rate = zeros(2,length(TFR.trial),...
    length(cfg.fBOSC.threshold.percentile));
false_alarm = zeros(2,length(TFR.trial),...
    length(cfg.fBOSC.threshold.percentile));

% Start of trial loop

for indTrial = 1:size(data_comb_osc,1)
    % get wavelet transform for single trial
    % tfr padding is removed to avoid edge artifacts from the wavelet
    % transform. Note that a padding fpr detection remains attached so that there
    % is no problems with too few sample points at the edges to
    % fulfill the duration criterion.
    TFR_ = TFR.trial{indTrial}(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample);
    disp(indTrial);
    cfg.tmp.trial = indTrial; % encode current trial for later
    cfg.tmp.channel = 1;
    
    for thresh = 1:length(cfg.fBOSC.threshold.percentile)
        
        %% fBOSC
        detected = zeros(size(TFR_));
        
        for f = 1:length(cfg.fBOSC.F)
            detected(f,:) = BOSC_detect(TFR_(f,:),pt_all(thresh,f),dt(f),...
                cfg.fBOSC.fsample);
        end; clear f
        
        % remove padding for detection (matrix with padding required for refinement)
        detected = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
        
        %% fBOSC
        detected_times = single(squeeze(nanmean(detected(...
            cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),1))>0);
        
        % Calculate hit rate
        hits = detected_times(theta_logical);
        hit_rate(1,indTrial,thresh) = length(find(hits==1))/length(theta_osc_array);
        
        % Calculate false alarm rate
        falses = detected_times(~theta_logical);
        false_alarm(1,indTrial,thresh) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
        
        clear hits falses detected_times detected
        
        %% eBOSC
        detected = zeros(size(TFR_));
        
        for f = 1:length(cfg.fBOSC.F)
            detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_eBOSC(thresh,f),dt(f),...
                cfg.fBOSC.fsample);
        end; clear f
        
        % remove padding for detection (matrix with padding required for refinement)
        detected = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
        
        detected_times = single(squeeze(nanmean(detected(...
            cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),1))>0);
        
        % Calculate hit rate
        hits = detected_times(theta_logical);
        hit_rate(2,indTrial,thresh) = length(find(hits==1))/length(theta_osc_array);
        
        % Calculate false alarm rate
        falses = detected_times(~theta_logical);
        false_alarm(2,indTrial,thresh) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
        
        clear hits falses detected_times detected
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here we are not applying any post-processing:
        % Something to revisit? However it probably won't affect
        % the comparison
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

figure; scatter(squeeze(mean(false_alarm(1,:,:),2)),...
    squeeze(mean(hit_rate(1,:,:),2))); hold on

scatter(squeeze(mean(false_alarm(2,:,:),2)),...
    squeeze(mean(hit_rate(2,:,:),2)));

ylabel('Hit Rate');
xlabel('False Alarm Rate');

%     %% Copy fBOSC to EBOSC for consistency
%     cfg.eBOSC = cfg.fBOSC;
%     cfg.eBOSC.channel = 1;
%
%     %% Step 4 (optional): create table of separate rhythmic episodes
%
%       cfg.tmp.inputTime = 0.002:0.002:40;
%       cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample);
%       cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
%
%     [fBOSC.episodes, detected_ep] = eBOSC_episode_create(cfg,TFR_,detected,fBOSC);
%
%     % remove padding for detection (already done for eBOSC.episodes)
%     fBOSC.detected_ep(1, indTrial,:,:) = detected_ep(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
%     clear detected_ep;
clear detected

%% Repeat for eBOSC
% The next section applies both the power and the duration
% threshold to detect individual rhythmic segments in the continuous signals.
detected = zeros(size(TFR_));
for f = 1:length(cfg.fBOSC.F)
    detected(f,:) = BOSC_detect(TFR_(f,:),pt_eBOSC(f),dt_eBOSC(f),cfg.fBOSC.fsample);
end; clear f

% remove padding for detection (matrix with padding required for refinement)
eBOSC.detected(1, indTrial,:,:) = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);

% encode pepisode of detected rhythms (optional)
eBOSC.pepisode(1, indTrial,:) = mean(eBOSC.detected(1, indTrial,:,:),4);

end



false_alarm = [];
hit_rate    = [];

% For theta
for indTrial = 1:size(data_comb_osc,1)
    % fBOSC
    detected_times = single(squeeze(nanmean(fBOSC.detected(1, ...
        indTrial,cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),3))>0);
    
    % Calculate hit rate
    hits = detected_times(theta_logical);
    hit_rate(indTrial,1) = length(find(hits==1))/length(theta_osc_array);
    
    falses = detected_times(~theta_logical);
    
    false_alarm(indTrial,1) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
    clear falses hits detected_times
    
    % eBOSC
    detected_times = single(squeeze(nanmean(eBOSC.detected(1, ...
        indTrial,cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),3))>0);
    
    % Calculate hit rate
    hits = detected_times(theta_logical);
    hit_rate(indTrial,2) = length(find(hits==1))/length(theta_osc_array);
    
    falses = detected_times(~theta_logical);
    
    false_alarm(indTrial,2) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
    clear falses hits detected_times
    
end

figure;boxplot(false_alarm,{'fBOSC','eBOSC'});












%%
h = figure('units','normalized','position',[.1 .1 .6 .3]);
hold on;
plot(origData_time,squeeze(origData), 'k');
tmpDetected = single(squeeze(nanmean(fBOSC.detected_ep(1, ...
    indTrial,cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),3))>0); ...
    tmpDetected(tmpDetected==0) = NaN;
plot(origData_time,squeeze(origData).*tmpDetected', 'r');
xlim([0 40])
xlabel('Time (s)'); ylabel('Power (a.u)');
[~, hobj, ~, ~] = legend({'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'northoutside'); legend('boxoff')
hl = findobj(hobj,'type','line');
ht = findobj(hobj,'type','text')
set(ht,'FontSize',20);
set(hl,'LineWidth',3);
















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
        
        
        
        
        
