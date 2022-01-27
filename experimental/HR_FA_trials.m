%%
% Investigate how the hit-rate and false alarm rate differ between BOSC,
% eBOSC and fBOSC
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
start_fBOSC
cd(fullfile(cd,'simulate'));

%% Specify conditions for simulations

ntrials                     = 200; % Number of trials
Fs                          = 500; % Sampling Rate
SNR                         = [1.5 4]; % This will make size of peaks the same
                                        % 5.1 14

% Simulate nonlinear 1/f
condition   = 'synaptic';

% Frequency of oscillations
freq = [4 10];

%% Array to hold results of simulations
HR_thresh_alpha = zeros(length(freq),3,ntrials);
FA_thresh_alpha = zeros(length(freq),3,ntrials);

% For each frequency
for fr = 1:length(freq)
    disp(['Frequency: ' num2str(freq(fr)) 'Hz']);
    %% Load the aperiodic data from python
    if strcmp(condition,'synaptic')
        load('synaptic.mat');
        aperiodic = synaptic;
    else
        load('linear.mat');
        aperiodic = linear;
    end
    
    
    % HP-Filter @ 1Hz, like in typical EEG/MEG studies
    for n = 1:ntrials
        aperiodic(n,:) = ft_preproc_highpassfilter(aperiodic(n,:), 500,...
            1, 3);
    end
    
    %% Simulate
    cfg                         = [];
    cfg.freq                    = freq(fr); % Simulated a 10Hz oscillation
    cfg.amplitude               = SNR(fr); % SNR of the simulated oscillation
    cfg.cycles                  = 6/(1/freq(fr)); % 6s bursts
    cfg.time                    = 20; % 20s of simulated data
    cfg.trial                   = ntrials; % How many 'trials'?
    [aperiodic_out_alpha, ...
        osc_alpha, ...
        data_alpha, ...
        alpha_osc_array]      = sim_fBOSC(cfg,Fs,aperiodic);
    
    
    % Combine all trials into one array
    data_comb_osc       = zeros(cfg.trial,(cfg.time*Fs));
    for i = 1:cfg.trial
        data_comb_osc(i,:) = horzcat(data_alpha.trial{i});
    end
    
    
    % Make time variable for combined data
    t = 1/Fs;
    time = [t:t:cfg.time];
    
    
    % Make logical arrays
    log_alpha = zeros(1,size(data_comb_osc,2));
    log_alpha(alpha_osc_array) = 1;
    log_alpha = logical(log_alpha);
    
    
    % Make Figure showing ground truth times of bursts
    ddd = data_comb_osc(1,:);
    figure; plot(time,ddd,'k'); hold on;
    ddd(~(log_alpha)) = NaN;
    plot(time,ddd,'r');
    title('Ground truth times of bursts');
    xlabel('Time (s)');
    
    
    %% Set up BOSC parameters
    cfg                     = [];
    cfg.fBOSC.F             = 2.^[1:.125:5.4];    % frequency sampling
    %cfg.fBOSC.F             = [2:1:40];
    cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
    cfg.fBOSC.fsample       = Fs;         % current sampling frequency of EEG data
    
    % padding
    cfg.fBOSC.pad.tfr_s         = 1;       % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
    cfg.fBOSC.pad.detection_s   = .5;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
    cfg.fBOSC.pad.background_s  = 1;       % padding of segments for BG (only avoiding edge artifacts)
    
    % fooof
    %         if strcmp(condition{c},'synaptic')
    %             cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
    %         else
    cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
    %         end
    
    % threshold settings
    cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
    
    % Use all trials to calculate 1/f
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
    
    % Plot for sanity
    figure; subplot(3,1,1); plot(time,data_comb_osc(1,:)); hold on;
    subplot(3,1,2);
    imagesc(time,cfg.fBOSC.F,log(TFR.trial{1}));
    set(gca,'YDir','normal')
    subplot(3,1,3);
    [val,find_freq]=min(abs(cfg.fBOSC.F-freq(fr)));
    plot(time,zscore((TFR.trial{1}(find_freq,:)))); hold on;
    set(gca,'FontSize',20);
    ylabel('Power (au)');
    xlabel('Time (s)');
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Using multiple thresholds for eBOSC and fBOSC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set threshold at .99
    cfg.fBOSC.threshold.percentile  = 0.99;
    
    % fBOSC
    [fBOSC, pt_all, dt, mean_pow, mp]  = fBOSC_getThresholds_multi_thresholds(cfg, TFR, []);
    
    figure;
    data = plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
    plot(log10(cfg.fBOSC.F), log10(mp), '--r','LineWidth',3);
    plot(log10(cfg.fBOSC.F), log10(pt_all), 'g','LineWidth',0.3);
    title('fBOSC');
    
    % Get thresholds using eBOSC
    cfg.eBOSC = cfg.fBOSC;
    % Remove freqs with oscillatory peaks
    cfg.eBOSC.F_fit = cfg.fBOSC.F(~(cfg.fBOSC.F >= freq(fr)-2 & cfg.fBOSC.F <= freq(fr)+2));
    [eBOSC, pt_all_eBOSC,pt_all_BOSC, dt]       = eBOSC_getThresholds_multi_thresholds(cfg, TFR, []);
    
    figure;
    data = plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
    plot(log10(cfg.fBOSC.F), log10(pt_all_BOSC), 'g','LineWidth',0.3);
    title('eBOSC');
    
    %% Start of trial loop
    for indTrial = 1:ntrials
        
        TFR_ = TFR.trial{indTrial};
        cfg.tmp.trial = indTrial; % encode current trial for later
        cfg.tmp.channel = 1;
        
        %% fBOSC
        detected = zeros(size(TFR_));
        
        for f = 1:length(cfg.fBOSC.F)
            detected(f,:) = BOSC_detect(TFR_(f,:),pt_all(1,f),dt(f),...
                cfg.fBOSC.fsample);
        end; clear f
        
        % ALPHA
        detected_times = single(squeeze(nanmean(detected(...
            cfg.fBOSC.F > freq(fr)-0.5 & cfg.fBOSC.F < freq(fr)+0.5,:),1))>0);
        
        %             % Plot for sanity
        %             ddd = data_comb_osc(1,:);
        %             figure; plot(time,ddd,'k'); hold on;
        %             ddd(~(logical(detected_times))) = NaN;
        %             plot(time,ddd,'r');
        %             xlabel('Time (s)');
        %             title('RED = Alpha Burst Detected');
        
        % Calculate hit rate
        hits = detected_times(log_alpha);
        HR_thresh_alpha(fr,3,indTrial) = ...
            length(find(hits==1))/length(alpha_osc_array);
        
        % Calculate false alarm rate
        falses = detected_times(~log_alpha);
        FA_thresh_alpha(fr,3,indTrial) = ...
            length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
        
        clear hits falses detected_times detected
        
        
        %% eBOSC
        detected = zeros(size(TFR_));
        
        for f = 1:length(cfg.fBOSC.F)
            detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_eBOSC(1,f),dt(f),...
                cfg.fBOSC.fsample);
        end; clear f
        
        % ALPHA
        detected_times = single(squeeze(nanmean(detected(...
            cfg.fBOSC.F > freq(fr)-0.5 & cfg.fBOSC.F < freq(fr)+0.5,:),1))>0);
        
        % Calculate hit rate
        hits = detected_times(log_alpha);
        HR_thresh_alpha(fr,2,indTrial) = ...
            length(find(hits==1))/length(alpha_osc_array);
        
        % Calculate false alarm rate
        falses = detected_times(~log_alpha);
        FA_thresh_alpha(fr,2,indTrial) = ...
            length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
        
        clear hits falses detected_times detected
        
        %% BOSC
        detected = zeros(size(TFR_));
        
        for f = 1:length(cfg.fBOSC.F)
            detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_BOSC(1,f),dt(f),...
                cfg.fBOSC.fsample);
        end; clear f
        
        
        % ALPHA
        detected_times = single(squeeze(nanmean(detected(...
            cfg.fBOSC.F > freq(fr)-0.5 & cfg.fBOSC.F < freq(fr)+0.5,:),1))>0);
        
        % Calculate hit rate
        hits = detected_times(log_alpha);
        HR_thresh_alpha(fr,1,indTrial) = ...
            length(find(hits==1))/length(alpha_osc_array);
        
        % Calculate false alarm rate
        falses = detected_times(~log_alpha);
        FA_thresh_alpha(fr,1,indTrial) = ...
            length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
        
        clear hits falses detected_times detected
    end
end

    
figure; 
for i = 1:3
    subplot(1,3,i);
    boxplot(squeeze(HR_thresh_alpha(2,i,:)) - ...
        squeeze(HR_thresh_alpha(1,i,:))); hold on;
end

figure; 
for i = 1:3
    subplot(1,3,i);
    boxplot([squeeze(HR_thresh_alpha(2,i,:))...
        squeeze(HR_thresh_alpha(1,i,:))]); hold on;
end


%% Export for plotting in python

freq_label = repmat(repmat({'theta';'alpha'},3,1),200,1);
condition = repmat({'BOSC';'BOSC';'eBOSC';'eBOSC';'fBOSC';'fBOSC'},200,1);

HR   = reshape(HR_thresh_alpha,[1200 1]);
FA   = reshape(FA_thresh_alpha,[1200 1]);
                          

t = table(HR,FA,condition,freq_label);
writetable(t,'HR_FA.csv')   

%% Export difference
diff_HR = squeeze(HR_thresh_alpha(2,:,:)) - ...
        squeeze(HR_thresh_alpha(1,:,:));

    diff_FA = squeeze(FA_thresh_alpha(2,:,:)) - ...
        squeeze(FA_thresh_alpha(1,:,:));
    
diff_HR   = reshape(diff_HR,[600 1]);
diff_FA   = reshape(diff_FA,[600 1]);

condition = repmat({'BOSC';'eBOSC';'fBOSC'},200,1);

t = table(diff_HR,diff_FA,condition);
writetable(t,'HR_FA_diff.csv')   




