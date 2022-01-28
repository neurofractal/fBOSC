%%
% Calculate error in 1/f fitting for BOSC, eBOSC and fBOSC under various
% different simulated conditions
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
root = fileparts(matlab.desktop.editor.getActiveFilename);
cd(fullfile(root,'..'));

start_fBOSC
cd(fullfile(root,'..','validation'));

%% Specify the different simulation options
condition     = {'synaptic','linear'};
oscillations  = {'theta','alpha',};
freq          = [4 10];
Fs            = 500; % Sampling Rate of Simulated Data

SNR = [5.1 14]; % SNRs chosen so amplitude of oscillations is equivalent in
                % log-log space is the same
ntrials = 200;

% Create empty array to hold root mean square error values
RMS_error = zeros(ntrials,3,length(condition),length(oscillations));

% For each condition
for c = 1:length(condition)
    for o = 1:length(oscillations)
        
        disp(['Simulating ' oscillations{o} ' + ' condition{c} ' 1/f']);
        
        % Load the relevent aperiodic data
        if strcmp(condition{c},'synaptic')
            load('synaptic.mat');
            aperiodic = synaptic;
        else
            disp('Using linear 1/f');
            load('linear.mat');
            aperiodic = linear;
        end
        
        % HP-Filter @ 1Hz, like in typical EEG/MEG studies
        for n = 1:ntrials
            aperiodic(n,:) = ft_preproc_highpassfilter(aperiodic(n,:), 500,...
                1, 3);
        end
        
        %% Simulate Data
        cfg                         = [];
        cfg.freq                    = freq(o);
        cfg.amplitude               = SNR(o); % SNR of the simulated oscillation
        cfg.cycles                  = 6/(1/freq(o)); % 6s bursts;
        cfg.time                    = 20; % 20s of simulated data
        cfg.trial                   = ntrials; % How many 'trials'?
        [aperiodic_out, ...
            osc, ...
            data_out, ...
            osc_array]              = sim_fBOSC(cfg,Fs,aperiodic);
        
        % Combine all trials into one array
        data_comb_osc       = zeros(cfg.trial,(cfg.time*Fs));
        for i = 1:cfg.trial
            data_comb_osc(i,:) = horzcat(data_out.trial{i});
        end
        
        % Make time variable for combined data
        t = 1/Fs;
        time = [t:t:cfg.time];
        
        % Make logical arrays
        log_alpha = zeros(1,size(data_comb_osc,2));
        log_alpha(osc_array) = 1;
        log_alpha = logical(log_alpha);
        
        % Make Figure showing ground truth times of bursts
        ddd = data_comb_osc(1,:);
        figure; plot(time,ddd,'k'); hold on;
        ddd(~(log_alpha)) = NaN;
        plot(time,ddd,'r');
        title('Ground truth times of bursts');
        xlabel('Time (s)');
        
        % Set up BOSC, eBOSC and fBOSC parameters
        cfg                     = [];
        cfg.fBOSC.F             = 2.^[1:.125:5.4];    % frequency sampling
        %cfg.fBOSC.F             = [2:0.5:50];
        cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
        cfg.fBOSC.fsample       = Fs;         % current sampling frequency of EEG data
        
        % padding
        cfg.fBOSC.pad.tfr_s         = 1;       % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
        cfg.fBOSC.pad.detection_s   = .5;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
        cfg.fBOSC.pad.background_s  = 1;       % padding of segments for BG (only avoiding edge artifacts)
        
        % fooof
        if strcmp(condition{c},'synaptic')
            cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
        else
            cfg.fBOSC.fooof.aperiodic_mode    = 'fixed';
        end
        
        % threshold settings
        cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
        cfg.fBOSC.threshold.percentile  = .99;                                      % percentile of background fit for power threshold
        
        cfg.fBOSC.trial_background = 1:ntrials;
        % calculate the sample points for paddding
        cfg.fBOSC.pad.tfr_sample = cfg.fBOSC.pad.tfr_s.*cfg.fBOSC.fsample;                          % automatic sample point calculation
        cfg.fBOSC.pad.detection_sample = cfg.fBOSC.pad.detection_s.*cfg.fBOSC.fsample;              % automatic sample point calculation
        cfg.fBOSC.pad.total_s = cfg.fBOSC.pad.tfr_s + cfg.fBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
        cfg.fBOSC.pad.total_sample = cfg.fBOSC.pad.tfr_sample + cfg.fBOSC.pad.detection_sample;
        cfg.fBOSC.pad.background_sample = cfg.fBOSC.pad.tfr_sample;
        cfg.tmp.channel = 1;
        
        % For every trial
        % Compute wavelet transform for the combined data
        TFR = [];
        for indTrial = 1:size(data_comb_osc,1)
            TFR.trial{indTrial} = BOSC_tf(data_comb_osc(indTrial,:),...
                cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
        end; clear indTrial
        
        % Get mean power
        BG = [TFR.trial{:}];
        mean_pow = 10.^(mean(log10(BG(:,:)),2))';
        clear BG
        
        %% Get thresholds + fit line
        % 1. Using fBOSC (FOOOF)
        [fBOSC, pt_fBOSC, dt_fBOSC] = fBOSC_getThresholds(cfg, TFR, []);
        
        % 2. Using eBOSC and BOSC
        cfg.eBOSC = cfg.fBOSC;
        % Remove freqs with oscillatory peaks
        if freq(o) == 4
            cfg.eBOSC.F_fit = cfg.fBOSC.F(~(cfg.fBOSC.F >= 2 ...
                & cfg.fBOSC.F <= 6));
        elseif  freq(o) == 10
            cfg.eBOSC.F_fit = cfg.eBOSC.F(~(cfg.eBOSC.F >= 8 ...
                & cfg.eBOSC.F <= 12));
        end
        [eBOSC, pt_eBOSC, dt_eBOSC] = eBOSC_getThresholds_multipeak(cfg, TFR, []);
        
        % For every trial
        % Compute wavelet transform for just the aperiodic data
        % without the simulated oscillations
        disp('Performing TFR on aperiodic data without oscillations');
        TFR_aperiodic = [];
        for indTrial = 1:size(data_comb_osc,1)
            
            TFR_aperiodic.trial{indTrial} = BOSC_tf(aperiodic_out(indTrial,:),...
                cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
        end
        
        % Concatenate TFRs into a 3D array
        TFR_pad = [];
        for tr = 1:length(TFR_aperiodic.trial)
            TFR_pad{tr} = TFR_aperiodic.trial{tr}(:,cfg.fBOSC.pad.background_sample+1:...
                end-cfg.fBOSC.pad.background_sample);
        end
        
        % Produce Figures showinng 1/f fit and threshold 
        % for BOSC, eBOSC and fBOSC
        figure;
        plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
        plot(log10(cfg.fBOSC.F), log10(fBOSC.static.mp), '--r','LineWidth',3);
        plot(log10(cfg.fBOSC.F), log10(fBOSC.static.pt), 'g','LineWidth',0.3);
        title('fBOSC');
        
        figure;
        data = plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
        plot(log10(cfg.fBOSC.F), log10(eBOSC.static.mp), '--r','LineWidth',3);
        plot(log10(cfg.fBOSC.F), log10(eBOSC.static.pt), 'g','LineWidth',0.3);        
        title('eBOSC');
        drawnow;
        
        %% Calculate RMSE for every trial
        figure;
        plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
        
        for indTrial = 1:ntrials
            % First BOSC
            trial_aperiodic_mean = mean(log10(TFR_aperiodic.trial{indTrial}(:,:)),2);
            p1 = plot(log10(cfg.fBOSC.F),trial_aperiodic_mean,'LineWidth',0.5);
            p1.Color(4) = 0.4;
            fit_to_use = log10(fBOSC.static.mp_old)';
            RMSE_trial = sqrt(mean((trial_aperiodic_mean-fit_to_use)).^2); % Root mean square error
            RMS_error(indTrial,1,c,o) = RMSE_trial;
            
            % Now eBOSC with frequencies excluded
            fit_to_use = log10(eBOSC.static.mp)';
            RMSE_trial = sqrt(mean((fit_to_use-trial_aperiodic_mean)).^2); % Root mean square error
            RMS_error(indTrial,2,c,o) = RMSE_trial;
            
            % Now fBOSC (using FOOOF)
            fit_to_use = log10(fBOSC.static.mp)';
            RMSE_trial = sqrt(mean((fit_to_use-trial_aperiodic_mean)).^2); % Root mean square error
            RMS_error(indTrial,3,c,o) = RMSE_trial;
            
        end
    end
end
 
%% Plot results roughly using MATLAB
for c = 1:length(condition)
    for o = 1:length(oscillations)
        ddd = squeeze(RMS_error(:,:,c,o));
        
        figure;boxplot(ddd);
        title([oscillations{o} ' + ' condition{c} ' 1/f']);
    end
end

%% Export to .csv file for fancier plotting in Python
condition                 = [];
[condition{1,1:200}]      = deal('BOSC');
[condition{1,101:400}]    = deal('eBOSC');
[condition{1,401:600}]    = deal('fBOSC');

nonlinear_theta = reshape(squeeze(RMS_error(:,:,1,1)),[600,1]);
nonlinear_alpha = reshape(squeeze(RMS_error(:,:,1,2)),[600,1]);
linear_theta    = reshape(squeeze(RMS_error(:,:,2,1)),[600,1]);
linear_alpha    = reshape(squeeze(RMS_error(:,:,2,2)),[600,1]);

condition                   = condition';

t = table(nonlinear_theta,nonlinear_alpha,linear_theta,linear_alpha,condition);
writetable(t,'error1f.csv');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's repeat the simulations for data where BOTH alpha and theta are present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condition     = {'synaptic','linear'};
freq          = [4 10];

ntrials = 200;

RMS_error = zeros(ntrials,3,length(condition));

% For each condition
for c = 1:length(condition)
    
    disp(['Simulating two oscillations + ' condition{c} ' 1/f']);
    
    % Load the aperiodic data
    if strcmp(condition{c},'synaptic')
        load('synaptic.mat');
        aperiodic = synaptic;
        clear synaptic
    else
        disp('Using linear 1/f');
        load('linear.mat');
        aperiodic = linear;
        clear linear
    end
    
    % HP-Filter @ 1Hz, like in typical EEG/MEG studies
    for n = 1:ntrials
        aperiodic(n,:) = ft_preproc_highpassfilter(aperiodic(n,:), 500,...
            1, 3);
    end
    
    %% Simulate Data
    % Theta
    cfg                         = [];
    cfg.freq                    = 4;
    cfg.amplitude               = 5.1; % SNR of the simulated oscillation
    cfg.cycles                  = 6/(1/4); % 6s bursts;
    cfg.time                    = 20; % 20s of simulated data
    cfg.trial                   = ntrials; % How many 'trials'?
    [aperiodic_out, ...
        osc, ...
        data_out1, ...
        osc_array]              = sim_fBOSC(cfg,Fs,aperiodic);
    
    % Alpha
    cfg                         = [];
    cfg.freq                    = 10;
    cfg.amplitude               = 14; % SNR of the simulated oscillation
    cfg.cycles                  = 6/(1/10); % 6s bursts;
    cfg.time                    = 20; % 20s of simulated data
    cfg.trial                   = ntrials; % How many 'trials'?
    [aperiodic_out, ...
        osc, ...
        data_out2, ...
        osc_array]              = sim_fBOSC(cfg,Fs,aperiodic);
    
    % Combine all trials into one array
    data_comb_osc       = zeros(cfg.trial,(2*cfg.time*Fs));
    for i = 1:cfg.trial
        data_comb_osc(i,:) = horzcat(data_out1.trial{i},...
            data_out2.trial{i});
    end
    
    % Make time variable for combined data
    t = 1/Fs;
    time = [t:t:2*cfg.time];
    
    % Make Figure showing ground first trial of data
    ddd = data_comb_osc(1,:);
    figure; plot(time,ddd,'k'); hold on;
    xlabel('Time (s)');
    
    %% fBOSC
    % Set up BOSC parameters
    cfg                     = [];
    cfg.fBOSC.F             = 2.^[1:.125:5.4];    % frequency sampling
    %cfg.fBOSC.F             = [2:0.5:50];
    cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
    cfg.fBOSC.fsample       = Fs;         % current sampling frequency of EEG data
    
    % padding
    cfg.fBOSC.pad.tfr_s         = 1;       % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
    cfg.fBOSC.pad.detection_s   = .5;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
    cfg.fBOSC.pad.background_s  = 1;       % padding of segments for BG (only avoiding edge artifacts)
    
    % fooof
    if strcmp(condition{c},'synaptic')
        cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
    else
        cfg.fBOSC.fooof.aperiodic_mode    = 'fixed';
    end
    
    % threshold settings
    cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.fBOSC.threshold.percentile  = .99;                                      % percentile of background fit for power threshold
    
    cfg.fBOSC.trial_background = 1:ntrials;
    % calculate the sample points for paddding
    cfg.fBOSC.pad.tfr_sample = cfg.fBOSC.pad.tfr_s.*cfg.fBOSC.fsample;                          % automatic sample point calculation
    cfg.fBOSC.pad.detection_sample = cfg.fBOSC.pad.detection_s.*cfg.fBOSC.fsample;              % automatic sample point calculation
    cfg.fBOSC.pad.total_s = cfg.fBOSC.pad.tfr_s + cfg.fBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
    cfg.fBOSC.pad.total_sample = cfg.fBOSC.pad.tfr_sample + cfg.fBOSC.pad.detection_sample;
    cfg.fBOSC.pad.background_sample = cfg.fBOSC.pad.tfr_sample;
    cfg.tmp.channel = 1;
    
    % For every trial
    % Compute wavelet transform for the combined data
    TFR = [];
    for indTrial = 1:size(data_comb_osc,1)
        TFR.trial{indTrial} = BOSC_tf(data_comb_osc(indTrial,:),...
            cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
    end; clear indTrial
    
    % Get mean power
    BG = [TFR.trial{:}];
    mean_pow = 10.^(mean(log10(BG(:,:)),2))';
    
    clear BG
    
    %% Get thresholds + fit line
    % 1. Using fBOSC (FOOOF)
    [fBOSC, pt_fBOSC, dt_fBOSC] = fBOSC_getThresholds(cfg, TFR, []);
    
    % 2. Using eBOSC and BOSC
    cfg.eBOSC = cfg.fBOSC;
    % Remove freqs with oscillatory peaks in theta and alpha
    cfg.eBOSC.F_fit = cfg.fBOSC.F(~(cfg.fBOSC.F >= 2 & cfg.fBOSC.F <= 6));
    cfg.eBOSC.F_fit = cfg.eBOSC.F_fit(~(cfg.eBOSC.F_fit >= 8 & cfg.eBOSC.F_fit <= 12));
    [eBOSC, pt_eBOSC, dt_eBOSC] = eBOSC_getThresholds_multipeak(cfg, TFR, []);
    
    % For every trial
    % Compute wavelet transform for just the aperiodic data
    disp('Performing TFR on aperiodic data without oscillations');
    TFR_aperiodic = [];
    for indTrial = 1:size(data_comb_osc,1)
        
        TFR_aperiodic.trial{indTrial} = BOSC_tf(aperiodic_out(indTrial,:),...
            cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
    end
    
    % Concatenate TFRs into a 3D array
    TFR_pad = [];
    for tr = 1:length(TFR_aperiodic.trial)
        TFR_pad{tr} = TFR_aperiodic.trial{tr}(:,cfg.fBOSC.pad.background_sample+1:...
            end-cfg.fBOSC.pad.background_sample);
    end
    
    % Produce Figure for fBOSC and eBOSC showing fit for the first
    % trial
    figure;
    plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
    plot(log10(cfg.fBOSC.F), log10(fBOSC.static.mp), '--r','LineWidth',3);
    plot(log10(cfg.fBOSC.F), log10(fBOSC.static.pt), 'g','LineWidth',3);
    set(gca,'FontSize',18);
    ylabel('Log Power (a.u)','FontSize',25);
    xlabel('Log Frequency (Hz)','FontSize',25);   
    xlim([0.25 1.7]);
    ylim([2.3 5]);
    title('fBOSC');
    %print('fBOSC','-dpng','-r400');
    
    figure;
    data = plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
    plot(log10(cfg.fBOSC.F), log10(eBOSC.static.mp), '--r','LineWidth',3);
    plot(log10(cfg.fBOSC.F), log10(eBOSC.static.pt), 'g','LineWidth',3);
    set(gca,'FontSize',18);
    ylabel('Log Power (a.u)','FontSize',25);
    xlabel('Log Frequency (Hz)','FontSize',25);
    title('eBOSC');
    xlim([0.25 1.7]);
    %ylim([2.3 5]);
    %print('eBOSC','-dpng','-r400');
    drawnow;
    
    figure;
    data = plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
    plot(log10(cfg.fBOSC.F), log10(fBOSC.static.mp_old), '--r','LineWidth',3);
    plot(log10(cfg.fBOSC.F), log10(fBOSC.static.pt_old), 'g','LineWidth',3);
    set(gca,'FontSize',18);
    ylabel('Log Power (a.u)','FontSize',25);
    xlabel('Log Frequency (Hz)','FontSize',25);
    xlim([0.25 1.7]);
    %ylim([2.3 5]);
    title('BOSC');
    drawnow;
    %print('BOSC','-dpng','-r400');


    %% Calculate RMSE for every trial
        
    figure;
    plot(log10(cfg.fBOSC.F), log10(mean_pow), 'black','LineWidth',4);  hold on;
    
    for indTrial = 1:ntrials
        % First BOSC
        trial_aperiodic_mean = mean(log10(TFR_aperiodic.trial{indTrial}(:,:)),2);
        p1 = plot(log10(cfg.fBOSC.F),trial_aperiodic_mean,'LineWidth',0.5);
        p1.Color(4) = 0.4;
        fit_to_use = log10(fBOSC.static.mp_old)';
        RMSE_trial = sqrt(mean((trial_aperiodic_mean-fit_to_use)).^2); % Root mean square error
        RMS_error(indTrial,1,c) = RMSE_trial;
        
        % Now eBOSC with frequencies excluded
        fit_to_use = log10(eBOSC.static.mp)';
        RMSE_trial = sqrt(mean((fit_to_use-trial_aperiodic_mean)).^2); % Root mean square error
        RMS_error(indTrial,2,c) = RMSE_trial;
        
        % Now fBOSC (using FOOOF)
        fit_to_use = log10(fBOSC.static.mp)';
        RMSE_trial = sqrt(mean((fit_to_use-trial_aperiodic_mean)).^2); % Root mean square error
        RMS_error(indTrial,3,c) = RMSE_trial;
        
    end
end

%% Export to .csv file for fancy plotting in python
condition                 = [];
[condition{1,1:200}]      = deal('BOSC');
[condition{1,101:400}]    = deal('eBOSC');
[condition{1,401:600}]    = deal('fBOSC');

nonlinear_theta_alpha     = reshape(squeeze(RMS_error(:,:,1)),[600,1]);
linear_theta_alpha        = reshape(squeeze(RMS_error(:,:,2)),[600,1]);

condition                   = condition';

t = table(nonlinear_theta_alpha,linear_theta_alpha,condition);
writetable(t,'2osc_error1f.csv');





