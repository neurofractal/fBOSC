%%
% An example analysis pipeline using fBOSC, and simulated data
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
root = fileparts(matlab.desktop.editor.getActiveFilename);
cd(fullfile(root,'..'));

start_fBOSC

time1   = [];
time2   = [];

trials  = horzcat([1:20:200],200);

for time_ = 1:length(trials)
    
    %% Make simulated data with 10hz oscillations embedded within a
    %  non-linear 1/f aperiodic signal
    
    % Load the previously computed non-linear 1/f aperiodic signal
    load(fullfile(root,'..','..','fBOSC_validation','synaptic.mat'));
    
    % Settings
    ntrials                     = trials(time_);   % Number of trials
    Fs                          = 500;  % Sampling Rate
    SNR                         = 40;   % SNR of oscillation
    len_of_trial                = 20;   % Length of each trial in s
    
    % Trim aperiodic data down to size of ntrials*len_of_trial
    synaptic                   = synaptic(1:ntrials,[1:(len_of_trial*Fs)]);
    
    % HP-Filter @ 1Hz, like in typical EEG/MEG studies
    for n = 1:ntrials
        synaptic(n,:) = ft_preproc_highpassfilter(synaptic(n,:), Fs,...
            1, 3);
    end
    
    % Simulate
    cfg                         = [];
    cfg.freq                    = 10; % Simulated a 10Hz oscillation
    cfg.amplitude               = SNR; % SNR of the simulated oscillation
    cfg.cycles                  = 4/(1/10); % How many cycles: For 6s
    cfg.time                    = 20; % 20s of simulated data
    cfg.trial                   = ntrials; % How many 'trials'?
    [aperiodic_out_alpha,...
        osc_alpha,...
        data_alpha,~]...
        = sim_fBOSC(cfg,Fs,synaptic);
    
    %% Plot the aperiodic and periodic data separately
    figure; subplot(3,1,1);plot(data_alpha.time{1},osc_alpha(1,:));
    title('Periodic Data');
    subplot(3,1,2);plot(data_alpha.time{1},aperiodic_out_alpha(1,:));
    title('Aperiodic Data');
    subplot(3,1,3);plot(data_alpha.time{1},data_alpha.trial{1});
    title('Combined Data');
    drawnow;
    
    %% Set-up fBOSC parameters
    
    % general setup
    cfg.fBOSC.F                 = 2.^[1:.125:5.4];
    %cfg.fBOSC.F                     = [2:0.5:40];
    cfg.fBOSC.wavenumber            = 6;
    cfg.fBOSC.fsample               = Fs;
    
    % padding
    cfg.fBOSC.pad.tfr_s             = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
    cfg.fBOSC.pad.detection_s       = 0.1;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
    cfg.fBOSC.pad.background_s      = 0.1;      % padding of segments for BG (only avoiding edge artifacts)
    
    % fooof parameters - fit with fixed line or allow a knee
    cfg.fBOSC.fooof.aperiodic_mode  = 'knee';
    cfg.fBOSC.fooof.version         = 'python';
    
    % threshold settings
    cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold
    
    % episode post-processing
    cfg.fBOSC.postproc.use          = 'no';        % Post-processing turned off for now
    
    % general processing settings
    cfg.fBOSC.channel               = [1]; % select posterior channels (default: all)
    cfg.fBOSC.trial                 = []; % select trials (default: all)
    cfg.fBOSC.trial_background      = []; % select trials for background (default: all)
    
    %% New version
    clear fBOSC
    tic
    [fBOSC, cfg] = fBOSC_wrapper2(cfg, data_alpha);
    time1(time_) = toc;
    
    %% Old version
    clear fBOSC
    tic
    [fBOSC, cfg] = fBOSC_wrapper(cfg, data_alpha);
    time2(time_) = toc;
    
    disp(time1);
    disp(time2);
    
    clear data_alpha
    
end


figure; plot(trials*20,time1,'LineWidth',2); hold on;
plot(trials*20,time2,'LineWidth',2);
set(gca, 'YScale', 'log');
set(gca,'FontSize',18);
ylabel({'fBOSC processing'; 'time (s)'},'FontSize',22);
xlabel('Data length (s)','FontSize',22);
legend({'speedy fBOSC','old fBOSC'},'Location','SouthWest');









