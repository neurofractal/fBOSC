%%
% Investigate how the hit-rate and false alarm rate differ between BOSC,
% eBOSC and fBOSC, when ALSO varying SNR
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
start_fBOSC
cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

%% Load the aperiodic data from python
% aperiodic1 is the synaptic activity
% aperiodic2 is the linear 1/f
load('aperiodic1.mat');
load('aperiodic2.mat');

%% Specify the different simulation options
condition     = {'synaptic'};

SNR = [0:2:24];

ntrials = 200;

error_SNR = zeros(ntrials,3,length(SNR));

% For each condition
for c = 1:length(condition)
    for s = 1:length(SNR)
    %% Get the correct data
    if strcmp(condition{c},'synaptic')
        aperiodic_for_this_loop = aperiodic1;
        
    else
        aperiodic_for_this_loop = aperiodic2;
    end
    
   
        disp(['Simulating data with ' condition{c}...
            ' 1/f and oscillations of ' num2str(SNR(s))...
            ' SNR']);
        
        %% Simulate Data
        % If 2 oscillations, simulate alpha, or else similate theta
        cfg                         = [];
        
        cfg.freq                    = 10; % Simulated a 10Hz oscillation
        cfg.amplitude               = SNR(s); % SNR of the simulated oscillation
        cfg.cycles                  = 64; % How many cycles?
        cfg.time                    = 20; % 20s of simulated data
        cfg.trial                   = ntrials; % How many 'trials'?
        [aperiodic_out_alpha, ...
            osc_alpha, ...
            data_alpha, ...
            alpha_osc_array]      = sim_fBOSC(cfg,500,aperiodic_for_this_loop);
        
        % Theta-band
        cfg                         = [];
        cfg.freq                    = 4; % Simulated a 10Hz oscillation
        cfg.amplitude               = SNR(s); % SNR of the simulated oscillation
        cfg.cycles                  = 26; % How many cycles?
        cfg.time                    = 20; % 20s of simulated data
        cfg.trial                   = ntrials; % How many 'trials'?
        [aperiodic_out_theta, ...
            osc_theta, ...
            data_theta,...
            theta_osc_array]        = sim_fBOSC(cfg,500,aperiodic_for_this_loop);
        
        % Plot the oscillations for sanity
        figure; subplot(2,1,1); plot(osc_alpha(1,:));
        subplot(2,1,2);plot(osc_theta(1,:));
        drawnow;
        
        % Combine all trials into one array
        data_comb_osc       = zeros(cfg.trial,(2*cfg.time*500));
        for i = 1:cfg.trial
            data_comb_osc(i,:) = horzcat(data_alpha.trial{i},data_theta.trial{i});
        end
        
        % Combine all aperiodic data into one array
        data_comb_aperiodic = horzcat(aperiodic_out_alpha(:,:),...
            aperiodic_out_theta(:,:));
        
        %% fBOSC
        % Set up BOSC parameters
        cfg                     = [];
        cfg.fBOSC.F             = 2.^[1:.125:5.4];    % frequency sampling
        cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
        cfg.fBOSC.fsample       = 500;         % current sampling frequency of EEG data
        
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
        cfg.fBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold
        
        
        cfg.fBOSC.trial_background = 1:size(data_comb_osc,1);
        % calculate the sample points for paddding
        cfg.fBOSC.pad.tfr_sample = cfg.fBOSC.pad.tfr_s.*cfg.fBOSC.fsample;                          % automatic sample point calculation
        cfg.fBOSC.pad.detection_sample = cfg.fBOSC.pad.detection_s.*cfg.fBOSC.fsample;              % automatic sample point calculation
        cfg.fBOSC.pad.total_s = cfg.fBOSC.pad.tfr_s + cfg.fBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
        cfg.fBOSC.pad.total_sample = cfg.fBOSC.pad.tfr_sample + cfg.fBOSC.pad.detection_sample;
        cfg.fBOSC.pad.background_sample = cfg.fBOSC.pad.tfr_sample;
        cfg.tmp.channel = 1;
        
        %% Compute Wavelet Transform for the combined data
        TFR = [];
        for indTrial = 1:size(data_comb_osc,1)
            TFR.trial{indTrial} = BOSC_tf(data_comb_osc(indTrial,:),...
                cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
        end; clear indTrial
       
        %% Get thresholds + fit line
        % 1. Using FOOOF
        [fBOSC, pt_fBOSC, dt_fBOSC] = fBOSC_getThresholds(cfg, TFR, []);
        
        % 2. Using eBOSC and BOSC
        cfg.eBOSC = cfg.fBOSC;
        % Remove freqs with oscillatory peaks
        cfg.eBOSC.F_fit = cfg.fBOSC.F(~(cfg.fBOSC.F >= 2 & cfg.fBOSC.F <= 6));
        cfg.eBOSC.F_fit = cfg.eBOSC.F_fit(~(cfg.eBOSC.F_fit >= 8 & cfg.eBOSC.F_fit <= 12));
        [eBOSC, pt_eBOSC, dt_eBOSC] = eBOSC_getThresholds_multipeak(cfg, TFR, []);
        
        %% Compute wavelet transform for just the aperiodic data
        disp('Performing TFR on aperiodic data without oscillations');
        TFR_aperiodic = [];
        for indTrial = 1:1:size(data_comb_osc,1)
            TFR_aperiodic.trial{indTrial} = BOSC_tf(data_comb_aperiodic(indTrial,:),...
                cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
        end; clear indTrial
        
        % Concatenate TFRs into a 3D array
        TFR_pad = [];
        for tr = 1:length(TFR_aperiodic.trial)
            TFR_pad{tr} = TFR_aperiodic.trial{tr}(:,cfg.fBOSC.pad.background_sample+1:...
                end-cfg.fBOSC.pad.background_sample);
        end
        
        BG = [TFR_pad{:}];
        mean_pow = 10.^(mean(log10(BG(:,:)),2))';
        
        
        %% Calculate RMSE
        RMSE_all    = zeros(ntrials,3);
        
        for i = 1:ntrials
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
        
        % Add to variable outside the loops
        error_SNR(:,:,s) = RMSE_all;
    end
end
  
mean_error = squeeze(mean(error_SNR,1));

figure;
plot(SNR,mean_error(1,:),'--','Color',...
    [1.0000    0.3608    0.0157],'LineWidth',3); hold on;
plot(SNR,mean_error(2,:),'--','Color',...
    [0.7843    0.3451    0.5569],'LineWidth',3); hold on;
plot(SNR,mean_error(3,:),'--','Color',...
    [0.0157    0.3608    1.0000],'LineWidth',3); hold on;
set(gca,'FontSize',20);
ylabel('RMSE','FontSize',30);
xlabel('SNR','FontSize',30);
legend({'BOSC','eBOSC','fBOSC'},'Location','NorthWest');
drawnow;
ylim([0 0.2]);
print('RMSE_SNR','-dpng','-r300');

figure;boxplot(RMSE_all,{'BOSC','eBOSC','fBOSC'});
drawnow;

        RMSE                        = reshape(RMSE_all,[1500,1]);
        fit_method                   = [];
        [fit_method{1,1:500}]        = deal('BOSC');
        [fit_method{1,501:1000}]     = deal('eBOSC');
        [fit_method{1,1001:1500}]    = deal('fBOSC');
        fit_method                   = fit_method';
        
        t = table(RMSE,fit_method);
        writetable(t,[fit_method{c},'_with_' num2str(number_of_osc(n)) 'osc.csv']);
    end
end





