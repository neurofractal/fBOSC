%%
start_fBOSC
cd(fullfile(cd,'simulate'));

%% Make simulated data with oscillations embedded within a 
%  non-linear 1/f aperiodic signal

% Load the previously computed non-linear 1/f aperiodic signal
load('aperiodic1.mat');


ntrials                     = 200;
Fs                          = 500; %sampling rate
SNR                         = [0:2:24];

hit_rate            = zeros(3,2,ntrials,length(SNR));
false_alarm         = zeros(3,2,ntrials,length(SNR));

hit_rate_alpha      = zeros(3,2,ntrials,length(SNR));
false_alarm_alpha   = zeros(3,2,ntrials,length(SNR));

condition   = {'synaptic'};

RMSE_SNR            = zeros(3,1,ntrials,length(SNR));

for c = 1:length(condition)
    for s = 1:length(SNR)
        disp(['Condition: ' condition{c} '; SNR: ' num2str(SNR(s))]);
        %% Load the aperiodic data from python
        % aperiodic3 is the synaptic activity
        % aperiodic4 is the linear 1/f
        
        if strcmp(condition{c},'synaptic')
            load('aperiodic1.mat');
            aperiodic = aperiodic1;
        else
            load('aperiodic2.mat');
            aperiodic = aperiodic2;
        end
        
        %%
        cfg                         = [];
        cfg.freq                    = 10; % Simulated a 10Hz oscillation
        cfg.amplitude               = SNR(s); % SNR of the simulated oscillation
        cfg.cycles                  = 64; % How many cycles?
        cfg.time                    = 20; % 20s of simulated data
        cfg.trial                   = ntrials; % How many 'trials'?
        [aperiodic_out_alpha, ...
            osc_alpha, ...
            data_alpha, ...
            alpha_osc_array]      = sim_fBOSC(cfg,Fs,aperiodic);
        
        cfg                         = [];
        cfg.freq                    = 4; % Simulated a 10Hz oscillation
        cfg.amplitude               = SNR(s); % SNR of the simulated oscillation
        cfg.cycles                  = 26; % How many cycles?
        cfg.time                    = 20; % 20s of simulated data
        cfg.trial                   = ntrials; % How many 'trials'?
        [aperiodic_out_theta, ...
            osc_theta, ...
            data_theta,...
            theta_osc_array]        = sim_fBOSC(cfg,Fs,aperiodic);
        
        % Plot the theta for sanity
        %     figure; subplot(2,1,1); plot(osc_alpha(1,:));
        %     subplot(2,1,2);plot(osc_theta(1,:));
        
        
        % Combine all trials into one array
        data_comb_osc       = zeros(cfg.trial,(2*cfg.time*Fs));
        for i = 1:cfg.trial
            data_comb_osc(i,:) = horzcat(data_alpha.trial{i},data_theta.trial{i});
        end
        
        % Combine all aperiodic data into one array
        data_comb_aperiodic = horzcat(aperiodic_out_alpha(:,:),...
            aperiodic_out_theta(:,:));
        
        %%
        % Set up BOSC parameters
        cfg                     = [];
        cfg.fBOSC.F             = 2.^[1:.125:5.4];    % frequency sampling
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
        
        %         % Plot for sanity
        %         figure; subplot(2,1,1); plot(data_comb_osc(2,:)); hold on;
        %         subplot(2,1,2);imagesc(TFR.trial{2});
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Using multiple thresholds for eBOSC and fBOSC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % For hit-rate vs false alarm rates we need to use multiple thresholds
        cfg.fBOSC.threshold.percentile  = [0.95];
        
        % fBOSC
        [fBOSC, pt_all, dt, mean_pow, mp]  = fBOSC_getThresholds_multi_thresholds(cfg, TFR, []);
        
        % Get thresholds using eBOSC
        cfg.eBOSC = cfg.fBOSC;
        % Remove freqs with oscillatory peaks
        cfg.eBOSC.F_fit = cfg.fBOSC.F(~(cfg.fBOSC.F >= 2 & cfg.fBOSC.F <= 6));
        cfg.eBOSC.F_fit = cfg.eBOSC.F_fit(~(cfg.eBOSC.F_fit >= 8 & cfg.eBOSC.F_fit <= 12));
        [eBOSC, pt_all_eBOSC,pt_all_BOSC, dt]       = eBOSC_getThresholds_multi_thresholds(cfg, TFR, []);
                
        %
        %Plot the differnece between the thresholds
        figure;
        plot(log10(cfg.fBOSC.F),log10(pt_all_BOSC),'--','Color',...
            [1.0000    0.3608    0.0157],'LineWidth',3); hold on;
        plot(log10(cfg.fBOSC.F),log10(pt_all_eBOSC),'--','Color',...
            [0.7843    0.3451    0.5569],'LineWidth',3); hold on;
        plot(log10(cfg.fBOSC.F),log10(pt_all),'--','Color',...
            [0.0157    0.3608    1.0000],'LineWidth',3); hold on;
        %plot(log10(cfg.fBOSC.F),log10(mp)); hold on;
        plot(log10(cfg.fBOSC.F),log10(mean_pow),'k','LineWidth',4);
        set(gca,'FontSize',18);
        ylabel('Power (a.u.)','FontSize',22);
        xlabel('Frequency (Hz)','FontSize',22);
        xlim([0.2 1.8]);
        ylim([2 5]);
        legend({'BOSC','eBOSC','fBOSC'});
        title(['SNR: ' num2str(SNR(s))]);
        drawnow;
        print(['SNR: ' num2str(SNR(s)) '_' condition{c}],'-dpng','-r300');
        
        %% Calculate RMSE for each trial
        
        % Combine all aperiodic data into one array
        data_comb_aperiodic = horzcat(aperiodic_out_alpha(:,:),...
            aperiodic_out_theta(:,:));
        
        disp('Performing TFR on aperiodic data without oscillations');
        TFR_aperiodic = [];
        for indTrial = 1:1:ntrials
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
        mean_pow2 = 10.^(mean(log10(BG(:,:)),2))';
  
        %% Calculate RMSE
        RMSE_all    = zeros(ntrials,3);
        
        for i = 1:ntrials
            % First BOSC
            
            figure; plot(log10(cfg.fBOSC.F),log10(mean_pow2),'k','LineWidth',5); hold on;
            plot(log10(cfg.fBOSC.F),log10(pt_all_BOSC),'--','Color',...
                [1.0000    0.3608    0.0157],'LineWidth',3); hold on;
            plot(log10(cfg.fBOSC.F),log10(pt_all_eBOSC),'--','Color',...
                [0.7843    0.3451    0.5569],'LineWidth',3); hold on;
            plot(log10(cfg.fBOSC.F),log10(pt_all),'--','Color',...
                [0.0157    0.3608    1.0000],'LineWidth',3); hold on;
            
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
        
        %% Get logical array for known times of alpha and theta oscillations
        alpha_osc_array;
        theta_osc_array = theta_osc_array + size(osc_alpha,2);
        
        osc_combined = [osc_alpha(1,:) osc_theta(1,:)];
        
        % Create logical array for theta
        theta_logical = zeros(1,size(osc_combined,2));
        theta_logical(theta_osc_array) = 1;
        %     figure; plot(osc_combined); hold on;
        %     plot(theta_logical);
        theta_logical = theta_logical(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample)
        theta_logical = theta_logical(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
        theta_logical = logical(theta_logical);
        
        % Create logical array for alpha
        alpha_logical = zeros(1,size(osc_combined,2));
        alpha_logical(alpha_osc_array) = 1;
        %     figure; plot(osc_combined); hold on;
        %plot(alpha_logical);
        alpha_logical = alpha_logical(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample)
        alpha_logical = alpha_logical(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
        alpha_logical = logical(alpha_logical);
        
        %% Start of trial loop
        
        for indTrial = 1:ntrials
            % get wavelet transform for single trial
            % tfr padding is removed to avoid edge artifacts from the wavelet
            % transform. Note that a padding fpr detection remains attached so that there
            % is no problems with too few sample points at the edges to
            % fulfill the duration criterion.
            TFR_ = TFR.trial{indTrial}(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample);
            disp(indTrial);
            cfg.tmp.trial = indTrial; % encode current trial for later
            cfg.tmp.channel = 1;
            
            %% fBOSC
            detected = zeros(size(TFR_));
            
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt_all(1,f),dt(f),...
                    cfg.fBOSC.fsample);
            end; clear f
            
            % remove padding for detection (matrix with padding required for refinement)
            detected = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
            
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),1))>0);
            
            %% THETA
            % Calculate hit rate
            hits = detected_times(theta_logical);
            hit_rate(3,c,indTrial,s) = length(find(hits==1))/length(theta_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~theta_logical);
            false_alarm(3,c,indTrial,s) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
            
            clear hits falses detected_times 
            
             %% ALPHA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > 9.5 & cfg.fBOSC.F < 10.5,:),1))>0);
                        
            % Calculate hit rate
            hits = detected_times(alpha_logical);
            hit_rate_alpha(3,c,indTrial,s) = length(find(hits==1))/length(alpha_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~alpha_logical);
            false_alarm_alpha(3,c,indTrial,s) = length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
            
            clear hits falses detected_times detected
            
            %% eBOSC
            detected = zeros(size(TFR_));
            
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_eBOSC(1,f),dt(f),...
                    cfg.fBOSC.fsample);
            end; clear f
            
            % remove padding for detection (matrix with padding required for refinement)
            detected = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
            
            %% THETA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),1))>0);
            
            % Calculate hit rate
            hits = detected_times(theta_logical);
            hit_rate(2,c,indTrial,s) = length(find(hits==1))/length(theta_osc_array);
            
            
            % Calculate false alarm rate
            falses = detected_times(~theta_logical);
            false_alarm(2,c,indTrial,s) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
            
            clear hits falses detected_times
            
            %% ALPHA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > 9.5 & cfg.fBOSC.F < 10.5,:),1))>0);
            
            % Calculate hit rate
            hits = detected_times(alpha_logical);
            hit_rate_alpha(2,c,indTrial,s) = length(find(hits==1))/length(alpha_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~alpha_logical);
            false_alarm_alpha(2,c,indTrial,s) = length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
            
            clear hits falses detected_times detected
            
             %% BOSC
            detected = zeros(size(TFR_));
            
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_BOSC(1,f),dt(f),...
                    cfg.fBOSC.fsample);
            end; clear f
            
            % remove padding for detection (matrix with padding required for refinement)
            detected = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
            
            %% THETA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > 3.5 & cfg.fBOSC.F < 4.5,:),1))>0);
            
            % Calculate hit rate
            hits = detected_times(theta_logical);
            hit_rate(1,c,indTrial,s) = length(find(hits==1))/length(theta_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~theta_logical);
            false_alarm(1,c,indTrial,s) = length(find(falses==1))/(length(detected_times)-length(theta_osc_array));
            
            clear hits falses detected_times
            
            %% ALPHA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > 9.5 & cfg.fBOSC.F < 10.5,:),1))>0);
            
            % Calculate hit rate
            hits = detected_times(alpha_logical);
            hit_rate_alpha(1,c,indTrial,s) = length(find(hits==1))/length(alpha_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~alpha_logical);
            false_alarm_alpha(1,c,indTrial,s) = length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
            
            clear hits falses detected_times detected
            
        end
    end
end

%% Plot Hit-Rate vs SNR
h1= figure;
set(h1,'Position',[10 10 2000 500])

cols = [0.0157    0.3608    1.0000;0.7843    0.3451    0.5569; 
    1.0000    0.3608    0.0157];

for i = 1:3
    subplot(1,3,i);
    plot(SNR,squeeze(mean(hit_rate_alpha(i,1,:,:))),'--s',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    plot(SNR,squeeze(mean(hit_rate(i,1,:,:))),'-o',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    set(gca,'FontSize',20);
    xlabel('SNR','FontSize',24);
    ylabel('Hit Rate','FontSize',24);
    xticks([0:4:24]);
    ylim([0 1]);
end
print('HR_SNR','-dpng','-r300');


%% Plot False Alarm Rate vs SNR

h1= figure;
set(h1,'Position',[10 10 2000 500])

cols = [0.0157    0.3608    1.0000;0.7843    0.3451    0.5569; 
    1.0000    0.3608    0.0157];

for i = 1:3
    subplot(1,3,i);
    plot(SNR,squeeze(mean(false_alarm_alpha(i,1,:,:))),'--s',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    plot(SNR,squeeze(mean(false_alarm(i,1,:,:))),'-o',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    set(gca,'FontSize',20);
    xlabel('SNR','FontSize',24);
    ylabel('False Alarm Rate','FontSize',24);
    %xticks([0:4:24]);
    %ylim([0 1]);
end

print('FA_SNR','-dpng','-r300');

%% d'
h1= figure;
set(h1,'Position',[10 10 2000 500])

cols = [0.0157    0.3608    1.0000;0.7843    0.3451    0.5569; 
    1.0000    0.3608    0.0157];

for i = 1:3
    subplot(1,3,i);
    
    % Theta
    FA = squeeze(mean(false_alarm(i,1,:,:)));
    HR = squeeze(mean(hit_rate(i,1,:,:)));
    
    d_prime = (norminv(HR)*-1)-(norminv(FA)*-1);
    
    % Alpha
    FA_alpha = squeeze(mean(false_alarm_alpha(i,1,:,:)));
    HR_alpha = squeeze(mean(hit_rate_alpha(i,1,:,:)));
    
    d_prime_alpha = (norminv(HR_alpha)*-1)-(norminv(FA_alpha)*-1);
    
    
    plot(SNR,d_prime,'--s',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    plot(SNR,d_prime_alpha,'-o',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    set(gca,'FontSize',20);
    xlabel('SNR','FontSize',24);
    ylabel('d prime','FontSize',24);
    %xticks([0:4:24]);
    %ylim([0 1]);
end






% Calculate d'

figure;
plot(SNR,squeeze(mean(false_alarm_alpha(2,1,:,:))),'--s',...
    'Color',[0.7843    0.3451    0.5569],'LineWidth',2); hold on;
plot(SNR,squeeze(mean(false_alarm(2,1,:,:))),'-o',...
    'Color',[0.7843    0.3451    0.5569],'LineWidth',2); hold on;
set(gca,'FontSize',20);
xlabel('SNR','FontSize',24);
ylabel('Hit Rate','FontSize',24);
subplot(1,3,3);
plot(SNR,squeeze(mean(hit_rate_alpha(3,1,:,:))),'--s',...
    'Color',[1.0000    0.3608    0.0157],'LineWidth',2); hold on;
plot(SNR,squeeze(mean(hit_rate(3,1,:,:))),'-o',...
    'Color',[1.0000    0.3608    0.0157],'LineWidth',2); hold on;
set(gca,'FontSize',20);
xlabel('SNR','FontSize',24);
ylabel('Hit Rate','FontSize',24);






[1.0000    0.3608    0.0157],'LineWidth',3); hold on;
%         plot(log10(cfg.fBOSC.F),log10(pt_all_eBOSC),'--','Color',...
%             [0.7843    0.3451    0.5569],'LineWidth',3); hold on;
%         plot(log10(cfg.fBOSC.F),log10(pt_all),'--','Color',...
%             [0.0157    0.3608    1.0000]
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'synaptic_HR_alpha.gif';

for s=fliplr(1:12)
    for i = 1:3
    subplot(3,1,i);
    scatter(squeeze(mean(hit_rate_alpha(i,1,:,s))),squeeze(mean(hit_rate(i,1,:,s)))); hold on;
    
%     legend({'BOSC Alpha','eBOSC Alpha','fBOSC Alpha',...
%         'BOSC','eBOSC','fBOSC'},'Location','EastOutside');
    % scatter(squeeze(false_alarm(1,1,:)),squeeze(hit_rate(1,1,:)),2,'k'); hold on;
    % scatter(squeeze(false_alarm(1,2,:)),squeeze(hit_rate(1,2,:)),2,'b'); hold on;
    % scatter(squeeze(false_alarm(2,1,:)),squeeze(hit_rate(2,1,:)),2,'r'); hold on;
    % scatter(squeeze(false_alarm(2,2,:)),squeeze(hit_rate(2,2,:)),2,'g'); hold on;
    set(gca,'FontSize',18);
    ylim([0 1]);
    xlim([0 0.15]);
    ylabel('Hit Rate','FontSize',20);
    xlabel({'False Alarm';'Rate'},'FontSize',20);
    title(['SNR: ' num2str(SNR(s))]);
    end
    drawnow;
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File 
      if s == 12 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
      
end

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'linear_HR.gif';

for s=fliplr(1:12)
    scatter(squeeze(mean(false_alarm(1,2,:,s))),squeeze(mean(hit_rate(1,2,:,s))),100,...
        [0.0157    0.3608    1.0000],'filled'); hold on;
    scatter(squeeze(mean(false_alarm(2,2,:,s))),squeeze(mean(hit_rate(2,2,:,s))),100,...
        [0.7843    0.3451    0.5569],'filled'); hold on;
    scatter(squeeze(mean(false_alarm(3,2,:,s))),squeeze(mean(hit_rate(3,2,:,s))),100,...
        [1.0000    0.3608    0.0157],'filled'); hold on;
    legend({'BOSC','eBOSC','fBOSC'},'Location','EastOutside');
    % scatter(squeeze(false_alarm(1,1,:)),squeeze(hit_rate(1,1,:)),2,'k'); hold on;
    % scatter(squeeze(false_alarm(1,2,:)),squeeze(hit_rate(1,2,:)),2,'b'); hold on;
    % scatter(squeeze(false_alarm(2,1,:)),squeeze(hit_rate(2,1,:)),2,'r'); hold on;
    % scatter(squeeze(false_alarm(2,2,:)),squeeze(hit_rate(2,2,:)),2,'g'); hold on;
    set(gca,'FontSize',18);
    ylim([0 1]);
    xlim([0.0 0.1]);
    ylabel('Hit Rate','FontSize',20);
    xlabel({'False Alarm';'Rate'},'FontSize',20);
    title(['SNR: ' num2str(SNR(s))]);
    drawnow;
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File 
      if s == 12 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end

%%
figure; 
scatter(squeeze(mean(false_alarm(1,1,:,:),3)),squeeze(mean(hit_rate(1,1,:,:),3))); hold on;
scatter(squeeze(mean(false_alarm(1,2,:,:),3)),squeeze(mean(hit_rate(1,2,:,:),3))); hold on;
scatter(squeeze(mean(false_alarm(2,1,:,:),3)),squeeze(mean(hit_rate(2,1,:,:),3))); hold on;
scatter(squeeze(mean(false_alarm(2,2,:,:),3)),squeeze(mean(hit_rate(2,2,:,:),3))); hold on;

figure;
plot(squeeze(mean(hit_rate(1,1,:,:),3)),SNR); hold on;
plot(squeeze(mean(false_alarm(2,1,:,:),3)),SNR); hold on;


%% Make an animated gif
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'synaptic.gif';
for n = 1:12
    % Draw plot for y = x.^n
    d = imread(['SNR: ' num2str(SNR(n)) '_synaptic.png']);
    % Capture the plot as an image
    
    [imind,cm] = rgb2ind(d,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'linear.gif';
for n = 1:12
    % Draw plot for y = x.^n
    d = imread(['SNR: ' num2str(SNR(n)) '_linear.png']);
    % Capture the plot as an image
    
    [imind,cm] = rgb2ind(d,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%% Let's export the data for 
% 
for s=1:12
    
    theta                             = squeeze(hit_rate(:,1,:,s));
    alpha                             = squeeze(hit_rate_alpha(:,1,:,s));
    
    theta_collapsed                 = reshape(theta,[1500,1]);
    alpha_collapsed                 = reshape(alpha,[1500,1]);
   
    condition                       = [];
    [condition{1,1:500}]            = deal('BOSC');
    [condition{1,501:1000}]         = deal('eBOSC');
    [condition{1,1001:1500}]        = deal('fBOSC');
    condition                       = condition';
    condition                       = repmat(condition,2,1);
    
    oscillation                       = [];
    [oscillation{1,1:1500}]           = deal('theta');
    [oscillation{1,1501:3000}]        = deal('alpha');
    oscillation                     = oscillation';
    
    HR = vertcat(theta_collapsed, alpha_collapsed);
    
    t = table(HR,condition, oscillation);
    writetable(t,['hr_SNR_' num2str(s) '.csv'])
end


scatter(squeeze(mean(false_alarm(1,1,:,:),3)),squeeze(mean(hit_rate(1,1,:,:),3))); hold on;







        
        
