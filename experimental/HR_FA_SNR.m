%%
start_fBOSC
cd(fullfile(cd,'simulate'));

%% Make simulated data with oscillations embedded within a 
%  non-linear 1/f aperiodic signal

ntrials                     = 100;
Fs                          = 500; %sampling rate
SNR                         = [5.1 13]; %1.5 4
            
condition   = {'synaptic'};

freq = [4 10];

%% Loop of thresholds
% For hit-rate vs false alarm rates we need to use multiple thresholds
thresholds      = horzcat([0.05:0.1:0.9],[0.91:0.01:0.99]);
HR_thresh_alpha = zeros(length(freq),3,length(thresholds),ntrials);
FA_thresh_alpha = zeros(length(freq),3,length(thresholds),ntrials);


for fr = 1:length(freq)
    %disp(['Condition: ' condition{c} '; SNR: ' num2str(SNR(s))]);
    disp(['Frequency: ' num2str(freq(fr)) 'Hz']);
    %% Load the aperiodic data from python
    % aperiodic3 is the synaptic activity
    % aperiodic4 is the linear 1/f
    
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
    cfg.cycles                  = 6/(1/freq(fr)); % How many cycles?
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
    
    %thresholds                      = horzcat([0.91:0.01:0.99],[0.991:0.001:0.999]);
    cfg.fBOSC.threshold.percentile  = thresholds;
    
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
    
    
    
    
    %% For each threshold use the three BOSC approaches to find the oscillation!
    for thresh = 1:length(thresholds)
        disp(['Threshold: ' num2str(thresholds(thresh))]);
        
        %% Start of trial loop
        for indTrial = 1:ntrials
            
            TFR_ = TFR.trial{indTrial};
            cfg.tmp.trial = indTrial; % encode current trial for later
            cfg.tmp.channel = 1;
            
            %% fBOSC
            detected = zeros(size(TFR_));
            
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt_all(thresh,f),dt(f),...
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
            HR_thresh_alpha(fr,3,thresh,indTrial) = length(find(hits==1))/length(alpha_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~log_alpha);
            FA_thresh_alpha(fr,3,thresh,indTrial) = length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
            
            clear hits falses detected_times detected
            
            
            %% eBOSC
            detected = zeros(size(TFR_));
            
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_eBOSC(thresh,f),dt(f),...
                    cfg.fBOSC.fsample);
            end; clear f
            
            % ALPHA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > freq(fr)-0.5 & cfg.fBOSC.F < freq(fr)+0.5,:),1))>0);
            
            % Calculate hit rate
            hits = detected_times(log_alpha);
            HR_thresh_alpha(fr,2,thresh,indTrial) = length(find(hits==1))/length(alpha_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~log_alpha);
            FA_thresh_alpha(fr,2,thresh,indTrial) = length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
            
            clear hits falses detected_times detected
            
            %% BOSC
            detected = zeros(size(TFR_));
            
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt_all_BOSC(thresh,f),dt(f),...
                    cfg.fBOSC.fsample);
            end; clear f
            
            
            % ALPHA
            detected_times = single(squeeze(nanmean(detected(...
                cfg.fBOSC.F > freq(fr)-0.5 & cfg.fBOSC.F < freq(fr)+0.5,:),1))>0);
            
            % Calculate hit rate
            hits = detected_times(log_alpha);
            HR_thresh_alpha(fr,1,thresh,indTrial) = length(find(hits==1))/length(alpha_osc_array);
            
            % Calculate false alarm rate
            falses = detected_times(~log_alpha);
            FA_thresh_alpha(fr,1,thresh,indTrial) = length(find(falses==1))/(length(detected_times)-length(alpha_osc_array));
            
            clear hits falses detected_times detected
        end
    end
end
    

% Plot ROC curve
figure;
ddd = {'BOSC','eBOSC','fBOSC'};
for i = 1:3
    subplot(1,3,i);
    for fr = 1:2
        plot([squeeze(mean(FA_thresh_alpha(fr,i,:,:),4))]',...
            [squeeze(mean(HR_thresh_alpha(fr,i,:,:),4))']','-o'); hold on;
    end
    
    xlim([0 1]); ylim([0 1]);
    title(ddd{i});
    ylabel('Hit Rate','FontSize',20);
    xlabel('FA Rate','FontSize',20);
    if i == 3
        legend({'alpha','theta'},'Location','SouthEastOutside');
    end
    
end

% Plot ROC curve, but lines between the same thresholds
figure;
ddd = {'BOSC','eBOSC','fBOSC'};
for i = 1:3
    subplot(1,3,i);
    for j = 1:18
        plot([squeeze(mean(FA_thresh_alpha(:,i,j,:),4))]',...
            [squeeze(mean(HR_thresh_alpha(:,i,j,:),4))']','-ko'); hold on;
    end
    
    xlim([0 1]); ylim([0 1]);
    title(ddd{i});
    ylabel('Hit Rate','FontSize',20);
    xlabel('FA Rate','FontSize',20);
    if i == 3
        legend({'alpha','theta'},'Location','SouthEastOutside');
    end
    
end
  
% FA Rate
cols = [1.0000    0.3608    0.0157;
    0.7843    0.3451    0.5569;
    0.0157    0.3608    1.0000];

plot_style = {'--s','-o'};

figure;
set(gcf,'Position',[1 1 1600 400]);
ddd = {'BOSC','eBOSC','fBOSC'};
for i = 1:3
    subplot(1,3,i);
    for fr = flip(1:length(freq))
        plot(thresholds,squeeze(mean(FA_thresh_alpha(fr,i,:,:),4))',...
            plot_style{fr},'Color',cols(i,:),'LineWidth',2); hold on;
    end
    xlim([0 1]); ylim([0 1]);
    set(gca,'FontSize',20);
    ylabel('False Alarm Rate','FontSize',25);
    xlabel('Thresholds','FontSize',25);
    %title(ddd{i});

%     if i == 3
%         legend({'alpha','theta'},'Location','SouthEastOutside');
%     end
end
drawnow;
print('FA_diff','-dpng','-r400');


% Hit Rate
figure;
set(gcf,'Position',[1 1 1600 400]);
ddd = {'BOSC','eBOSC','fBOSC'};
for i = 1:3
    subplot(1,3,i);
    for fr = flip(1:length(freq))
        plot(thresholds,squeeze(mean(HR_thresh_alpha(fr,i,:,:),4))',...
            plot_style{fr},'Color',cols(i,:),'LineWidth',2); hold on;
    end
    xlim([0 1]); ylim([0.7 1]);
    set(gca,'FontSize',20);
    ylabel('Hit Rate','FontSize',25);
    xlabel('Thresholds','FontSize',25);
    %title(ddd{i});

%     if i == 3
%         legend({'alpha','theta'},'Location','SouthEastOutside');
%     end
end
drawnow;
print('HR_diff','-dpng','-r400');


%% Compare HR and FA rate per 
HR_diff = [];
for i = 1:3        
     diff = squeeze(mean(HR_thresh_alpha(1,i,:,:),4))-squeeze(mean(HR_thresh_alpha(2,i,:,:),4));
     HR_diff = [HR_diff diff];
end

figure; 
boxplot(abs(HR_diff));
title('HR Differences');

FA_diff = [];
for i = 1:3        
     diff = squeeze(mean(FA_thresh_alpha(1,i,:,:),4))-squeeze(mean(FA_thresh_alpha(2,i,:,:),4));
     FA_diff = [FA_diff diff];
end

figure; 
boxplot(abs(FA_diff));
title('FA Differences');

%% Export to csv, so we can plot in Python-land

condition                   = [];
[condition{1,1:size(FA_diff,1)}]        = deal('BOSC');
[condition{1,size(FA_diff,1)+1:size(FA_diff,1)*2}]     = deal('eBOSC');
[condition{1,(size(FA_diff,1)*2)+1:size(FA_diff,1)*3}]    = deal('fBOSC');
condition                   = condition';

FA_diff                      = abs(reshape(FA_diff,[size(FA_diff,1)...
                                *size(FA_diff,2),1]));
                            
HR_diff                      = abs(reshape(HR_diff,[size(HR_diff,1)...
                                *size(HR_diff,2),1]));

t = table(FA_diff,HR_diff,condition);
writetable(t,'diff_thresh.csv')                           
                            
                            
                            
                            
                            
%%                            
                            
% Calculate AUC
eee = [];

area_alpha = trapz([1 mean(FA_thresh_alpha,2)' 0],...
    [1 mean(HR_thresh_alpha,2)' 0]);
area_theta = trapz([1 mean(FA_thresh,2)' 0],...
    [1 mean(HR_thresh,2)' 0]);



% D Prime
d_prime = [];
for thresh = 1:length(thresholds)
    for fr = 1:2
        for i = 1:3
            HR = squeeze(mean(HR_thresh_alpha(fr,i,thresh,:),4))';
            FA = squeeze(mean(FA_thresh_alpha(fr,i,thresh,:),4))';
            
            d_prime(fr,i,thresh) = norminv(HR)-norminv(FA);
        end
    end
end

figure;
ddd = {'BOSC','eBOSC','fBOSC'};
for i = 1:3
    subplot(1,3,i);
    for fr = flip(1:length(freq))
        plot(thresholds,...
            squeeze(d_prime(fr,i,:)),'-o'); hold on;
    end
    %xlim([0 1]); ylim([0 1]);
    title(ddd{i});
    %ylabel('Hit Rate','FontSize',20);
    %xlabel('FA Rate','FontSize',20);
    if i == 3
        legend({'alpha','theta'},'Location','SouthEastOutside');
    end
end



%% Plot Hit-Rate vs SNR
h1= figure;
set(h1,'Position',[10 10 2000 500])

cols = [0.7843    0.3451    0.5569; 
    1.0000    0.3608    0.0157;0.0157    0.3608    1.0000];

for i = 1:3
%     subplot(1,3,i);
%     plot(SNR,squeeze(mean(hit_rate_alpha(i,1,:,:))),'--s',...
%         'Color',cols(i,:),'LineWidth',2,...
%         'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
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

cols = [0.7843    0.3451    0.5569; 
    1.0000    0.3608    0.0157;0.0157    0.3608    1.0000];

for i = 1:3
%     subplot(1,3,i);
%     plot(SNR,squeeze(mean(false_alarm_alpha(i,1,:,:))),'--s',...
%         'Color',cols(i,:),'LineWidth',2,...
%         'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    plot(SNR,squeeze(mean(false_alarm(i,1,:,:))),'-o',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    set(gca,'FontSize',20);
    xlabel('SNR','FontSize',24);
    ylabel('False Alarm Rate','FontSize',24);
    %xticks([0:4:24]);
    ylim([0 0.2]);
end

print('FA_SNR','-dpng','-r300');

%% d'
h1= figure;
set(h1,'Position',[10 10 2000 500])

cols = [0.7843    0.3451    0.5569; 
    1.0000    0.3608    0.0157;0.0157    0.3608    1.0000];

for i = 1:3
    %subplot(1,3,i);
    % Theta
    FA = squeeze(mean(false_alarm(i,1,:,:)));
    HR = squeeze(mean(hit_rate(i,1,:,:)));
    
    d_prime = (norminv(HR))-(norminv(FA));
    
    % Alpha
    FA_alpha = squeeze(mean(false_alarm_alpha(i,1,:,:)));
    HR_alpha = squeeze(mean(hit_rate_alpha(i,1,:,:)));
    
    d_prime_alpha = (norminv(HR_alpha))-(norminv(FA_alpha));
    
    
    plot(SNR,d_prime,'--s',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    plot(SNR,d_prime_alpha,'-o',...
        'Color',cols(i,:),'LineWidth',2,...
        'MarkerSize',10,'MarkerFaceColor',cols(i,:)); hold on;
    set(gca,'FontSize',20);
    xlabel('SNR','FontSize',24);
    ylabel('d prime','FontSize',24);
    xticks([0:4:24]);
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







        
        
