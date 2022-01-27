%% Paths (RS)
fieldtripDir    = '/Users/rseymoue/Documents/scripts/fieldtrip-20191213';
script_dir      = '/Users/rseymoue/Documents/GitHub/analyse_OPMEG';
denoise_dir     = '/Users/rseymoue/Documents/scripts/NoiseTools';
scannercast_dir = '/Users/rseymoue/Documents/GitHub/scannercast/examples/NA';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));

% Add NoiseTools to path
addpath(denoise_dir)

%% Load the data
data_path     = '/Users/rseymoue';
cd(data_path);
load('VE.mat');
VE.fsample = 600;

%% Try fooof in Fieldtrip
cfg               = [];
cfg.channel       = VE.label{30};
cfg.foilim        = [1 40];
cfg.pad           = 10;
cfg.tapsmofrq     = 2;
cfg.fooof.aperiodic_mode = 'knee';
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, VE);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, VE);

% subtract the fractal component from the power spectrum
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'log10(x2)-log10(x1)';
oscillatory = ft_math(cfg, fractal, original);

% display the spectra on a log-log scale
figure();
hold on;
plot(log10(original.freq), log10(original.powspctrm),'k');
plot(log10(fractal.freq), log10(fractal.powspctrm));
figure;
plot(log10(fractal.freq), oscillatory.powspctrm);
xlabel('log-freq'); ylabel('log-power');
legend({'original','fractal','oscillatory'},'location','southwest');


%% Load fBOSC
cd('/Users/rseymoue/Documents/GitHub/fBOSC');
start_fBOSC;

%% Reformat the VE to be one *long* continuous segment of data!
ddd = cat(2,VE.trial{:});
eee = cat(2,VE.time{:});


VE_combined                 = [];
VE_combined.fsample         = VE.fsample;
VE_combined.time{1}         = [0:(1/VE.fsample):(length(eee)/VE.fsample)];
VE_combined.time{1}(end)    = [];
VE_combined.trial{1}        = ddd;
VE_combined.label           = VE.label;

%% fBOSC parameters

% general setup
cfg.fBOSC.F             = [3:1:40];    % frequency sampling
cfg.fBOSC.wavenumber	= 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample       = 600;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

% fooof
cfg.fBOSC.fooof.aperiodic_mode    = 'knee';

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .99;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.fBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.fBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.fBOSC.postproc.edgeOnly = 'no';        % Deconvolution only at on- and offsets of fBOSC.episodes? (default = 'yes')
cfg.fBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.fBOSC.channel = [6 8 28 30]; % select posterior channels (default: all)
cfg.fBOSC.trial = []; % select trials (default: all)
cfg.fBOSC.trial_background = []; % select trials for background (default: all)

%% run fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, VE_combined);

%%
eee2 = eee(cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);    

%%
cfg.log_freqs = 1;
cfg.plot_old = 0;
fBOSC_fooof_plot(cfg,fBOSC)

%%
rty = [];
for t = 1:4

tmpDetected = single(squeeze(nanmean(fBOSC.detected(t, ...
        1,cfg.fBOSC.F > 8 & cfg.fBOSC.F < 13,:),3))>0);

   
[indx] = find(eee2==0);    



for i = 1:length(indx)
    ggg(i,:) = tmpDetected(indx(i)-900:indx(i)+3000);
end

rty(t,:) = mean(ggg,1);
end


figure;plot([-1.5:(1/600):5],rty(1:2,:),'LineWidth',2);
ylabel('Probability of mu burst','FontSize',14);
xlabel('Time (s)','FontSize',14);
legend({'L Motor';'L Premotor'});
set(gca,'FontSize',12);

sss = rty./(rty(:,1:900));
    



%%

indChan = 1; indTrial = 1; % Here we select the first trial and first channel we encoded (see cfg.fBSOSC.channel).

disp(['Results are for trial ', num2str(cfg.fBSOSC.trial(indTrial)),...
    ' at channel ', VE.label{cfg.fBSOSC.channel(indChan)}])

%% get original time series for plotting
origData = VE.trial{indTrial}(cfg.fBSOSC.channel(indChan),...
    cfg.fBSOSC.pad.total_sample+1:end-cfg.fBSOSC.pad.total_sample);

origData_time = VE.time{indTrial}(...
    cfg.fBSOSC.pad.total_sample+1:end-cfg.fBSOSC.pad.total_sample);

%%
% Supplementary Figure: plot estimated background + power threshold
figure; hold on;
plot((cfg.fBSOSC.F), log10(eBOSC.static.mp(indChan,:)), 'k--','LineWidth', 1.5); 
plot((cfg.fBSOSC.F), log10(eBOSC.static.pt(indChan,:)), 'k-', 'LineWidth', 1.5)
plot((cfg.fBSOSC.F), eBOSC.static.bg_log10_pow(indChan,:), 'r-', 'LineWidth', 2)

% plot(log10(cfg.fBSOSC.F), log10(eBOSC.static.mp(indChan,:)), 'k--','LineWidth', 1.5); 
% plot(log10(cfg.fBSOSC.F), log10(eBOSC.static.pt(indChan,:)), 'k-', 'LineWidth', 1.5)
% plot(log10(cfg.fBSOSC.F), eBOSC.static.bg_log10_pow(indChan,:), 'r-', 'LineWidth', 2)
set(gca,'fontsize',14)
xlabel('Frequency (log10 Hz)','FontSize',20); 
ylabel('Power (log 10 a.u.)','FontSize',20);
legend({'Aperiodic fit', 'Statistical power threshold', 'Avg. spectrum'}, ...
    'orientation', 'vertical', 'location', 'northoutside'); legend('boxoff');
%xlim([.3, 1.75])

%
figure();
hold on;

for i = 1:10
    
    origData = VE.trial{i}(cfg.fBSOSC.channel(indChan),...
        cfg.fBSOSC.pad.total_sample+1:end-cfg.fBSOSC.pad.total_sample);
    
    origData_time = VE.time{i}(...
        cfg.fBSOSC.pad.total_sample+1:end-cfg.fBSOSC.pad.total_sample);
    subplot(10,1,i);
    plot(origData_time,squeeze(origData), 'k');
    tmpDetected = single(squeeze(nanmean(eBOSC.detected(indChan, ...
        i,cfg.fBSOSC.F > 8 & cfg.fBSOSC.F < 13,:),3))>0); ...
        tmpDetected(tmpDetected==0) = NaN; hold on;
    plot(origData_time,squeeze(origData).*tmpDetected', 'r');
    %xlim([12 15])
    xlabel('Time (s)'); %ylabel('Power (a.u)');
    % [~, hobj, ~, ~] = legend({'Original signal'; 'Rhythmic signal'}, ...
    %     'orientation', 'horizontal', 'location', 'northoutside'); legend('boxoff')
    % hl = findobj(hobj,'type','line');
    % ht = findobj(hobj,'type','text')
    % set(ht,'FontSize',20);
    % set(hl,'LineWidth',3);
    %set(findall(gcf,'-property','FontSize'),'FontSize',14)
end

%%
tmpDetected = single(squeeze(nanmean(eBOSC.detected(indChan, ...
        i,cfg.fBSOSC.F > 8 & cfg.fBSOSC.F < 13,:),3))>0); ...
        tmpDetected(tmpDetected==0) = NaN; hold on;
    
ttt = single(squeeze(nanmean(eBOSC.detected(indChan, ...
        :,cfg.fBSOSC.F > 8 & cfg.fBSOSC.F < 13,:),3))>0); ...
        tmpDetected(tmpDetected==0) = NaN; hold on;
    
origData_time = VE.time{1}(...
        cfg.fBSOSC.pad.total_sample+1:end-cfg.fBSOSC.pad.total_sample);    
    
figure;plot(origData_time,mean(ttt,1));
ylabel({'Beta Burst Probability'; '(0-1)'},'FontSize',26);
xlabel('Time (s)','FontSize',26);
ax = gca;
ax.XAxis.FontSize = 18 %for y-axis 
ax.YAxis.FontSize = 18 %for y-axis 


%%
figure; 
subplot(3,2,1); histogram(eBOSC.episodes.SNRMean); title('SNR distribution')
subplot(3,2,2); histogram(log10(eBOSC.episodes.SNRMean)); title('SNR distribution(log10)')
subplot(3,2,3); histogram(eBOSC.episodes.DurationC); title('Duration distribution')
subplot(3,2,4); histogram(log10(eBOSC.episodes.DurationC)); title('Duration distribution(log10)')
subplot(3,2,5); histogram(eBOSC.episodes.FrequencyMean); title('Frequency distribution')
subplot(3,2,6); hold on; plot(squeeze(eBOSC.pepisode(indChan, indTrial,:))); plot(squeeze(eBOSC.abundance_ep(indChan, indTrial,:))); title('Pepisode, abundance')

%%
idx_alpha = find(eBOSC.episodes.Trial == indTrial & eBOSC.episodes.Channel == cfg.fBSOSC.channel(indChan) &...
    eBOSC.episodes.FrequencyMean > 8 & eBOSC.episodes.FrequencyMean <13);
% filter for alpha by mean frequency (!) of episode
idx_onset = []; idx_onsetTime = [];
alphaDetected = NaN(1,numel(origData));
for indEp = 1:numel(idx_alpha)
    % These are two alternative ways to extract the onset timepoint from the table
    idx_onsetTime(indEp) = find(cfg.tmp.finalTime>= eBOSC.episodes.Onset(idx_alpha(indEp)), 1, 'first');
    idx_onset(indEp) = eBOSC.episodes.ColID{idx_alpha(indEp)}(1);
    % Mark all periods with episodes falling into the alpha range
    alphaDetected(eBOSC.episodes.ColID{idx_alpha(indEp)}(1):eBOSC.episodes.ColID{idx_alpha(indEp)}(end)) = 1;
end; clear idx_alpha;
h = figure('units','normalized','position',[.1 .1 .7 .3]); hold on; 
scatter(idx_onset, repmat(0.3,1,numel(idx_onset)), 75, [.5 .5 .5], 'filled')
OnsetLine = squeeze(origData);
OnsetLine(idx_onset) = 0.3; clear idx_onset idx_onsetTime;
plot(OnsetLine, 'Color', [.5 .5 .5]); clear OnsetLine;
[orig]=plot(origData_time,squeeze(origData), 'k');
[rhythm]=plot(origData_time,squeeze(origData).*alphaDetected, 'r');
xlim([0 20])
xlabel('Time (s)'); ylabel('Power [µV]');
legend([orig, rhythm], {'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'northoutside'); legend('boxoff')
set(findall(gcf,'-property','FontSize'),'FontSize',26)




%%















