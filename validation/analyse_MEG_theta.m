%%
% Analyse MEG theta
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
root = fileparts(matlab.desktop.editor.getActiveFilename);
cd(fullfile(root,'..'));

start_fBOSC
cd(fullfile(root,'..','validation'));

% %% SPM
% % addpath /path/to/spm
% addpath('/Users/rseymoue/Documents/scripts/spm12');
% 
% spm('defaults','eeg')

%% Data directory
data_dir = '/Volumes/Robert T5/Navigation MEG Data/Controls/';
cd(data_dir);

%% Load all the Behavioural Data
load('allBehaviouralData.mat');

ft_hastoolbox('spm8',1);
ft_hastoolbox('spm12',1);

%% Loop Over Subjects and Runs

for sub = 1:length(ctrNos)
    for run = 1:3
        try
            D = spm_eeg_load(fullfile(data_dir,['MdpSubject' num2str(ctrNos(sub)) ...
                '_Run' num2str(run) '.mat']));
        catch
            ft_warning(['Could not load data from subject ' num2str(ctrNos(sub))...
                ' Run: ' num2str(run)]);
        end
        data = spm2fieldtrip(D);
        
        % Get fiducial information
        fids = D.fiducials;
        fids.fid.pos = fids.fid.pnt;
        
        % Set sampling freq
        data.fsample = D.fsample;
        data.hdr.Fs = D.fsample;
        data.hdr.label = data.label;
        clear D
        
               
%         cfg          = [];
%         cfg.viewmode = 'vertical';
%         cfg.channel  = 'UADC004';
%         ft_databrowser(cfg,data);
        
        %% Segment based on trigger
        cfg                         = [];
        cfg.rawData                 = data;
        cfg.trialdef.trigchan       = 'UADC004';
        cfg.correct_time            = [];
        cfg.trialdef.prestim        = 4.0;        % pre-stimulus interval
        cfg.trialdef.poststim       = 4.0;        % post-stimulus interval
        cfg.trialfun                = 'OPM_trialfun_usemat';
        banana                      = ft_definetrial(cfg);
        
        time_of_trig = (banana.trl(:,2)/data.hdr.Fs)-4.0;
        
        % Select the trigger closest to timing from event file
        time_of_event = ctrRespData{sub,run}(9:end,1);
        
        trl2 = [];
        
        for i = 1:length(time_of_event)
            [c index] = min(abs(time_of_trig-time_of_event(i)));
            
            % Warn if greater than 0.2 difference
            if time_of_trig(index)-time_of_event(i) > 0.2
                ft_warning(['Trial ' num2str(i) ' = ' num2str(time_of_trig(index)-time_of_event(i)) 's difference']);
                
            else
                % Adjust by 3s (onset of the cue)
                banana.trl(index,1) = banana.trl(index,1)+(data.hdr.Fs*3);
                banana.trl(index,2) = banana.trl(index,2)+(data.hdr.Fs*3);
                
                % Add to trl2
                trl2(i,:) = banana.trl(index,:);
            end
        end
        
        % Replace the trl
        banana2     = banana;
        banana2.trl = trl2;
        
        %% Filter
        % HP Filter
        disp('Filtering...');
        cfg                 = [];
        cfg.hpfilter        = 'yes';
        cfg.hpfreq          = 1;
        data                = ft_preprocessing(cfg,data);
        
        % LP Filter
        cfg                 = [];
        cfg.lpfilter        = 'yes';
        cfg.lpfreq          = 40;
        data                = ft_preprocessing(cfg,data);
        
        %% Epoch
        % Redefines the data
        data_clean  = ft_redefinetrial(banana2,data);
        
        %% Get only MEG data
        cfg             = [];
        cfg.channel     = 'MEG';
        data_clean      = ft_selectdata(cfg,data_clean);
        
        % cfg          = [];
        % cfg.viewmode = 'vertical';
        % cfg.channel  = 'MEG';
        % ft_databrowser(cfg,data_clean);
        
        
        %% Transform the template sourcemodel
        [ftver, ftpath] = ft_version;
        load(fullfile(ftpath,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
        sourcemodel_headshape = sourcemodel; clear sourcemodel;
        sourcemodel_headshape = ft_convert_units(sourcemodel_headshape,'mm');
        
        load(fullfile(ftpath,'template','headmodel','standard_singleshell.mat'));
        headmodel  = vol;
        headmodel  = ft_convert_units(headmodel,'mm');
        
        cfg         = [];
        cfg.method 	= 'fids';
        cfg.verbose = 'yes';
        [sourcemodel_warped, M1] = param_12_affine(cfg,fids,sourcemodel_headshape);
        hold on;
        
        % Transform the headmodel
        headmodel.bnd.pos = ft_warp_apply(M1,headmodel.bnd.pos);
        
        %% Source analysis
        % Now let's process the data
        disp('Computing covariance');
        cfg                     = [];
        cfg.channel             = 'MEG';
        cfg.covariance          = 'yes';
        cfg.vartrllength        = 2;
        cfg.keeptrials          = 'no';
        cfg.covariancewindow    = 'all';
        avg_bpf                 = ft_timelockanalysis(cfg,data_clean);
        
        figure; imagesc(avg_bpf.cov);
        
        figure; hold on;
        ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
        alpha 0.5; camlight;
        ft_plot_mesh(sourcemodel_warped.pos(sourcemodel_warped.inside,:));
        ft_plot_sens(data_clean.grad, 'style', 'r*');
        view([0,0]);
        
        %% Do beamforming
        disp('Beamforming...');
        cfg                 = [];
        cfg.channel         = 'MEG';
        cfg.method          = 'lcmv';
        cfg.grid            = sourcemodel_warped;
        cfg.headmodel       = headmodel;
        cfg.keeptrials      = 'no';
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.fixedori   = 'yes';
        cfg.lcmv.lambda     = '5%';
        sourceavg           = ft_sourceanalysis(cfg, avg_bpf);
        
        sourceavg           = rmfield(sourceavg,'cfg');
        
        %% Load atlas
        load('/Users/rseymoue/Documents/GitHub/HMM/atlas_HCPMMP.mat');
        
        % Get the path and version of Fieldtrip
        [~, r] = ft_version;
        
        % Load 8mm sourcemodel
        load(fullfile(r,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
        template_grid   = sourcemodel;
        clear sourcemodel
        template_grid     = ft_convert_units(template_grid,'mm');
        
        % Create VE
        [VE] = atlas2VE(atlas_HCPMMP,template_grid,...
            atlas_HCPMMP.tissuelabel, data_clean,...
            sourceavg, avg_bpf);
        
        cfg             = [];
        cfg.channel     = {'RH_Anterior_Cingulate_and_Medial_Prefrontal_Cortex',...
            'LH_Anterior_Cingulate_and_Medial_Prefrontal_Cortex'};
        VE_frontal      = ft_selectdata(cfg,VE);
        
        save(['VE_frontal_Subject' num2str(ctrNos(sub)) '_Run' num2str(run) '.mat'],...
            'VE_frontal');
        
        % cfg          = [];
        % %cfg.viewmode = 'vertical';
        % cfg.channel  = 'all';
        % ft_databrowser(cfg,VE_frontal);
        
        %% Calculate TFRs using a hanning taper
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.pad          = 'nextpow2';
        cfg.foi          = 1:1:40;
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
        cfg.toi          = -4:0.05:4;
        TFR              = ft_freqanalysis(cfg, VE_frontal);
        
        for chan = 1:length(VE_frontal.label)
            cfg                 = [];
            cfg.channel         = VE_frontal.label{chan};
            %cfg.colormap        = colormap123;
            cfg.ylim            = [1 40];
            %cfg.xlim            = [-3 3];
            cfg.baseline        = [-3 0];
            cfg.baselinetype    = 'db';
            cfg.zlim            = 'maxabs';
            cfg.comment         = 'no';
            figure;ft_singleplotTFR(cfg,TFR); hold on;
            title(['VE Frontal Subject' num2str(ctrNos(sub)) ' Run' num2str(run)]);
            drawnow;
        end
        
    end
end




%% Let's try and fBOSC this
fff = [];
for run = 1:3
    fff{run} = load(['VE_frontal_Subject15_Run' num2str(run) '.mat']);
end

cfg                 = [];
cfg.keepsampleinfo  = 'no';
data_comb           = ft_appenddata(cfg,fff{1}.VE_frontal,fff{2}.VE_frontal,...
                        fff{3}.VE_frontal);


%% Set-up fBOSC parameters

% general setup
cfg.fBOSC.F                 = 2.^[1:.125:5];
%cfg.fBOSC.F                 = [2:0.5:40];
cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample           = 250;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

% fooof parameters - fit with fixed line or allow a knee
cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
cfg.fBOSC.fooof.version           = 'python';

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold

% episode post-processing
cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now

% general processing settings
cfg.fBOSC.channel           = [1]; % select posterior channels (default: all)
cfg.fBOSC.trial             = []; % select trials (default: all)
cfg.fBOSC.trial_background  = []; % select trials for background (default: all)

%% Run fBOSC
clear fBOSC
[fBOSC, cfg] = fBOSC_wrapper(cfg, data_comb);

%% Plot the Results of the 1/f fit
cfg.log_freqs = 1;
cfg.plot_old = 0;
fBOSC_fooof_plot(cfg,fBOSC)


%% Plot the time-series in black alongside detected bursts in red

figure;
set(gcf,'Position',[100 100 1200 600]);

for indTrial = 1:10
    
    % This is the crucial bit of code that gets the detected time(s) of
    % oscillatory bursts from each trial
    tmpDetected = single(squeeze(nanmean(fBOSC.detected(1, ...
        indTrial,cfg.fBOSC.F > 3 & cfg.fBOSC.F < 8,:),3))>0); ...
        tmpDetected(tmpDetected==0) = NaN;
    
    origData = data_comb.trial{indTrial}(cfg.fBOSC.channel(1),...
        cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
    
    origData_time = data_comb.time{indTrial}(...
        cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
    
    % Plot
    subplot(10,1,indTrial); hold on;
    plot(origData_time,squeeze(origData), 'k');
    plot(origData_time,squeeze(origData).*tmpDetected', 'r');
    ylabel({'Power';'(a.u.)'});
    set(get(gca,'YLabel'),'Rotation',45)
    if indTrial == length(data_comb.trial)
        xlabel('Time (s)');
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    
end


%% Let's plot the same thing but using TFRs

% Compute TFRs
TFR = [];
for indTrial = 1:length(data_comb.trial)
    TFR.trial{indTrial} = BOSC_tf(data_comb.trial{indTrial}(1,:),...
        cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
end; clear indTrial


% Plot the figure
figure;
set(gcf,'Position',[100 100 1200 900]);

for indTrial = 1:10
    subplot(10,1,indTrial)
    TFR_ = TFR.trial{indTrial}; % Account for padding
    
    % This is the crucial bit of code 
    tmpDetected = logical(squeeze(fBOSC.detected(1,indTrial,:,:)));

    % Plot
    imagesc('YData',log10(cfg.fBOSC.F),'XData',data_comb.time{1},'CData',...
        TFR_,'AlphaData', tmpDetected); hold on;
    imagesc('YData',log10(cfg.fBOSC.F),'XData',data_comb.time{1},'CData',...
        TFR_,'AlphaData', 0.2);
    set(gca,'YDir','normal');
    xlim([-3 3]);
    %ylim([0.4 1.4]);
    set(gca,'FontSize',20);
    ylabel({'Power';'(a.u.)'},'FontSize',10);
    set(get(gca,'YLabel'),'Rotation',45);

    if indTrial==length(data_comb.trial)
        xlabel('Time (s)');
    end
    
end







