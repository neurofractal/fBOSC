%%
% MEG Resting State (Theta)
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
root = fileparts(matlab.desktop.editor.getActiveFilename);
cd(fullfile(root,'..'));

start_fBOSC
cd(fullfile(root,'..','validation'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let's load some source localised MEG data.
% 
% This is from HCP MEG Resting-State Dataset
% See: https://www.humanconnectome.org/study/hcp-young-adult/project-protocol/resting-state-meg
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir    = '/Volumes/Robert T5/RS_HCP_VE';

subject = {'100307','102816','105923','106521','108323','109123',...
    '111514','112920','113922','116524','116726','133019','140117',...
    '146129','149741','153732','154532','158136',...
    '162026','162935','164636','166438','169040','172029','174841',...
    '175237','175540','177746','179245','181232',...
    '187547','189349','191033','191437','191841','192641','195041',...
    '198653','204521','205119','212318','212823','214524','221319',...
    '223929','233326','248339','250427','255639','257845'};

for sub = 8:20
    disp(['Subject: ' num2str(subject{sub})]);
    
    cd(fullfile(data_dir,subject{sub}));
    load('VE_for_HMM.mat');
    
%     %%%%%%%%%%
%     % fBOSC
%     %%%%%%%%%%
%     %% Set-up fBOSC parameters
%     
%     % general setup
%     cfg.fBOSC.F                 = 2.^[1:.125:5];
%     %cfg.fBOSC.F                 = [2:0.5:40];
%     cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
%     cfg.fBOSC.fsample           = 250;         % current sampling frequency of EEG data
%     
%     % padding
%     cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
%     cfg.fBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
%     cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)
%     
%     % fooof parameters - fit with fixed line or allow a knee
%     cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
%     cfg.fBOSC.fooof.version           = 'python';
%     
%     % threshold settings
%     cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
%     cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold
%     
%     % episode post-processing
%     cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now
%     
%     % general processing settings
%     cfg.fBOSC.channel           = []; % select posterior channels (default: all)
%     cfg.fBOSC.trial             = []; % select trials (default: all)
%     cfg.fBOSC.trial_background  = []; % select trials for background (default: all)
%     
%     %% Run fBOSC
%     clear fBOSC
%     [fBOSC, cfg] = fBOSC_wrapper(cfg, VE_for_HMM);
%     
%     %% Get mean SNR of theta bursts
%     SNR_mean = [];
%     for chan = 1:44
%         T = fBOSC.episodes(fBOSC.episodes.Channel == chan,:);
%         T = T(T.FrequencyMean >=3 & T.FrequencyMean < 8,:);
%         SNR_mean(chan) = mean(T.SNRMean);
%         clear T
%     end
%     
%     save SNR_mean SNR_mean
%     
%     %% Get abundance in theta
%     freqs       = 2.^[1:.125:5];
%     theta_freqs = find(freqs < 8 & freqs > 3);
% 
%     theta_abundance = mean(squeeze(fBOSC.abundance_ep(:,:,theta_freqs)),2);
%     
%     save theta_abundance theta_abundance
%     
%     %% Clear variables
%     clear theta_abundance SNR_mean T fBOSC cfg
%     
    %%%%%%%%%%%%%%%%%
    % Now BOSC
    %%%%%%%%%%%%%%%%%
    %% Set-up BOSC parameters
    
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
    cfg.fBOSC.fooof.aperiodic_mode    = 'old';
    cfg.fBOSC.fooof.version           = 'python';
    
    % threshold settings
    cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold
    
    % episode post-processing
    cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now
    
    % general processing settings
    cfg.fBOSC.channel           = []; % select posterior channels (default: all)
    cfg.fBOSC.trial             = []; % select trials (default: all)
    cfg.fBOSC.trial_background  = []; % select trials for background (default: all)
    
    %% Run fBOSC
    
    clear fBOSC
    [fBOSC, cfg] = fBOSC_wrapper(cfg, VE_for_HMM);
    
    
    
%     %% Plot the Results of the 1/f fit
%     cfg.log_freqs = 1;
%     cfg.plot_old = 0;
%     fBOSC_fooof_plot(cfg,fBOSC)
    
    %% Get mean SNR of theta bursts
    SNR_mean = [];
    for chan = 1:44
        T = fBOSC.episodes(fBOSC.episodes.Channel == chan,:);
        T = T(T.FrequencyMean >=3 & T.FrequencyMean < 8,:);
        SNR_mean(chan) = mean(T.SNRMean);
        clear T
    end
    
    save SNR_mean_BOSC SNR_mean
    
    %% Get abundance in theta
    freqs       = 2.^[1:.125:5];
    theta_freqs = find(freqs < 8 & freqs > 3);

    theta_abundance = mean(squeeze(fBOSC.abundance_ep(:,:,theta_freqs)),2);
    
    save theta_abundance_BOSC theta_abundance
    
    %% Clear variables
    clear theta_abundance SNR_mean T fBOSC cfg
    
end


%% Load all data
subject = {'100307','102816','105923','106521','108323','109123',...
    '111514','112920','113922','116524','116726','133019','140117',...
    '146129','149741','153732','154532','158136',...
    '162026','162935'};

SNR_mean_all = [];
theta_abundance_all = [];

for sub=1:length(subject)
    disp(['Subject: ' num2str(subject{sub})]);
    
    cd(fullfile(data_dir,subject{sub}));
    load('SNR_mean_BOSC.mat');
    SNR_mean_all(:,sub) = SNR_mean;
    load('theta_abundance_BOSC.mat');
    theta_abundance_all(:,sub) = theta_abundance;
    clear SNR_mean theta_abundance
end

SNR_mean_all            = mean(SNR_mean_all,2);
theta_abundance_all     = mean(theta_abundance_all,2);

%% Load Template Grid
% Get the path and version of Fieldtrip
[t, r] = ft_version;

load(fullfile(r,'template','sourcemodel','standard_sourcemodel3d8mm.mat'));
template_grid   = sourcemodel;
clear sourcemodel
template_grid     = ft_convert_units(template_grid,'mm');

%%
HMM_dir       = '/Users/rseymoue/Documents/GitHub/HMM';

% Load hippocampal atlas
atlas_hippocampus = ft_read_atlas(fullfile(HMM_dir,'Melb_Subcortical_Atlas',...
    'melb_sub.mat'));

% Load hippocampal atlas
atlas_subcortical = ft_read_atlas(fullfile(HMM_dir,'Melb_Subcortical_Atlas',...
    'melb_sub.mat'));

% Load atlas_HCPMMP
load(fullfile(HMM_dir,'atlas_HCPMMP.mat'));

% and call ft_sourceinterpolate:
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'tissue';
sourcemodel2                = ft_sourceinterpolate(cfg, atlas_HCPMMP,template_grid);
cfg.parameter               = 'tissue';
sourcemodel_subcortical     = ft_sourceinterpolate(cfg, atlas_hippocampus,template_grid);

sourceavg = template_grid;

% Load the SPM Brain
spm_brain = ft_read_mri(fullfile(r,'template/anatomy/single_subj_T1.nii'));

% in order to do this we segment the brain from the anatomical template
spm_brain.coordsys = 'mni';
cfg         = [];
cfg.output  = {'brain'};
spm_brain_seg = ft_volumesegment(cfg,spm_brain);

% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

% Create empty power array
pow  = nan(size(template_grid.inside,1),1);
pow(template_grid.inside) = 0;

% Find the atlas labels used in our (cortical) parcellation
[a b] = ismember(VE_for_HMM.label,sourcemodel2.tissuelabel);
atlas_labels_to_use = b(a);

% For each region make power = theta abundance
for i=1:length(atlas_labels_to_use)
    disp([num2str(atlas_labels_to_use(i))...
        '. ' sourcemodel2.tissuelabel{atlas_labels_to_use(i)}]);
    
    d = find(sourcemodel2.tissue==atlas_labels_to_use(i));
    pow(d) = theta_abundance_all(i);
end

% For hippocampus - use hardcoded
atlas_labels_to_use = [1 9];

for i=1:2
    disp(sourcemodel_subcortical.tissuelabel{atlas_labels_to_use(i)});
    
    d = find(sourcemodel_subcortical.tissue==atlas_labels_to_use(i));
    pow(d) = theta_abundance_all(i);
end

%%
sourceavg.pow = pow;

% and call ft_sourceinterpolate:
cfg                 = [];
cfg.interpmethod    = 'nearest';
cfg.parameter       = 'pow';
sourceI             = ft_sourceinterpolate(cfg, sourceavg, spm_brain);
sourceI.coordsys    = 'mni';

% Mask bits outside the brain
sourceI.anat_mask = spm_brain_seg.brain .* double(sourceI.anatomy);

cfg                 = [];
cfg.funparameter    = 'pow';
cfg.funcolormap        = cmap;
%cfg.funcolorlim     = 'maxabs';
cfg.maskparameter   = 'anat_mask';
cfg.atlas           = atlas_HCPMMP;
ft_sourceplot(cfg,sourceI);

%% Export to nifti formt and use your favourite MRI software to visualise
cd('/Users/rseymoue/Documents/');
cfg = [];
cfg.filetype = 'nifti';
cfg.filename = ['theta_abundance_BOSC'];
cfg.parameter = 'pow';
ft_sourcewrite(cfg,sourceI);

%% Export to connectome workbench (specfic to my computer)
    try
        system(['/Applications/workbench/bin_macosx64/wb_command -volume-to-surface-mapping /Users/rseymoue/Documents/' ...
            'theta_abundance_BOSC.nii /Applications/workbench/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii /Users/rseymoue/Documents/' 'theta_abundance_BOSC_LEFT.shape.gii -trilinear'])
                system(['/Applications/workbench/bin_macosx64/wb_command -volume-to-surface-mapping /Users/rseymoue/Documents/' ...
            'theta_abundance_BOSC.nii /Applications/workbench/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.surf.gii /Users/rseymoue/Documents/' 'theta_abundance_BOSC_RIGHT.shape.gii -trilinear'])
    catch
        disp('Could not convert to gifti format');
    end

saveas(gcf, ['state' num2str(k) '.png'],'png');
drawnow;

%%

%% Load all data
subject = {'100307','102816','105923','106521','108323','109123',...
    '111514','112920','113922','116524','116726','133019','140117',...
    '146129','149741','153732','154532','158136',...
    '162026','162935'};

SNR_mean_all = [];
theta_abundance_all = [];
SNR_mean_all_BOSC = [];
theta_abundance_all_BOSC = [];

for sub=1:length(subject)
    disp(['Subject: ' num2str(subject{sub})]);
    
    cd(fullfile(data_dir,subject{sub}));
    % fBOSC
    load('SNR_mean.mat');
    SNR_mean_all(:,sub) = SNR_mean;
    load('theta_abundance.mat');
    theta_abundance_all(:,sub) = max(theta_abundance);
    clear SNR_mean theta_abundance
    
    % BOSC
    load('SNR_mean_BOSC.mat');
    SNR_mean_all_BOSC(:,sub) = SNR_mean;
    load('theta_abundance_BOSC.mat');
    theta_abundance_all_BOSC(:,sub) = max(theta_abundance);
    clear SNR_mean theta_abundance
end


figure; boxplot([theta_abundance_all',theta_abundance_all_BOSC']);


figure;boxplot([error_fixed error_knee],{'fixed','knee'});

error_both                  = vertcat(error_fixed, error_knee);
condition                   = vertcat(repmat({'fixed'},length(error_fixed),1),...
    repmat({'knee'},length(error_knee),1));

t = table(error_both,condition);
cd(fullfile(root,'..','validation'));
writetable(t,'MEG_source_localised.csv');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let's load some Resting State ECOG Data from 10 patients
%
% This dataset is fixation_PAC.zip, which can be downloaded from:
% https://searchworks.stanford.edu/view/zk881ps0522
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear error_fixed error_knee

data_dir = '/Volumes/Robert T5/fixation_PAC';

subject_list = {'bp','cc','hl','jc','jm','jp','ug','wc','wm','zt'};
Fs       = 1000;

error_fixed     = [];
error_knee      = [];

for sub = 1:length(subject_list)
    disp(subject_list{sub});
    cd(fullfile(data_dir,'data',subject_list{sub}));
    load([subject_list{sub} '_base.mat']);
    
    % Calculate power spectra with Welch's method
    [psds, freqs] = pwelch(data, 500, [], [], Fs);
    
    figure;
    title(subject_list{sub});
    semilogy(log10(freqs),psds);
    xlim([0.1 2.5]);
    xlabel('Log Frequency (Hz','FontSize',22);
    ylabel('Log Power (au.)','FontSize',22);
    set(gca,'FontSize',20)
    
    % Transpose, to make inputs row vectors
    freqs = freqs';
    
    % FOOOF settings
    settings = struct();
    settings.aperiodic_mode = 'fixed';
    settings.verbose = 0;
    f_range = [1, 40];
    
    % Run FOOOF across a group of power spectra
    fooof_results = fooof_group(freqs, psds, f_range, settings);
    
    error_fixed = vertcat(error_fixed,[fooof_results.error]');
    
    % Now repeat with 'knee'
    clear fooof_results
    settings.aperiodic_mode = 'knee';
    
    % Run FOOOF across a group of power spectra
    fooof_results = fooof_group(freqs, psds, f_range, settings);
    
    error_knee = vertcat(error_knee,[fooof_results.error]');
    
    clear data psds freqs
end

figure;boxplot([error_fixed error_knee],{'fixed','knee'});

% Export csv file
cd(fullfile(root,'..','validation'));

error_both                  = vertcat(error_fixed, error_knee)
condition                   = vertcat(repmat({'fixed'},length(error_fixed),1),...
    repmat({'knee'},length(error_fixed),1));

t = table(error_both,condition);
writetable(t,'ECOG.csv')














