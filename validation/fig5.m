%%
% Does modelling 'knee' vs 'fixed' reduce the model fit error for MEG and
% ECOG data betwen 1-40Hz?
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

error_fixed     = [];
error_knee      = [];

for sub = 1:length(subject)
    disp(['Subject: ' num2str(subject{sub})]);
    
    cd(fullfile(data_dir,subject{sub}));
    load('VE_for_HMM.mat');
    
    % Calculate power spectra with Welch's method
    [psds, freqs] = pwelch(VE_for_HMM.trial{1}', 500, [], [], 250);
    
%     figure;
%     semilogy(log10(freqs),psds);
%     xlim([0.2 1.6]);
    
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
    
    clear VE_for_HMM psds freqs
end


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














