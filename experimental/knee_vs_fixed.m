%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some code to test the prevalence of linear
% versus non-linear aperiodic in MEG data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_fBOSC
cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let's load some source localised MEG data
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
    disp(subject{sub});
    
    cd(fullfile(data_dir,subject{sub}));
    load('VE_for_HMM.mat');
    
    % Calculate power spectra with Welch's method
    [psds, freqs] = pwelch(VE_for_HMM.trial{1}', 500, [], [], 250);
    
%     figure;
%     semilogy(freqs,psds);
    
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

cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

figure;boxplot([error_fixed error_knee],{'fixed','knee'});

error_both                  = vertcat(error_fixed, error_knee)
condition                   = vertcat(repmat({'fixed'},2200,1),...
    repmat({'knee'},2200,1));

t = table(error_both,condition);
writetable(t,'MEG_source_localised.csv')



%% Now let's run on some source-localised OPM-MEG data
data_path     = '/Users/rseymoue';
cd(data_path);
load('VE.mat');
VE.fsample = 600;

error_fixed     = [];
error_knee      = [];

% Calculate power spectra with Welch's method
[psds, freqs] = pwelch(VE.trial{1}', 500, [], [], 600);

figure;
semilogy(freqs,psds);

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


cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

figure;boxplot([error_fixed error_knee],{'fixed','knee'});

error_both                  = vertcat(error_fixed, error_knee)
condition                   = vertcat(repmat({'fixed'},44,1),...
    repmat({'knee'},44,1));

t = table(error_both,condition);
writetable(t,'OPMMEG.csv')


%% Let's try it out on some ECOG data
data_dir = '/Users/rseymoue/Downloads/fixation_PAC'

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
    semilogy(log10(freqs),psds);
    xlim([1 40]);
    
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

error_knee = vertcat(error_knee,[fooof_results.error]');


cd('/Users/rseymoue/Documents/GitHub/fBOSC/experimental');

figure;boxplot([error_fixed error_knee],{'fixed','knee'});

error_both                  = vertcat(error_fixed, error_knee)
condition                   = vertcat(repmat({'fixed'},535,1),...
    repmat({'knee'},535,1));

t = table(error_both,condition);
writetable(t,'ECOG.csv')














