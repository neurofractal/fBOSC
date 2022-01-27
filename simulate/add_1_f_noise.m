% This script generates the scale-free background that is used for the
% simulations of standard BOSC and eBOSC.

% Note that the function 'f_alpha_gaussian' from the CNOISE toolbox 
% has to be added to the path. The toolbox is available at
% https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/cnoise.html

% For reproducibility, we use a fixed seed.


disp('Simulating linear 1/f');

for i = 1:200
    
    linear(i,:) = f_alpha_gaussian(10000,1,1)/24; 
end



TFR = zeros(length(cfg.fBOSC.F),size(linear,2),size(linear,1));
for indTrial = 1:size(linear,1)
    TFR(:,:,indTrial) = BOSC_tf(linear(indTrial,:),...
        cfg.fBOSC.F,cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
end; clear indTrial

% Get mean power
BG = reshape(TFR,[36,10000*200]);
mean_pow = 10.^(mean(log10(BG(:,:)),2));

figure; plot(log10(cfg.fBOSC.F),log10(mean_pow),'LineWidth',3); hold on;


%% Calculate RMSE for every trial
RMSE_trial   = [];
for indTrial = 1:ntrials
    % First BOSC
    trial_aperiodic_mean = mean(log10(squeeze(TFR(:,:,indTrial))),2);
    RMSE_trial(indTrial) = sqrt(immse(log10(mean_pow), trial_aperiodic_mean));
    p1 = plot(log10(cfg.fBOSC.F),trial_aperiodic_mean);
    p1.Color(4) = 0.1;
end

figure; boxplot(RMSE_trial);
save linear linear