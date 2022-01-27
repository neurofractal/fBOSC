% This script generates the scale-free background that is used for the
% simulations

% Note that the function 'f_alpha_gaussian' from the CNOISE toolbox 
% has to be added to the path. The toolbox is available at
% https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/cnoise.html

disp('Simulating linear 1/f');

linear = [];

for i = 1:200
    
    linear(i,:) = f_alpha_gaussian(10000,1,1); 
end


% Perform TFR
TFR = zeros(50,size(linear,2),size(linear,1));
for indTrial = 1:size(linear,1)
    TFR(:,:,indTrial) = BOSC_tf(linear(indTrial,:),...
        1:1:50,500,6);
end; clear indTrial

% Get mean power
BG = reshape(TFR,[50,10000*200]);
mean_pow = 10.^(mean(log10(BG(:,:)),2));

figure; plot(log10(1:1:50),log10(mean_pow),'LineWidth',3); hold on;

%% Calculate RMSE for every trial
RMSE_trial   = [];
for indTrial = 1:ntrials
    % First BOSC
    trial_aperiodic_mean = mean(log10(squeeze(TFR(:,:,indTrial))),2);
    RMSE_trial(indTrial) = sqrt(immse(log10(mean_pow), trial_aperiodic_mean));
    p1 = plot(log10(1:1:50),trial_aperiodic_mean);
    p1.Color(4) = 0.1;
end

figure; boxplot(RMSE_trial);
disp('Saving Linear 1/f');
save linear linear