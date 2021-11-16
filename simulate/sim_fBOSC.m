function [aperiodic_out, osc_out, combined,AlphaPlace] = sim_fBOSC(cfg,Fs,aperiodic)
% Function to combine simulated aperiodic data with a simulated
% oscillatory burst of known duration and amplitude.
%
%   INPUTS:
%   cfg.amplitude   = SNR of oscillatory burst versus aperiodic signal
%   cfg.cycles      = Number of cycles for oscillation
%   cfg.time        = Length of each 'trial'
%   cfg.trial       = Number of trials to simulate
%   cfg.freq        = Frequency of oscillation
%   Fs              = Sampling Rate
%   aperiodic       = aperiodic data (to be changed)
%
% % EXAMPLE USEAGE:   [aperiodic_out, osc_out, combined] = sim_fBOSC(cfg,Fs,aperiodic)
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Trust Centre for Neuroimaging
%
% Author: Robert Seymour      (rob.seymour@ucl.ac.uk);
%
% Modifed from Kosciessa et al., (2020) NI
% https://github.com/jkosciessa/eBOSC_resources_NI2020/blob/master/...
% simulations/B_simulation_standardBOSC.m
%__________________________________________________________________________
%

if ~isfield(cfg, 'amplitude')
    cfg.amplitude = 24;
end

if ~isfield(cfg, 'cycles')
    cfg.cycles = 64;
end

if ~isfield(cfg, 'time')
    cfg.time = 20;
end

if ~isfield(cfg, 'trial')
    cfg.trial = 20;
end

if ~isfield(cfg, 'freq')
    cfg.freq = 10;
end

% time vector
t = 1/Fs;
time = [t:t:cfg.time];

combined            = [];
combined.label      = {'chan1'};
combined.time       = [];
combined.trial      = [];
combined.fsample    = Fs;

osc_out         = zeros(cfg.trial,length(time));
aperiodic_out   = zeros(cfg.trial,length(time));

%% Start of trial loop
for k = 1:cfg.trial
    disp(k);
    alphaCycles = cfg.cycles;
    alphaFreq = cfg.freq;
    
    alphaTime = round((alphaCycles/alphaFreq),3);

    % make alpha symmetrical around the middle
    timeNew = round(alphaTime/t,0);
    if mod(timeNew,2) ~= 0
        timeNew = timeNew + 1;
        alphaTime = timeNew.*t;
    else alphaTime = timeNew.*t;
    end
    timeAlpha = [t:t:alphaTime];
    AlphaPlace = (numel(time)/2)-(numel(timeAlpha)/2)+1:(numel(time)/2)+(numel(timeAlpha)/2);
    
%%    % filter entire signal    
try
    [tmp_bpsignal] = ft_preproc_bandpassfilter(aperiodic(k,:), 500,...
        [cfg.freq-2 cfg.freq+2], 5);
catch
    [tmp_bpsignal] = ft_preproc_bandpassfilter(aperiodic(k,:), 500,...
        [cfg.freq-1 cfg.freq+1], 4);
    if k == 1
        warning('Using BP filter of cfg.freq +-1Hz');
    end
end
    
    % Generate signal with specified SNR
    VarBG = var(tmp_bpsignal);
    
    targetPower = cfg.amplitude*VarBG;
    amplitudeFromRMS = (sqrt(targetPower)*sqrt(2));
    alpha_sim = sin(timeAlpha*2*pi*alphaFreq)*amplitudeFromRMS;
    osc = zeros(1,numel(time));
    osc(1,AlphaPlace) = alpha_sim;
    
    
    combined.time{k}  = time;
    combined.trial{k} = aperiodic(k,:) + osc;
    osc_out(k,:)      = osc;
    aperiodic_out(k,:)      = aperiodic(k,:);
    
    clear osc VarBG targetPower alpha_sim amplitudeFromRMS

end




