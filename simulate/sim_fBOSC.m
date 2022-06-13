function [aperiodic_out, osc_out, combined,AlphaPlace] = sim_fBOSC(cfg,Fs,aperiodic)
% Function to combine simulated aperiodic data with a simulated
% oscillatory burst of known duration and amplitude.
%
%   INPUTS:
%   cfg.time         = Length of each 'trial'
%   cfg.trial        = Number of trials to simulate
%   cfg.amplitude    = SNR of burst with regard to aperiodic signal
%   cfg.maxBurstTime = Max total burst time (in seconds, [min max])
%   cfg.minGap       = Minimum gap between bursts (in seconds)
%   cfg.cycles       = Number of cycles for the bursts [min max]
%   cfg.freq         = Frequency of burst
%   Fs               = Sampling Rate
%   aperiodic        = aperiodic data (to be changed)
%
% % EXAMPLE USEAGE:   [aperiodic_out, osc_out, combined] = sim_fBOSC(cfg,Fs,aperiodic)
%__________________________________________________________________________
% Copyright (C) 2021-22 Wellcome Trust Centre for Neuroimaging
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
    cfg.cycles = [4 8];
end

if ~isfield(cfg, 'time')
    cfg.time = 20;
end

if ~isfield(cfg, 'maxBurstTime')
    cfg.maxBurstTime = [1 2];
end

if ~isfield(cfg, 'min_gap')
    cfg.min_gap = 0.5;
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
    %disp(k);
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
        [cfg.freq-1 cfg.freq+1], 3);
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
    
    % We need to make distinct 'bursts'
    time_so_far  = 0;     
    bursts       = {};
    count        = 1;
    cumul_length = 0;
    
    % Make bursts of variable cycle numbers up to cfg.maxBurstTime
    while time_so_far < cfg.maxBurstTime
        
        val       = randi([cfg.cycles(1) cfg.cycles(2)]);
        cycle_val = round((val/alphaFreq),3);

        timeAlpha = [t:t:cycle_val];
        
        bursts{count} = sin(timeAlpha*2*pi*alphaFreq)*amplitudeFromRMS;
        cumul_length = cumul_length+length(bursts{count});
        
        time_so_far = time_so_far+(length(bursts{count})/Fs);
        count = count+1;
    end
    
    
    % Create random gaps (0s) equal to time remaining in the trial
    % Make sure gaps are longer than cfg.min_gap
    min_time = 1;
    while min_time
        randint = randi([cfg.min_gap*Fs numel(time)-cumul_length],numel(bursts),1);
        randint = [0 randint' numel(time)-cumul_length];
        
        zero_lengths = (diff(sort(randint)));
        
        if min(zero_lengths) > cfg.min_gap*Fs
            min_time = 0;
        end
        
    end
    
    % Add each element to trial and create logical array for when a burst
    % was simulated in the trial
    osc = [];
    alpha_sim = zeros(length(aperiodic),1);
    
    for b = 1:length(zero_lengths)
        osc_zero = zeros(1,zero_lengths(b));
        
        if b <= length(bursts)
            osc_zero = horzcat(osc_zero,bursts{b});
            alpha_sim((length(osc)+zero_lengths(b)):(length(osc)+length(osc_zero))) = 1;
        end
        
        osc = horzcat(osc,osc_zero);
        
    end
    
    % Figure for debugging
    % figure; plot(osc); hold on; plot(alpha_sim);
    
    % Variables to export
    AlphaPlace(k,:) = alpha_sim;
    combined.time{k}  = time;
    combined.trial{k} = aperiodic(k,:) + osc;
    osc_out(k,:)      = osc;
    aperiodic_out(k,:)      = aperiodic(k,:);
    
    clear osc VarBG targetPower alpha_sim amplitudeFromRMS

end



