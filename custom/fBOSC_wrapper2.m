%%
% fBOSC_wrapper() - Wrapper for fBOSC
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________
%    
%    This file is part of the fBOSC library.
%    License: The GNU General Public License v3.0
%
%    Built on work by Kosciessa and colleagues
%    https://github.com/jkosciessa/eBOSC
%__________________________________________________________________________

function [fBOSC, cfg] = fBOSC_wrapper2(cfg, data)
% Main fBOSC wrapper function. Executes fBOSC subfunctions.
%
% Use as
%   [fBOSC, cfg] = fBOSC_wrapper(cfg, data)
% where "data" is a Fieldtrip data structure with trial, time and channel
% fields.

% The configuration should contain:
%
%       GENERAL OPTIONS:
%   cfg.fBOSC.channel               = which channels to use (default: [] = all)
%   cfg.fBOSC.trial                 = which trials to use (default: [] = all)
%   cfg.fBOSC.trial_background      = use subset of trials for background fit? 
%                                     (default: [] = all)
%   cfg.fBOSC.F                     = frequency sampling
%   cfg.fBOSC.wavenumber            = wavelet number parameter (TFR tradeoff)
%   cfg.fBOSC.fsample               = sampling frequency of the data
%
%       PADDING OPTIONS:
%   cfg.fBOSC.pad.tfr_s             = padding following wavelet transform to 
%                                     avoid edge artifacts in seconds (bi-lateral)
%   cfg.fBOSC.pad.detection_s       = padding following rhythm detection in 
%                                     seconds (bi-lateral); 'shoulder' for 
%                                     BOSC eBOSC.detected matrix to account 
%                                     for duration threshold
%   cfg.fBOSC.pad.total_s           = complete padding (WL + shoulder)
%   cfg.fBOSC.pad.background_s      = padding of segments for background fit

%       FOOOF OPTIONS:
%   cfg.fBOSC.fooof.version         = 'matlab' or 'python' (default = 'python')
%   cfg.fBOSC.fooof.aperiodic_mode  = which approach to take for fitting the 
%                                     aperiodic component ('fixed' or'knee')
%   cfg.fBOSC.fooof.max_peaks       = max number of peaks allowed for modelling 
%                                     the power spectrum (default = 3)
%   cfg.fBOSC.fooof.min_peak_height = minimum peak height when modelling the
%                                     the power spectrum (default = 0.1)
%   cfg.fBOSC.fooof.fit_function    = fit the gaussian using 'fooof' or 
%                                     'matlab' gaussian fit function 
%                                     (default = 'fooof')
%   cfg.fBOSC.fooof.verbose         = fooof verbosity mode. If True, 
%                                     prints out warnings and general 
%                                     status updates (default = 0)
%
%       BURST DETECTION OPTIONS:
%   cfg.fBOSC.threshold.duration    = vector of duration thresholds at each frequency
%   cfg.fBOSC.threshold.percentile  = percentile of background fit for
%                                     power threshold (default = 0.95)
%       
%       POST-PROCESSING
%   cfg.fBOSC.postproc.use          = Post-processing of rhythmic 
%                                     fBOSC.episodes, i.e., wavelet 
%                                     'deconvolution' (default = 'no')
%   cfg.fBOSC.postproc.method       = Deconvolution method (default = 
%                                     'MaxBias', FWHM: 'FWHM')
%   cfg.fBOSC.postproc.edgeOnly     = Deconvolution only at on- and offsets 
%                                     of fBOSC.episodes? (default = 'yes')
%   cfg.fBOSC.postproc.effSignal	= Power deconvolution on whole signal 
%                                     or signal above power threshold? 
%                                     (default = 'PT')
%
% Outputs: 
%           fBOSC | main fBOSC output structure
%               fBOSC.episodes | table of individual rhythmic episodes (see fBOSC_episode_create)
%               fBOSC.detected | binary matrix of detected time-frequency points (prior to episode creation)
%               fBOSC.pepisode | temporal average of detected rhythms (prior to episode creation)
%               fBOSC.detected_ep | binary matrix of detected time-frequency points (following episode creation)
%               fBOSC.abundance_ep | temporal average of detected rhythms (following episode creation)
%           cfg | config structure

    fBOSC           = [];
    fBOSC.label     = [];
    
    % set some defaults for included channels and trials, if not specified
    if isempty(cfg.fBOSC.channel)
        cfg.fBOSC.channel = 1:numel(data.label);
    end
    if isempty(cfg.fBOSC.trial)
        cfg.fBOSC.trial = 1:numel(data.trial);
    end
    if isempty(cfg.fBOSC.trial_background)
        cfg.fBOSC.trial_background = 1:numel(data.trial);
    end
    
    % Default for % threshold at 95%
    if isempty(cfg.fBOSC.threshold.percentile)
        cfg.fBOSC.trial_background = 0.95;
    end
    
    % Some defaults for fooof:
    % Use python implementation by default (for now)
    if isempty(cfg.fBOSC.fooof.version)
        cfg.fBOSC.fooof.version = 'python';
    end  
    
    % Use fixed aperiodic mode as default
    if ~isfield(cfg.fBOSC.fooof,'aperiodic_mode')
        cfg.fBOSC.fooof.aperiodic_mode = 'fixed';
        disp('FOOOF will fit without a knee parameter');
    end
    
    % Verbose for fooof
    if ~isfield(cfg.fBOSC.fooof,'verbose')
        cfg.fBOSC.fooof.verbose = 0;
    end

    % fit_function for fooof
    if ~isfield(cfg.fBOSC.fooof,'fit_function')
        cfg.fBOSC.fooof.fit_function = 'fooof';
    end
    
    % max_peaks for fooof
    if ~isfield(cfg.fBOSC.fooof,'max_peaks')
        cfg.fBOSC.fooof.max_peaks =  3;
    end
    
    % min_peak_height for fooof
    if ~isfield(cfg.fBOSC.fooof,'min_peak_height')
        cfg.fBOSC.fooof.min_peak_height =  0.1;
    end
    
    % calculate the sample points for paddding
    cfg.fBOSC.pad.tfr_sample = cfg.fBOSC.pad.tfr_s.*cfg.fBOSC.fsample;                          % automatic sample point calculation
    cfg.fBOSC.pad.detection_sample = cfg.fBOSC.pad.detection_s.*cfg.fBOSC.fsample;              % automatic sample point calculation
    cfg.fBOSC.pad.total_s = cfg.fBOSC.pad.tfr_s + cfg.fBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
    cfg.fBOSC.pad.total_sample = cfg.fBOSC.pad.tfr_sample + cfg.fBOSC.pad.detection_sample;
    cfg.fBOSC.pad.background_sample = cfg.fBOSC.pad.tfr_sample;
    
    % Set up empty array for memory efficiency
    fBOSC.detected  = zeros(numel(cfg.fBOSC.channel),numel(cfg.fBOSC.trial),...
        length(cfg.fBOSC.F),length(data.trial{1})-...
        (cfg.fBOSC.pad.detection_sample*2)-(cfg.fBOSC.pad.detection_sample*2));
    
    fBOSC.detected_ep  = zeros(numel(cfg.fBOSC.channel),numel(cfg.fBOSC.trial),...
        length(cfg.fBOSC.F),length(data.trial{1})-...
        (cfg.fBOSC.pad.detection_sample*2)-(cfg.fBOSC.pad.detection_sample*2));
    
    fBOSC.pepisode  = zeros(numel(cfg.fBOSC.channel),numel(cfg.fBOSC.trial),...
        length(cfg.fBOSC.F));
    
    fBOSC.abundance_ep  = zeros(numel(cfg.fBOSC.channel),numel(cfg.fBOSC.trial),...
        length(cfg.fBOSC.F));

    
    % Start of Channel Loop:
    for indChan = 1: numel(cfg.fBOSC.channel)
    
        display(['Channel ',num2str(indChan), '/', ...
            num2str(numel(cfg.fBOSC.channel)),': ',...
            data.label{cfg.fBOSC.channel(indChan)}]);
        
        % Add data channel label to fBOSC
        fBOSC.label{1,indChan} = data.label{cfg.fBOSC.channel(indChan)};
        
        cfg.tmp.channel = indChan; % encode current channel for later

        %% Step 1: time-frequency wavelet decomposition for whole signal to prepare background fit
            
        % Preallocate memory - 8.4s
        TFR = zeros(numel(cfg.fBOSC.trial),length(cfg.fBOSC.F),...
            length(data.trial{1}));

        for indTrial = 1:numel(cfg.fBOSC.trial)
            TFR(indTrial,:,:) = BOSC_tf(data.trial{indTrial}...
                (cfg.fBOSC.channel(indChan),:),cfg.fBOSC.F,...
                cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
        end; clear indTrial
        
        %% Step 2: robust background power fit
        [fBOSC, pt, dt] = fBOSC_getThresholds2(cfg, TFR, fBOSC);

        %% Application of thresholds to single trials

        for indTrial = 1:numel(cfg.fBOSC.trial)

            cfg.tmp.trial = cfg.fBOSC.trial(indTrial); % encode current trial for later

            % get wavelet transform for single trial
            % tfr padding is removed to avoid edge artifacts from the wavelet
            % transform. Note that a padding fpr detection remains attached so that there
            % is no problems with too few sample points at the edges to
            % fulfill the duration criterion.
            TFR_ = squeeze(TFR(indTrial,:,...
                cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample));
            
            %% Step 3: detect rhythms and calculate Pepisode

            % The next section applies both the power and the duration
            % threshold to detect individual rhythmic segments in the continuous signals.
            detected = zeros(size(TFR_));
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect...
                    (TFR_(f,:),pt(f),dt(f),cfg.fBOSC.fsample);
            end; clear f
            
            % remove padding for detection (matrix with padding required for refinement)
            fBOSC.detected(indChan, indTrial,:,:) = ...
                detected(:,cfg.fBOSC.pad.detection_sample...
                +1:end-cfg.fBOSC.pad.detection_sample);

            % encode pepisode of detected rhythms (optional)
            fBOSC.pepisode(indChan, indTrial,:) = ...
                mean(fBOSC.detected(indChan, indTrial,:,:),4);
            
            %% Copy fBOSC to EBOSC for consistency
            cfg.eBOSC = cfg.fBOSC;
           
            %% Step 4 (optional): create table of separate rhythmic episodes
            
            cfg.tmp.inputTime = data.time{cfg.tmp.trial};
            cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample);
            cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
            
            
            [fBOSC.episodes, detected_ep] = eBOSC_episode_create(cfg,TFR_,detected,fBOSC);
            
            % remove padding for detection (already done for fBOSC.episodes)
            
            fBOSC.detected_ep(indChan, indTrial,:,:) = detected_ep(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
            clear detected_ep;
            

            % encode abundance of fBOSC.episodes (optional)
            fBOSC.abundance_ep(indChan, indTrial,:) = mean(squeeze(fBOSC.detected_ep(indChan, indTrial,:,:)),2);
            
            clear TFR_
        end; clear indTrial; % trial loop

    end % channel loop
end