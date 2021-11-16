function [fBOSC, cfg] = fBOSC_wrapper(cfg, data)
% Main fBOSC wrapper function. Executes fBOSC subfunctions.
%
% Inputs: 
%           cfg | config structure with cfg.eBOSC field:
%                     cfg.fBOSC.F                     | frequency sampling
%                     cfg.fBOSC.wavenumber            | wavelet family parameter (time-frequency tradeoff)
%                     cfg.fBOSC.fsample               | current sampling frequency of EEG data
%                     cfg.fBOSC.pad.tfr_s             | padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
%                     cfg.fBOSC.pad.detection_s       | padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
%                     cfg.fBOSC.pad.total_s           | complete padding (WL + shoulder)
%                     cfg.fBOSC.pad.background_s      | padding of segments for BG (only avoiding edge artifacts)
%                     cfg.fBOSC.fooof.aperiodic_mode  | which approach to take for fitting the aperiodic component ('fixed' or'knee')
%                     cfg.fBOSC.fooof.verbose         | fooof verbosity mode. If True, prints out warnings and general status updates (default = 0)
%                     cfg.fBOSC.threshold.duration    | vector of duration thresholds at each frequency (previously: ncyc)
%                     cfg.fBOSC.threshold.percentile  | percentile of background fit for power threshold
%                     cfg.fBOSC.postproc.use          | Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
%                     cfg.fBOSC.postproc.method       | Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
%                     cfg.fBOSC.postproc.edgeOnly     | Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
%                     cfg.fBOSC.postproc.effSignal	  | Power deconvolution on whole signal or signal above power threshold? (default = 'PT')
%                     cfg.fBOSC.channel               | Subset of channels? (default: [] = all)
%                     cfg.fBOSC.trial                 | Subset of trials? (default: [] = all)
%                     cfg.fBOSC.trial_background      | Subset of trials for background? (default: [] = all)
%           data | input time series data in FieldTrip format with:
%                | .trial field: {trial}(channel x time)
%                | .time field: {trial}(channel x time)
%                | .label field: {channelName}
%
% Outputs: 
%           fBOSC | main eBOSC output structure
%               fBOSC.episodes | table of individual rhythmic episodes (see eBOSC_episode_create)
%               fBOSC.detected | binary matrix of detected time-frequency points (prior to episode creation)
%               fBOSC.pepisode | temporal average of detected rhythms (prior to episode creation)
%               fBOSC.detected_ep | binary matrix of detected time-frequency points (following episode creation)
%               fBOSC.abundance_ep | temporal average of detected rhythms (following episode creation)
%           cfg | config structure

    fBOSC       = [];
    fBOSC.label = [];
    
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
    
    % Some defaults for fooof
    if ~isfield(cfg.fBOSC.fooof,'aperiodic_mode')
        cfg.fBOSC.fooof.aperiodic_mode = 'fixed';
        disp('FOOOF will fit without a knee parameter');
    end
    
    if ~isfield(cfg.fBOSC.fooof,'verbose')
        cfg.fBOSC.fooof.verbose = 0;
    end   
    
    % calculate the sample points for paddding
    cfg.fBOSC.pad.tfr_sample = cfg.fBOSC.pad.tfr_s.*cfg.fBOSC.fsample;                          % automatic sample point calculation
    cfg.fBOSC.pad.detection_sample = cfg.fBOSC.pad.detection_s.*cfg.fBOSC.fsample;              % automatic sample point calculation
    cfg.fBOSC.pad.total_s = cfg.fBOSC.pad.tfr_s + cfg.fBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
    cfg.fBOSC.pad.total_sample = cfg.fBOSC.pad.tfr_sample + cfg.fBOSC.pad.detection_sample;
    cfg.fBOSC.pad.background_sample = cfg.fBOSC.pad.tfr_sample;

    for indChan = 1: numel(cfg.fBOSC.channel)
    
        display(['Channel ',num2str(indChan), '/', ...
            num2str(numel(cfg.fBOSC.channel)),': ',...
            data.label{cfg.fBOSC.channel(indChan)}]);
        
        % Add data channel label to fBOSC
        fBOSC.label{1,indChan} = data.label{cfg.fBOSC.channel(indChan)};
        
        cfg.tmp.channel = indChan; % encode current channel for later

        %% Step 1: time-frequency wavelet decomposition for whole signal to prepare background fit

        TFR = [];
        for indTrial = 1:numel(cfg.fBOSC.trial)
            TFR.trial{indTrial} = BOSC_tf(data.trial{indTrial}...
                (cfg.fBOSC.channel(indChan),:),cfg.fBOSC.F,...
                cfg.fBOSC.fsample,cfg.fBOSC.wavenumber);
        end; clear indTrial

        %% Step 2: robust background power fit (see 2020 NeuroImage paper)
        [fBOSC, pt, dt] = fBOSC_getThresholds(cfg, TFR, fBOSC);

        %% Application of thresholds to single trials

        for indTrial = 1:numel(cfg.fBOSC.trial)

            cfg.tmp.trial = cfg.fBOSC.trial(indTrial); % encode current trial for later

            % get wavelet transform for single trial
            % tfr padding is removed to avoid edge artifacts from the wavelet
            % transform. Note that a padding fpr detection remains attached so that there
            % is no problems with too few sample points at the edges to
            % fulfill the duration criterion.
            TFR_ = TFR.trial{indTrial}(:,cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample);

            %% Step 3: detect rhythms and calculate Pepisode

            % The next section applies both the power and the duration
            % threshold to detect individual rhythmic segments in the continuous signals.
            detected = zeros(size(TFR_));
            for f = 1:length(cfg.fBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt(f),dt(f),cfg.fBOSC.fsample);
            end; clear f

            % remove padding for detection (matrix with padding required for refinement)
            fBOSC.detected(indChan, indTrial,:,:) = detected(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);

            % encode pepisode of detected rhythms (optional)
            fBOSC.pepisode(indChan, indTrial,:) = mean(fBOSC.detected(indChan, indTrial,:,:),4);
            
%             %% Copy fBOSC to EBOSC for consistency
%             cfg.eBOSC = cfg.fBOSC;
%             
%             %% Step 4 (optional): create table of separate rhythmic episodes
% 
%             cfg.tmp.inputTime = data.time{cfg.tmp.trial};
%             cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.fBOSC.pad.tfr_sample+1:end-cfg.fBOSC.pad.tfr_sample);
%             cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
%             
%             [fBOSC.episodes, detected_ep] = eBOSC_episode_create(cfg,TFR_,detected,fBOSC);
%             
%             % remove padding for detection (already done for fBOSC.episodes)
%             fBOSC.detected_ep(indChan, indTrial,:,:) = detected_ep(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample);
%             clear detected_ep;
% 
%             % encode abundance of fBOSC.episodes (optional)
%             fBOSC.abundance_ep(indChan, indTrial,:) = mean(squeeze(fBOSC.detected_ep(indChan, indTrial,:,:)),2);

%           % Supplementary Plot: original eBOSC.detected vs. sparse episode power
%           figure; 
%           subplot(121); imagesc(squeeze(eBOSC.detected).*TFR_(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample));
%           subplot(122); imagesc(squeeze(eBOSC.detected_ep).*TFR_(:,cfg.fBOSC.pad.detection_sample+1:end-cfg.fBOSC.pad.detection_sample));

        end; clear indTrial; % trial loop

    end % channel loop
end