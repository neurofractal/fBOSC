%%
% fooof_plot() - Plot a FOOOF model.
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


function fBOSC_fooof_plot(cfg,fBOSC)

if ~isfield(cfg,'log_freqs')
    cfg.log_freqs = 1;
end

if ~isfield(cfg,'plot_old')
    cfg.plot_old = 0;
end

%% Set Up

if cfg.log_freqs
    plt_freqs = log10(cfg.fBOSC.F);
else
    plt_freqs = cfg.fBOSC.F;
end

% Plot settings
lw = 2.5;

%% Create the plots

for lab = 1:length(fBOSC.label)
    
    figure();
    set(gcf,'Position',[100 100 1200 700])
    hold on
    
    % Plot the original data
    data = plot(plt_freqs, fBOSC.static.bg_log10_pow(lab,:), 'black');
    
    if cfg.plot_old
        % Plot the aperiodic fit
        ap_fit = plot(plt_freqs, log10(fBOSC.static.mp_old(lab,:)), 'b--');
        
        % Plot the Statistical Threshold
        pt_fit= plot(plt_freqs, log10(fBOSC.static.pt_old(lab,:)), 'g--');
    else
        % Plot the full model fit
        model = plot(plt_freqs, fBOSC.static.fooofed_spectrum(lab,:), 'red');
        
        % Plot the aperiodic fit
        ap_fit = plot(plt_freqs, fBOSC.static.ap_fit(lab,:), 'b--');
        
        % Plot the Statistical Threshold
        pt_fit = plot(plt_freqs, log10(fBOSC.static.pt(lab,:)), 'g--');
    end
    
    %% Plot Settings
    
    % Apply general plot settings
    if cfg.plot_old
        for plt = [data, ap_fit pt_fit]
            set(plt, 'LineWidth', lw);
        end
    else
        for plt = [data, model, ap_fit pt_fit]
            set(plt, 'LineWidth', lw);
        end
    end
    
    % Set alpha value for model - in a wonky way, because Matlab
    %   Note: the '4' is magical and mysterious. No idea.
    model.Color(4) = 0.5;
    
    grid on
    legend('Original Spectrum', 'Full Model Fit', 'Aperiodic Fit')
    
    set(gca,'Fontsize',20)
    xlabel('log10(Frequency) (Hz)','FontSize',30);
    ylabel('log10(Power) (a.u.)','FontSize',30);
    if cfg.plot_old
        legend({'Original Spectrum', ...
            'Aperiodic Fit','BOSC Threshold'}, ...
            'orientation', 'vertical', 'location', 'northeastoutside','FontSize',20);
    else
        legend({'Original Spectrum', 'FOOOF Model Fit', ...
            'Aperiodic Fit','BOSC Threshold'}, ...
            'orientation', 'vertical', 'location', 'northeastoutside','FontSize',20);
    end
    legend('boxoff');
    
    if cfg.plot_old
        title(['Channel = ' fBOSC.label{lab}],'interpreter','none');
    else
        title({['Channel = ' fBOSC.label{lab}]; ...
            ['Goodness of Fit = ' num2str(fBOSC.static.r_squared(lab))]},...
            'interpreter','none');
    end
    
    hold off
    drawnow;
end

end