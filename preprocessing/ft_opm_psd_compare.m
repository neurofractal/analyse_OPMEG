function [shield,freq] = ft_opm_psd_compare(cfg,rawData1,rawData2)
% Function to compute relative PSD of two OPM datasets 
% (for checking shielding factors). Designed for (OPMEG) data acquired 
% from the UCL Wellcome Centre for Neuroimaging.
%
% EXAMPLE USEAGE:   [shield,freq] = ft_opm_psd_compare(cfg,rawData1,rawData2)
% ...where, cfg is the input structure and rawData1 and rawData2 are  
% raw OPM data organised in Fieldtrip structures.
% 
%   cfg. method         = 'fieldtrip' to use ft_freqanalysis or 'tim' for
%                       custom PSD code from Tim Tierney
%   cfg.foi             = frequencies of interest in form [X Y].
%                       Default = [1 100]
%   cfg.trial_length    = length of segments (in seconds). Default = 1.
%   cfg.channel         = 'all', 'MEG', 'RAD', 'TAN'. Default = 'all'.
%   cfg.plot            = 'yes' or 'no'
%   cfg.dB              = do you want to plot dB (default = 'yes',
%                       highly recommended)

%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   N.B. If the inputs have different numbers of channels, please specify
%   exactly which channels to use, via cfg.channel
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_rpsd (Tim Tierney)

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%__________________________________________________________________________

%% Function housekeeping
if ~isfield(cfg, 'channel')
    cfg.channel = 'all';
end

if ~isfield(cfg, 'trial_length')
    cfg.trial_length = 1;
end

if ~isfield(cfg, 'method')
    cfg.method = 'tim';
end

if ~isfield(cfg, 'foi')
    cfg.foi = [1 100];
end

if ~isfield(cfg, 'plot')
    cfg.plot = 'yes';
end

if ~isfield(cfg, 'dB')
    cfg.dB = 'yes';
end

% Warn the user if the input isn't as expected
if ~strcmp(cfg.method,'tim')
    ft_warning('Currently only cfg.method = ''tim'' is supported');
    %pause(1.0);
end

if ~strcmp(cfg.dB,'yes')
    ft_warning(['Consider using cfg.dB = ''yes''. Subtracting PSDs',...
        ' from one another will result in odd-looking plots']);
    pause(1.0);
end


% Calculate the PSDs
cfg2 = cfg;
cfg2.plot = 'no';
[pow1, ~, label1] = ft_opm_psd(cfg2,rawData1);
[pow2, freq, ~] = ft_opm_psd(cfg2,rawData2);

% Check the size of the arrays is the same
if size(pow1,2) ~= size(pow2,2)
    error('Possible mismatch between the number of channels'); 
end

if size(pow1,1) ~= size(pow2,1)
    error('Possible mismatch between sampling rate (Hz)');
end

% Calculate Shielding Factor
if strcmp(cfg.method,'tim')
    pow1 = median(pow1(:,:,:),3);
    pow2 = median(pow2(:,:,:),3);
end

if strcmp(cfg.dB,'yes')
    shield = 20*log10(pow1./pow2);
else
    shield = pow1-pow2;
end

% Plot
if strcmp(cfg.plot,'yes')
    % Calculate median out of all epochs
    figure()
    set(gcf,'Position',[100 100 1200 800]);
    % Plot all channels
    if strcmp(cfg.dB,'yes')
        plot(freq,shield,'LineWidth',1); hold on;
        % Plot the mean in black
        if strcmp(cfg.method,'matlab')
            plot(freq,squeeze(mean(shield,1)),'-k','LineWidth',2);
        else
            plot(freq,squeeze(mean(shield,2)),'-k','LineWidth',2);
        end
    else
        semilogy(freq,shield,'LineWidth',1); hold on;
        % Plot the mean in black
        semilogy(freq,squeeze(mean(shield,2)),'-k','LineWidth',2);
    end
    hold on;
    grid on;
    ax = gca; % current axes
    ax.FontSize = 20;
    ax.TickLength = [0.02 0.02];
    fig= gcf;
    fig.Color=[1,1,1];
    xlabel('Frequency (Hz)','FontSize',20)
    if strcmp(cfg.dB,'yes')
        labY = ['Shielding Factor (dB)'];
    else
        labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
    end
    ylabel(labY,'interpreter','latex','FontSize',20)
    
    % Adjust limits based on cfg.foi
    xlim([cfg.foi(1), cfg.foi(end)]);
    %ylim([1, 1000]);
    % Legend
    [~, hobj, ~, ~] = legend(vertcat(label1, 'mean'),'location','eastoutside');
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',4);
    ht = findobj(hobj,'type','text');
    set(ht,'FontSize',12);

else
    disp('NOT PLOTTING');
end


end





