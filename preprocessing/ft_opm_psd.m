function [pow, freq, label] = ft_opm_psd(cfg,rawData)
% Function to calculate PSD on optically-pumped magnetencephalography
% (OPMEG) data acquired from the UCL Wellcome Centre for Neuroimaging.
%
% EXAMPLE USEAGE:   data = ft_opm_psd(cfg,rawData)
% ...where, cfg is the input structure and rawData is the raw OPM data
% loaded using ft_opm_create
%
%   cfg. method         = 'fieldtrip' to use ft_freqanalysis or 'tim' for
%                       custom PSD code from Tim Tierney
%   cfg.foi             = frequencies of interest in form [X Y].
%                       Default = [1 100]
%   cfg.trial_length    = length of segments (in seconds). Default = 1.
%   cfg.channel        = 'all', 'MEG', 'RAD', 'TAN'. Default = 'all'.
%   cfg.plot            = 'yes' or 'no'
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_create (Tim Tierney)

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%           Nicholas Alexander  (n.alexander@ucl.ac.uk)
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

%% Select the data based on cfg.channel option
if strcmp(cfg.channel,'all')
    disp('Calculating PSD for ALL channels');
else
    try
        chan = cfg.channel;
        cfg2 = [];
        cfg2.channel = ft_channelselection_opm(chan,rawData,...
            'quspin_g2');
        rawData = ft_selectdata(cfg2,rawData);
        
    catch
        ft_warning(['Could not select ' chan]);
    end
end

label = rawData.label;

%%
method_for_fft = cfg.method;

switch method_for_fft
    case 'fieldtrip'
        ft_warning('NOT TESTED DO NOT USE');
        foi = [cfg.foi(1):0.05:cfg.foi(end)];
        
        disp(['Calculating PSD from ' num2str(cfg.foi(1)) 'Hz to ' ...
            num2str(cfg.foi(end)) 'Hz']);
        
        trial_length    = cfg.trial_length;
        
        cfg             = [];
        cfg.length      = trial_length;
        cfg.overlap     = 0;
        rawData         = ft_redefinetrial(cfg, rawData);
        
        cfg2 = [];
        cfg2.output      = 'pow';          % Return PSD
        cfg2.channel     = 'all';
        cfg2.method      = 'mtmfft';
        %cfg.pad         = 'nextpow2';
        cfg2.taper       = 'boxcar';      % Hann window as taper
        cfg2.foi         = foi;         % Frequency range
        
        psd_hann = ft_freqanalysis(cfg2, rawData);
        pow = psd_hann.powspctrm;

        figure;
        set(gcf,'Position',[100 100 1200 800]);
        semilogy(psd_hann.freq,psd_hann.powspctrm);
        set(gca,'FontSize',14);
        xlabel('Frequency (Hz)');
        ylabel('Power/Frequency (au^2/Hz)')
        hold on
        semilogy(psd_hann.freq,mean(log(psd_hann.powspctrm)),'-k','LineWidth',2);
        lgd = legend(vertcat(rawData.label, 'mean'));
        set(lgd,'Location','BestOutside');
        
        % Another method from Tim T
    case 'tim'
        
        % Split data into epochs
        nsamps = (cfg.trial_length)*rawData.fsample;
        beg = 1:nsamps:size(rawData.trial{1},2);
        endsamp =  beg+(nsamps-1);
        inRange = ~(beg>size(rawData.trial{1},2)|endsamp>size(rawData.trial{1},2));
        chans = 1:1:(size(rawData.trial{1},1));
        eD = zeros(length(chans),nsamps,sum(inRange));
        for i =1:length(inRange)
            if(inRange(i))
                eD(:,:,i)=rawData.trial{1}(chans,beg(i):endsamp(i),:);
            end
        end
        
        % Get variables for calculating PSD
        fs = rawData.fsample;
        N = size(eD,2);
        Nf = ceil((N+1)/2);
        nepochs = size(eD,3);
        pow = zeros(Nf,size(eD,1),nepochs);
        wind  = window(@flattopwin,size(eD,2));
        wind = repmat(wind,1,size(eD,1));
        
        % Calculate PSD for each epoch
        for j = 1:nepochs
            Btemp=eD(:,:,j)';
            Btemp = Btemp.*wind;
            mu=mean(Btemp);
            zf = bsxfun(@minus,Btemp,mu);
            fzf = zf;
            N= length(fzf);
            xdft = fft(fzf);
            xdft=xdft(1:floor(N/2+1),:);
            psdx = abs(xdft)./sqrt(N*fs);
            freq = 0:fs/size(fzf,1):fs/2;
            odd=mod(size(fzf,1),2)==1;
            if(odd)
                psdx(2:end) = sqrt(2)*psdx(2:end);
            else
                psdx(2:end-1) = sqrt(2)*psdx(2:end-1);
            end
            pow(:,:,j) =psdx;
        end
        
        % Plot
        if strcmp(cfg.plot,'yes')
            % Calculate median out of all epochs
            po = median(pow(:,:,:),3);
            figure()
            set(gcf,'Position',[100 100 1200 800]);
            % Plot all channels
            semilogy(freq,po,'LineWidth',1); hold on;
            % Plot the mean in black
            semilogy(freq,squeeze(mean(po,2)),'-k','LineWidth',2);
            hold on
            xp2 =0:round(freq(end));
            yp2=ones(1,round(freq(end))+1)*15;
            p2 =plot(xp2,yp2,'--k');
            p2.LineWidth=2;
            xlabel('Frequency (Hz)')
            labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
            ylabel(labY,'interpreter','latex')
            grid on
            ax = gca; % current axes
            ax.FontSize = 13;
            ax.TickLength = [0.02 0.02];
            fig= gcf;
            fig.Color=[1,1,1];
            
            % Adjust limits based on cfg.foi
            xlim([cfg.foi(1), cfg.foi(end)]);
            %ylim([1, 1000]);
            lgd = legend(vertcat(label, 'mean'));
        else
            disp('NOT PLOTTING');
        end

        
end


end




