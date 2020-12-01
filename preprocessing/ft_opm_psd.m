function [pow, freq, label] = ft_opm_psd(cfg,rawData)
% Function to calculate PSD on optically-pumped magnetencephalography
% (OPMEG) data acquired from the UCL Wellcome Centre for Neuroimaging.
%
% EXAMPLE USEAGE:   data = ft_opm_psd(cfg,rawData)
% ...where, cfg is the input structure and rawData is the raw OPM data
% loaded using ft_opm_create
%
%   cfg. method         = 'tim', 'fieldtrip' or 'matlab' (default = 'tim')
%   cfg.foi             = frequencies of interest in form [X Y].
%                       (default = [1 100]).
%   cfg.trial_length    = length of segments in seconds (default = 1).
%                       Set to [] to make width of trial.
%   cfg.channel         = 'all', 'MEG', 'RAD', 'TAN' (default = 'all').
%   cfg.plot            = 'yes' or 'no' (default = 'yes')
%   cfg.plot_legend     = 'yes' or 'no' (default = 'yes')
%   cfg.average         = how to average over blocks: 'median' or 
%                         'mean' (default = 'median')
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_create (Tim Tierney)

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%           Nicholas Alexander  (n.alexander@ucl.ac.uk)
%           George O'Neill      
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

if ~isfield(cfg, 'plot_legend')
    cfg.plot_legend = 'yes';
end

if ~isfield(cfg, 'average')
    cfg.average = 'median';
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
    % 1. MATLAB FFT Method, with code adapted from ft_icabrowser
    case 'matlab'
        smo = 50;
        steps = 10;
        Fs = rawData.fsample;
        N = floor(size(rawData.trial{1},2));
        chan = size(rawData.trial{1},1);
        
        % Perhaps there is a better way to do this than per channel in a
        % loop?
        
        po = [];
        
        ft_progress('init', 'text', 'Please wait...')
        
        for c = 1:chan
            
            % Display the progress of the function
            ft_progress(c/chan,'Calculating PSD for channel: %s',label{c});
            
            xdft = fft(rawData.trial{1}(c,:));
            xdft = xdft(1:floor(N/2)+1);
            psdx = (1/(Fs*N)).*abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            
            smoothed = conv(psdx, ones(smo,1), 'valid');
            smoothed = smoothed(1:steps:fix(length(psdx) - smo))/smo;
            
            freq = linspace(0,Fs/2,size(smoothed,2));
            
            po(c,:) = smoothed;
        end
        
        strt = find(freq > cfg.foi(1),1,'first');
        stp  = find(freq < cfg.foi(2),1,'last');
        
        pow = po;
        
        ft_progress('close');
        
        % Funky colorscheme, looks OK...
        colormap123     = linspecer(length(label));
        
        if strcmp(cfg.plot,'yes')
            % Make a Figure
            figure()
            set(gcf,'Position',[100 100 1200 800]);
            
            h = plot(freq(strt:stp),log10(po(:,strt:stp)),...
                'LineWidth',1);
            set(h, {'color'},num2cell(colormap123,2));
            hold on;
            plot(freq(strt:stp),log10(mean(po(:,strt:stp),1)),'-k','LineWidth',2);
            grid on
            ax = gca; % current axes
            ax.FontSize = 20;
            ax.TickLength = [0.02 0.02];
            ylabel('PSD (T^2/Hz)','FontSize',30);
            xlabel('Frequency (Hz)','FontSize',30);
            % Legend
            [~, hobj, ~, ~] = legend(vertcat(label, 'mean'),'location','eastoutside');
            hl = findobj(hobj,'type','line');
            set(hl,'LineWidth',4);
            ht = findobj(hobj,'type','text');
            set(ht,'FontSize',12);
        end
        
        % 2. RS Homebrewed Script Using ft_freqanalysis. Not quite right...
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
        
        % 3. The method from Tim T using data split up into different trial
        % lengths
        
    case 'tim'
        
        % if cfg.trial_length is an empty array we'll get it to guess the
        % length of the trial
        if isempty(cfg.trial_length)
            cfg.trial_length = range(rawData.time{1});
        end
        
        if cfg.trial_length > range(rawData.time{1})
            error('You are specifying PSD windows longer than the trial length!')
        end
        
        % Split data into epochs
        nsamps = (cfg.trial_length)*rawData.fsample;
        ntrials = numel(rawData.trial);
        
        for ii = 1:ntrials
            beg = 1:nsamps:size(rawData.trial{ii},2);
            endsamp =  beg+(nsamps-1);
            inRange = ~(beg>=size(rawData.trial{ii},2)|endsamp>=size(rawData.trial{ii},2));
            chans = 1:1:(size(rawData.trial{ii},1));
            if ii == 1
                eD = zeros(length(chans),nsamps,sum(inRange),ntrials);
            end
            for jj = 1:length(inRange)
                if(inRange(jj))
                    eD(:,:,jj,ii)=rawData.trial{ii}(chans,beg(jj):endsamp(jj),:);
                end
            end
        end
        
        % Reshape eD so its 3D
        eD = reshape(eD,size(eD,1),size(eD,2),[]);
        
        % Get variables for calculating PSD
        fs = rawData.fsample;
        N = size(eD,2);
        Nf = ceil((N+1)/2);
        nepochs = size(eD,3);
        pow = zeros(Nf,size(eD,1),nepochs);
        wind  = window(@hanning,size(eD,2));
        coFac = max(wind)/mean(wind);
        wind = repmat(wind,1,size(eD,1));
        
        % Calculate PSD for each epoch
        for j = 1:nepochs
            Btemp=eD(:,:,j)';
            Btemp = Btemp.*wind*coFac;
            Btemp = Btemp;
            
            mu=median(Btemp);
            zf = bsxfun(@minus,Btemp,mu);
            fzf = zf;
            N= length(fzf);
            xdft = fft(fzf);
            xdft=xdft(1:floor(N/2+1),:);
            psdx = abs(xdft)./sqrt(N*fs);
            freq = 0:fs/size(fzf,1):fs/2;
            odd=mod(size(fzf,1),2)==1;
            if(odd)
                %psdx(2:end) = sqrt(2)*psdx(2:end);
                psdx(2:end) = psdx(2:end);
            else
                %psdx(2:end-1) = sqrt(2)*psdx(2:end-1);
                psdx(2:end-1) = psdx(2:end-1);
            end
            pow(:,:,j) =psdx;
        end
        
        % Plot
        if strcmp(cfg.plot,'yes')
            
            % Try loading a fancy colorscheme
            try
                colormap123     = linspecer(length(label));
            catch
                disp('Using default colorscheme')
            end
            
            % Calculate median or mean out of all epochs
            if strcmp(cfg.average,'median')
                po = median(pow(:,:,:),3);
            else
                po = mean(pow(:,:,:),3);
            end
            
            figure()
            set(gcf,'Position',[100 100 1200 800]);
            fig= gcf;
            fig.Color=[1,1,1];
            % Plot all channels
            h = semilogy(freq,po,'LineWidth',1);
            try
                set(h, {'color'},num2cell(colormap123,2));
            catch
            end
            hold on;
            % Plot the mean in black
            semilogy(freq,squeeze(mean(po,2)),'-k','LineWidth',2);
            hold on
            xp2 =0:round(freq(end));
            yp2=ones(1,round(freq(end))+1)*15;
            p2 =plot(xp2,yp2,'--k');
            p2.LineWidth=2;
            grid on
            ax = gca; % current axes
            ax.FontSize = 20;
            ax.TickLength = [0.02 0.02];
            
            xlabel('Frequency (Hz)','FontSize',25)
            labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
            ylabel(labY,'interpreter','latex','FontSize',25)
            
            % Adjust limits based on cfg.foi
            xlim([cfg.foi(1), cfg.foi(end)]);
            
            % Plot legend
            if strcmp(cfg.plot_legend,'yes')
                if length(label) > 24
                    [rrr,object_h] = ...
                        columnlegend(2, vertcat(label, 'mean'),...
                        'Location','northeastoutside','FontSize',14);
                    
                    % Fix Me: change FontSize?
                    
                    hl = findobj(object_h,'type','line');
                    h2 = findobj(object_h,'type','text');
                    set(hl,'LineWidth',5);
                    
                else
                    % Legend
                    [~, hobj, ~, ~] = legend(vertcat(label, 'mean'),...
                        'location','eastoutside');
                    hl = findobj(hobj,'type','line');
                    set(hl,'LineWidth',4);
                    ht = findobj(hobj,'type','text');
                    set(ht,'FontSize',12);
                end
            end
        end
        
        %4. MATLAB's inbuilt pwelch
    case 'pwelch'
        
        if isempty(cfg.trial_length)
            numsamp = length(x);
        else 
            numsamp = rawData.fsample*cfg.trial_length;
        end
        
        x = rawData.trial{1};
        [pow,freq] = pwelch(x',numsamp,[],[],rawData.fsample);
        
        % Is this correct?
        pow = pow';
        
        strt = find(freq > cfg.foi(1),1,'first');
        stp  = find(freq < cfg.foi(2),1,'last');

        % Funky colorscheme, looks OK...
        colormap123     = linspecer(length(label));

         if strcmp(cfg.plot,'yes')
            % Make a Figure
            figure()
            set(gcf,'Position',[100 100 1200 800]);
            
            h = plot(freq(strt:stp),log10(pow(:,strt:stp)),...
                'LineWidth',1);
            set(h, {'color'},num2cell(colormap123,2));
            hold on;
            plot(freq(strt:stp),log10(mean(pow(:,strt:stp),1)),'-k','LineWidth',2);
            grid on
            ax = gca; % current axes
            ax.FontSize = 20;
            ax.TickLength = [0.02 0.02];
            ylabel('PSD (dB/Hz)','FontSize',30);
            xlabel('Frequency (Hz)','FontSize',30);
            % Legend
            [~, hobj, ~, ~] = legend(vertcat(label, 'mean'),'location','eastoutside');
            hl = findobj(hobj,'type','line');
            set(hl,'LineWidth',4);
            ht = findobj(hobj,'type','text');
            set(ht,'FontSize',12);
         end
        
end


end




