function [comp2remove] = ft_RS_ICviewer(cfg,comp)
% Function to display ICs with a GUI, with options for KEEP/REMOVE.
% Output is an array of the comps to be rejected.
%
% Currently tested on MATLAB 2021b, Fieldtrip v20210606. Will work best for
% continuous datasets, rather than datasets already broken up into 
% shorter epochs.
%
% EXAMPLE USEAGE:   [comp2remove] = ft_RS_ICviewer(cfg,comp)
% ...where, cfg is the input structure and comp is the output from
% ft_componentanalysis
%   cfg.layout          = layout file created by ft_prepare_layout
%   cfg.foi             = frequencies of interest in form [X Y] for PSD
%                       (default = [1 100]).
%   cfg.trial_length    = size of chunks used for PSD calculation
%                       (default = 10)
%   cfg.sort            = statistical property to sort presentation of ICs
%                         (['no'],'kurtosis','autocorrelation')
%   cfg.save            = save a .png of removed comps? (default = 'no')
%   
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)   
%__________________________________________________________________________

%% Function housekeeping
if ~isfield(cfg, 'foi')
    cfg.foi = [1 100];
end

if ~isfield(cfg, 'trial_length')
    cfg.trial_length = 10;
end

if ~isfield(cfg, 'save')
    cfg.save = 'no';
end

if ~isfield(cfg, 'layout')
    warning('Will not work without a layout file')
end

if ~isfield(cfg, 'sort')
    cfg.sort = 'no';
end

%% Start of function proper
% Calculate PSD

cfg2                 = [];
cfg2.channel         = 'all';
cfg2.trial_length    = cfg.trial_length;
cfg2.foi             = [cfg.foi(1) cfg.foi(2)];
cfg2.plot            = 'no';
[pow, freq]          = calc_PSD(cfg2,comp);

% Make sure there are no nans (and take mean)
po                  = nanmean(pow(:,:,:),3);

%% Sort by various statistical properties of the ICs
if strcmp(cfg.sort,'kurtosis')
    disp('Sorting componenets by kurtosis value')
    % Kurtosis
kurt = kurtosis(comp.trial{1},[],2);
kurt = kurt(:);
metrics.kurtosis.raw = kurt-3;
% metrics.kurtosis.log = log(kurt-3);
% metrics.kurtosis.abs = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis

[~,sorted_comps] = sort(metrics.kurtosis.raw,'descend');

elseif strcmp(cfg.sort,'autocorrelation')
    disp('Sorting components by autocorrelation with cardiac signal')
    bpm_range = [50 100]./60;
    maxlags = 2*comp.fsample;
    [~,lags] = xcorr(comp.trial{1}(1,:),maxlags,'coeff');
    lags = lags./comp.fsample;
    lags_range = lags>bpm_range(1) & lags<bpm_range(2);
    ac_max = zeros(1,length(comp.label));
    for ic = 1:length(comp.label)
        [tc_bp] =  ft_preproc_lowpassfilter(comp.trial{1}(ic,:), comp.fsample, 48);
        ac = xcorr(tc_bp,maxlags,'coeff');
        ac_max(ic) = max(ac(lags_range));
    end
    metrics.cardiac_autocorrelation.value = ac_max(:);
    [~,sorted_comps] = sort(metrics.cardiac_autocorrelation.value,'descend');

else
    sorted_comps = [1:length(comp.label)];
end

% Use Brewermap :colors RdBu
ft_hastoolbox('brewermap',1);
colormap123 = colormap(flipud(brewermap(64,'RdBu')));

comp2remove = [];
c = 1;
cont = 1;

% Create Figure
S.f = figure;

while cont

    set(gcf, 'Position',  [300, 100, 1100, 900]);

    %create two pushbttons
    S.pb = uicontrol(S.f,'style','togglebutton',...
        'units','pix',...
        'position',[850 300 200 40],...
        'backgroundcol',[1 0 0],...
        'fontsize',14,...
        'Tag','REMOVE',...
        'string','REMOVE',...
        'callback',@pb_remove);

    S.pb = uicontrol(S.f,'style','push',...
        'units','pix',...
        'position',[650 300 200 40],...
        'backgroundcol',[0.2824 1 0],...
        'fontsize',14,...
        'Tag','flip_button',...
        'string','KEEP',...
        'callback',@pb_keep);

        S.pb = uicontrol(S.f,'style','push',...
        'units','pix',...
        'position',[20 300 100 40],...
        'backgroundcol',[1 1 1],...
        'fontsize',14,...
        'Tag','EXIT',...
        'string','EXIT',...
        'callback',@pb_exit);

    % Find limits of y-axis
    [minDistance, indexOfMin] = min(abs(freq-cfg.foi(1)));
    [minDistance, indexOfMax] = min(abs(freq-cfg.foi(2)));
    max_lim = max(po(indexOfMin:indexOfMax,sorted_comps(c)))*1.1;
    min_lim = min(po(indexOfMin:indexOfMax,sorted_comps(c)))*1.1;

    % Plot topoplot
    cfg2 = [];
    cfg2.component = sorted_comps(c);       % specify the component(s) that should be plotted
    cfg2.layout    = cfg.layout; % specify the layout file that should be used for plotting
    cfg2.comment   = 'no';
    cfg2.figure    = S.f;
    cfg2.marker    = 'off';
    cfg2.colorbar  = 'EastOutside';
    cfg2.zlim      = 'maxabs';
    cfg2.highlight = 'off';
    cfg2.highlightfontsize = 1;
    cfg2.colormap  = colormap123;
    subplot(8,4,[3:4 7:8 11:12 15:16]);
    ft_topoplotIC(cfg2, comp)
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18

    % Plot the raw data from 1-11s
    subplot(8,4,[28:32]);
    plot(comp.time{1},comp.trial{1}(sorted_comps(c),1:end),'linewidth',2);
    xlim([1 11]);
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    xlabel('Time(s)','FontSize',23);
    ylabel({'Magnetic';'Field'},'FontSize',23);

    % Plot PSD
    subplot(8,4,[1:2 5:6 9:10 13:14]);
    semilogy(freq,po(:,sorted_comps(c)),'-k','LineWidth',2);
    xlim([cfg.foi(1) cfg.foi(2)]);
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    xlabel('Frequency (Hz)','FontSize',23)
    labY = ['$$PSD (a.u.) $$'];
    ylabel(labY,'interpreter','latex','FontSize',23)
    ylim([min_lim max_lim]);
    
    annotation('textbox', [0.21,0.332222223248747,...
        0.294090900597247,0.045555554529031],...
        'string', ['Plotted Component: ' num2str(c) ' of ' ...
        num2str(length(comp.label))],'FontSize',14,'EdgeColor','None');


    % Wait for user input
    uiwait(S.f)
    
    % Perform various actions based on the button pressed
    if strcmp(S.f.UserData,'EXIT')
        cont = 0;
    elseif strcmp(S.f.UserData,'REMOVE')
        comp2remove = vertcat(comp2remove,sorted_comps(c));
        disp('');
        disp(['Component ' num2str(sorted_comps(c)) ' = BAD']);
        disp('');

        % Save a picture if required
        if strcmp(cfg.save, 'yes')
            print(['component_' num2str(sorted_comps(c))],'-dpng','-r100');
        end

    end

%     % Clear the Figure for the next IC
    clf(S.f);
    
    % Don't continue if max components reached
    c = c+1;
    if c > length(comp.label)
        cont = 0;
    end
    
    % Close Figure if cont == 0
    if ~cont
        disp('Closing Figure...')
        close(S.f);
    end
   

end

%% Sub-functions
function [pow, freq] = calc_PSD(cfg,rawData)
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
end

function pb_keep(hObject,~)
    hObject.Parent.UserData = 'moveon';
    uiresume;
    return
    end

function pb_remove(hObject,~)
    disp('');
    %At the end of the callback function:
    hObject.Parent.UserData = 'REMOVE';
    uiresume;
    return
end

function pb_exit(hObject,~)
    disp('EXIT');
    hObject.Parent.UserData = 'EXIT';
    uiresume;
    return
end

end
