function [comp2remove] = ft_RS_ICviewer(cfg,comp,data_ICA)
% Function to display ICs and select which to remove
%
% EXAMPLE USEAGE:   [comp2remove] = ft_RS_ICviewer(cfg,comp)
% ...where, cfg is the input structure and comp is the output from
% ft_componentanalysis
%   cfg.layout          = lay;
%   cfg.foi             = frequencies of interest in form [X Y].
%                       (default = [1 100]).
%   cfg.save            = save a .png of removed comps? (default = 'no')
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)   
%__________________________________________________________________________

%% Function housekeeping
if ~isfield(cfg, 'foi')
    cfg.foi = [1 100];
end

if ~isfield(cfg, 'save')
    cfg.save = 'no';
end

if ~isfield(cfg, 'layout')
    warning('Will not work without a layout file')
end

%% Start of function proper
% Calculate PSD
cfg2                 = [];
cfg2.channel         = 'all';
cfg2.trial_length    = 10;
cfg2.method          = 'tim';
cfg2.foi             = [cfg.foi(1) cfg.foi(2)];
cfg2.plot            = 'no';
[pow freq]          = ft_opm_psd(cfg2,data_ICA);
po                  = nanmean(pow(:,:,:),3);

% Use Brewermap :colors RdBu
ft_hastoolbox('brewermap',1);
colormap123 = colormap(flipud(brewermap(64,'RdBu')));

comp2remove = [];
c = 1;
cont = 1;

% Create Figure
S.f = figure;

while cont

    set(gcf, 'Position',  [300, 0, 1100, 900]);
    %create two pushbttons
    S.pb = uicontrol(S.f,'style','togglebutton',...
        'units','pix',...
        'position',[850 300 200 40],...
        'backgroundcol',[1 0 0],...
        'fontsize',14,...
        'Tag','REMOVE',...
        'string','REMOVE',...
        'callback',@pb_call2);

    S.pb = uicontrol(S.f,'style','push',...
        'units','pix',...
        'position',[650 300 200 40],...
        'backgroundcol',[0.2824 1 0],...
        'fontsize',14,...
        'Tag','flip_button',...
        'string','KEEP',...
        'callback',@pb_call);

        S.pb = uicontrol(S.f,'style','push',...
        'units','pix',...
        'position',[0 300 150 40],...
        'backgroundcol',[1 1 1],...
        'fontsize',14,...
        'Tag','EXIT',...
        'string','EXIT',...
        'callback',@pb_call3);

    % Find limits of y-axis
    [minDistance, indexOfMin] = min(abs(freq-cfg.foi(1)));
    [minDistance, indexOfMax] = min(abs(freq-cfg.foi(2)));
    max_lim = max(po(indexOfMin:indexOfMax,c))*1.1;
    min_lim = min(po(indexOfMin:indexOfMax,c))*1.1;

    % Plot topoplot
    cfg2 = [];
    cfg2.component = c;       % specify the component(s) that should be plotted
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
    plot(comp.time{1},comp.trial{1}(c,1:end),'linewidth',2);
    xlim([1 11]);
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    xlabel('Time(s)','FontSize',23);
    ylabel({'Magnetic';'Field'},'FontSize',23);

        % Plot PSD
    subplot(8,4,[1:2 5:6 9:10 13:14]);
    semilogy(freq,po(:,c),'-k','LineWidth',2);
    xlim([1 80]);
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    xlabel('Frequency (Hz)','FontSize',23)
    labY = ['$$PSD (a.u.) $$)'];
    ylabel(labY,'interpreter','latex','FontSize',23)
    ylim([min_lim max_lim]);
    %print(['component_PSD' num2str(c)],'-dpng','-r300');



    % Wait for user input
    uiwait(S.f)
    
    if strcmp(S.f.UserData,'EXIT')
        cont = 0;
    elseif strcmp(S.f.UserData,'REMOVE')
        comp2remove = vertcat(comp2remove,c);
        disp('');
        disp(['Component ' num2str(c) ' = BAD']);
        disp('');

         % Save a picture
    if strcmp(cfg.save, 'yes')
        print(['component_' num2str(c)],'-dpng','-r300');
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