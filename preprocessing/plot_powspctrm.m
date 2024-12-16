function plot_powspctrm(po,cfg,label,freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the PSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Try loading a fancy colorscheme
try
    colormap123     = linspecer(length(label));
catch
    disp('Using default colorscheme')
end

% Calculate STD
stdev_po = std(po');

% Calculate mean
mean_po = squeeze(mean(po',1));

figure()
set(gcf,'Position',[100 100 1200 800]);
fig= gcf;
fig.Color=[1,1,1];

if strcmp(cfg.plot_chans,'yes')
    % Plot all channels
    h = plot(freq,po','LineWidth',0.5);
    set(h, {'color'},num2cell(colormap123,2));
    
    % Change the transparency of the lines
    if cfg.transparency ~= 1
        for i = 1:size(po,2)
            h(i).Color(4) = cfg.transparency;
        end
    end
end

hold on;
% Plot the mean in black
if strcmp(cfg.plot_mean,'yes')
    plot(freq,mean_po,'-k','LineWidth',2);
    hold on
end

xp2 =0:round(freq(end));
yp2=ones(1,round(freq(end))+1)*15;
p2 =plot(xp2,yp2,'--k');
p2.LineWidth=2;
grid on
ax = gca; % current axes
ax.FontSize = 20;
ax.TickLength = [0.02 0.02];


% Plot the CI in gray
if strcmp(cfg.plot_ci,'yes')
    plot(freq,mean_po+stdev_po,'--','LineWidth',1,...
        'Color',[0.5 0.5 0.5]); hold on;
    plot(freq,mean_po-stdev_po,'--','LineWidth',1,...
        'Color',[0.5 0.5 0.5]); hold on;
end

%             % Plot Confidence Interval
%             ciplot(mean_po-stdev_po,mean_po+stdev_po,freq,'k',0.9);

set(gca, 'YScale', 'log'); hold on;

xlabel('Frequency (Hz)','FontSize',25)
labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex','FontSize',25)

% Adjust limits based on cfg.foi
xlim([cfg.foi(1), cfg.foi(end)]);

% Plot legend
if strcmp(cfg.plot_chans,'yes')
    if strcmp(cfg.plot_legend,'yes')
        if length(label) > 24
            [rrr,object_h] = ...
                columnlegend(3, vertcat(label, 'mean'),...
                'Location','northeastoutside','FontSize',10);
            
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
    
    if strcmp(cfg.interactive,'yes')
        % Assign unique labels to each line
        for i = 1:length(h)
            h(i).UserData = label{i}; % Customize as needed
        end

        disp(h(1).UserData)

        % Enable data cursor mode
        dcm = datacursormode(gcf);
        dcm.Enable = 'on';
        dcm.UpdateFcn = @(~, event) customLabel(event);
    end
end
end

function txt = customLabel(event)
    disp('CustomLabel called'); % Debug message
    targetLine = event.Target;

    if isprop(targetLine, 'UserData') && ~isempty(targetLine.UserData)
        label = targetLine.UserData; % Retrieve the label from UserData
    else
        label = 'No Label'; % Default if UserData is empty
    end

    % Adjust Y value for logarithmic scale
    x = event.Position(1); % X-coordinate
    y = event.Position(2); % Y-coordinate (already in log scale for display)

    % Format the displayed text
    txt = sprintf('Label: %s\nX: %.2f Hz\nY: %.2f (log)', label, x, y);
end

