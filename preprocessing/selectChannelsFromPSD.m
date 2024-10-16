function [selectedLabel, selectedIdx] = selectChannelsFromPSD(cfg, pow, freq, label)
    % Reduce the third dimension (e.g., by taking the mean across time)
	if length(size(pow)) > 2
		if strcmp(cfg.average,'median')
			pow = median(pow, 3); %
		elseif strcmp(cfg.average,'mean')
			pow = mean(pow, 3); %
		else
			disp("Default to median average");
			pow = mean(pow, 3); %
		end
		pow = mean(pow, 3); %
	end
	
    % Calculate median across channels
    medianPo = median(pow, 2); % Initial mean of all channels

    % Try loading a fancy colour scheme
    try
        colormap = linspecer(length(label)); % Custom colour scheme
    catch
        disp('Using default colour scheme');
        colormap = lines(length(label)); % Default colours
    end

    % Make a figure
    fig = figure();
    fig.Color = [1, 1, 1];

    % Initialize the plot handles array
    numChannels = size(pow, 2);
    h = gobjects(numChannels, 1); % Preallocate graphics object array
    selectedIdx = []; % Array to store selected line indices

    % Plot all channels
    h = plot(freq, pow, 'LineWidth', 0.5);
    set(h, {'color'}, num2cell(colormap, 2));
    hold on;

	for i = 1:numel(h)
		set(h(i), 'UserData', false); % Initialize UserData to false
	end

    % Plot the initial mean in black
    medianLine = plot(freq, medianPo, '-k', 'LineWidth', 2);

    % Plot a threshold line at 15
    xp2 = 0:round(freq(end));
    yp2 = ones(1, round(freq(end)) + 1) * 15;
    p2 = plot(xp2, yp2, '--k');
    p2.LineWidth = 2;

    % Set axis properties
    grid on;
    ax = gca;
    ax.FontSize = 20;
    ax.TickLength = [0.02, 0.02];
    set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)', 'FontSize', 25);
    ylabel('PSD (fT/sqrt(Hz))', 'interpreter', 'latex', 'FontSize', 25);

    % Adjust x-axis limits based on cfg.foi
    xlim([cfg.foi(1), cfg.foi(end)]);

    % Set up selection mechanism
    set(fig, 'WindowButtonDownFcn', @startSelection);

    % Wait for user to interact and close the figure
    waitfor(fig);

    % Get labels of selected channels
    if ~isempty(selectedIdx)
        selectedLabel = label(selectedIdx);
    else
        selectedLabel = {};
        selectedIdx = [];
    end

    % Nested function to start selection
    function startSelection(~, ~)
        % Get initial point of the selection
        point1 = get(ax, 'CurrentPoint');
        rbbox; % Draw selection box (rubberband box)
        point2 = get(ax, 'CurrentPoint');

        % Get the x and y limits of the selection box
        xLimits = [min(point1(1,1), point2(1,1)), max(point1(1,1), point2(1,1))];
        yLimits = [min(point1(1,2), point2(1,2)), max(point1(1,2), point2(1,2))];

        % Check which lines are inside the box
        for i = 1:numChannels
            xData = get(h(i), 'XData');
            yData = get(h(i), 'YData');
            
            % Find any points within the x and y limits of the selection box
            inX = (xData >= xLimits(1)) & (xData <= xLimits(2));
            inY = (yData >= yLimits(1)) & (yData <= yLimits(2));

            % If any point of the line is within the selection box, toggle selection
            if any(inX & inY)
                if get(h(i), 'UserData') == false
                    set(h(i), 'UserData', true);  % Mark as selected
					h(i).Color(4) = 0;
                    selectedIdx = [selectedIdx, i]; % Add to selected list
% 				else
%                     set(h(i), 'UserData', false); % Mark as not selected
% 					h(i).Color(4) = 1;
%                     selectedIdx(selectedIdx == i) = []; % Remove from selected list
                end
                % Recalculate the mean excluding selected lines
                updateMeanLine();
            end
        end
    end

    % Nested function to update the mean line
    function updateMeanLine()
        % Exclude selected channels from the mean calculation
        unselectedIdx = setdiff(1:numChannels, selectedIdx);
        newMean = mean(pow(:, unselectedIdx), 2); % Recalculate mean
        set(medianLine, 'YData', newMean); % Update the mean line
    end
end
