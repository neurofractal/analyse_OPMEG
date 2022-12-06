function lay = addOutlineToLayout(slotLayout,scale,plot)
% Draws a circle as the head boundary around an existing layout and adds a 
% nose and ears, adds a nose and ears. Scale should be ~1 and is used to
% shrink or inflate the outline as required.

% Clear any outlines/masks added
slotLayout.outline = [];
slotLayout.mask = [];

% Get centre position
xMax = max(slotLayout.pos(:,1));
xMin = min(slotLayout.pos(:,1));
yMax = max(slotLayout.pos(:,2));
yMin = min(slotLayout.pos(:,2));
% Make an outline/mask circle for the head bounds.
radius	= scale*(max([yMax-yMin, xMax-xMin])/2);
samples = 128;
theta	= linspace(0,2*pi,samples);
x		= radius*cos(theta);
y		= radius*sin(theta);
slotLayout.outline{1} = [x; y]';


% Add a nose to the outline
noseScale = radius*1.1;
noseHalfWidth = 0.07;
[~, index] = min(abs(theta-noseHalfWidth)); % Index for theta closest to 0.1 radians
slotLayout.outline{1,2}(1,1:2) = slotLayout.outline{1}((samples/4) - (index - 1),:); % left nostril
slotLayout.outline{1,2}(2,1:2) = [0, noseScale];
slotLayout.outline{1,2}(3,1:2) = slotLayout.outline{1,2}(1,1:2) .* [-1 1]; % right nostril
% slotLayout.outline{1,2}(4,1:2) = slotLayout.outline{1,2}(1,1:2); % complete

% Add left ear
earHalfWidth = 0.16;
[~, index] = min(abs(theta-earHalfWidth)); % Index for theta closest to 0.1 radians
slotLayout.outline{1,3}(1,1:2) = slotLayout.outline{1}((samples/2) - (index - 1),:); % top start
slotLayout.outline{1,3}(2,1:2) = [-radius*1.037, radius*0.235];
slotLayout.outline{1,3}(3,1:2) = [-radius*1.07, radius*0.2];
slotLayout.outline{1,3}(4,1:2) = [-radius*1.08, radius*0.165];

slotLayout.outline{1,3}(5,1:2) = [-radius*1.05, radius*0.02]; % left middle
slotLayout.outline{1,3}(6,1:2) = [-radius*1.09, -radius*0.15];

slotLayout.outline{1,3}(7,1:2) = [-radius*1.075, -radius*0.23];
slotLayout.outline{1,3}(7,1:2) = [-radius*1.065, -radius*0.21];
slotLayout.outline{1,3}(8,1:2) = slotLayout.outline{1,3}(1,1:2) .* [1 -1]; % bottom end

% slotLayout.outline{1,3}(8,1:2) = slotLayout.outline{1,3}(1,1:2); % complete

% Add right ear
slotLayout.outline{1,4} = slotLayout.outline{1,3} .* [-1 1];

% Apply offsets
xMean = mean(slotLayout.pos(:,1));
yMean = mean(slotLayout.pos(:,2));

slotLayout.outline{1} = slotLayout.outline{1} + [xMean, yMean];
slotLayout.outline{2} = slotLayout.outline{2} + [xMean, yMean];
slotLayout.outline{3} = slotLayout.outline{3} + [xMean, yMean];
slotLayout.outline{4} = slotLayout.outline{4} + [xMean, yMean];

slotLayout.mask{1} = slotLayout.outline{1};

if plot
	ft_plot_layout(slotLayout)
end

lay = slotLayout;
