function lay = makeTriAxialLayout(slotLayout,slot2sens,plot)
% Applies sensor labels to a layout file setup labelled with slot numbers.
% The output is a structure with four layouts: one with all sensors and
% then each of the three possible axes labels independent of each other. 
%
% Example usage:
% cfg			= [];
% cfg.image		= 'D:\data\..BIDS..\layoutImage.png';
% slotLayout	= ft_prepare_layout(cfg);
%
% % Follow fieldtrip procedure of labelling parts of the image, drawing the
% % outline etc.
% 
% save('slotLayout.mat','slotLayout');
% 
% % At this point you have a layout file with 1-N channel positions and
% % labels saying 1-N. Use slot2sens to work that out
% slot2sens = readtable('D:\data\FilburyOPM\P01\sub-OP00057\slot2sens.csv');
% plot		= false;
% lay		= makeTriAxialLayout(slotLayout,slot2sens,plot);
% % save('sensorLayout.mat','sensorLayout');
% % ft_plot_layout(lay.Z) % Check Z axis
% 
% Author:		Nicholas Alexander	(n.alexander@ucl.ac.uk)

%% Create a layout for each axis, and one with all of them. Get the outline
%  etc.
lay					= [];
lay.all.mask		= slotLayout.mask;
lay.all.outline		= slotLayout.outline;
lay.X				= lay.all;
lay.Y				= lay.all;
lay.Z				= lay.all;

% Count for loop
labCount.all	= 0;
labCount.X		= 0;
labCount.Y		= 0;
labCount.Z		= 0;

% Go through the unlabelled positions and asign them to the sensors in the
% slots. 
for labIdx = 1:length(slot2sens.sensor)
	% Asign axes appropriately for sensor type
	if strcmp(slot2sens.sensor{labIdx}(1:2), 'G3')
		oris = {'-Y','-Z','-X'};
	elseif strcmp(slot2sens.sensor{labIdx}(1:2), 'G2')
		oris = {'-Y','-Z'};
	else
		error("Somethings up!");
	end

	% Go through the orientations one by one
	for oriIdx = 1:length(oris)
		% Top up the local indexing
		labCount.all = labCount.all + 1;
		labCount.(oris{oriIdx}(2))		= labCount.(oris{oriIdx}(2)) + 1;
		
		% Asign the right values for all axes
		lay.all.label(labCount.all, 1)	= strcat(slot2sens.sensor{labIdx},oris(oriIdx));
		lay.all.width(labCount.all,1)	= slotLayout.width(labIdx,1);
		lay.all.height(labCount.all,1)	= slotLayout.height(labIdx,1);
		lay.all.pos(labCount.all, 1:2)	= slotLayout.pos(labIdx,1:2);

		% And for the specific axis
		lay.(oris{oriIdx}(2)).label(labCount.(oris{oriIdx}(2)), 1)	= strcat(slot2sens.sensor{labIdx},oris(oriIdx));
		lay.(oris{oriIdx}(2)).width(labCount.(oris{oriIdx}(2)),1)	= slotLayout.width(labIdx,1);
		lay.(oris{oriIdx}(2)).height(labCount.(oris{oriIdx}(2)),1)	= slotLayout.height(labIdx,1);
		lay.(oris{oriIdx}(2)).pos(labCount.(oris{oriIdx}(2)),1:2)	= slotLayout.pos(labIdx,1:2);
	end
end

% Plot to check
if (plot)
	subplot(1,3,1);
	ft_plot_layout(lay.X)
	subplot(1,3,2);
	ft_plot_layout(lay.Y)
	subplot(1,3,3);
	ft_plot_layout(lay.Z)
end