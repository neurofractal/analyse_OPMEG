function [filename] = bidsFilename(bids,description,type,derivative,detailed)
% Function to create a BIDS style filename. It is useful for saving and
% loading data. Use as:
% Usual bids info as structure containing strings
% bids.sub		= 'pilot01';
% bids.ses		= '002';
% bids.run		= '003';
% bids.task		= 'oddball';
% bids.directory	= 'D:\data\MyBidsStudy'; % Folder must exist already
% description		= 'test'; % Describe the file. e.g. 'freqData'
% type			= '.mat'; % Usually a .mat
% derivative		= true; % Place in derivative folder or main folder?
% detailed		= true; % Include ses,run,task info?
%
% Author: Nicholas Alexander, n.alexander@ucl.ac.uk

% Check if the initial directory exists. 
if isfolder(bids.directory)
	% Check the backslash at the end
	if bids.directory(end) ~= '\'
		bids.directory = strcat(bids.directory,'\');
	end
	% If using the derivative folder, add that on, including sub
	if derivative
		bids.directory = sprintf('%1$sderivatives\\sub-%2$s\\',bids.directory,bids.sub);

		% Check that directory exists with the derivative folder
		if ~isfolder(bids.directory)
			warning("Creating new derivatives directory")
			mkdir(bids.directory);
		end
	else
		% Otherwise just add the subject folder
		bids.directory = sprintf('%1$ssub-%2$s\\',bids.directory,bids.sub);

		if ~isfolder(bids.directory)
			error("Specified BIDS participant folder does not exist");
		end

	end
else
	error("Specified BIDS directory does not exist");
end

% Now create the filename with required detail
if detailed
	filename = sprintf('%1$ssub-%2$s_ses-%3$s_task-%4$s_run-%5$s_%6$s%7$s',bids.directory,bids.sub,bids.ses,bids.task,bids.run,description,type);
else
	filename = sprintf('%1$ssub-%2$s_%3$s%4$s',bids.directory,bids.sub,description,type);
end