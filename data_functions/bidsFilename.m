function [filename, folder] = bidsFilename(cfg,bids)
% Function to create a BIDS style filename. It is useful for saving and
% loading data. Use as:
% Usual bids info as structure containing strings:
% bids.sub		= 'pilot01';
% bids.ses		= '002';
% bids.run		= '003';
% bids.task		= 'oddball';
% bids.directory	= 'D:\data\MyBidsStudy'; % Folder must exist already
%
% Config for variations on BIDS format:
% cfg.category		= 'meg'; % BIDS category. If derivative, this is the current pipeline.
% cfg.description	= 'test'; % Describe the file. e.g. 'freqData' or 'meg'.
% cfg.type			= '.mat'; % File type or '.bin'. Default '.mat'.
% cfg.derivative	= true; % or false. Reference derivative folder?
% cfg.detailed		= true; % or false. Include ses/task/run info?
% If cfg is empty, or not enough info is provided only a folder will be output
% 
% Author: Nicholas Alexander, n.alexander@ucl.ac.uk

%% Check input
% cfg
folderOnly = false;
if ~isstruct(cfg)
	disp('cfg is empty. Only providing folder output.');
	folderOnly = true;
end
if folderOnly
	cfg.category = 'folderOnly';
elseif	(~isfield(cfg,'category') || isempty(cfg.category))
	disp('cfg.category not provided. Only returning folder output');
	cfg.category = 'folderOnly';
	folderOnly = true;
end
if folderOnly
	cfg.description = 'folderOnly';
elseif (~isfield(cfg,'description') || isempty(cfg.description))
	disp('cfg.description is empty. Only returning folder output.')
	cfg.description = 'folderOnly';
	folderOnly = true;
end
if folderOnly
	cfg.type = '.mat';
elseif (~isfield(cfg,'type') || isempty(cfg.type))
	disp('Defaulting cfg.type to .mat')
	cfg.type = '.mat';
end
if (~isfield(cfg,'derivative') || isempty(cfg.derivative))
	disp('Defaulting cfg.derivative to true');
	cfg.derivative = true;
end
if (~folderOnly && (~isfield(cfg,'detailed') || isempty(cfg.detailed)))
	disp('Defaulting cfg.detailed to false');
	cfg.detailed = false;
end

% bids
if ~isstruct(bids)
	error('bids is empty. Please provide at least bids.directory and bids.sub');
elseif ((~isfield(bids,'directory') || isempty(bids.directory)) || (~isfield(bids,'sub') || isempty(bids.sub)))
	error('Please provide both bids.directory and bids.sub');
end

if (~isfield(bids,'ses') || isempty(bids.ses))
	cfg.detailed = false;
	if ~cfg.derivative
		error("Must provide bids.ses if cfg.derivative = false.")
	end
end
if (~isfield(bids,'task') || isempty(bids.task))
	cfg.detailed = false;
end
if (~isfield(bids,'run') || isempty(bids.run))
	cfg.detailed = false;
end

%% Produce output filename
% Check if the initial directory exists. 
if isfolder(bids.directory)
	% Check the backslash at the end
	if bids.directory(end) ~= '\'
		bids.directory = strcat(bids.directory,'\');
	end
	% If using the derivative folder, add that on, including sub
	if cfg.derivative
		if ~folderOnly
			bids.directory = sprintf('%1$sderivatives\\%2$s\\sub-%3$s\\',bids.directory,cfg.category,bids.sub);
		else
			bids.directory = sprintf('%1$sderivatives\\sub-%3$s\\',bids.directory,bids.sub);
		end

		% Check that directory exists with the derivative folder
		if ~isfolder(bids.directory)
			warning("Creating new derivatives directory")
			mkdir(bids.directory);
		end

	% Otherwise just add the subject folder
	else
		if ~folderOnly
			bids.directory = sprintf('%1$ssub-%2$s\\ses-%3$s\\%4$s\\',bids.directory,bids.sub,bids.ses,cfg.category);
		else
			bids.directory = sprintf('%1$ssub-%2$s\\ses-%3$s\\',bids.directory,bids.sub,bids.ses);
		end

		% This function should not interfere with the original data folder.
		if ~isfolder(bids.directory)
			error("Specified BIDS folder does not exist.");
		end
	end
else
	error("Specified BIDS directory does not exist");
end

% Now create the filename with required detail
if cfg.detailed
	filename = sprintf('%1$ssub-%2$s_ses-%3$s_task-%4$s_run-%5$s_%6$s%7$s',bids.directory,bids.sub,bids.ses,bids.task,bids.run,cfg.description,cfg.type);
else
	filename = sprintf('%1$ssub-%2$s_%3$s%4$s',bids.directory,bids.sub,cfg.description,cfg.type);
end

% Get the folder separately
[folder,~,~] = fileparts(filename);

% If not enough information was provided in the cfg, set filename to folder.
if folderOnly
	filename = folder;
end



















