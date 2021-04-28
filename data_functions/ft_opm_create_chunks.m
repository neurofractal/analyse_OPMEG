function [rawData] = ft_opm_create_chunks(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WARNING: USE THIS AT YOUR OWN RISK - NOT FULLY TESTED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to read optically-pumped magnetencephalography (OPMEG) data
% acquired from the UCL Wellcome Centre for Neuroimaging in chunks of Ns
% and save to disk
%
% Hopefully this will be useful for very long OPM datasets.
%
% EXAMPLE USEAGE:   data = ft_opm_create_chunks(cfg)
% ...where, cfg is the input structure
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
% cfg.chunksize       = number in seconds for the length of the chunk
%                       (default = 300s)
% cfg.downsample      = 'yes' or 'no' to 1000Hz (default = 'no')
%
%%%%%%%%%%%
% Option 1
%%%%%%%%%%%
%   cfg.data          = path to raw .bin file
%                       (requires .json and _channels.tsv in same folder
%                       with same name as the .bin file)
%%%%%%%%%%%
% Option 2
%%%%%%%%%%%
%   cfg.folder          = path to folder containing data
%                       organised in BIDS format
%   cfg.bids.task       = 'AEF'
%   cfg.bids.sub        = '001'
%   cfg.bids.ses        = '003'
%   cfg.bids.run        = '002'
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_create (Tim Tierney)

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%__________________________________________________________________________

%% Set default values
if ~isfield(cfg, 'fs')
    cfg.fs   = 6000;
end
if ~isfield(cfg, 'nSamples')
    cfg.nSamples   = 1000;
end

if ~isfield(cfg, 'precision')
    cfg.precision = 'double';
end

if ~isfield(cfg, 'chunksize')
    cfg.chunksize = 300;
end

if ~isfield(cfg, 'downsample')
    cfg.downsample = 'no';
end

%% Determine whether to use cfg.data or cfg.folder
if ~isfield(cfg, 'data')
    use_bids = 1;
else
    path_to_bin_file = cfg.data;
    use_bids = 0;
end

% If the user has specified both cfg.data and cfg.folder default to
% cfg.data, but warn the user
if isfield(cfg, 'data') && isfield(cfg, 'folder')
    ft_warning(['Both cfg.data and cfg.folder have been supplied.'...
        ' Defaulting to cfg.data']);
end

%% If using option 2 = BIDS!
if use_bids
    try
        file_name_bids = ['sub-' cfg.bids.sub '_ses-' cfg.bids.ses ...
            '_task-' cfg.bids.task '_run-' cfg.bids.run '_meg.bin'];
    catch
        error('Did you specify all the required cfg.bids information?')
    end
    path_to_bin_file = fullfile(cfg.folder,['sub-' cfg.bids.sub],...
        ['ses-' cfg.bids.ses],'meg',file_name_bids);
end


%% Read Binary File
try % to read data
    [direc, dataFile] = fileparts(path_to_bin_file);
    dat         = fopen(path_to_bin_file);
    %    data_raw    = fread(dat,Inf,cfg.precision,0,'b');
    %fclose(dat);
    binData     = 1;
catch % if not readable check if it is numeric
    if ~isa(path_to_bin_file,'numeric') % if not numeric throw error
        error(['Cannot read the file: ' path_to_bin_file])
    end
    binData     = 0;
    direc       = pwd();
    dataFile    = cfg.fname;
end

%%
% Identify potential BIDS Files
base        = strsplit(dataFile,'meg');
chanFile    = fullfile(direc,[base{1},'channels.tsv']);
megFile     = fullfile(direc,[base{1},'meg.json']);
posFile     = fullfile(direc,[base{1},'positions.tsv']);
coordFile   = fullfile(direc,[base{1},'coordsystem.json']);

%% Load a channels file.
try
    channels    = ft_read_chan_tsv(chanFile);
catch
    error(['Cannot read the file: ' chanFile]);
end

% Reformat data according to channel info.
numChan     = size(channels.name,1);

%% Get header information from json files
% Check for MEG Info
try
    meg             = ft_read_json(megFile);
catch
    error(['ERROR. Attempted to read: ' megFile]);
end

try
    coordsystem     = ft_read_json(coordFile);
catch
    ft_warning('Not loaded any coordinate system');
    coordsystem = [];
end

% Position File check
try % to load a channels file
    posOri          = ft_read_pos_tsv(posFile);
    positions       = 1;
catch
    try % to load a BIDS channel file
        posOri       = ft_read_pos_tsv(posFile);
        positions    = 1;
    catch
        positions    = 0;
    end
end

%%
time_for_loop       = 0;
sampleinfo_for_loop = 0;
count               = 1;

%% Chunk it up
% compute the number of samples and chunks

% Read file in chunks to save memory
Chunk = (numChan*meg.SamplingFrequency*cfg.chunksize)+numChan;
Count = Chunk;    % Dummy value to satisfy WHILE condition

while Count == Chunk
    [data_raw, Count] = fread(dat,Chunk,cfg.precision,0,'b');
    
    
    %% Load in channels and reorganise the raw data
    
    if binData
        data    = reshape(data_raw,numChan,numel(data_raw)/numChan);
    elseif numChan ~= size(data_raw,1)
        error('number of channels in cfg.data different to cfg.channels')
    else
        data    = data_raw;
    end
    
    clear data_raw
    
    %% Make the main fieldtrip structure
    % Main structure
    rawData                     = [];
    rawData.label               = channels.name;
    rawData.fsample             = meg.SamplingFrequency;
    rawData.time{1}             = ((0:1:(size(data,2)-1))./rawData.fsample)+time_for_loop;
    rawData.trial{1}            = data;
    rawData.sampleinfo          = [sampleinfo_for_loop sampleinfo_for_loop+length(data(1,:))];
    
    % Save for next loop
%    time_for_loop       = rawData.time{1}(end)+(1/meg.SamplingFrequency);
    time_for_loop       = rawData.time{1}(end);
    sampleinfo_for_loop = rawData.sampleinfo(2)+1;
    
    % Construct a cfg file
    %%%%%%%%%
    % TO DO
    %%%%%%%%%%
    %rawData.cfg                 = cfg;
    
    % Now do the same for the header information
    rawData.hdr                 = [];
    rawData.hdr.label           = channels.name;
    rawData.hdr.chantype        = channels.type;
    rawData.hdr.chanunit        = channels.units;
    rawData.hdr.dimord          = 'chan_time';
    rawData.hdr.Fs              = meg.SamplingFrequency;
    rawData.hdr.nSamples        = length(data(1,:));
    rawData.hdr.nTrials         = 1;
    rawData.hdr.nSamplesPre     = 0;
    rawData.hdr.fieldori        = channels.fieldori;
    
    if positions
        rawData.grad.chanpos    = posOri.chanpos;
        rawData.grad.chanori    = posOri.chanori;
        rawData.grad.label      = posOri.label;
        rawData.grad.coilpos    = posOri.chanpos;
        rawData.grad.coilori    = posOri.chanori;
        
        % Find the position of the channels in overall channel list
        find_chan = contains(channels.name,posOri.label);
        rawData.grad.chanunit   = channels.units(find_chan);
        rawData.grad.tra = diag(ones(1,length(rawData.grad.label)));
        
    end
    
%     % Resample
%     sampleinfo = rawData.sampleinfo;
    
    if strcmp(cfg.downsample,'yes')
        %% Resample to 1000Hz
        disp(['Resampling Chunk... ' num2str(count)]);
        cfg2                 = [];
        cfg2.resamplefs      = 250;
        [rawData]            = ft_resampledata(cfg2, rawData);
        %rawData.sampleinfo   = sampleinfo./10;
    end
    
    %% Save
    disp(['Saving Chunk... ' num2str(count)]);
    save(['rawData_' num2str(count)],'rawData');
    count = count+1;
    
end

end







