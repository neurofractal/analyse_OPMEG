function [rawData] = ft_opm_create(cfg)
% Function to read optically-pumped magnetencephalography (OPMEG) data
% acquired from the UCL Wellcome Centre for Neuroimaging.
%
% FORMAT data = ft_opm_create(cfg)
%   cfg               - input structure
% Optional fields of cfg:
% File information
%   cfg.folder          = 'string', BIDS folder.
% Provide either the filename or BIDS format: task, sub, ses and run.
%   cfg.filename        = 'MEG.bin'
%   cfg.filename.task   = 'AEF'
%   cfg.filename.sub    = '001'
%   cfg.filename.ses    = '003'
%   cfg.filename.run    = '002'
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_create (Tim Tierney)

% Authors:  Nicholas Alexander  (n.alexander@ucl.ac.uk)
%           Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Set default values
if ~isfield(cfg, 'voltype')
    cfg.voltype = 'Single Shell';
end
if ~isfield(cfg, 'meshres')
    cfg.meshres = 1;
end
if ~isfield(cfg, 'scalp')
    cfg.scalp = [];
end
if ~isfield(cfg, 'cortex')
    cfg.cortex = [];
end
if ~isfield(cfg, 'iskull')
    cfg.iskull = [];
end
if ~isfield(cfg, 'oskull')
    cfg.oskull = [];
end
if ~isfield(cfg, 'lead')
    cfg.lead = 0;
end
if ~isfield(cfg, 'fs')
    cfg.fs   = 1000;
end
if ~isfield(cfg, 'nSamples')
    cfg.nSamples   = 1000;
end
if ~isfield(cfg, 'nDens')
    cfg.nDens   = 40;
end
if ~isfield(cfg, 'space')
    cfg.space  = 30;
end
if ~isfield(cfg, 'offset')
    cfg.offset  = 6.5;
end
if ~isfield(cfg, 'data')
    cfg.data = zeros(1,cfg.nSamples);
end
if ~isfield(cfg, 'wholehead')
    cfg.wholehead = 1;
end
if ~isfield(cfg, 'fname')
    cfg.fname = 'sim_opm';
end
if ~isfield(cfg, 'precision')
    cfg.precision = 'single';
end

%% Read Binary File
try % to read data 
    [direc, dataFile] = fileparts(cfg.data);
    dat         = fopen(cfg.data);
    data_raw    = fread(dat,Inf,cfg.precision,0,'b');
    fclose(dat);
    binData     = 1;
catch % if not readable check if it is numeric 
    if ~isa(cfg.data,'numeric') % if not numeric throw error
        error('A valid dataest or file was not supplied')
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

%% Load in channels and reorganise the raw data
% Load a channels file.
channels    = ft_read_chan_tsv(chanFile);

% Reformat data according to channel info.
numChan     = size(channels.name,1);

if binData
    data    = reshape(data_raw,numChan,numel(data_raw)/numChan);
elseif numChan ~= size(data_raw,1)
    error('number of channels in cfg.data different to cfg.channels')
else
    data    = data_raw;
end

clear data_raw

%% Get header information from json files
% Check for MEG Info
meg             = ft_read_json(megFile);
coordsystem     = ft_read_json(coordFile);

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

%% Make the main fieldtrip structure
% Main structure
rawData                     = [];
rawData.label               = channels.name;
rawData.fsample             = meg.SamplingFrequency;
rawData.time{1}             = (0:1:(size(data,2)-1))./rawData.fsample;
rawData.trial{1}            = data;
rawData.sampleinfo          = [1 length(data(1,:))];

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

if positions
    rawData.grad.chanpos    = posOri.chanpos;
    rawData.grad.chanori    = posOri.chanori;
    indx = find(contains(channels.type,'megmag'));
    rawData.grad.chanunit   = channels.units((indx));
    rawData.grad.label      = channels.name;
end

% % Forward model Check
% subjectSource  = (positions|isfield(cfg,'space')) & isfield(cfg,'sMRI');
% subjectSensor = ~subjectSource;
% 
% if subjectSource
%     forward     = 1;
%     template    = 0;
% elseif subjectSensor
%     forward 	= 0;
%     template    = 0;
% else
%     forward     = 1;
%     template    = 1;
end







