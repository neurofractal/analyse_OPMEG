function [rawData] = ft_opm_create(cfg)
% Read or simulate magnetometer data and optionally set up forward model
% FORMAT D = spm_opm_create(S)
%   cfg               - input structure
% Optional fields of cfg:
% SENSOR LEVEL INFO
%   cfg.data          - filepath to .bin file        - Default: Simulates data 
%   cfg.channels      - channels.tsv file            - Default: REQUIRED  
%   cfg.fs            - Sampling frequency (Hz)      - Default: REQUIRED if cfg.meg is empty
%   cfg.meg           - meg.json file                - Default: REQUIRED if cfg.fs is empty
%   cfg.precision     - 'single' or 'double'         - Default: 'single'
% SIMULATION
%   cfg.wholehead     - whole head coverage flag     - Deafult: 0
%   cfg.space         - space between sensors(mm)    - Default: 25
%   cfg.offset        - scalp to sensor distance(mm) - Default: 6.5
%   cfg.nSamples      - number of samples            - Default: 1000
%   cfg.Dens          - number of density checks     - Default: 40

% SOURCE LEVEL INFO
%   cfg.coordsystem   - coordsystem.json file        - Default: 
%   cfg.positions     - positions.tsv file           - Default:
%   cfg.sMRI          - Filepath to  MRI file        - Default: uses template
%   cfg.cortex        - Custom cortical mesh         - Default: Use inverse normalised cortical mesh
%   cfg.scalp         - Custom scalp mesh            - Default: Use inverse normalised scalp mesh
%   cfg.oskull        - Custom outer skull mesh      - Default: Use inverse normalised outer skull mesh
%   cfg.iskull        - Custom inner skull mesh      - Default: Use inverse normalised inner skull mesh
%   cfg.voltype       - Volume conducter Model type  - Default: 'Single Shell'
%   cfg.meshres       - mesh resolution(1,2,3)       - Default: 1
%   cfg.lead          - flag to compute lead field   - Default: 0
% Output:
%  D           - MEEG object (also written to disk)
%  L           - Lead field (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% Adapted from spm_opm_create (Tim Tierney) by Nicholas Alexander

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
    cfg.data    = fread(dat,Inf,cfg.precision,0,'b');
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

%% identify potential BIDS Files
base        = strsplit(dataFile,'meg');
chanFile    = fullfile(direc,[base{1},'channels.tsv']);
megFile     = fullfile(direc,[base{1},'meg.json']);
posFile     = fullfile(direc,[base{1},'positions.tsv']);
coordFile   = fullfile(direc,[base{1},'coordsystem.json']);

% Load a channels file.
channels    = ft_read_chan_tsv(chanFile);

% Reformat data according to channel info.
numChan     = size(channels.label,1);

if binData
    data    = reshape(cfg.data,numChan,numel(cfg.data)/numChan);
elseif numChan ~= size(cfg.data,1)
    error('number of channels in cfg.data different to cfg.channels')
else
    data    = cfg.data;
end

% Check for MEG Info
meg             = ft_read_json(cfg.meg);
coordsystem     = ft_read_json(coordFile);

% Position File check 
try % to load a channels file
    posOri          = ft_read_pos_tsv(cfg.positions);
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
rawData         = [];
rawData.label       = cellstr(channels.label);
rawData.chantype    = cellstr(channels.chantype);
rawData.chanunit    = cellstr(channels.unit);
rawData.trial{1}    = data;
rawData.time{1}     = 0:(1/meg.SamplingFrequency):(length(data(1,:))/meg.SamplingFrequency);
rawData.dimord      = 'chan_time';
rawData.cfg         = cfg;
rawData.sampleinfo  = [1 length(data(1,:))];
rawData.hdr.Fs          = meg.SamplingFrequency;
rawData.hdr.nChans      = length(channels.label);
rawData.hdr.nSamples    = length(data(1,:));
rawData.hdr.nTrials     = 1;
rawData.hdr.label       = channels.label;
rawData.hdr.chantype    = channels.chantype;
rawData.hdr.chanunit    = channels.unit;
rawData.grad.chanpos    = posOri.chanpos;
rawData.grad.chanori    = posOri.chanori;
rawData.grad.chanunit   = channels.unit;
rawData.grad.label      = channels.label;

% Forward model Check
subjectSource  = (positions|isfield(cfg,'space')) & isfield(cfg,'sMRI');
subjectSensor = ~subjectSource;

if subjectSource
    forward     = 1;
    template    = 0;
elseif subjectSensor
    forward 	= 0;
    template    = 0;
else
    forward     = 1;
    template    = 1;
end







