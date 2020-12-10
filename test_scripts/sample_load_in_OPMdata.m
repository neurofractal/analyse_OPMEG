%% Paths (RS)
fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
data_dir        = 'D:\data\20201112_motor';
save_dir        = 'C:\Users\rseymour\Dropbox\Research\Projects\2020\opm_benchmarking\motor';
denoise_dir     = 'D:\scripts\NoiseTools';
NR4M_dir        = 'D:\Github\NR4M';
scannercast_dir = 'D:\Github\scannercast\examples\NA\';

% Add Fieldtrip to path
disp('Adding Fieldtrip, NoiseTools, NR4M and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));

% Add NoiseTools to path
addpath(denoise_dir)
% Add NR4M to path
addpath(genpath(NR4M_dir));
% cd to save dir
cd(save_dir)

%% (2) Start preprocessing.
% Read in the raw data using BIDS
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'motor_right';
cfg.bids.sub    = 'NA';
cfg.bids.ses    = '001';
cfg.bids.run    = '001';
rawData         = ft_opm_create(cfg);

%% Resample to 600Hz
cfg                 = [];
cfg.resamplefs      = 600;
[rawData]           = ft_resampledata(cfg, rawData);

%% Plot using ft_databrowser
cfg             = [];
cfg.viewmode    = 'butterfly';
cfg.blocksize   = 20;
ft_databrowser([],rawData);

%% Select only OPM channels
cfg             = [];
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
rawData_MEG     = ft_selectdata(cfg,rawData);

%% PSD
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [1 150];
cfg.plot            = 'yes';
[pow freq]          = ft_opm_psd(cfg,rawData_MEG);
ylim([1 1e4])

%% Mean Field Correction
% Please ask George O'N for ft_denoise_mfc.m
[rawData_MEG_mfc, M, chan_inds] = ft_denoise_mfc(rawData_MEG);

% Plot PSD
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [1 150];
cfg.plot            = 'yes';
cfg.plot_legend      = 'no';
[pow freq]          = ft_opm_psd(cfg,rawData_MEG_mfc);
ylim([1 1e4])

%% DSSP
% Load the sourcemodel (3D mesh) and headmodel
load(fullfile(scannercast_dir,...
    'sourcemodel_10mm.mat'));
load(fullfile(scannercast_dir,'headmodel.mat'));

disp('Preparing Leadfield');
cfg             = [];
cfg.method      = 'lcmv';
cfg.channel     = rawData.label;
cfg.grid        = sourcemodel;
cfg.headmodel   = headmodel;
cfg.grad        = rawData.grad;
cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
cfg.normalize       = 'no' ; %Normalise Leadfield: 'yes' for beamformer
%cfg.normalizeparam  = 1;
lf_for_DSSP = ft_prepare_leadfield(cfg);

% make a figure of the single subject{i} headmodel, and grid positions
figure; hold on;
ft_plot_headmodel(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; camlight;
ft_plot_mesh(lf_for_DSSP.pos(lf_for_DSSP.inside,:));
ft_plot_sens(rawData.grad, 'style', 'r*','orientation','true'); view([0,0]);

% DSSP
cfg                     = [];
cfg.sourcemodel         = lf_for_DSSP;
cfg.dssp.n_space        = length(rawData_MEG.label)-2;
cfg.dssp.n_in           = length(rawData_MEG.label)-2;
cfg.dssp.n_out          = length(rawData_MEG.label)-2;
cfg.dssp.n_intersect    = 5;
cfg.winsize             = 10;
rawData_DSSP            = ft_dssp_window(cfg, rawData_MEG);

% Plot the PSD
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [1 150];
cfg.plot            = 'yes';
cfg.plot_legend      = 'no';
pow                 = ft_opm_psd(cfg,rawData_DSSP);
ylim([1 1e4])





