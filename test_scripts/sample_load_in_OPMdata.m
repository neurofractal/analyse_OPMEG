%% Paths (RS)
fieldtripDir    = '/Users/rseymoue/Documents/scripts/fieldtrip-20191213';
script_dir      = '/Users/rseymoue/Documents/GitHub/analyse_OPMEG';
data_dir        = '/Volumes/Robert T5/20201209';
save_dir        = '/Users/rseymoue/Dropbox/Research/Projects/2020/opm_benchmarking/optitrack';
denoise_dir     = '/Users/rseymoue/Documents/scripts/NoiseTools';
NR4M_dir        = '/Users/rseymoue/Documents/GitHub/NR4M';
scannercast_dir = '/Users/rseymoue/Documents/GitHub/scannercast/examples/NA';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));

% Add NoiseTools to path
addpath(denoise_dir)
addpath(genpath(NR4M_dir));
% cd to save dir
cd(save_dir)

%% (2) Start preprocessing.
% Read in the raw data using BIDS
disp('Loading Data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'noise';
cfg.bids.sub    = 'noise';
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

%% Work in progress code for mean field correction
cfg                     = [];
cfg.path_to_SPM         = '/Users/rseymoue/Documents/scripts/spm12';
cfg.path_to_OPM_repo    = '/Users/rseymoue/Documents/GitHub/OPM';
cfg.correctgrad         = 'yes';
[rawData_MEG_mfc]       = ft_wrapper_spm_opm_mfc(cfg,rawData_MEG);

%% PSD
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [1 150];
cfg.plot            = 'yes';
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





