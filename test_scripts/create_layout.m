%% Create layouts
fieldtripDir    = '/Users/rseymoue/Documents/scripts/fieldtrip-20191213';
script_dir      = '/Users/rseymoue/Documents/GitHub/analyse_OPMEG';
data_dir        = '/Volumes/Robert T5/OPM_data/benchmarking_26_02/';
save_dir        = '/Users/rseymoue/Documents/GitHub/opm_benchmarking_feb_2020/';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));

% cd to save dir
cd(save_dir)

%% (2) Start preprocessing.
% Read in the raw data using BIDS
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'faces';
cfg.bids.sub    = '002';
cfg.bids.ses    = '001';
cfg.bids.run    = '001';
rawData         = ft_opm_create(cfg);

%% Load in the MRI
mri = ft_read_mri(fullfile(data_dir,'mmsMQ0484_orig.img'));
mri.coordsys = 'neuromag';

%% Extract Scalp Surface from the MRI and create a mesh
cfg                     = [];
cfg.output              = 'scalp';
cfg.scalpsmooth         = 5;
cfg.scalpthreshold      = 0.08;
scalp                   = ft_volumesegment(cfg, mri);

% Create mesh
cfg                     = [];
cfg.method              = 'isosurface';
cfg.numvertices         = 10000;
mesh                    = ft_prepare_mesh(cfg,scalp);
mesh                    = ft_convert_units(mesh,'mm');

% Plot Mesh
figure;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.8); camlight; drawnow;
ft_plot_mesh(gareth_head_trans); camlight;

%% Script to create layout
which_ori = 'BOTH';
paradigm  = 'faces_M170';

% TANS
if strcmp(which_ori,'TAN') || strcmp(which_ori,'BOTH')
    cfg             = [];
    cfg.output      = ['lay_tan_' paradigm '.mat'];
    cfg.grad        = rawData.grad;
    cfg.headshape   = mesh;
    cfg.projection  = 'orthographic';
    cfg.channel     = vertcat(ft_channelselection_opm('TAN',rawData),...
        '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
    % Change this depending on where the sensors are!
    cfg.viewpoint   = 'posterior';
    lay_tan         = ft_prepare_layout(cfg);
    
    figure; ft_plot_layout(lay_tan); title('TAN');
end

% Now RADS
if strcmp(which_ori,'RAD') || strcmp(which_ori,'BOTH')
    cfg             = [];
    cfg.output      = ['lay_rad_' paradigm '.mat'];
    cfg.grad        = rawData.grad;
    cfg.headshape   = mesh;
    cfg.projection  = 'orthographic';
    cfg.channel     = vertcat(ft_channelselection_opm('RAD',rawData),...
        '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
    % Change this depending on where the sensors are!
    cfg.viewpoint   = 'posterior';
    lay_rad         = ft_prepare_layout(cfg);
    
    figure; ft_plot_layout(lay_rad); title('RAD');
end







