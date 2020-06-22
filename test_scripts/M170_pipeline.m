%% Paths (RS)
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

% % Plot using ft_databrowser
% ft_databrowser([],rawData);

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData),...
     '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [1 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,rawData);
ylim([1 1e4])

%% Resample to 600Hz to speed things up 
cfg                 = [];
cfg.resamplefs      = 600;
[rawData]           = ft_resampledata(cfg, rawData);

%% Filter
cfg = [];
cfg.hpfilter        = 'yes';
cfg.hpfreq      	= 2;
rawData             = ft_preprocessing(cfg,rawData);

cfg = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 100;
rawData         = ft_preprocessing(cfg,rawData);

cfg = [];
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [115 125];
rawData         = ft_preprocessing(cfg,rawData);

cfg = [];
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [97 103];
rawData         = ft_preprocessing(cfg,rawData);

%% Perform synthetic gradiometry
cfg = [];
cfg.channel = vertcat(ft_channelselection_opm('MEG',rawData),...
     '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN','-N0-RAD','-N4-RAD',...
     '-N3-RAD','-MV-RAD');
cfg.refchannel = 'MEGREF';
cfg.filter_ref = [2 30; 40 60; 60 80];
cfg.derivative = 'yes';
cfg.return_all = 'no';
cfg.winsize    = 50;
[rawData_meg_window] = ft_opm_synth_gradiometer_window(cfg,rawData);

% Plot the PSD
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [0.2 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,rawData_meg_window);

% Plot the Shielding Factor
cfg                 = [];
cfg.channel         = rawData_meg_window.label;
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [0.2 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd_compare(cfg,rawData,rawData_meg_window);

% Low-pass filter at 30
cfg             = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 30;
rawData_meg     = ft_preprocessing(cfg,rawData_meg_window);

% % Plot using ft_databrowser
ft_databrowser([],rawData_meg_window);

%% Trial Definition
% Using 'OPM_TrialFun_RS'
cd([data_dir 'sub-002/ses-001/meg']);

cfg = [];
cfg.dataset                 = 'sub-002_ses-001_task-faces_run-001_meg.bin';
cfg.trialdef.trigchan       = 'FluxZ-B';
cfg.trialdef.downsample     = 600;
cfg.correct_time            = 0.1;
cfg.trialdef.prestim        = 0.5        % pre-stimulus interval
cfg.trialdef.poststim       = 0.7;        % post-stimulus interval
cfg.trialfun                = 'OPM_TrialFun_RS';
banana                      = ft_definetrial(cfg);

% Redefines the filtered data
cfg = [];
cfg.detrend = 'yes';
data = ft_redefinetrial(banana,rawData_meg_window);

% Plot the data
cfg             = [];
cfg.viewmode    = 'vertical';
cfg.colormode   = 'allblack';
ft_databrowser(cfg,data_detrend);

%% Load the csv file to know which stimuli was on screen
csvfile = readtable([data_dir 'M170_RS_run_1_20200226_1141.csv']);
%csvfile = readtable([data_dir 'M170_RS_run_2_20200226_1158.csv']);

% I made a mistake.. so you need to extract the identity of the stimuli
% from the csv file
stimfile = csvfile.Stim_Path;
ddd = [];

for s = 1:length(stimfile)
    xxx = split(stimfile{s},'\');
    if strcmp(xxx{2}(1),'f')
        ddd{s,1} = 'Famous';
    elseif strcmp(xxx{2}(1),'s')
        ddd{s,1} = 'Scrambled';
    elseif strcmp(xxx{2}(1),'u')
        ddd{s,1} = 'Unfamiliar';
    end
end

data.condition = ddd;



%% Select the data for each condition
cfg                 = [];
cfg.trials          = contains(data.condition,'Famous');
famous_faces        = ft_selectdata(cfg,data);
cfg.trials          = contains(data.condition,'Unfamiliar');
unfamiliar_faces    = ft_selectdata(cfg,data);
cfg.trials          = contains(data.condition,'Scrambled');
scrambled_faces     = ft_selectdata(cfg,data);

% Append the face data together
faces_all = ft_appenddata([], famous_faces, unfamiliar_faces);

% Perform timelockanalysis
cfg             = [];
avg_all         = ft_timelockanalysis([],data);
avg_famous      = ft_timelockanalysis([],famous_faces);
avg_unfamiliar  = ft_timelockanalysis([],unfamiliar_faces);
avg_scrambled   = ft_timelockanalysis([],scrambled_faces);
avg_faces_all   = ft_timelockanalysis([],faces_all);

% %% Plot each sensor in turn
% for i = 1:length(avg_famous.label)
%     cfg = [];
%     cfg.channel = avg_famous.label{i}
%     cfg.parameter = 'avg';
%     cfg.title = avg_famous.label{i};
%     cfg.baseline = [-0.5 0];
%     %cfg.baselinetype = 'relative';
%     %cfg.showlegend    = 'yes';
%     %cfg.ylim = [-250 250];
%     figure; ft_singleplotER(cfg,avg_famous, avg_unfamiliar,...
%         avg_scrambled);
% end
% 
% %% Now plot a nice sample channel
% color_scheme = [0.0157         0    1.0000;
%      0.1529    0.6000    0.0039;
%      0.8000    0.0667         0];
%     
% cfg = [];
% cfg.channel = 'DG-TAN';
% cfg.parameter = 'avg';
% cfg.title = 'DG-TAN';
% cfg.baseline = [-0.5 0];
% cfg.linecolor = color_scheme;
% cfg.linewidth = 2;
% %cfg.showlegend    = 'yes';
% cfg.ylim = [-300 300];
% figure; ft_singleplotER(cfg,avg_famous, avg_unfamiliar,...
%     avg_scrambled);legend({'Famous Face'; 'Unfamiliar Face'; ...
%     'Scrambled Face'},'Location','southwest');
% xlabel('Time (s)');
% ylabel('fT');
% set(gca,'FontSize',14);
% print('MT_faces_M170_TAN','-dpng','-r300');

%% Plot all the channels avg in one plot

cfg = [];
cfg.baseline = [-0.7 0];
[avg_all] = ft_timelockbaseline(cfg, avg_all)

cfg = [];
%cfg.ylim = [-455 455];
cfg.parameter = 'avg';
ft_databrowser(cfg,avg_all);

% cfg = [];
% cfg.channel = {'all','-MX-TAN'};
% data = ft_selectdata(cfg,data);

%% Plot sensor-level topoplot
cfg = [];
cfg.xlim = [0.1 0.2];
cfg.colorbar = 'yes';
cfg.baseline = [-0.5 0];
%cfg.layout = lay;
figure;ft_topoplotER(cfg,avg_faces_all);

% %%
% cfg = [];
% cfg.channel = {'all','-MX-TAN','-MX-RAD'};
% data = ft_selectdata(cfg,data);

%% Source Analysis
load([data_dir 'headmodel.mat']);
load([data_dir 'sourcemodel3d.mat']);

%% Prepare Leadfield

disp('Preparing Leadfield');
cfg = [];
cfg.method='lcmv';
cfg.channel = data.label;
cfg.grid = sourcemodel3d;
cfg.headmodel = headmodel;
cfg.grad = rawData.grad;
%cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
%cfg.normalizeparam  = 1;
lf = ft_prepare_leadfield(cfg);

% make a figure of the single subject{i} headmodel, and grid positions
figure; hold on;
ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; camlight;
%ft_plot_mesh(lf.pos(lf.inside,:));
ft_plot_sens(rawData.grad, 'style', 'r*','orientation','true'); view([0,0]);
%print('lf_headmodel_sens','-dpng','-r100');

%% Compute covariance
cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.covariancewindow = [-0.3 0.3];
avg                  = ft_timelockanalysis(cfg,data);

% %Make a dummy variable with covariance matrices averaged
% avg_combined        = avg_deviant;
% avg_combined.cov    = (avg_deviant.cov+avg_predeviant.cov)./2;

%% Source Analysis Proper (using LCMV)
cfg                    = [];
cfg.channel            = data.label;
cfg.grad               = rawData.grad;
cfg.method             = 'lcmv';
cfg.grid               = lf;
cfg.headmodel          = headmodel;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.fixedori      = 'yes';
cfg.lcmv.projectnoise  = 'yes';
cfg.lcmv.weightnorm    = 'nai';
cfg.lcmv.lambda        = '5%';
sourceall              = ft_sourceanalysis(cfg, avg);

% Replace .pos field with template_grid.pos

[t, r] = ft_version;

load(fullfile(r,'template/sourcemodel/standard_sourcemodel3d5mm.mat'));
template_grid = sourcemodel;
clear sourcemodel

sourceall.pos = template_grid.pos;

% Remove cfg field to save memory
sourceall = rmfield(sourceall,'cfg');

addpath(genpath('/Users/rseymoue/Documents/Github/MQ_MEG_Scripts'));

source_pow_post = get_source_pow(data,sourceall,[0.15 0.2]);
source_pow_pre = get_source_pow(data,sourceall,[-0.2 -0.15]);

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'pow';
source_R = ft_math(cfg,source_pow_post,source_pow_pre);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.location = 'max';
ft_sourceplot(cfg,source_R);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap


%% Interpolate onto SPM brain
spm_brain = ft_read_mri(fullfile(r,'template/anatomy/single_subj_T1.nii'));       
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourceI  = ft_sourceinterpolate(cfg, source_R, spm_brain);

% Plot
cfg = [];
cfg.funparameter = 'pow';
cfg.location = 'max';
ft_sourceplot(cfg,sourceI);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%saveas(gcf,'sourceI_new3.png'); drawnow;

% 
% %%
% cfg = [];
% cfg.method          = 'surface';
% cfg.funparameter    = 'pow';
% cfg.colorbar        = 'yes';
% %cfg.funcolorlim     = [-5.1e-14 5.1e-14];
% %cfg.maskparameter = 'mask';
% %cfg.funcolorlim    = [0 20];
% %cfg.downsample     = 6;
% cfg.projmethod     = 'nearest';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% %cfg.surfdownsample = 10
% %cfg.projthresh     = 0.6;
% cfg.camlight       = 'no';
% ft_sourceplot(cfg, sourceI);
% try
%     ft_hastoolbox('brewermap', 1);
%     colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
% catch
%     warning('default colormap used instead of red/blue one');
% end
% light ('Position',[-180 0 0])
% light ('Position',[180 0 0])
% material dull;
% title('');
% view([180 0]);
% set(gca,'FontSize',14);
% drawnow;
% 
% %% Do virtual electrode analysis
% %%
% % Now we compute nonlinear warping between MNI and
% mri = ft_read_mri([data_dir 'mmsMQ0484_orig.img']);
% mri.coordsys = 'neuromag';
% 
% cfg             = [];
% cfg.template    = mri;
% cfg.nonlinear   = 'yes';
% norm            = ft_volumenormalise([],mri);
% 
% pos = [-28 -72 8];
% 
% % Now we warp the MNI coordinates using the nonlinear warping parameters
% posback         = ft_warp_apply(norm.params,pos,'sn2individual');
% % xyz positions in individual coordinates
% pos_grid        = ft_warp_apply(pinv(norm.initial),posback);
%     
% %% Prepare Leadfield
% cfg = [];
% cfg.method='lcmv';
% cfg.channel = data.label;
% cfg.grid.pos = pos_grid;
% cfg.grid.unit = 'mm';
% cfg.headmodel = headmodel;
% cfg.grad = rawData.grad;
% %cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
% cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
% %cfg.normalizeparam  = 1;
% lf_2 = ft_prepare_leadfield(cfg);
% 
% % make a figure of the single subject{i} headmodel, and grid positions
% figure; hold on;
% ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
% alpha 0.5; camlight;
% ft_plot_mesh(lf_2.pos(lf_2.inside,:),'vertexsize',20,'vertexcolor','r');
% ft_plot_sens(rawData.grad, 'style', 'r*'); view([0,0]);
% 
% 
% cfg                    = [];
% cfg.channel            = data.label;
% cfg.grad               = rawData.grad;
% cfg.method             = 'lcmv';
% cfg.grid               = lf_2;
% cfg.headmodel          = headmodel;
% cfg.lcmv.keepfilter    = 'yes';
% cfg.lcmv.fixedori      = 'yes';
% cfg.lcmv.projectnoise  = 'yes';
% %cfg.lcmv.weightnorm    = 'nai';
% cfg.lcmv.lambda        = '5%';
% sourceall              = ft_sourceanalysis(cfg, avg);
% 
% % Find filter from max point
% filter = sourceall.avg.filter{1,1};
% 
% VE = [];
% VE.label = {'max'};
% VE.trialinfo = data.trialinfo;
% VE.condition = data.condition;
% for subs=1:(length(data.trialinfo))
%     % note that this is the non-filtered "raw" data
%     VE.time{subs}       = data.time{subs};
%     VE.trial{subs}(1,:) = filter(1,:)*data.trial{subs}(:,:);
% end
% 
% % Select data based on conditions
% cfg                 = [];
% cfg.trials          = contains(VE.condition,'Famous');
% famous_faces        = ft_selectdata(cfg,VE);
% cfg.trials          = contains(VE.condition,'Unfamiliar');
% unfamiliar_faces    = ft_selectdata(cfg,VE);
% cfg.trials          = contains(VE.condition,'Scrambled');
% scrambled_faces     = ft_selectdata(cfg,VE);
% 
% 
% % Append the face data together
% faces_all = ft_appenddata([], famous_faces, unfamiliar_faces);
% 
% % Perform timelockanalysis
% cfg             = [];
% avg_famous      = ft_timelockanalysis([],famous_faces);
% avg_unfamiliar  = ft_timelockanalysis([],unfamiliar_faces);
% avg_scrambled   = ft_timelockanalysis([],scrambled_faces);
% avg_faces_all   = ft_timelockanalysis([],faces_all);
% 
% %% Plot each sensor in turn
% cfg = [];
% cfg.parameter = 'avg';
% cfg.baseline = [-0.5 0];
% cfg.showlegend    = 'yes';
% cfg.linewidth = 2;
% %cfg.ylim = [-250 250];
% figure; ft_singleplotER(cfg,avg_famous, avg_unfamiliar,...
%     avg_scrambled);
% legend({'Famous';'Unfamiliar';'Scrambled'});
% set(gca,'FontSize',20);
% xlabel('Time (s)');
% title('');
% 
% 
