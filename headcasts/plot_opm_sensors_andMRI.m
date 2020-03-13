%% Load in Gareth's MRI
cd('D:\data\Gareth_MRI');

mri = ft_read_mri('mmsMQ0484_orig.img');

cfg = [];

ft_sourceplot(cfg,mri);
mri.coordsys = 'neuromag';

%% Segment the brain
cfg           = [];
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri);

%% Create singleshell headmodel
cfg = [];
cfg.tissue = 'brain';
cfg.method='singleshell';
headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented);

% Create Figure
figure;ft_plot_vol(headmodel_singleshell);
alpha 0.3; view([0,0]);

%% Load the sensor locations
cd('D:\Github\analyse_OPMEG\headcasts');
sens_loc = readtable('Headcast_allAxes.txt');
pos = horzcat(sens_loc.Position_x,sens_loc.Position_y,sens_loc.Position_z);

%% Rotate sensors

ttt = cos(90*(pi/180));
rrr = sin(90*(pi/180));

trans = [ttt 0 rrr 0;
    0 1 0 0;
    -rrr 0 ttt 0;
    0 0 0 1];

ttt = cos(90*(pi/180));
rrr = sin(90*(pi/180));

trans2 = [1 0 0 0;
    0 ttt -rrr 0;
    0 rrr ttt 0;
    0 0 0 1];
    
allsens = ft_warp_apply(trans,pos);
allsens = ft_warp_apply(trans2,allsens);

ft_determine_coordsys(mri, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_mesh(allsens,'vertexcolor','r','vertexsize',40);

%% Extract Scalp Surface
cfg = [];
cfg.output    = 'scalp';
cfg.scalpsmooth = 5;
cfg.scalpthreshold = 0.08;
scalp  = ft_volumesegment(cfg, mri);

%% Create mesh out of scalp surface
cfg = [];
cfg.method = 'isosurface';
cfg.numvertices = 10000;
mesh = ft_prepare_mesh(cfg,scalp);
mesh = ft_convert_units(mesh,'mm');

figure;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.8); drawnow;

%% Now try ICP
[R, t, err] = icp(mesh.pos', allsens', 50, ...
    'Minimize', 'plane', 'Extrapolation', true,...
    'WorstRejection', 0.05);

%% Create figure to display how the ICP algorithm reduces error
clear plot;
figure; plot([1:1:51]',err,'LineWidth',8);
ylabel('Error'); xlabel('Iteration');
title('Error*Iteration');
set(gca,'FontSize',25);

%% Create transformation matrix
trans_matrix = inv([real(R) real(t);0 0 0 1]);

%% Create figure to assess accuracy of coregistration
pos_spare = allsens;
pos_spare = ft_warp_apply(inv(trans_matrix), pos_spare);

ft_determine_coordsys(mri, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_mesh(pos_spare,'vertexcolor','r','vertexsize',40); hold on;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.5); drawnow;

%%
% This where we define our MNI coordinates of interest
% visual cortex
% pos = [10 -90 0];
% FFA
pos = [44 -76 -14];

% Now we compute nonlinear warping between MNI and
cfg             = [];
cfg.template    = mri;
cfg.nonlinear   = 'yes';
norm            = ft_volumenormalise([],mri);

% Now we warp the MNI coordinates using the nonlinear warping parameters
posback         = ft_warp_apply(norm.params,pos,'sn2individual');
% xyz positions in individual coordinates
pos_grid        = ft_warp_apply(pinv(norm.initial),posback);

%% Figure with the works!
ft_determine_coordsys(mri, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_mesh(pos_spare,'vertexcolor','r','vertexsize',40); hold on;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.5); drawnow;
ft_plot_mesh(pos_grid,'vertexcolor','b','vertexsize',40); hold on;
for i = 1:length(pos_spare)
    text(pos_spare(i,1),pos_spare(i,2),pos_spare(i,3),num2str(i),...
        'FontSize',14)
end
ft_plot_vol(headmodel_singleshell);
alpha 0.3; view([0,0]);


%% Segment the brain
cfg                 = [];
cfg.output          = 'brain';
mri_segmented       = ft_volumesegment(cfg, mri);

%% Create singleshell headmodel
cfg                     = [];
cfg.tissue              = 'brain';
cfg.method              = 'singleshell';
headmodel               = ft_prepare_headmodel(cfg, mri_segmented);

%%
ft_determine_coordsys(mri, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_mesh(pos_spare,'vertexcolor','r','vertexsize',40); hold on;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none',...
    'facealpha',0.5); drawnow;
ft_plot_mesh(pos_grid,'vertexcolor','b','vertexsize',40); hold on;
for i = 1:length(pos_spare)
    text(pos_spare(i,1),pos_spare(i,2),pos_spare(i,3),num2str(i),...
        'FontSize',14)
end
ft_plot_vol(headmodel_singleshell);
alpha 0.3; view([0,0]);




