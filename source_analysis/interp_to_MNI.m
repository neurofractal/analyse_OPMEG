function [sourceI] = interp_to_MNI(source,normalise,spm_brain)

% First load the 5mm template_grid (currented hard-coded.. could change)
[t, r] = ft_version;

load(fullfile(r,'template/sourcemodel/standard_sourcemodel3d5mm.mat'));
template_grid = sourcemodel;
clear sourcemodel
template_grid = ft_convert_units(template_grid,'mm');

% Do the transform(s)
source.pos = ft_warp_apply(normalise.initial,source.pos);
source.pos = ft_warp_apply(normalise.params,source.pos, 'individual2sn');

% Plot for anity check
ft_determine_coordsys(spm_brain, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
hold on;
ft_plot_mesh(source.pos,'vertexcolor','r','vertexalpha',0.2); hold on;
ft_plot_mesh(template_grid.pos(template_grid.inside,:),'vertexsize',1);
drawnow;

disp('Interpolating...');
% Interpolate onto SPM brain
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
% Interpolate from individual to template grid
sourceI  = ft_sourceinterpolate(cfg, source, template_grid);
% Interpolate onto spm brain
sourceI  = ft_sourceinterpolate(cfg, source, spm_brain);

