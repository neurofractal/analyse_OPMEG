function sourceI = interpolate_parcels(template_grid,atlas,VE,map)
%__________________________________________________________________________
% Interpolate a map (e.g. ROI x power) onto SPM canonical brain
% 
% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

% Load 5mm sourcemodel
[t, r] = ft_version;

template_grid   = ft_convert_units(template_grid,'mm');
atlas    = ft_convert_units(atlas,'mm');

if ~isfield(atlas,'tissue')
    atlas.tissue = atlas.parcellation;
    atlas.tissuelabel = atlas.parcellationlabel;
end


% and call ft_sourceinterpolate:
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'tissue';
sourcemodel2                = ft_sourceinterpolate(cfg, atlas,template_grid);

%% Plotting stuff
% Load the SPM Brain
spm_brain = ft_read_mri(fullfile(r,'template/anatomy/single_subj_T1.nii'));

% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

% Create empty power array
pow  = nan(size(template_grid.inside,1),1);
pow(template_grid.inside) = 0;

% Find the atlas labels used in our (cortical) parcellation
[a b] = ismember(VE.label,sourcemodel2.tissuelabel);
atlas_labels_to_use = b(a);

% For each region make power = map
for i=1:length(atlas_labels_to_use)
    disp([num2str(atlas_labels_to_use(i))...
        '. ' sourcemodel2.tissuelabel{atlas_labels_to_use(i)}]);

    d = find(sourcemodel2.tissue==atlas_labels_to_use(i));
    pow(d) = map(i);
end

%% Interpolate
template_grid.pow = pow;

% and call ft_sourceinterpolate:
cfg                 = [];
cfg.interpmethod    = 'nearest';
cfg.parameter       = 'pow';
sourceI             = ft_sourceinterpolate(cfg, template_grid, spm_brain);
sourceI.coordsys    = 'mni';
sourceI = rmfield(sourceI,"cfg");

