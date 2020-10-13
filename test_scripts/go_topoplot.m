function h = go_topoplot(cfg,dat)

if ~isstruct(cfg.layout)
    load(cfg.layout);
    cfg.layout = layout;
    clear layout
end

if ft_platform_supports('griddata-v4')
  default_interpmethod = 'v4';
else
  % Octave does not support 'v4', and 'cubic' not yet implemented
  default_interpmethod = 'linear';
end

cfg.contournum        = ft_getopt(cfg, 'contournum',        6);
cfg.gridscale         = ft_getopt(cfg, 'gridscale',         67);
cfg.interplimits      = ft_getopt(cfg, 'interplimits',     'head');
cfg.interpolation     = ft_getopt(cfg, 'interpolation',     default_interpmethod);
cfg.shading           = ft_getopt(cfg, 'shading',          'flat');
cfg.style             = ft_getopt(cfg, 'style',            'both');

% Set ft_plot_topo specific options
if strcmp(cfg.interplimits, 'head')
  interplimits = 'mask';
else
  interplimits = cfg.interplimits;
end

if strcmp(cfg.style, 'both');            style = 'surfiso';     end
if strcmp(cfg.style, 'straight');        style = 'surf';        end
if strcmp(cfg.style, 'contour');         style = 'iso';         end
if strcmp(cfg.style, 'fill');            style = 'isofill';     end
if strcmp(cfg.style, 'straight_imsat');  style = 'imsat';       end
if strcmp(cfg.style, 'both_imsat');      style = 'imsatiso';    end

msk = [];

% Select the channels in the data that match with the layout:
seldat = find(contains(cfg.layout.label,cfg.channels));
if isempty(seldat)
  ft_error('labels in data and labels in layout do not match');
end

dat = dat(seldat);
if ~isempty(msk)
  msk = msk(seldat);
end

% Select x and y coordinates and labels of the channels in the data
chanX = cfg.layout.pos(seldat,1);
chanY = cfg.layout.pos(seldat,2);

opt = {'interpmethod', cfg.interpolation, ...
    'interplim',    interplimits, ...
    'gridscale',    cfg.gridscale, ...
    'outline',      cfg.layout.outline, ...
    'shading',      cfg.shading, ...
    'isolines',     cfg.contournum, ...
    'mask',         cfg.layout.mask, ...
    'style',        style, ...
    'datmask',      msk};

% figure
ft_plot_topo(chanX, chanY, dat, opt{:});
hold on
scatter(chanX,chanY,'k.')
hold off
axis equal
axis off

