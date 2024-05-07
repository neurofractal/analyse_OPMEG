function [VE] = atlas2VE(cfg,atlas, data_clean,sourceall, avg_cov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mq_atlas2VE: a function to create a virtual electrode using an atlas!
%
% EXAMPLE USEAGE:   [VE] = atlas2VE(cfg, atlas, data_clean, sourceall, avg_cov)
% ...where, cfg is a configuration structure.
%
% - atlas           = atlas you wish to use for VE computation. Should have
%                   been read into MATLAB using ft_read_atlas
% - data_clean      = your clean sensor-level data
% - source_all      = result of ft_sourceanalysis with the filters you wish
%                   to use for your VE computation
% - avg_cov         = covariance matrix used for source analysis
%                   (output of ft_timelockanalysis with .cov field)
% - cfg.channel     = name of the ROIs/chans in your atlas 
%                     (default = 'all')
% - cfg.template_grid   = template_grid used for source analysis (read into
%                   MATLAB from .../fieldtrip-XXXXX/template/sourcemodel/
%                   standard_sourcemodel3dXmm.mat
% - cfg.vis         = visualisation option ('no','yes','fancy')
%
%%%%%%%%%
% Output:
%%%%%%%%%
% - VE              = virtual electrode with label, time, trial and
%                   trialinfo fields
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)   
%__________________________________________________________________________

if ~isfield(cfg, 'channel')
    cfg.channel = 'all';
end

if ~isfield(cfg, 'vis')
    cfg.vis = 'fancy';
end

% If 'all' specified use all the tissuelabels
if strcmp(cfg.channel,'all')
    cfg.channel = atlas.tissuelabel;
end

% If atlas has a parcellation field but no tissue field, fix this
if ~isfield(atlas,'tissue') && isfield(atlas,'parcellation')
    atlas.tissue = atlas.parcellation;
    try
        atlas.tissuelabel = atlas.parcellationlabel;
    catch
    end
end

% Make atlas.color field if not present
if ~isfield(atlas, 'color')
    try
        atlas.color = linspecer(length(atlas.tissuelabel));
    catch
        atlas.color = repmat(0,length(atlas.tissuelabel),3);
    end
end
    
% Interpolate atlas
cfg2 = [];
cfg2.interpmethod = 'nearest';
cfg2.parameter = 'tissue';
sourcemodel2 = ft_sourceinterpolate(cfg2, atlas, cfg.template_grid);

% Fancy data viz
if strcmp(cfg.vis,'fancy')
    DT = delaunayTriangulation(cfg.template_grid.pos(cfg.template_grid.inside,:));

    DT_ = [];
    DT_.vertices = DT.Points;
    [K,v] = convexHull(DT);
    DT_.faces = K;


    figure; ft_plot_mesh(DT_,'facealpha',0.1,'edgecolor','none');
    view([0 74]);
    camlight;
    hold on;

elseif strcmp(cfg.vis,'yes')
    % Simple viz
    figure; ft_plot_mesh(cfg.template_grid.pos(cfg.template_grid.inside,:),'vertexsize',1);
    view([-56 45]);
    hold on;
end

%% Create VE
VE = [];
try
    VE.trialinfo = data_clean.trialinfo;
catch
end

count = 1;

for lab = 1:length(cfg.channel)

    fprintf('Computing VE for %10s\n',cfg.channel{lab});

    % Find atlas points
    try
        atlas_points = find(sourcemodel2.tissue==...
            find(contains(sourcemodel2.tissuelabel,cfg.channel{lab})));

        % Fancy data viz
        if strcmp(cfg.vis,'fancy')
            DT = delaunayTriangulation(cfg.template_grid.pos(atlas_points,:));

            DT_ = [];
            DT_.vertices = DT.Points;
            [K,~] = convexHull(DT);
            DT_.faces = K;

            ft_plot_mesh(DT_,'facealpha',0.3,'edgecolor','none','facecolor',[atlas.color(lab,:)]);
            drawnow; hold on;
            
        elseif strcmp(cfg.vis,'yes')
            % Simple data viz
            ft_plot_mesh(cfg.template_grid.pos(atlas_points,:),'vertexsize',10,...
                'vertexcolor',[atlas.color(lab,:)]); drawnow; hold on;
        end

        % Concatenate spatial filter

        try
            F    = cat(1,sourceall.avg.filter{atlas_points});
        catch
            fprintf(['Cannot find the avg.filter ... did you specfify ',...
                'cfg.keepfilter = yes \nwhen performing ft_sourceanalysis?\n']);
        end

        % Do SVD
        [u,s,v] = svd(F*avg_cov.cov*F');
        filter = u'*F;

        % Create VE
        for subs=1:(length(data_clean.trial))
            % note that this is the non-filtered "raw" data
            VE.time{subs}       = data_clean.time{subs};
            VE.trial{subs}(count,:) = filter(1,:)*data_clean.trial{subs}(:,:);
        end

        % Add label
        VE.label{count,1} = cfg.channel{lab};
        count = count+1;

    catch
        warning(['Could not find any atlas points for ROI: ' cfg.channel{lab}]);
    end

end







end