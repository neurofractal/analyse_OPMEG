function [VE] = atlas2VE(atlas,template_grid,labels, data_clean,...
    sourceall, avg_cov, headmodel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mq_atlas2VE: a function to create a virtual electrode using an atlas!
%
% Author: Robert Seymour June 2019 (robert.seymour@mq.edu.au)
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
% - atlas           = atlas you wish to use for VE computation. Should have
%                   been read into MATLAB using ft_read_atlas
% - template_grid   = template_grid used for source analysis (read into
%                   MATLAB from .../fieldtrip-XXXXX/template/sourcemodel/
%                   standard_sourcemodel3dXmm.mat
% - labels          = name of the ROIs in your atlas
% - data_clean      = your clean sensor-level data
% - source_all      = result of ft_sourceanalysis with the filters you wish
%                   to use for your VE computation
% - avg_cov         = covariance matrix used for source analysis
%                   (output of ft_timelockanalysis with .cov field)
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
% - VE              = virtual electrode with label, time, trial and
%                   trialinfo fields
%
% Example: [VE] = mq_atlas2VE(atlas,template_grid,labels,data_clean,...
%                 sourceall,avg_cov);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If 'all' specified use all the tissuelabels
if strcmp(labels,'all')
    labels = atlas.tissuelabel;
end

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodel2 = ft_sourceinterpolate(cfg, atlas, template_grid);

figure; ft_plot_mesh(template_grid.pos(template_grid.inside,:),'vertexsize',1);
view([-56 45]);
hold on;


% Create VE
VE = [];
try
    VE.trialinfo = data_clean.trialinfo;
catch
end

count = 1;

for lab = 1:length(labels)

    fprintf('Computing VE for %10s\n',labels{lab});

    % Find atlas points
    try
        atlas_points = find(sourcemodel2.tissue==...
            find(contains(sourcemodel2.tissuelabel,labels{lab})));

        ft_plot_mesh(template_grid.pos(atlas_points,:),'vertexsize',10,...
            'vertexcolor','r'); drawnow; hold on;

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
        VE.label{count,1} = labels{lab};
        count = count+1;

    catch
        warning(['Could not find any atlas points for ROI: ' labels{lab}]);
    end

end







end