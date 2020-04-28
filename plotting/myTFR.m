function myTFR(cfg,freqData)
% TFR plot function which smoothes the image and allows for log scale.
% Inputs should include:
% cfg.logscale  = 'yes' or 'no' (default = 'yes')
% cfg.latency   = [startTime endTime] or 'all' (default = 'all')
% cfg.frequency = [list of frequencies to be plotted] or 'all' (default = 'all')
% cfg.channel   = channel or list of channels to be plotted/averaged and
%                   plotted. e.g. {'CPz','Cz'}
% cfg.parameter = 'parameterName', default = 'powspctrm'
% cfg.zlim      = [min max];

% Select on the input parameters
cfg.avgoverchan     = 'yes';
selectedData        = ft_selectdata(cfg,freqData);

parameterData       = getfield(selectedData,cfg.parameter);

numInt              = length(cfg.colourmap);

%% setup contour thresholds and z axis.
if strcmp(cfg.zlim,'maxabs')
    minZ                = min(min(min(parameterData)));
    maxZ                = max(max(max(parameterData)));
    maxAbsZ             = max(abs([minZ maxZ]));
    colourLimits        = -maxAbsZ:((maxAbsZ*2)/numInt):maxAbsZ;
elseif strcmp(cfg.zlim,'maxmin')
    minZ                = min(min(min(parameterData)));
    maxZ                = max(max(max(parameterData)));
    lowerWidth          = abs(minZ/(numInt/2));
    upperWidth          = maxZ/(numInt/2);
    lowerLimits         = minZ:lowerWidth:(0-(lowerWidth/2));
    upperLimits         = (0+(upperWidth/2)):upperWidth:maxZ;
    colourLimits        = [lowerLimits,upperLimits];
elseif ~isstring(cfg.zlim)
    minZ                = cfg.zlim(1);
    maxZ                = cfg.zlim(2);
    lowerWidth          = abs(minZ/(numInt/2));
    upperWidth          = maxZ/(numInt/2);
    lowerLimits         = minZ:lowerWidth:(0-(lowerWidth/2));
    upperLimits         = (0+(upperWidth/2)):upperWidth:maxZ;
    colourLimits        = [lowerLimits,upperLimits];

else
    maxAbsZ             = max(cfg.zlim);
    colourLimits        = -maxAbsZ:((maxAbsZ*2)/numInt):maxAbsZ;
end
contourThresholds        = colourLimits;

% contourf(selectedData.time, selectedData.freq, squeeze(parameterData),contourThresholds);
% hold on
% surf(selectedData.time, selectedData.freq, zeros(size(squeeze(parameterData))), squeeze(parameterData), 'EdgeColor', 'none');
% contour(selectedData.time, selectedData.freq, squeeze(parameterData),cfg.contours);

contourfcmap(selectedData.time, selectedData.freq, squeeze(parameterData),contourThresholds,cfg.colourmap, ...
     'lo', [.8 .8 .8], ...
     'hi', [.2 .2 .2], ...
     'cbarloc', 'eastoutside', ...
     'method', 'calccontour' ,...
     'evencb',true);

if strcmp(cfg.logscale,'yes')
    set(gca,'yscale','log','ytick',[2 4 8 16 32 64],'TickDir','out');
    colormap(cfg.colourmap)
end

% if strcmp(cfg.zlim,'maxabs')
%     caxis([-maxAbsZ maxAbsZ])
% elseif strcmp(cfg.zlim,'maxmin')    caxis([minZ maxZ])
% else
%     caxis([-maxAbsZ maxAbsZ])
% end
% 
% if strcmp(cfg.colourbar,'yes')
%     colorbar
% end



