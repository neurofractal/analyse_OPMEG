function [niftiF, smoothNiftiF, leftGiftiF, rightGiftiF,minMaxRange] = exportFtStats(statStruct, outputFolder, outputFilename, cwDir, lSurf, rSurf)
% Checks fieldtrip stat function and makes an unthresholded and masked
% version. Exports both as nifti images and smoothes them. For The masked
% image, smoothing is only applied within the mask. i.e., the extents of
% the mask are not changed. The smooth nifti is then projected onto a mesh
% (gifti) with some useful palette settings pre-applied.
%
% Example use:
% numSubjects = length(data1);
% cfg                     = [];
% cfg.dim                 = data1.dim;
% cfg.method              = 'analytic';
% cfg.alpha				  = 0.001;
% cfg.correctm			  = 'no';
% cfg.statistic           = 'depsamplesT';
% cfg.parameter           = 'pow';
% cfg.correcttail		  = 'alpha';
% cfg.tail                = 0;
% cfg.design(1,:)         = [1:numSubjects, 1:numSubjects];
% cfg.design(2,:)         = [ones(1,numSubjects), ones(1,numSubjects)*2];
% cfg.uvar                = 1;
% cfg.ivar                = 2; 
% cfg.computecritval	  = 'yes';
% sourceStats			  = ft_sourcestatistics(cfg, data1{:}, data2{:});
% wbDir					  = 'C:\workbench-windows64-v1.5.0\workbench\bin_windows64';
% lSurf					  = 'C:\YOUR\FOLDER\Conte69_atlas-v2.LR.32k_fs_LR.wb\Conte69.L.white.32k_fs_LR.surf.gii';
% rSurf					  = 'C:\YOUR\FOLDER\Conte69_atlas-v2.LR.32k_fs_LR.wb\Conte69.R.white.32k_fs_LR.surf.gii';
% outputFolder			  = 'C:\results';
% filename				  = 'theta_0_4';
% [~, smoothNiftiF, ~, ~,minMaxRange]	= exportFtStats(sourceStats, outputFolder, filename, cwDir, lSurf, rSurf);
% 
% The output references to filenames and ranges etc are useful for feeding
% into volumetric plotting functions. e.g. slicer:
% slicer({2, smoothNiftiF.Masked},...
% 		'limits', {[], [-minMaxRange minMaxRange]},...
% 		'minClusterSize', {0,0},....
% 		'alpha', {1, 0.66},...
% 		'colormaps', {1, 49},...
% 		'cbLocation', 'best',...
% 		'title', filename,...
% 		'output', [outputFolder,'\',filename, '_mixed_ax'],...
% 		'resolution', 600,...
% 		'slices', 10:4:80,...
% 		'view', 'ax')
%
% Author: Nicholas Alexander (n.alexander@ucl.ac.uk)

% Useful to assign this to a field
stat = [];
stat.Unthresholded = statStruct;

% Some work needed to produce a mask for cluster stats
if isfield(stat.Unthresholded,'posclusters') || isfield(stat.Unthresholded,'negclusters')
	% Create a mask for significant clusters
	mask = zeros(size(stat.Unthresholded.mask));
	
	% There is not always a positive/negative cluster struct
	if isfield(stat.Unthresholded,'posclusters') && ~isempty(stat.Unthresholded.posclusters)
		posClusterIdx = find([stat.Unthresholded.posclusters.prob] <= stat.Unthresholded.cfg.clusteralpha);
	else
		posClusterIdx = [];
	end
	if isfield(stat.Unthresholded,'negclusters') && ~isempty(stat.Unthresholded.negclusters)
		negClusterIdx = find([stat.Unthresholded.negclusters.prob] <= stat.Unthresholded.cfg.clusteralpha);
	else
		negClusterIdx = [];
	end
	if ~isempty(posClusterIdx)
		posMask = zeros(length(mask), length(posClusterIdx));
		for i = 1:length(posClusterIdx)
			posMask(:,i) = stat.Unthresholded.posclusterslabelmat == posClusterIdx(i);
		end
	else
		posMask = zeros(length(mask),1);
	end
	if ~isempty(negClusterIdx)
		negMask = zeros(length(mask), length(negClusterIdx));
		for i = 1:length(negClusterIdx)
			negMask(:,i) = stat.Unthresholded.negclusterslabelmat == negClusterIdx(i);
		end
	else
		negMask = zeros(length(mask),1);
	end

	% The masks can be combined
	mask = sum([posMask,negMask],2);

	% These fields cause an error with interpolation, so remove them
	try
		stat.Unthresholded = rmfield(stat.Unthresholded,'posclusters');
	catch
	end
	try
		stat.Unthresholded = rmfield(stat.Unthresholded,'negclusters');
	catch 
	end
	stat.Masked = stat.Unthresholded;
	stat.Masked.stat = stat.Unthresholded.stat .* mask;
else
	% Masked
	stat.Masked = stat.Unthresholded;
	stat.Masked.stat = stat.Unthresholded.stat .* stat.Unthresholded.mask;
end

% Interpolate to template
fieldTripPath = fileparts(which('ft_defaults'));
templateMri = ft_read_mri([fieldTripPath,'\template\anatomy\single_subj_T1.nii']);

statI = [];
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = {'stat'};
cfg.interpmethod = 'nearest';
statI.Masked = ft_sourceinterpolate(cfg, stat.Masked, templateMri);
statI.Unthresholded = ft_sourceinterpolate(cfg, stat.Unthresholded, templateMri);

statFields = fieldnames(statI);
for fieldIdx = 1:length(statFields)
	niftiF.(statFields{fieldIdx}) = fullfile(outputFolder, [statFields{fieldIdx},'_',outputFilename,'.nii']);

	cfg = [];
	cfg.filetype = 'nifti';
	cfg.filename = niftiF.(statFields{fieldIdx});
	cfg.parameter = 'stat';
	ft_sourcewrite(cfg,statI.(statFields{fieldIdx}));

	% Smooth
	smoothNiftiF.(statFields{fieldIdx}) = fullfile(outputFolder, [statFields{fieldIdx},outputFilename,'_smooth_','.nii']);

	commandLine = strcat([cwDir, '\\wb_command -volume-smoothing  ', ...
				niftiF.(statFields{fieldIdx}), ...
				' 4 ', ...
				smoothNiftiF.(statFields{fieldIdx}),...
				' -fix-zeros -fwhm']);
	system(commandLine)

	% Left hemisphere
	leftGiftiF.(statFields{fieldIdx}) = fullfile(outputFolder, [statFields{fieldIdx},outputFilename,'_smooth_left','.shape.gii']);
	commandLine = [cwDir, '\\wb_command -volume-to-surface-mapping ', ...
				niftiF.(statFields{fieldIdx}), ' ',...
				lSurf, ' ', ...
				leftGiftiF.(statFields{fieldIdx}) ' -trilinear'];
	system(commandLine)

	% Right hemisphere
	rightGiftiF.(statFields{fieldIdx}) = fullfile(outputFolder, [statFields{fieldIdx},outputFilename,'_smooth_right','.shape.gii']);
	commandLine = strcat([cwDir, '\\wb_command -volume-to-surface-mapping ', ...
				niftiF.(statFields{fieldIdx}), ' ',...
				rSurf, ' ', ...
				rightGiftiF.(statFields{fieldIdx}) ' -trilinear']);
	system(commandLine)

	% Sort out colourmap usding a scaled back maxabs value due to smoothing
	originalMinMaxRange = max(max(max(abs(statI.(statFields{fieldIdx}).stat))));
	smoothedData = ft_read_mri(niftiF.(statFields{fieldIdx}));
	minMaxRange = max(max(max(abs(smoothedData.anatomy)))) * 0.9;
	originalToSmoothRatio = minMaxRange/originalMinMaxRange;
	criticalVal = originalToSmoothRatio * max(abs(stat.(statFields{fieldIdx}).critval));

	if ~isnan(criticalVal) || criticalVal == 0 || minMaxRange == 0 || ~isnan(minMaxRange)
		commandLine = strcat([cwDir, '\\wb_command -metric-palette ', ...
					leftGiftiF.(statFields{fieldIdx}), ...
					' MODE_USER_SCALE -palette-name spectral ' ...
					'-neg-user 0 ', num2str(-minMaxRange), ' ', ...
					'-pos-user 0 ', num2str(minMaxRange), ' ',...
				 	'-thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE ', ...
				 	num2str(-criticalVal), ' ', num2str(criticalVal)]);
		system(commandLine)

		commandLine = strcat([cwDir, '\\wb_command -metric-palette ', ...
					rightGiftiF.(statFields{fieldIdx}), ...
					' MODE_USER_SCALE -palette-name spectral ' ...
					'-neg-user 0 ', num2str(-minMaxRange), ' ', ...
					'-pos-user 0 ', num2str(minMaxRange), ' ',...
				 	'-thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE ', ...
				 	num2str(-criticalVal), ' ', num2str(criticalVal)]);
		system(commandLine)
	end
end













