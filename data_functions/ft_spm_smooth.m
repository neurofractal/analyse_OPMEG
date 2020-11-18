function [source_smooth] = ft_spm_smooth(source,path_to_mesh,fwhm,spm_path)
%
% Smooth functional data on mesh using SPM
%
% EXAMPLE USEAGE:  [source_smooth] = ft_spm_smooth(source,...
%                                    path_to_mesh,fwhm,spm_path)
%
% INPUTS:
%
%   sourceall      = output of ft_sourceanalysis with a .pow field
%   path_to_mesh   = path to .gii used during source analysis
%   fwhm           = in mm
%   spm_path       = path to SPM12
%
% OUTPUTS:
%
%   source_smooth  = your smoothed
%
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from code by Tim Tierney
%
% Author: Robert Seymour      (rob.seymour@ucl.ac.uk)
%__________________________________________________________________________
%%

% Load mesh
brain = gifti(path_to_mesh);

% Get source power. Replace Nan with 0
r = source.pow;
r(isnan(r))=0;

% Add SPM
addpath(spm_path);

% distance to nearest neighbours for each vertex
[~,Di]=spm_mesh_neighbours(brain,1);
muD= mean(mean(Di));

% Tim made this step up so..... I'm not sure it's online anywhere
n= round((fwhm/muD)^2);
disp(n);

% smooth r on brain with n iterations
col=spm_mesh_smooth(brain,r,n);

source_smooth = source;
source_smooth.pow = col;

%%%%% 
% Dunno if I should do this?
%%%%% 

% Replace NaNs
r = source.pow;
source_smooth.pow(isnan(r))=NaN;

rmpath(spm_path);
