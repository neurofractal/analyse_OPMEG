function [data_detrend] = ft_robust_detrend(cfg,data)
% Function to perform robust detrending (via the NoiseTools) on data
% arranged in a Fieldtrip structure
%
% EXAMPLE USEAGE:   data_detrend = ft_opm_create(cfg,data)
% ...where, cfg is the input structure
% 
%   cfg.poly_num          = order of polynomials (default = 3)
%   cfg.thresh            = threshold for outliers (default = 3)
%   cfg.niter             = number of iterations performed (default = 3)
%   cfg.wsize             = size of window to perform detrending 
%                          (default = [], i.e. the whole trial) 
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging

% Authors: Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

% Check if NoiseTools is on your path 
try
    nt_greetings;
catch
    error(['You need NoiseTools on your MATLAB path. '...
        'Please visit http://audition.ens.fr/adc/NoiseTools']);
end

%% Set default values
if ~isfield(cfg, 'poly_num')
    cfg.poly_num = 3;
end

if ~isfield(cfg, 'thresh')
    cfg.thresh = 3;
end

if ~isfield(cfg, 'niter')
    cfg.niter = 3;
end

if ~isfield(cfg, 'wsize')
    cfg.wsize = [];
end

%%
data_detrend = data;

% For every trial
disp('Detrending...');
for t = 1:length(data.trial)
    % Get data and conjugate
    data_4_detrend = data.trial{t}';
    
    % Detrend
    [data_out] = nt_detrend(data_4_detrend,cfg.poly_num,[],...
        'polynomials',cfg.thresh,cfg.niter,cfg.wsize);

    data_detrend.trial{t} = transpose(data_out);
end
    
end
