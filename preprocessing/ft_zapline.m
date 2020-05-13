function [data_zapline] = ft_zapline(cfg,data)
% Function to perform robust detrending (via the NoiseTools) on data
% arranged in a Fieldtrip structure
%
% EXAMPLE USEAGE:   data_zapped = ft_zapline(cfg,data)
% ...where, cfg is the input structure
% 
%   cfg.ln_freq           = Line noise (50 or 60Hz)
%   cfg.n_remove          = number of components to remove. Increase until
%                           good result is achieved
%   cfg.truncate_PC       = number of principal components for truncating
%                           (default = 1-number_of_chans)
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
if ~isfield(cfg, 'ln_freq')
    cfg.ln_freq = 50;
end

if ~isfield(cfg, 'n_remove')
    cfg.n_remove = 1;
end

if ~isfield(cfg, 'truncate_PC')
    truncate_PC = length(data)-1;
end

%%
data_zapline    = data;
line_noise      = cfg.ln_freq/data.fsample;
p=[];
p.nkeep=cfg.truncate_PC; 

% For every trial
disp('Zaplining...');
for t = 1:length(data.trial)
    data_4_zapline = data.trial{t}';
    data_4_zapline = nt_demean(data_4_zapline);
    data_out = nt_zapline(data_4_zapline,line_noise,cfg.n_remove);
    data_zapline.trial{t} = transpose(data_out);
end
    
end
