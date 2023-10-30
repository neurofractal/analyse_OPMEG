function [data_out_AMM] = ft_wrapper_spm_amm(cfg,data)
% Fieldtrip wrapper for spm_opm_amm.m
%
%   cfg.path_to_SPM         = Path to SPM12
%   cfg.correctgrad         = Correct the forward model ('yes' or no';
%                             default = 'no'). This is experimental and
%                             should be used with caution.
%	cfg.corrLim				= default 1, adjust to 0.98 or 0.9 as necessary.
%	cfg.li					= internal harmonic order (default = 9)
%	cfg.le					= external harmonic order (default = 2)
%	cfg.window				= temporal window size in seconds (default = 10)
% 
% Copyright (C) 2023 Wellcome Trust Centre for Neuroimaging
%
% Author for this wrapper:      Robert Seymour (rob.seymour@ucl.ac.uk)
% Author of original AMM code:  Tim Tierney
%
% CITATION TO FOLLOW
%__________________________________________________________________________

%% Set default values
if ~isfield(cfg, 'correctgrad')
    cfg.correctgrad = 'yes';
end

<<<<<<< HEAD
if ~isfield(cfg, 'corrLim')
	cfg.corrLim = 1;
end

if ~isfield(cfg, 'li')
	cfg.li = 9;
end

if ~isfield(cfg, 'le')
	cfg.le = 2;
end

if ~isfield(cfg, 'window')
	cfg.window = 10;
end

=======
if ~isfield(cfg, 'li')
    cfg.li = 9;
end

if ~isfield(cfg, 'le')
    cfg.le = 2;
end

if ~isfield(cfg, 'corrLim')
    cfg.corrLim = 0.98;
end


>>>>>>> 166b81613d6d57bc8f5b51dbd317d811d4393b51
%% Add paths
% Check if fieldtrip is currently in the path
try 
	[~, ftpath] = ft_version;
	rmpath(ftpath);
catch
end

addpath(cfg.path_to_SPM)
warning('off','all');
spm('defaults', 'eeg')

%% Remove sensors from grad structure not in data (i.e. usually bad chans)
try
    data = fixsens(data);
catch
    warning('on','all');
    warning('Could not remove sensors from grad structure not present in data');
end

%% Convert data from FT to SPM
disp('Converting data to SPM format...')
[data_SPM] = spm_eeg_ft2spm_OPM(data,'data');
warning('on','all');

%% AMM
S               = [];
S.D             = data_SPM;
<<<<<<< HEAD
S.corrLim       = cfg.corrLim;
S.li			= cfg.li;
S.le			= cfg.le;
S.window		= cfg.window;
try
	data_SPM_AMM    = spm_opm_amm(S); % Contact Tim Tierney for these scripts
catch
	error('Please check that you have the latest version of SPM which includes spm_opm_amm.');
end
=======
S.le            = cfg.le;
S.li            = cfg.li;
S.corrLim       = cfg.corrLim;
data_SPM_AMM    = spm_opm_amm(S); % Contact Tim Tierney for these scripts
>>>>>>> 166b81613d6d57bc8f5b51dbd317d811d4393b51

%% To Fieldtrip
data_out = data_SPM_AMM.ftraw;

%% Return the data
% SPM seems to lose some of the original Fieldtrip information, so the
% safest thing to do is replace the trial (and grad) fields + return
data_out_AMM = data;
for i = 1:length(data.trial)
    data_out_AMM.trial{i} = data_out.trial{i};
end

if strcmp(cfg.correctgrad,'yes')
    disp('Replacing grad structure');
    data_out_AMM.grad = data_out.grad;
end

%% Remove SPM and delete
close all force
warning('off','all');
rmpath(genpath(cfg.path_to_SPM));
warning('on','all');

% Delete SPM data
try
    delete data.dat
    delete data.mat
    delete mdata.dat
    delete mdata.mat
catch
end


% Add fieldtrip back
try
	addpath(ftpath);
	ft_defaults;
catch
end


end


% Subfunction to remove sensors from the grad structure not in the data
function dataout = fixsens(data1)
    chans2keep = ismember(data1.grad.label,data1.label);

    grad_mod         = data1.grad;

    % For doubles
    pot_grad_fields = {'chanori','chanpos','coilori','coilpos'};

    for p = 1:length(pot_grad_fields)
        if isfield(grad_mod,pot_grad_fields{p})
            % Define the field name to update
			fieldName = pot_grad_fields{p};
			
			% Check if the field exists in grad_mod
			if isfield(grad_mod, fieldName)
    			% Access and update the field directly
    			grad_mod.(fieldName) = grad_mod.(fieldName)(chans2keep, :);
			else
    			% Handle the case when the field does not exist
    			disp(['Field "', fieldName, '" does not exist in grad_mod.']);
			end
        end
    end

    % For cells
    pot_grad_fields = {'chantype','chanunit','label'};

    for p = 1:length(pot_grad_fields)
        if isfield(grad_mod,pot_grad_fields{p})
			fieldName = pot_grad_fields{p};
			
			if isfield(grad_mod, fieldName)
    			grad_mod.(fieldName) = grad_mod.(fieldName)(chans2keep);
			else
    			disp(['Field "', fieldName, '" does not exist in grad_mod.']);
			end
        end
    end

    % For tra matrix
    if isfield(grad_mod,'tra')
        grad_mod.tra = grad_mod.tra(chans2keep,chans2keep);
    end

    % Replace
    dataout=data1;
    dataout.grad = grad_mod;
end
