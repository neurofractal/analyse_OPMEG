function [data_out_AMM] = ft_wrapper_spm_amm(cfg,data)
% Fieldtrip wrapper for spm_opm_amm.m
%
%   cfg.path_to_SPM         = Path to SPM12
%   cfg.correctgrad         = Correct the forward model ('yes' or no';
%                             default = 'no'). This is experimental and
%                             should be used with caution.
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

%% Add paths
addpath('D:\scripts\tim_amm'); % FIXME Hardcoded for Rob at the moment
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
% FIXME Currently hard-coded
S               = [];
S.D             = data_SPM;
S.le            = 1;
S.corrLim       = 0.98;
data_SPM_AMM    = spm_opm_amm(S); % Contact Tim Tierney for these scripts

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
% disp('Adding Fieldtrip back your path...');
% addpath(ft_path);
%ft_defaults

% Delete SPM data
try
    delete data.dat
    delete data.mat
    delete mdata.dat
    delete mdata.mat
catch
end

% Subfunction to remove sensors from the grad structure not in the data
function dataout = fixsens(data1)
    chans2keep = ismember(data1.grad.label,data1.label);

    grad_mod         = data1.grad;

    % For doubles
    pot_grad_fields = {'chanori','chanpos','coilori','coilpos'};

    for p = 1:length(pot_grad_fields)
        if isfield(grad_mod,pot_grad_fields{p})
            eval(['grad_mod.' pot_grad_fields{p} ' = grad_mod.' pot_grad_fields{p} '(chans2keep,:)']);
        end
    end

    % For cells
    pot_grad_fields = {'chantype','chanunit','label'};

    for p = 1:length(pot_grad_fields)
        if isfield(grad_mod,pot_grad_fields{p})
            eval(['grad_mod.' pot_grad_fields{p} ...
                ' = grad_mod.' pot_grad_fields{p} '(chans2keep)']);
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

end
