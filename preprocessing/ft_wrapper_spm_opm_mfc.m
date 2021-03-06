function [data_out_mfc] = ft_wrapper_spm_opm_mfc(cfg,data)
% Fieldtrip wrapper for spm_opm_mfc.m
% Please see https://github.com/tierneytim/OPM/blob/master/spm_opm_mfc.m
%
%   cfg.path_to_SPM         = Path to SPM12
%   cfg.path_to_OPM_repo    = Path to OPM repository which can be
%                           downloaded from:
%                           https://github.com/tierneytim/OPM
%   cfg.correctgrad         = Correct the forward model ('yes' or no';
%                             default = 'no'). This is experimental and
%                             should be used with caution.
%
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
%
% Author for this wrapper:      Robert Seymour (rob.seymour@ucl.ac.uk)
% Author of original SPM code:  Tim Tierney
%
% Citation: Tierney, T. M., Alexander, N., Mellor, S., Holmes, N.,
% Seymour, R., O'Neill, G. C., ... & Barnes, G. R. (2020). Modelling
% optically pumped magnetometer interference as a mean (magnetic) field.
% bioRxiv.
%__________________________________________________________________________

%% Set default values
if ~isfield(cfg, 'correctgrad')
    cfg.correctgrad = 'no';
end

if ~isfield(cfg, 'downsample')
    cfg.downsample = 'no';
end

%% Check the user has removed all non-MEG channels
if length(data.label) ~= length(data.grad.label)
    error('Different number of data channels between trial and grad structures');
end

%% Get Fieldtrip path and remove from your path
[~, ft_path] = ft_version;
% Little hack to turn off warnings
id = 'MATLAB:rmpath:DirNotFound';
warning('off',id);
% Remove Fieldtrip path and external/spm12 and external/spm8
rmpath(ft_path);
rmpath(fullfile(ft_path,'external','spm12'));
rmpath(fullfile(ft_path,'external','spm8'));

%% Add SPM
addpath(cfg.path_to_SPM);
spm('defaults', 'eeg');

%% Some jiggery pokery to make it work
% Change chantype of megmag
for i = 1:length(data.grad.chantype)
    data.grad.chantype{i,1} = 'megmag';
end

% Change chanunit to 'T' ??? Not sure why this works?
for ii = 1:length(data.grad.label)
    data.grad.chanunit{ii} = 'T';
end

%% Convert data from FT to SPM
[data_SPM] = spm_eeg_ft2spm_OPM(data,'data');

data_ft = data_SPM.ftraw;

% Turn warnings back on
warning('on',id);

%% Downsample
if isnumeric(cfg.downsample)
    disp('Downsampling using spm_eeg_downsample...');
    warning('Results might be suboptimal using this method...');
    S        = [];
    S.D      = data_SPM;
    S.fsample_new = cfg.downsample;
    data_SPM = spm_eeg_downsample(S)
end

%% Check the conversion has happened properly 
% (Ignore when the data has been downsampled)
if ~isnumeric(cfg.downsample)
    data_ft = data_SPM.ftraw;
    x = data.trial{1};
    y = data_ft.trial{1};
    
    if ~isequal(x,y)
        warning('FT to SPM conversion has not worked... Look at the data type conversion')
    end
    clear x y
end
%% Add Tim's OPM repo
try
    addpath(cfg.path_to_OPM_repo);
catch
    warning('Did you specify the OPM repo correctly?');
end

%% Perform MFC
disp('Performing MFC...');
S           = [];
S.D         = data_SPM;
if strcmp(cfg.correctgrad,'yes')
    S.balance   = 1;
else
    S.balance   = 0;
end
[data_SPM_mfc,Yinds] = spm_opm_mfc(S);

%% Convert from SPM to Fieldtrip
data_out = data_SPM_mfc.ftraw;

%% Undoing the jiggery pokery
% Change chanunit back to 'fT'
for ii = 1:length(data.grad.label)
    data_out.grad.chanunit{ii} = 'fT';
end

%% Return the data
% SPM seems to lose some of the original Fieldtrip information, so the
% safest thing to do is replace the trial (and grad) fields + return
data_out_mfc = data;
for i = 1:length(data.trial)
    data_out_mfc.trial{i} = data_out.trial{i};
end

if strcmp(cfg.correctgrad,'yes')
    disp('Replacing grad structure corrected for MFC...');
    data_out_mfc.grad = data_out.grad;
end

% Replace these values if the data has been downsampled
if isnumeric(cfg.downsample)
     data_out_mfc.time{i} = data_out.time{i};
     data_out_mfc.fsample = cfg.downsample;
end
    
%% Remove SPM, Tim's OPM repo and re-add Fieldtrip
disp('Removing SPM, OPM from your path...');
close all force
warning('off',id);
rmpath(genpath(cfg.path_to_SPM));
rmpath(cfg.path_to_OPM_repo);
warning('on',id);
disp('Adding Fieldtrip back your path...');
addpath(ft_path);
ft_defaults

% Delete SPM data
try
    if isnumeric(cfg.downsample)
        delete data.dat
        delete data.mat
        delete ddata.dat
        delete ddata.mat
        delete MF_ddata.dat
        delete MF_ddata.mat   
    else
        delete data.dat
        delete data.mat
        delete MF_data.dat
        delete MF_data.mat
    end
    
catch
    warning('Could not delete SPM files...');
end

end
