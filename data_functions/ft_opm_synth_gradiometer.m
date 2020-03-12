function [data_out] = ft_opm_synth_gradiometer(cfg,data)
% Function to regress reference signals from the optically-pumped
% magnetencephalography (OPMEG) data acquired from the 
% UCL Wellcome Centre for Neuroimaging.
%
% EXAMPLE USEAGE:   [data_out] = ft_opm_synth_gradiometer(cfg,data)
% ...where, cfg is the input structure
%
%   cfg.channel         = the channels to be denoised (default = 'MEG')
%   cfg.refchannel      = the channels used as reference signal 
%                       (default = 'MEGREF')
%   cfg.filter_ref      = filters to apply to the reference data. Enter in
%                       the form [0 30; 40 60; 110 130], where values 
%                       indicate frequency in Hz (0-30Hz, 40-60Hz, 
%                       110-130Hz). (Default = []).
%   cfg.derivative      = add the derivative of the reference signal to the
%                       regressors (default = 'no');
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_synth_gradiometer (Tim Tierney)

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Set default values
if ~isfield(cfg, 'channel')
    cfg.channel = 'MEG';
end

if ~isfield(cfg, 'refchannel')
    cfg.channel = 'MEGREF';
end

if ~isfield(cfg, 'filter_ref')
    cfg.filter_ref = [];
end

if ~isfield(cfg, 'derivative')
    cfg.derivative = 'no';
end

filter_ref  = cfg.filter_ref;
channel     = cfg.channel;
refchannel  = cfg.refchannel;

%% Select MEG data and reference data
disp('Selecting data channels');
if strcmp(channel,'all')
    
    ft_warning('Selecting MEG channels');
    
    cfg2 = [];
    cfg2.channel = ft_channelselection_opm('MEG',data);
    meg_data = ft_selectdata(cfg2,data);
    
else
    cfg2 = [];
    cfg2.channel = ft_channelselection_opm(channel,data)
    meg_data = ft_selectdata(cfg2,data);
end

disp('Selecting reference channels');
cfg2 = [];
cfg2.channel = ft_channelselection_opm(refchannel,data)
ref_data = ft_selectdata(cfg2,data);

%% 

megind = meg_data.trial{1};
megres = zeros(size(megind));
ref = ref_data.trial{1};

ref_size = size(ref,1);
winSize = size(megind,2);

% add a mean column to the reference regressors
intercept = ones(winSize,1);

%%
% If the user wants to filter the reference data...
if ~isempty(filter_ref)
    
    % Create array of zeros to hold the data
    reference = zeros(ref_size*size(filter_ref,1),winSize);
    
    % Bit of indexing (probably inefficient)
    indx = reshape(1:size(reference,1),ref_size,size(filter_ref,1));
    
    % For every pair of frequencies...
    for filt = 1:size(filter_ref,1)
        
        fprintf('Filtering reference data: %3dHz - %3dHz ... \n',...
            filter_ref(filt,1),filter_ref(filt,2));
        
        data_filt = ft_preproc_highpassfilter(ref,data.fsample,...
            filter_ref(filt,1),5,'but','twopass','reduce');
            
        data_filt = ft_preproc_lowpassfilter(data_filt,data.fsample,...
            filter_ref(filt,2),5,'but','twopass','reduce');
        
        reference(indx(1,filt):indx(end,filt),:) = data_filt;
    end
    
    reference = reference';
    
    
else

    reference = ref';
end

%% 
if strcmp(cfg.derivative,'yes')
    drefer = diff(reference);
    drefer = [drefer(1,:); drefer];
    reference = [drefer reference];
end

%%
% Here we are actually doing the regression
disp('Regressin''...');
reference= [reference ones(size(reference,1),1)];
% reference is column major so transpose sensors
beta = pinv(reference)*megind';
regressed_data = (megind'- reference*beta)';

%% Return the data!
data_out = meg_data;
data_out.trial{1} = regressed_data;













