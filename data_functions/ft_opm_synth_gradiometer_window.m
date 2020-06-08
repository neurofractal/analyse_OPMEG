function [data_out] = ft_opm_synth_gradiometer_window(cfg,data_in)
% Function to regress reference signals from the optically-pumped
% magnetoencephalography (OPMEG) data acquired from the
% UCL Wellcome Centre for Neuroimaging. 
%
% This function applies the regression on over-lapping windows, which has
% been shown to improve performance in some instances.
%
% EXAMPLE USEAGE:   [data_out] = ft_opm_synth_gradiometer_window(cfg,data)
% ...where, cfg is the input structure
%
%   cfg.channel         = the channels to be denoised (default = 'MEG')
%   cfg.refchannel      = the channels used as reference signal
%                       (default = 'MEGREF')
%   cfg.filter_ref      = filters to apply to the reference data. Enter in
%                       the form [40 60], where values
%                       indicate frequency in Hz (e.g.40-60Hz, 
%                       default = []).
%   cfg.derivative      = add the derivative of the reference signal to the
%                       regressors (default = 'no').
%   cfg.return_all      = 'yes' or 'no' - do you want to return all
%                       channels, not just those specified with 
%                       cfg.channel? (default = 'no').
%   cfg.winsize          = size of the window you wish to apply the DSSP
%                          (in seconds, default = 10)
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
    cfg.refchannel = 'MEGREF';
end

if ~isfield(cfg, 'filter_ref')
    cfg.filter_ref = [];
end

if ~isfield(cfg, 'derivative')
    cfg.derivative = 'no';
end

if ~isfield(cfg, 'return_all')
    cfg.return_all = 'no';
end

if ~isfield(cfg, 'winsize')
    cfg.winsize = 10;
end

%% Calculate window size in terms of number of data points
wsize = cfg.winsize*data_in.fsample;
fprintf('Using %d data-points per window\n',wsize);

%% Do I need to do this?
filter_ref          = cfg.filter_ref;
channel             = cfg.channel;
refchannel          = cfg.refchannel;

%% Select MEG data and reference data
% MEG data
if strcmp(channel,'all')
    
    ft_warning('Selecting MEG channels');
    
    cfg2            = [];
    cfg2.channel    = ft_channelselection_opm('MEG',data_in);
    meg_data        = ft_selectdata(cfg2,data_in);
    
else
    cfg2            = [];
    cfg2.channel    = ft_channelselection_opm(channel,data_in);
    meg_data        = ft_selectdata(cfg2,data_in);
end

fprintf('Selected %2d data channels\n',length(meg_data.label));

% Reference data
cfg2                = [];
cfg2.channel        = ft_channelselection_opm(refchannel,data_in);
ref_data            = ft_selectdata(cfg2,data_in);
fprintf('Selected %2d reference channels\n',length(ref_data.label));

% Save data for later
if strcmp(cfg.return_all,'yes')
    data_out        = data_in;
    ch_indx         = contains(data_in.label,meg_data.label);
else
    data_out        = meg_data;
end


%% Start of the trial loop
for tr = 1:numel(ref_data.trial)
    
    % Get MEG data
    megind          = meg_data.trial{tr};
    % Get reference data
    ref             = ref_data.trial{tr};
    
    ref_size        = size(ref,1);
    len_of_data     = size(megind,2);
    
    
    %% If the user wants to filter the reference data...
    if ~isempty(filter_ref)
        
        % Create array of zeros to hold the data
        reference = zeros(ref_size*size(filter_ref,1),len_of_data);
        
        % Bit of indexing (probably inefficient)
        indx = reshape(1:size(reference,1),ref_size,size(filter_ref,1));
        
        % For every pair of frequencies...tr
        for filt = 1:size(filter_ref,1)
            if tr == 1
                fprintf('Filtering reference data: %3dHz - %3dHz ... \n',...
                    filter_ref(filt,1),filter_ref(filt,2));
            end
            
            % High-pass filter, except where the user has specified 0
            if filter_ref(filt,1) ~= 0
            data_filt = ft_preproc_highpassfilter(ref,data_in.fsample,...
                filter_ref(filt,1),5,'but','twopass','reduce');
            else
                disp('NOT high-pass filtering');
                data_filt = ref;
            end
            
            % Low-pass filter, except where the user has specified 0
            if filter_ref(filt,2) ~= 0
                data_filt = ft_preproc_lowpassfilter(data_filt,data_in.fsample,...
                    filter_ref(filt,2),5,'but','twopass','reduce');
            else
                disp('Not low-pass filtering');
            end
            
            reference(indx(1,filt):indx(end,filt),:) = data_filt;
        end
        
        reference = reference';
        
    else
        
        reference = ref';
    end
    
    %% Derivative
    if strcmp(cfg.derivative,'yes')
        drefer      = diff(reference);
        drefer      = [drefer(1,:); drefer];
        reference   = [drefer reference];
    end
    
    %% Here we are actually doing the regression
    if tr == 1
        disp('Regressin''...');
    end
    
    reference       = [reference ones(size(reference,1),1)];

    %% Start of Window-Loop
    
    % Create array of zeros for data
    meg_ind_synth_grad  = zeros(size(megind));
    % Create array of zeros for triangular weighting
    a                   = zeros(1,size(megind,2));
    
    % Weighting? I could probably take this out later
    w=ones(size(megind));
    if size(w,1)==1; w=repmat(w,1,size(megind,1)); end
        
    % Start at 0
    offset=0;
    ft_progress('init', 'etf', 'Regressing...')
    while true
        ft_progress(offset/size(megind,2))        
        
        % Calculate start and stop points for this loop 
        start=offset+1;
        stop=min(size(megind,2),offset+wsize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These values are grown by a factor of 5
        % Is this needed??
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        counter=0;
        while any (sum(min(w(start:stop),1))) <wsize
            if counter <= 0 ; break; end
            start=max(1,start-wsize/2);
            stop=min(size(megind,2),stop+wsize/2);
            counter=counter-1;
        end
        if rem(stop-start+1,2)==1; stop=stop-1; end
        wsize2=stop-start+1;
        
       % Do gradiometry
        beta    = pinv(reference(start:stop,:))*megind(:,start:stop)';
        yy      = (megind(:,start:stop)' - reference(start:stop,:)*beta)';
        
        % triangular weighting (specified via b variable)
        if start==1
            b=[ones(1,wsize2/2)*wsize2/2, wsize2/2:-1:1];
        elseif stop==size(megind,2)
            b=[1:wsize2/2, ones(1,wsize2/2)*wsize2/2];
        else
            b=[1:wsize2/2, wsize2/2:-1:1];
        end
                
        % Add to meg_ind_synth_grad variable outside the loop, weighted by b
        meg_ind_synth_grad(:,start:stop)=meg_ind_synth_grad(:,start:stop)+bsxfun(@times,yy,b);
        
        % Add triangular weighting to variable outside the loop
        a(1,start:stop)=a(start:stop)+b;
        
        % Adjust offset parameter by window size divided by 5
        offset=offset+wsize/5;
        
        % If we have reached the end of the data BREAK 
        if offset>size(megind,2)-wsize/5; break; end
    end
    ft_progress('close')
    
    % Adjust for triangular weighting
    meg_ind_synth_grad=bsxfun(@times,meg_ind_synth_grad,1./a); 
    % Find any NaN values and convert to 0
    meg_ind_synth_grad(isnan(meg_ind_synth_grad))=0;

    %% Return the data!
    if strcmp(cfg.return_all,'yes')
        data_out.trial{tr}(ch_indx,:) = meg_ind_synth_grad;
    else
        data_out.trial{tr} = meg_ind_synth_grad;
    end
    
    
end












