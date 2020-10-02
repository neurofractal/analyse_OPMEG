function [data_zapline] = ft_zapline_window(cfg,data)
% Function to apply the Zapline algorithm (via NoiseTools) on overlapping
% windows of electrophysiological data arranged in a Fieldtrip structure
% (Please cite the original paper: 
% https://doi.org/10.1016/j.neuroimage.2019.116356)
%
% EXAMPLE USEAGE:   data_zapped = ft_zapline_window(cfg,data)
% ...where, cfg is the input structure
% 
%   cfg.winsize           = Length of window for applying the Zapline
%                           algorithm
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
if ~isfield(cfg, 'winsize')
    cfg.winsize = 10;
end

if ~isfield(cfg, 'ln_freq')
    cfg.ln_freq = 50;
end

if ~isfield(cfg, 'n_remove')
    cfg.n_remove = 1;
end

if ~isfield(cfg, 'truncate_PC')
  cfg.truncate_PC = size(data.trial{1},2)-1;
end

%%
% Calculate window size in terms of number of data points
wsize = cfg.winsize*data.fsample;
fprintf('Using %d data-points per window\n',wsize);

% Make copy of the data
data_zapline    = data;

% Get line_noise parameter for zapline
line_noise      = cfg.ln_freq/data.fsample;

% Get truncate principle components parameter for zapline
p=[];
p.nkeep=cfg.truncate_PC;

% For every trial
for trial = 1:length(data.trial)
    % Get data for this trial
    x = data.trial{trial}';
    % Demean
    x = nt_demean(x);
    
    % Create array of zeros for data
    y=zeros(size(x));
    % Create array of zeros for triangular weighting
    a=zeros(size(x,1),1);
    
    % Weighting? I could probably take this out later
    w=ones(size(x));
    if size(w,2)==1; w=repmat(w,1,size(x,2)); end
    
    
    % Start at 0
    offset=0;
    
    % Display progress using ft_progress
    ft_progress('init', 'etf', 'Zapping...')

    while true
        ft_progress(offset/size(x,1))        
        
        % Calculate start and stop times
        start=offset+1;
        stop=min(size(x,1),offset+wsize);
        
        if rem(stop-start+1,2)==1; stop=stop-1; end
        wsize2=stop-start+1;
        
        % Apply Zapline
        yy = nt_zapline(x(start:stop,:),line_noise,cfg.n_remove,p);
        
        % triangular weighting
        if start==1
            b=[ones(1,wsize2/2)*wsize2/2, wsize2/2:-1:1]';
        elseif stop==size(x,1)
            b=[1:wsize2/2, ones(1,wsize2/2)*wsize2/2]';
        else
            b=[1:wsize2/2, wsize2/2:-1:1]';
        end
        
        % overlap-add to output
        y(start:stop,:)=y(start:stop,:)+bsxfun(@times,yy,b);
        
        % Add triangular weighting to output
        a(start:stop,1)=a(start:stop)+b;
        
        % Adjust offset parameter by window size divided by 5
        offset=offset+wsize/5;
        
        % If we have reached the end of the data BREAK 
        if offset>size(x,1)-wsize/5; break; end
    end
    ft_progress('close')

    % Adjust triangular weighting
    y=bsxfun(@times,y,1./a); 
    % Find any NaN values and convert to 0
    y(isnan(y))=0;
    
    % Put the data back into Fieldtrip structure
    data_zapline.trial{trial} = transpose(y);
end
end



