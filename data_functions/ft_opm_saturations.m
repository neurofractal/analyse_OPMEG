function [sat] = ft_opm_saturations(cfg,data)
% Function to detect times/channels with periods of railing / saturations
% from optically-pumped magnetencephalography (OPMEG) data
% acquired from the UCL Wellcome Centre for Neuroimaging.
%
% EXAMPLE USEAGE:   sat = ft_opm_saturations(cfg,data)
% ...where, cfg is the input structure and rawData is the raw OPM data
% loaded using ft_opm_create
% 
%   cfg.channel        = 'all', 'MEG', 'RAD', 'TAN'. Default = 'MEG'.
%   cfg.plot            = 'yes' or 'no'
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% Adapted from spm_opm_create (Tim Tierney)

% Authors:  Nicholas Alexander  (n.alexander@ucl.ac.uk)
%           Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Set default values
if ~isfield(cfg, 'channel')
    cfg.channel = 'MEG';
end

if ~isfield(cfg, 'plot')
    cfg.plot = 'yes';
end

%% Select the data based on cfg.channel option
if strcmp(cfg.channel,'all')
    disp('Detecting Saturations for ALL channels');
else
    try
        chan = cfg.channel;
        cfg2 = [];
        cfg2.channel = ft_channelselection_opm(chan,data,...
            'quspin_g2');
        data = ft_selectdata(cfg2,data);
        
    catch
        ft_warning(['Could not select ' chan]);
    end
end

%%
method      = cfg.method;
min_length  = cfg.min_length;
plot        = cfg.plot;

%% Start of the code 

sat = [];
count = 1;

% Use ft_progress to track the progress!
disp('Searching the data channel by channel for signal saturation...');
ft_progress('init', 'text', 'Please wait...')

% Split the data into 0.1s segments (unsure of optimal time?)
nsamps = 0.1*data.fsample;
beg = 1:nsamps:size(data.trial{1},2);
endsamp =  beg+(nsamps-1);
inRange = ~(beg>size(data.trial{1},2)|endsamp>size(data.trial{1},2));

% For every channel
for chan = 1:length(data.label)
    
    % Display the progress of the function
    ft_progress(chan/length(data.label),...
        'Searching for saturations: %s',data.label{chan});
    
    % Get data
    ttt = data.trial{1}(chan,:);
    
    % Calculate the range of the data for all segments
    rng = zeros(length(beg),1);
    
    for t = 1:length(beg)
        if inRange(t)
            seg = ttt(beg(t):endsamp(t));
        else
            seg = ttt(beg(t):end);
        end
        
        rng(t) = max(seg)-min(seg);
    end
    
    %figure; plot(rng);
    
    % Find values below 5000 (again this is kind of arbitary but works well
    % so far...
    find_vals = find(rng>5e3);
    
    % Create an array of ones
    time_sat = ones(length(data.time{1}),1)*1;
    
    % Change values to 0 when there are saturations
    for vals = 1:length(find_vals)
        time_sat(beg(find_vals(vals)):endsamp(find_vals(vals))) = 0;
    end
    
    % Make sure it's the right length (final sample sometimes cut off)
    time_sat = time_sat(1:length(data.time{1}));
    % Make it logical (like a vulcan)
    time_sat = logical(time_sat);
    
    % Get the saturated tikes
    time_sat2 = data.time{1}(time_sat);
    
    % Add this to the sat array
    if ~isempty(time_sat2)
        sat.label{count,1}      = data.label{chan};
        sat.time{count,1}       = time_sat2;
        sat.time_log{count,1}   = time_sat;
        count = count+1;
    end
end

ft_progress('close');




if ~isempty(sat)
    
    %% Find the overall times when the data is saturated (on any channel)
    all_sat = zeros(length(sat.label),length(data.trial{1,1}));
    
    for i = 1:length(sat.label)
        all_sat(i,1:length(sat.time{i,1})) = sat.time{i,1};
    end
    
    unique_all_sat = unique(all_sat);
    unique_all_sat(1) = [];
    
    sat.alltime = unique_all_sat';
    
    clear all_sat;
    
    % Get only data from saturated channels
    cfg = [];
    cfg.channel = sat.label;
    data_saturations = ft_selectdata(cfg,data);
    
    %% Plot how much of the data is saturated
    
    if strcmp(plot,'yes')
        
        time_saturated = [];
        
        for r = 1:length(sat.label)
            time_saturated(r) = length(sat.time{r})./data.fsample;
        end
        
        % Now add the TOTAL time saturated over any channel
        time_saturated(length(time_saturated)+1) = length(sat.alltime)./...
            data.fsample;
        
        figure;
        set(gcf,'Position',[100 100 900 800]);
        stem(time_saturated,'r','LineWidth',2);
        set(gca,'xtick',[1:length(time_saturated)],'xticklabel',...
            vertcat(sat.label, 'TOTAL'))
        ylabel('Time Saturated (s)','FontSize',20);
        ax = gca;
        
        if length(sat.label) > 50
            ax.XAxis.FontSize = 7;
        else
            ax.XAxis.FontSize = 10;
        end
        
        ax.YAxis.FontSize = 16;
        view(90,90)
    end

end

