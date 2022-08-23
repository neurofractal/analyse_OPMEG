function trl = OPM_trialfun_usemat(cfg)
% OPM_trialfun_usemat is a custom trial function for data loaded into
% memory
%
% Copyright (C) 2021-22, Robert Seymour
% Wellcome Centre for Human Neuroimaging, UCL

% set the defaults
cfg.threshold               = ft_getopt(cfg, 'threshold', 3);
cfg.detectflank    =  ft_getopt(cfg, 'detectflank', 'up');

disp(cfg.detectflank)

% Find channel with trigger info
pos_of_trig = contains(cfg.rawData.label,cfg.trialdef.trigchan);


% % Read header
% hdr              = ft_read_header(strcat(cfg.dataset(1:end-3),'mat'));
% 
% % Read continuous event channel.
% cfgEvent            = [];
% cfgEvent.continuous = 'yes';
% cfgEvent.datafile   = cfg.dataset;
% cfgEvent.headerfile = strcat(strcat(cfg.dataset(1:end-3),'mat'));
% cfgEvent.channel    = cfg.trialdef.trigchan;
% cfgEvent.dataformat = 'spmeeg_mat';
% eventContStruct     = ft_preprocessing(cfgEvent);

% % Downsample from 6000Hz
% cfgDS                     = [];
% cfgDS.resamplefs          = cfg.trialdef.downsample;
% cfgDS.detrend             = 'no';
% eventContStruct         = ft_resampledata(cfgDS,eventContStruct);

eventCont           = cfg.rawData.trial{1}(pos_of_trig,:);

time                = cfg.rawData.time{1};
sample              = 1:length(time);
% convert to binary
eventDisc           = (eventCont > cfg.threshold);

% find first sample. first find differences
tmp                 = eventDisc - [0,eventDisc(1:end-1)];

% Trigger based on flank up or down
if strcmp(cfg.detectflank,'up')
    tmp2                = (tmp == 1);
elseif strcmp(cfg.detectflank,'down')
    tmp2                = (tmp == -1);
end
events              = sample(tmp2)';
events              = round(events);

if ~isempty(cfg.correct_time)
    disp(['Correcting timing by ' num2str(cfg.correct_time) 's']);
    time_to_correct = cfg.rawData.fsample*cfg.correct_time;    
    events = events+time_to_correct;
end


if isfield(cfg.trialdef,'downsample')
    events = round((events/(cfg.rawData.fsample/cfg.trialdef.downsample)));
end  

% create trl structure based upon the events
trl = [];

for j = 1:length(events)
    if isfield(cfg.trialdef,'downsample')
        trlbegin = (events(j) - (cfg.trialdef.prestim*cfg.trialdef.downsample));
        trlend   = (events(j) + (cfg.trialdef.poststim*cfg.trialdef.downsample));
        offset        = (-cfg.trialdef.prestim*cfg.trialdef.downsample);
        trl(end+1, :) = ([trlbegin trlend offset j]);
    else
        trlbegin = events(j) - cfg.trialdef.prestim*cfg.rawData.hdr.Fs;
        trlend   = events(j) + cfg.trialdef.poststim*cfg.rawData.hdr.Fs;
        offset        = -cfg.trialdef.prestim*cfg.rawData.hdr.Fs;
        trl(end+1, :) = ([trlbegin trlend offset j]);
    end
end
end