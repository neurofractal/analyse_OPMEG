function trl = OPM_TrialFun(cfg)

% Read header
header              = ft_read_header(strcat(cfg.dataset(1:end-3),'mat'));

% Read continuous event channel.
cfgEvent            = [];
cfgEvent.continuous = 'yes';
cfgEvent.datafile   = cfg.dataset;
cfgEvent.headerfile = strcat(strcat(cfg.dataset(1:end-3),'mat'));
cfgEvent.channel    = cfg.trialdef.trigchan;
cfgEvent.dataformat = 'spmeeg_mat';
eventContStruct     = ft_preprocessing(cfgEvent);

% Downsample from 6000Hz
cfgDS                     = [];
cfgDS.resamplefs          = cfg.trialdef.downsample;
cfgDS.detrend             = 'no';
eventContStruct         = ft_resampledata(cfgDS,eventContStruct);



eventCont           = eventContStruct.trial{1}(1,:);

time                = eventContStruct.time{1};
sample              = 1:length(time);
% convert to binary
eventDisc           = (eventCont > 1000000);

% find first sample. first find differences
tmp                 = eventDisc - [0,eventDisc(1:end-1)];
tmp2                = (tmp == 1);
events              = sample(tmp2)';
begsample           = events - cfg.trialdef.prestim/(1/cfg.trialdef.downsample);
endsample           = events + cfg.trialdef.poststim/(1/cfg.trialdef.downsample) - 1;

trl                 = [begsample endsample ones(size(begsample)).*(-cfg.trialdef.prestim/(1/cfg.trialdef.downsample))];