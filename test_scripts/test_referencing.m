cd('D:\data\benchmarking_26_02\sub-002\ses-001\meg');

cfg = [];
cfg.data = 'sub-002_ses-001_task-faces_run-001_meg.bin';
cfg.precision  = 'single';
[data] = ft_opm_create(cfg)

cfg = [];
cfg.resamplefs = 1000;
[data] = ft_resampledata(cfg, data);


cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [98 102];
cfg.bsfiltord = 5;
data        = ft_preprocessing(cfg,data);

% cfg = [];
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 3;
% %cfg.hpfiltord = 5;
% data        = ft_preprocessing(cfg,data);

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 100;
%cfg.lpfiltord = 5;
data        = ft_preprocessing(cfg,data);


cfg                 = [];
cfg.channel             = vertcat(ft_channelselection_opm('TAN',data),...
     '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [1 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,data);
ylim([1 1e3]);


cfg = [];
cfg.channel = vertcat(ft_channelselection_opm('TAN',data),...
     '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
cfg.refchannel = 'MEGREF';
cfg.filter_ref = [3 20; 35 40; 40 60; 80 100];
cfg.derivative = 'yes';
[data_out] = ft_opm_synth_gradiometer(cfg,data);

cfg                 = [];
cfg.channel         = 'all'
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [2 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,data_out);
ylim([1 1e3]);

ft_databrowser([],data_out);


%% Test Zapline
addpath('D:\scripts\NoiseTools');


cd('D:\data\benchmarking_26_02\sub-002\ses-001\meg');

cfg = [];
cfg.data = 'sub-002_ses-001_task-faces_run-001_meg.bin';
cfg.precision  = 'single';
[data] = ft_opm_create(cfg);



cfg                 = [];
cfg.channel             = vertcat(ft_channelselection_opm('TAN',data),...
     '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
data = ft_selectdata(cfg,data);

cfg = [];
cfg.resamplefs = 1000;
[data] = ft_resampledata(cfg, data);

x = data.trial{1}';
sr = data.fsample;

x=nt_demean(x);
x=nt_detrend(x,1);


% parameters
FLINE=50/sr; % line frequency
NREMOVE=1; % number of components to remove

tic; ttt = nt_zapline(x,FLINE,NREMOVE); toc;


data_out = data;
data_out.trial{1} = x';

ft_databrowser([],data)

cfg                 = [];
cfg.channel         = 'all'
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [1 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,data_out);
ylim([1 1e3]);



