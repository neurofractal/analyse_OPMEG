cd('D:\data\benchmarking_26_02\sub-002\ses-001\meg');

cfg = [];
cfg.data = 'sub-002_ses-001_task-faces_run-001_meg.bin';
cfg.precision  = 'single';
[data] = ft_opm_create(cfg)

% 
% cfg = [];
% cfg.bsfilter = 'yes';
% cfg.bsfreq = [98 102];
% cfg.bsfiltord = 5;
% data        = ft_preprocessing(cfg,data);

% cfg = [];
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 3;
% cfg.hpfiltord = 5;
% data        = ft_preprocessing(cfg,data);

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 100;
cfg.lpfiltord = 5;
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
cfg.filter_ref = [1 20; 35 40; 40 60; 70 90];
cfg.derivative = 'yes';
[data_out] = ft_opm_synth_gradiometer(cfg,data);

cfg                 = [];
cfg.channel         = 'all'
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [1 100];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,rawData_meg);
ylim([1 1e3]);

ft_databrowser([],data_out);



% Referencing
%Select Reference Channels
cfg                     = [];
cfg.channel          = {'N0-RAD','N4-RAD','N3-RAD','MV-RAD',...
    'N0-TAN','N4-TAN','N3-TAN','MV-TAN'};
refs                    = ft_selectdata(cfg,data);




cfg = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 3;
cfg.hpfiltord = 5;
refs1            = ft_preprocessing(cfg,refs);

cfg = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 20;
cfg.lpfiltord = 5;
refs1            = ft_preprocessing(cfg,refs1);





cfg = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 40;
cfg.hpfiltord = 5;
refs2            = ft_preprocessing(cfg,refs);

cfg = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 60;
cfg.lpfiltord = 5;
refs2            = ft_preprocessing(cfg,refs2);


cfg = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 35;
cfg.hpfiltord = 5;
refs3            = ft_preprocessing(cfg,refs);

cfg = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 40;
cfg.lpfiltord = 5;
refs3            = ft_preprocessing(cfg,refs3);



cfg = [];
cfg.channel = 

cfg = [];
cfg


ddd = ft_appenddata([],refs1,refs2,refs3);

ref_to_use = [];
ref_to_use.fsample = ddd.fsample;
ref_to_use.trial{1} = vertcat(ddd.trial{1,1},ddd.trial{1,2},ddd.trial{1,3});
ref_to_use.time{1} = vertcat(ddd.time{1,1},ddd.time{1,2},ddd.time{1,3});
ref_to_use.label = cellstr(string(1:24))';




% Denoise
cfg                     = [];
cfg.refchannel          = ref_to_use.label;
cfg.channel             = vertcat(ft_channelselection_opm('TAN',data),...
    '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');
cfg.truncate            = 'no';
cfg.zscore              = 'yes';
cfg.updatesens          = 'no';
rawData_meg             = ft_denoise_pca(cfg,data,ref_to_use);

cfg                     = [];
cfg.refchannel          = {'N0-RAD','N4-RAD','N3-RAD','MV-RAD',...
    'N0-TAN','N4-TAN','N3-TAN','MV-TAN'};
cfg.channel             = vertcat(ft_channelselection_opm('TAN',rawData),...
    '-N0-TAN','-N4-TAN','-N3-TAN','-MV-TAN');

rawData_meg        = ft_denoise_tsr(cfg,rawData,refs);



























