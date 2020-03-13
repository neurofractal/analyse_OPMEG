% Read in the raw data using BIDS
cfg             = [];
cfg.folder      = 'D:\data\HMD_Tests';
cfg.bids.task   = 'noise';
cfg.bids.sub    = 'noise_14_02';
cfg.bids.ses    = '003';
cfg.bids.run    = '005';
cfg.precision   = 'single';
data         = ft_opm_create(cfg);


cfg = [];
cfg.channel = 'MEG';
%cfg.min_length = 0.05;
sat = ft_opm_saturations(cfg,data);




cfg                 = [];
cfg.channel         = 'MEG';
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [1 200];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,data);

[data2] = rm_sat_data(sat,data)


cfg                 = [];
cfg.channel         = 'MEG';
cfg.trial_length    = 3;
cfg.method          = 'tim';
cfg.foi             = [1 200];
cfg.plot            = 'yes';
pow                 = ft_opm_psd(cfg,data2);
























