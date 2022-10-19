data_dir        = 'D:\data\auditory_moving_ERF_BIDS'

%% Load the OPM data
% Read in the raw data using BIDS
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'aef';
cfg.bids.sub    = '001';
cfg.bids.ses    = '001';
cfg.bids.run    = '001';
rawData         = ft_opm_create(cfg);

%% Resample to 1000Hz
cfg                 = [];
cfg.resamplefs      = 1000;
[rawData]           = ft_resampledata(cfg, rawData);

%%
ft_channelselection('meg',rawData,'QZFM_Gen2')
ft_channelselection('megref',rawData,'QZFM_Gen2')

ft_channelselection('MEGTAN',rawData,'QZFM_Gen2')



senstype = ft_senstype(rawData);
regexp(datachannel,'.*TAN$')
datachannel = rawData.label;

labelmeg(~cellfun(@isempty, regexp(labelmeg, '.*TAN$')))

endsWith(datachannel,'TAN')


datachantype = ft_chantype(rawData.hdr);


senstype = ft_senstype(rawData);

rawData.grad.type = 'quspin_g2'

