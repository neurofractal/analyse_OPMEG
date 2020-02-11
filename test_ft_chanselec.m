


[channel] = ft_channelselection_opm('MEG', rawData, 'quspin_g2')

[channel2] = ft_channelselection_opm('RAD', rawData, 'quspin_g2')


cfg = [];
cfg.channel = {'all','-NI-TRIG'};
fff = ft_selectdata(cfg,rawData);

cfg = [];
cfg.channel = ft_channelselection('megref', rawData, 'quspin_g2');
ref = ft_selectdata(cfg,fff);

cfg = [];
cfg.channel = ft_channelselection('megref', rawData, 'quspin_g2');
rad = ft_selectdata(cfg,rawData);


ft_channelselection('megref',fff,'quspin_g2')


  senstype = ft_senstype(rawData.label);
