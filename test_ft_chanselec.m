


[channel] = ft_channelselection_opm('MEG', rawData)

cfg = [];
cfg.channel = ft_channelselection_opm('MEG', rawData);
MEG = ft_selectdata(cfg,rawData);

cfg = [];
cfg.channel = ft_channelselection_opm('RAD', MEG);
MEG2 = ft_selectdata(cfg,MEG);

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
