function ft_FFT_OPM(cfg,rawData)

if ~isfield(cfg, 'voltype')
    cfg.voltype = 'Single Shell';
end




cfg = [];

cfg = [];
cfg.output      = 'pow';          % Return PSD
cfg.channel     = 'all'  % Calculate for MEG and EEG
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';      % Hann window as taper
cfg.foilim      = [0.5 500];         % Frequency range

psd_hann = ft_freqanalysis(cfg, rawData);



figure; plot(psd_hann.freq,log(psd_hann.powspctrm));




cfg = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306mag'; % Layout for MEG magnetometers
cfg.showlabels      = 'yes';
cfg.xlim            = [3 35];           % Frequencies to plot

figure;
ft_singleplotER(cfg, psd_hann);


