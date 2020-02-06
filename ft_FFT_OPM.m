function ft_FFT_OPM(cfg,rawData)


if ~isfield(cfg, 'remove_trig')
    cfg.remove_trig = 'no';
end

if ~isfield(cfg, 'foi')
    cfg.foi = [0:0.01:200];
end

disp(['Calculating PSD from ' num2str(cfg.foi(1)) 'Hz to ' ...
    num2str(cfg.foi(end)) 'Hz']);

if strcmp(cfg.remove_trig,'yes')
    try
        indx = ~contains(rawData.hdr.chantype,'trigger');
        
        cfg = [];
        cfg.channel = rawData.label(indx);
        rawData = ft_selectdata(cfg,rawData);
        
    catch
        ft_warning('Could not remove trigger channel');
    end
end

foi = cfg.foi;

cfg2 = [];
cfg2.output      = 'pow';          % Return PSD
cfg2.channel     = 'all';
cfg2.method      = 'mtmfft';
%cfg.pad         = 'nextpow2';
cfg2.taper       = 'boxcar';      % Hann window as taper
cfg2.foi         = foi;         % Frequency range

psd_hann = ft_freqanalysis(cfg2, rawData);

figure; 
set(gcf,'Position',[100 100 1200 800]);
plot(psd_hann.freq,log(psd_hann.powspctrm));
set(gca,'FontSize',14);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (au^2/Hz)')
hold on
plot(psd_hann.freq,mean(log(psd_hann.powspctrm)),'-k','LineWidth',2);
lgd = legend(vertcat(rawData.label, 'mean'));
set(lgd,'Location','BestOutside');
end




