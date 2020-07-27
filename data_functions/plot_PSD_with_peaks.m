function plot_PSD_with_peaks(cfg,pow,freq)

if ~isfield(cfg, 'foi')
    cfg.foi = [1 100];
end

if ~isfield(cfg, 'ylim')
    cfg.ylim = [];
end

if ~isfield(cfg, 'transparency')
    cfg.transparency = 0.5;
end

if ~isfield(cfg, 'MinPeakProminence')
    cfg.MinPeakProminence = 30;
end

if ~isfield(cfg, 'MinPeakDistance')
    cfg.MinPeakDistance = 4;
end

if ~isfield(cfg, 'plot_kind')
    cfg.plot_kind = 'log';
end


freq_min = cfg.foi(1);
freq_max = cfg.foi(2);

pow = median(pow,3);
pow = mean(pow,2);

pos_min = interp1(freq,1:length(freq),freq_min,'nearest');
pos_max = interp1(freq,1:length(freq),freq_max,'nearest');

pow = pow(pos_min:pos_max,1);
freq = freq(1,pos_min:pos_max);

[p,l] = findpeaks(pow,'MinPeakDistance',4,'MinPeakProminence',...
    cfg.MinPeakProminence);


figure; 
if strcmp(cfg.plot_kind,'log')

    t = semilogy(freq,pow,'LineWidth',2);
else
    t = plot(freq,pow,'LineWidth',2);
end

t.Color(4)=cfg.transparency;
hold on;
plot(freq(l),p,'o')
xlabel('Frequency (Hz)','FontSize',14)
labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex','FontSize',14)
if ~isempty(cfg.ylim)
   ylim([cfg.ylim(1) cfg.ylim(2)]);
end
for f = 1:length(freq(l))
    freq_list = freq(l);
    h=text(freq_list(f)+0.5,p(f)+50,num2str(freq_list(f)),'FontSize',8);
    set(h,'Rotation',45);
end