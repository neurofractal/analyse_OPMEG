function plot_mean_multiple_PSD(pow,freq,cols)

figure()
fig= gcf;
fig.Color=[1,1,1];
set(fig, 'DefaultFigureRenderer', 'painters');
set(gcf,'Position',[200 200 1200 500]);
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2)+0.01 pos(3) pos(4)-0.01]);

lims = [0.01 5; 5 60];

for d = 1:2
    subplot(1,2,d);
    for i = 1:length(pow)
        po = mean(pow{i}(:,:,:),3);
        mean_po = squeeze(mean(po',1));
            h = plot((freq),mean_po','color',cols{i},'LineWidth',3);
        h.Color(4) = 0.8;
        hold on;
        xp2 =-3:round((freq(end)));
        yp2=ones(1,length(xp2))*15;
        p2 =plot(xp2,yp2,'--k');
        p2.LineWidth=2;
        %     area(log2(freq),mean_po','FaceColor',cols{i},'FaceAlpha',0.7,...
        %         'EdgeAlpha',0);
        %grid on
        ax = gca; % current axes
        ax.FontSize = 20;
        ax.TickLength = [0.02 0.02];
        set(gca, 'YScale', 'log'); hold on;
        %set(gca, 'XScale', 'log2'); hold on;
        xlabel('Frequency (Hz)','FontSize',25)
        if d==1
            labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
            ylabel(labY,'interpreter','latex','FontSize',25)
            set(gca,'YLim',[7 2e6],'YTick',[10,100,1000,10000,100000,1e6]);
        else
            set(gca,'YLim',[7 3e3],'YTick',[10,100,1000]);
        end

    end
    xlim([lims(d,:)])
end



