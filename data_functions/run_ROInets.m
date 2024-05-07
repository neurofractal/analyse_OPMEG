function [mat,t] = run_ROInets(out,Fs,BPF)

% Bp-filter
[out_bpf] = ft_preproc_bandpassfilter(out, Fs, BPF, 4);

figure; plot(out_bpf(1,:))

% Orthoganalise data to remove source leakage
% fprintf(' Orthoganalising the data using "symmetric" algorithm". May take a while\n');
data_orthog = (ROInets.remove_source_leakage(out_bpf, 'closest'));
%data_orthog = out_bpf;

figure; plot(data_orthog(1,:))

% Create Figure to show the effect of orthoganalisation
figure; subplot(2,1,1);
plot(out_bpf(1,1:1000),'LineWidth',3); hold on;
plot(data_orthog(1,1:1000),'LineWidth',3); hold on;
set(gca,'FontSize',20);
subplot(2,1,2);
plot(out_bpf(2,1:1000),'LineWidth',3); hold on;
plot(data_orthog(2,1:1000),'LineWidth',3);
set(gca,'FontSize',20);
title('closest')
drawnow;

% Take Hilbert transform
[dat] = ft_preproc_hilbert(data_orthog, 'abs');

% figure;  plot(t,dat(1,:),'LineWidth',2); hold on;
%          plot(t,data_orthog(1,:),'LineWidth',0.3);

% % Make dummy time value
t = (1:size(dat,2))/Fs;
% 
% Envelope the data using a 2 second window
[envelopedData,newFs] = downsample_envelope(dat,2,t);

disp(['New sampling frequency:' num2str(newFs) 'Hz']);
% 
%figure; plot(envelopedData); hold on;

% Run Correlations
Settings.Regularize.do            = false;

mat = ROInets.run_correlation_analysis(data_orthog,     ...
    envelopedData, Settings.Regularize);

end
