function [data] = do_SPM_referencing(folder,filename)

cd(folder);

restoredefaultpath
addpath('D:\scripts\spm12');
spm('defaults', 'eeg')
addpath('D:\Github\OPM');

S =[];
S.data = [filename '_meg.bin'];
S.channels = [filename '_channels.tsv'];
S.meg= [filename '_meg.json'];
S.precision = 'single';
D = spm_opm_create(S);

S=[];
S.D=D;
S.triallength=3000;
S.bc=1;
S.channels='MEG';
S.plot=1;
spm_opm_psd(S);
xlim([1, 100])
ylim([1, 1000])

S=[]
S.D=D;
S.confounds={'REF'};
D = spm_opm_synth_gradiometer(S);

D.save()

data_file = fullfile(folder,['d' filename '_meg.mat']);

D           = spm_eeg_load(data_file);

restoredefaultpath

addpath('D:\scripts\fieldtrip-master');
ft_defaults

cd('D:\scripts\fieldtrip-master\external\spm8');

[data123]   = spm2fieldtrip(D);

cd(folder);
end








