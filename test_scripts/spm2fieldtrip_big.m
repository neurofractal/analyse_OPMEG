%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

addpath('D:\scripts\spm12');
addpath('D:\Github\OPM');

% Set up 
spm('defaults','eeg')
spm('createintwin')

%%
cd('D:\data\tutorial_OPM_data\sub-001\ses-002\meg');

S =[];
S.data = 'sub-001_ses-002_task-aef_run-001_meg.bin';
S.channels='sub-001_ses-002_task-aef_run-001_channels.tsv';
S.meg='sub-001_ses-002_task-aef_run-001_meg.json';
S.positions='sub-001_ses-002_task-aef_run-001_positions.tsv';
S.precision='single';
S.sMRI = 'D:\data\tutorial_OPM_data\sub-001\ses-002\anat\001.nii';
D = spm_opm_create(S);

%%
S=[];
S.D=D;
S.balance = 1;
[hfD, Yinds]=spm_opm_hfc(S);

%%
ftraw = hfD.ftraw;
ft_channelselection({'all','-XXX'},ftraw.label)

grad = hfD.sensors('MEG');


cfg = [];
cfg.channel = {'all','-G2-DQ-TAN'}
ftraw2 = ft_selectdata(cfg,ftraw);


mfD.inv{1}





