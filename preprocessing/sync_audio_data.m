%% Paths (RS)
fieldtripDir    = 'C:\Users\vr-admin\Documents\fieldtrip-20210418';
script_dir      = 'C:\Users\vr-admin\Documents\analyse_OPMEG';
data_dir        = 'C:\Users\vr-admin\Documents\microphone_test';
save_dir        = 'C:\Users\vr-admin\Documents';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));

cd(save_dir)

%% (2) Start preprocessing.
% Read in the raw data using BIDS
disp('Loading Data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'test';
cfg.bids.sub    = '001';
cfg.bids.ses    = '001';
cfg.bids.run    = '002';
rawData         = ft_opm_create(cfg);

%% Load Audio Data
[y,Fs] = audioread('test_data2.wav');
t = ([0:1:length(y)]/Fs)';
t(end) = [];

% Filter the data to capture 441Hz tone
[filt] = ft_preproc_bandpassfilter(y', Fs, [435 445], 2);

% Find first value that exceeds 0.5
find_filt = find(filt'>0.5); %Change the value 0.5 if too high or low
aud_2s = y(find_filt(1)-Fs:find_filt(1)+Fs);
% Look for onset of the tone in 2s of original data
find_aud_2s = find(aud_2s>0.5)
trig = find(y==aud_2s(find_aud_2s(1)));

% Clip the audio based on this trigger
y_clipped = y(trig:end,:);
t = ([0:1:length(y_clipped)]/Fs)';
t(end) = [];

% Play the first 2s of the clipped audio
sound(y_clipped(1:Fs*2,:),Fs);

% Plot the clipped audio for sanity
figure; plot(t,y_clipped);

%% Get audio trigger from OPMs
trig_chan = find(contains(rawData.label,'FluxZ-B'));

tChan=rawData.trial{1}(trig_chan,:);               % Get ith trigger
thresh=std(tChan)*2;                               % Threshold it

evSamples=find(diff(tChan>thresh)==1)+1;           % Find 1st sample exceeding threshold

% Get audio onset time and end time
opm_audio_onset = rawData.time{1}(evSamples);
opm_audio_end   = rawData.time{1}(end);

%% Pad the audio to match OPM data
pad_data = zeros([round(Fs*opm_audio_onset,0) 2]);

y_clipped2 = vertcat(pad_data,y_clipped);

t = ([0:1:length(y_clipped2)]/Fs)';
t(end) = [];

% Plot Figure to check they are lined up correctly
tChan_normalised = 2*((tChan-min(tChan))/(max(tChan)-min(tChan)))-1;
figure;plot(t,y_clipped2); hold on;
plot(rawData.time{1},tChan_normalised);

%% Trim audio to end of OPM recording 
if t(end) < opm_audio_end
    disp('Audio data shorter than the OPM data.. NOT SUPPORTED YET');
    % Here you would need to pad the end of the audio recording
else
    ind = interp1(t,1:length(t),opm_audio_end,'nearest');
    y_clipped3 = y_clipped2(1:ind,:);
    t = ([0:1:length(y_clipped3)]/Fs)';
    t(end) = [];
end

%% Output audio
audiowrite('audio_out.wav',y_clipped3,Fs)

% Plot Figure to check they are lined up correctly
tChan_normalised = 2*((tChan-min(tChan))/(max(tChan)-min(tChan)))-1;
figure;plot(t,y_clipped3); hold on;
plot(rawData.time{1},tChan_normalised);

%% Trigger audio based on NI-TRIG
trig_chan = find(contains(rawData.label,'NI-TRIG'));

tChan=rawData.trial{1}(trig_chan,:);               % Get ith trigger

% Upsample
tinit = ([0:1:length(tChan)]/rawData.fsample);
tinit(end) = [];
tfinal = ([0:1:length(tChan)]/Fs);
tfinal(end) = [];
tChan_upsamp = interp1(tinit, tChan, t, 'linear', 'extrap');

% Find threshold
thresh=std(tChan_upsamp)*2;                               % Threshold it
evSamples=find(diff(tChan_upsamp>thresh)==1)+1;           % Find 1st sample exceeding threshold

% Epoch and export audio data for 5s after each trigger
num_sec = 5;

for trial = 1:length(evSamples)
    
    y_epoched = y_clipped3(evSamples(trial):(evSamples(trial)+(Fs*num_sec)));
    
    audiowrite(['audio_trial_' num2str(trial) '.wav'],y_epoched,Fs)
    disp(['Written trial ' num2str(trial)]);
end





