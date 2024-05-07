function [data_out] = reset_time(data)
Fs = data.fsample;

% Reset the time index
for trl = 1:length(data.trial)
    data.time{trl} = ([1:1:length(data.trial{trl})])./Fs;
end

data_out = data;
end