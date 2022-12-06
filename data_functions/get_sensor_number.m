function [num_sens, chan_list] = get_sensor_number(rawData)
% Function to count the number of unique sensors in a Fieldtrip
% FIL-OPM data structure

labels = [];
count = 1;

for chan = 1:length(rawData.label)
    try
        label_regx = regexp(rawData.label{chan},...
            'G*-.*.-','match');
        labels{count} = label_regx{1};
        count = count+1;
    catch
        warning(['Unusual channel found: ' rawData.label{chan}]);
    end
end

chan_list = unique(labels);
num_sens =  length(chan_list);


