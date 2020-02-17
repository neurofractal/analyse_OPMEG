function [data2] = rm_sat_data(sat,data)

[tf] = ~ismember(data.time{1},sat.alltime);
data2 = data;
data2.trial{1} = data2.trial{1}(:,tf)
data2.time{1} = (0:1:(size(data2.trial{1},2)-1))./data.fsample;

end