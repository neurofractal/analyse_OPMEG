function [events] = ft_opm_make_events(data, trigger_chan,name)
%__________________________________________________________________________
% Make a table of events from an FIL OPM dataset.
%
% Inputs:
%   - data:         raw data loaded into Fieldtrip
%   - trigger_chan: full name of the trigger channels (e.g.
%   {'NI-TRIG-1','NI-TRIG-2'}
%   - name:         name/description corresponding to the trigger (e.g.
%   {'VR_onset','walking_onset'}
%
% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%__________________________________________________________________________

% Check if length of trigger_chan and name is equal
if length(trigger_chan) ~= length(name)
    warning(['Please make sure the length of trigger_name and name are the same']);
end

% Create empty arrays to hold the trigger info
onset       = [];
duration    = [];
sample      = [];
trial_type  = [];

% For each trigger channel...
for i = 1:length(trigger_chan)
    try
        cfg = [];
        cfg.rawData                 = data;
        cfg.trialdef.trigchan       = trigger_chan{i};
        %cfg.trialdef.downsample     = 500;
        cfg.correct_time            = [];
        cfg.trialdef.prestim        = 0;        % pre-stimulus interval
        cfg.trialdef.poststim       = 0;        % post-stimulus interval
        cfg.trialfun                = 'OPM_trialfun_Filbury';
        banana                      = ft_definetrial(cfg);
    catch
        ft_warning(['ft_definetrial failed for ' trigger_chan{i}]);
    end

    if ~isempty(banana.trl)
        % Onset
        onset = vertcat(onset,data.time{1}(banana.trl(:,1))');

        % Duration
        dur1 = data.time{1}(banana.trl(:,2))'-data.time{1}(banana.trl(:,1))';
        duration = vertcat(duration,dur1);
        clear dur1

        % Sample
        sample = vertcat(sample,banana.trl(:,1));

        % Name
        tt    = cell(size(banana.trl,1),1);
        tt(:) = {name{i}};
        trial_type = vertcat(trial_type,tt);
        clear tt

    else
        ft_warning(['Cannot find any valid triggers for ' trigger_chan{i}]);
    end
end

% Make table
t = table(onset,duration,sample,trial_type);

% Sort by onset
events = sortrows(t,'onset');

end



