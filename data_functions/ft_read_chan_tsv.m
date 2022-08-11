function channels = ft_read_chan_tsv(filename)

% Script for importing data from the following text file:
%
%    filename: D:\MATLAB\Analysis\Data\testData\channels.tsv
%
% Auto-generated by MATLAB on 28-Jan-2020 10:07:39

%% Setup the Import Options
opts        = detectImportOptions(filename,'filetype','text');
tsv         = readtable(filename, opts); 
tsv         = table2struct(tsv);

%% Reformat into FT-like structure with cell array
channels = [];

fn = fieldnames(tsv);
for i = 1:length(fn)
    if ischar(tsv(1).(fn{i}))
        channels.(fn{i}) = {tsv(:).(fn{i})}';
    end
end



%% Try to determine the orientation of the sensors
for i = 1:length(channels.name)
    try
        chan_end  = channels.name{i}(end-2:end);
        chan_end2 = channels.name{i}(end);
        if strcmp(chan_end,'TAN')
            channels.fieldori{i,1} = chan_end;
        elseif strcmp(chan_end,'RAD')
            channels.fieldori{i,1} = chan_end;
        elseif strcmp(chan_end2,'X')
            channels.fieldori{i,1} = chan_end2;
        elseif strcmp(chan_end2,'Y')
            channels.fieldori{i,1} = chan_end2;
        elseif strcmp(chan_end2,'Z')
            channels.fieldori{i,1} = chan_end2;
        else
            channels.fieldori{i,1} = 'UNKNOWN';
        end
    catch
        channels.fieldori{i,1} = 'UNKNOWN';
    end
end

%% Replace G2 prefix
% For older data from UCL OPM lab (where chan names ended -RAD or -TAN),
% analyse_OPMEG removed the G2 prefix. Don't do this for newer data, but
% do for older data:
if any(contains(channels.name,"RAD"))
    % Remove the G2 part of the channel name
    try
        % There is probably a more efficient way to do this...
        indx = find(contains(channels.name,'G2'));

        for i = 1:length(indx)
            to_replace = char(channels.name(indx(i)));
            to_replace = to_replace(4:end);
            channels.name(indx(i)) = {to_replace};
        end
    catch
    end
end

% For backwards compatibility, when channel labels contained
% just 'MEG'. Can be removed in future:
% Replace 'MEG' with 'megmag'
try
    indx = find(ismember(channels.type,'MEG'));
    channels.type(indx) = {'megmag'};
catch
    ft_warning('Cannot replace MEGMAG with megmag');
end

% Replace 'MEMAG' with 'megmag' (typo in some recordings.. can be removed
% later)
try
    indx = find(ismember(channels.type,'MEMAG'));
    channels.type(indx) = {'megmag'};
catch
    ft_warning('Cannot replace MEGMAG with megmag');
end

% Replace 'MEGMAG' with 'megmag'
try
    indx = find(ismember(channels.type,'MEGMAG'));
    channels.type(indx) = {'megmag'};
catch
    ft_warning('Cannot replace MEGMAG with megmag');
end

% Replace 'MEGREFMAG' with 'refmag'
try
    indx = find(ismember(channels.type,'MEGREFMAG'));
    channels.type(indx) = {'refmag'};
catch
    ft_warning('Cannot replace REF with refmag');
end

% Replace 'TRIG' with 'trigger
try
    indx = find(contains(channels.type,'TRIG'));
    channels.type(indx) = {'trigger'};
catch
    ft_warning('Cannot replace TRIG with trigger');
end

%%%%%%%%%
% TO FIX
%%%%%%%%%
% If the channel name contains flux, change the channel type
try
    indx = contains(channels.name,'Flux');
    channels.type(indx) = {'flux'};
catch
end



