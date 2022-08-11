function pos = ft_read_pos_tsv(filename)

% Script for importing data from the following text file:
%
%    filename: D:\MATLAB\Analysis\Data\testData\channels.tsv
%
% Auto-generated by MATLAB on 28-Jan-2020 10:07:39

%% Setup the Import Options

try
    opts                    = delimitedTextImportOptions('NumVariables', 4);
    
    % Specify range and delimiter
    opts.DataLines          = [2, Inf];
    opts.Delimiter          = '\t';
    
    % Specify column names and types
    opts.VariableNames      = ['label', 'Px', 'Py', 'Pz', 'Ox', 'Oy', 'Oz'];
    opts.VariableTypes      = ['string', 'double', 'double', 'double', 'double', 'double', 'double'];
    opts                    = setvaropts(opts, 1, 'WhitespaceRule', 'preserve');
    opts                    = setvaropts(opts, [1, 2, 3], 'EmptyFieldRule', 'auto');
    opts.ExtraColumnsRule   = 'ignore';
    opts.EmptyLineRule      = 'read';
    
    % Import the data
    posAll                  = table2struct(readtable(filename, opts),'ToScalar',true);
    
    % Where an earlier version of MATLAB is used - try the tsvread in /external
catch
    posAll = tsvread(filename);
    posAll.label = posAll.name;
end

pos.label               = posAll.label;
pos.chanpos             = [posAll.Px posAll.Py posAll.Pz];
pos.chanori             = [posAll.Ox posAll.Oy posAll.Oz];

pos.label = cellstr(pos.label);

% For older data from UCL OPM lab (where chan names ended -RAD or -TAN),
% analyse_OPMEG removed the G2 prefix. Don't do this for newer data, but
% do for older data:
if any(contains(pos.label,"RAD"))
    % Remove the G2 part of the channel name
    try
        % There is probably a more efficient way to do this...
        indx = find(contains(pos.label,'G2'));

        for i = 1:length(indx)
            to_replace = char(pos.label(indx(i)));
            to_replace = to_replace(4:end);
            pos.label(indx(i)) = {to_replace};
        end
    catch
    end
end

end


