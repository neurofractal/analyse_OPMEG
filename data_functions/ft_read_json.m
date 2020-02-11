function meg    = ft_read_json(filename)
%% Setup the Import Options
opts                    = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines          = [1, Inf];
opts.Delimiter          = "";
opts.VariableTypes      = "string";
opts                    = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts                    = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule   = "ignore";
opts.EmptyLineRule      = "read";

% Import the data
meg         = readtable(filename, opts);

%% Convert to output type
meg         = table2array(meg);
meg         = jsondecode(meg);