function [rawData_out] = ft_opm_reorder(rawData)
% Function to reorder the channels of data based on their spatial layout
% (front-to-back).
%
% EXAMPLE USEAGE:   rawData_out = ft_opm_reorder(cfg,rawData)
% ...where, cfg is the input structure and rawData is the raw OPM data
% loaded using ft_opm_create
%
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk)
%__________________________________________________________________________

%% Function housekeeping
if ~isfield(rawData, 'grad')
    error('ERROR: No grad structure present');
end

%% Start of function
label_data = rawData.label;
label_grad = rawData.grad.label;

% Sort from front to back
[a b] = sort(rawData.grad.chanpos(:,2),'descend');

% Sort the label_grad based on this
label_grad = label_grad(b,:);

% Find indices of label_grad for label_data
[sel1, sel2] = match_str(label_grad,label_data);

% Find indices of channels not in the grad
sel3 = find(~ismember(1:length(rawData.label),sel2));

comb = vertcat(sel2,sel3');

% Reorder the data
rawData_out = rawData;
rawData_out.label = rawData.label(comb);
rawData_out.trial{1} = rawData_out.trial{1}(comb,:);

if length(rawData_out.hdr.label) == length(label_data)
    
    rawData_out.hdr.label = rawData_out.hdr.label(comb);
    rawData_out.hdr.chantype = rawData_out.hdr.chantype(comb);
    rawData_out.hdr.chanunit = rawData_out.hdr.chanunit(comb);
    %rawData_out.hdr.fieldori = rawData_out.hdr.fieldori(comb);
else
    warning('Removing hdr information');
    rawData_out = rmfield(rawData_out,'hdr');
end
end

