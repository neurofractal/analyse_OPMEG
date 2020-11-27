function [data, M, chan_inds] = ft_denoise_mfc(data)

% Regression of ambient interference through the modelling of a mean 
% magnetic field. Tierney et al, https://doi.org/10.1101/2020.11.25.397778 

% This is primarily designed for use on raw magnetometer data which hasn't
% been preprocessed in any way. 

% NOTE: This is developmental code, no results are guarenteed!
% NOTE: This codebase is not offically supported by WCHN
% George O'Neill (g.o'neill [at] ucl.ac.uk)

% generate projector
N = data.grad.chanori;
M = eye(length(N)) - N*pinv(N);

% Update the tra to account for this (its probably identity already,
% but no harm in a belt and braces approach). Essential to correct the lead
% fields going forward.
data.grad.tra = M*data.grad.tra;

% Identify which channels need projecting through M
coil_labels = data.grad.label;
chan_labels = data.label;

chan_inds = zeros(1,numel(coil_labels));
for ii = 1:numel(coil_labels)
    try
        chan_inds(ii) = match_str(chan_labels,coil_labels{ii});
    catch
        chan_inds(ii) = NaN;
    end
end

denoised_data = cell(data.hdr.nTrials,1);

for ii = 1:data.hdr.nTrials
    
    trial = data.trial{ii};
    Y = trial(chan_inds,:);
    
    Yc = M*Y;
    trial(chan_inds,:) = Yc;
    denoised_data{ii} = trial;
    
end

data.trial = denoised_data;
