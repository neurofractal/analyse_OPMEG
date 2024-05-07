function [freq_bc] = ft_bc_freq(freq1,freq2,method)
%__________________________________________________________________________
% Baseline correct a Fieldtrip TFR structure (freq1) using a second
% structure of a different time (freq2)
% 
% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

if nargin == 2
    method = 'dB';
end

pow_baseline = freq2.powspctrm;

meanVals = nan+zeros(size(freq1.powspctrm));
for k = 1:size(pow_baseline,2)
    meanVals(:,k,:) = repmat(nanmean(pow_baseline(:,k,:), 3), [1 1 size(freq1.powspctrm, 3)]);
end

freq_bc = freq1;

switch method
    case 'db'
        freq_bc.powspctrm = 10*log10(freq_bc.powspctrm ./ meanVals);

    case 'perc_change'
        freq_bc.powspctrm = 100.*((freq_bc.powspctrm-meanVals)./ meanVals);
end