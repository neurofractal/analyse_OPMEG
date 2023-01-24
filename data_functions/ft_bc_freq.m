function [freq_bc] = ft_bc_freq(freq1,freq2)
%__________________________________________________________________________
% Baseline correct a Fieldtrip TFR structure (freq1) using a second
% structure of a different time (freq2)
% 
% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

pow_baseline = freq2.powspctrm;

meanVals = nan+zeros(size(freq1.powspctrm));
for k = 1:size(pow_baseline,2)
    meanVals(:,k,:) = repmat(nanmean(pow_baseline(:,k,:), 3), [1 1 size(freq1.powspctrm, 3)]);
end

freq_bc = freq1;
freq_bc.powspctrm = 10*log10(freq_bc.powspctrm ./ meanVals);
