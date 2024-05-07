function freq_out = median_BC(freq_end)
% Some jiggery pokery to get median across trials
%rrr2 = squeeze(trimmean(freq_end.powspctrm,50,1));

rrr2 = squeeze(nanmedian(freq_end.powspctrm,1));

freq_end.powspctrm = rrr2;
freq_end.dimord = 'chan_freq_time';
freq_out = freq_end;
end