function path = go_bids2data(cfg)

cfg.bids.modality = ft_getopt(cfg.bids,'modality','meg');
cfg.bids.extention = ft_getopt(cfg.bids,'extention','bin');

try
    file_name_bids = ['sub-' cfg.bids.sub '_ses-' cfg.bids.ses ...
        '_task-' cfg.bids.task '_run-' cfg.bids.run '_' cfg.bids.modality '.' cfg.bids.extention];
catch
    error('Did you specify all the required cfg.bids information?')
end
path = fullfile(cfg.folder,['sub-' cfg.bids.sub],...
    ['ses-' cfg.bids.ses],cfg.bids.modality,file_name_bids);