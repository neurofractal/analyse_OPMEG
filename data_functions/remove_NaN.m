function [data_out] = remove_NaN(data)

 trial2keep_r   = [];
    trial2reject   = [];
    count          = 1;
    count2         = 1;

    for t = 1:length(data.trial)
        result = sum(isnan(data.trial{t}(:)));
        if ~result
            trial2keep_r(count) = t;
            count=count+1;
        else
            trial2reject(count2) = t;
            count2       = count2+1;
        end
    end

    disp(['Keeping ' num2str(length(trial2keep_r)) ' trials' ...
        '. Rejecting ' num2str(length(trial2reject)) ' trials']);

    % Select the data
    cfg            = [];
    cfg.trials     = trial2keep_r;
    data_out   = ft_selectdata(cfg,data);
end