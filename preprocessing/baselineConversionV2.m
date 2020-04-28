function baseConvData = baselineConversionV2(baselineFreqData,activityFreqData,cfg)
% This function allows for comparison between a baseline period not
% time-locked to the same event as an activity period. As far as I can
% tell, the FieldTrip toolbox does not have this functionality.
% 
% Both baseline conversion at the invididual level and trial level are
% calculated.
% 
% The following cfg settings should be used.
% 
%   cfg.baselinetype            = 'db','relative','relchange','normchange',
%                                   'absolute' or 'perchange')
%   cfg.baseline.toi            = [start end], eg. [-0.4 -0.2].
%   cfg.baseline.foi            = [start end], eg. [3 35].
% 
% NA 15/08/2019.

% Grab the cfg inputs
foi                     = cfg.baseline.foi;
toi                     = cfg.baseline.toi;
baselinetype            = cfg.baselinetype;

% Find the baseline data for the input activity. 
activityTrials          = activityFreqData.trialinfo(:,1);
baselineTrials          = baselineFreqData.trialinfo(:,1);

[~,baselineTrialIdx,~]  = intersect(baselineTrials,activityTrials);

clear activityTrials baselineTrials 

% Select the frequencies, time and trials of the baseline data. Average
% over time, because I don't need that for the baseline conv function. 
cfg                     = [];
cfg.trials              = baselineTrialIdx;
cfg.latency             = toi;
cfg.avgovertime         = 'yes';
cfg.frequency           = foi;
cfg.nanmean             = 'yes';
baselineFreqData        = ft_selectdata(cfg,baselineFreqData);

clear baselineTrialIndex

% Select the requested frequencies of the activity data. 
cfg                     = [];
cfg.frequency           = foi;
activityFreqData        = ft_selectdata(cfg,activityFreqData);

clear cfgSelectData cfg

% Find the average power within individual.
cfg                     = [];
cfg.avgoverrpt          = 'yes';
cfg.nanmean             = 'yes';
baselineAvgFreqData     = ft_selectdata(cfg,baselineFreqData);
activityAvgFreqData     = ft_selectdata(cfg,activityFreqData);

% Run baseline conversion of the activity data against the baseline. The
% dimensions of the baselineConversion matrix needs to be trial x channel x
% freq x time.
baseConvTrl         = nan(size(activityFreqData.powspctrm));
baseConvInd    = nan(size(activityFreqData.powspctrm));
for trial = 1:length(activityFreqData.trialinfo(:,1))
    for channel = 1:length(activityFreqData.label)
        for freq = 1:length(activityFreqData.freq)
            % Select the baseline power.
            baselineTrl     = squeeze(baselineFreqData.powspctrm(trial,channel,freq));
            baselineInd     = squeeze(baselineAvgFreqData.powspctrm(channel,freq));
            activityTrl     = squeeze(activityFreqData.powspctrm(trial,channel,freq,:));
            activityInd     = squeeze(activityAvgFreqData.powspctrm(channel,freq,:));
            switch baselinetype
                case 'db'
                    baseConvTrl(trial,channel,freq,:)...
                        = 10*log10(activityTrl ./ baselineTrl);

                    baseConvInd(trial,channel,freq,:)...
                        = 10*log10(activityInd ./ baselineInd);
                case 'relative'
                    baseConvTrl(trial,channel,freq,:)...
                        = activityTrl ./ baselineTrl;

                    baseConvInd(trial,channel,freq,:)...
                        = activityInd ./ baselineInd;
                case 'relchange'
                    baseConvTrl(trial,channel,freq,:)...
                        = (activityTrl - baselineTrl) ./ baselineTrl;

                    baseConvInd(trial,channel,freq,:)...
                        = (activityInd - baselineInd) ./ baselineInd;
                case 'normchange'
                    baseConvTrl(trial,channel,freq,:)...
                        = (activityTrl - baselineTrl) ./ (activityTrl + baselineTrl);

                    baseConvInd(trial,channel,freq,:)...
                        = (activityInd - baselineInd) ./ (activityInd + baselineInd);
                case 'absolute'
                    baseConvTrl(trial,channel,freq,:)...
                        = activityTrl - baselineTrl;

                    baseConvInd(trial,channel,freq,:)...
                        = activityInd - baselineInd;
                case 'perchange'
                    baseConvTrl(trial,channel,freq,:)...
                        = 100*(activityTrl - baselineTrl) ./ baselineTrl;

                    baseConvInd(trial,channel,freq,:)...
                        = 100*(activityInd - baselineInd) ./ baselineInd;
                case 'dbraw'
                    baseConvTrl(trial,channel,freq,:)...
                        = 10*log10(activityTrl ./ baselineTrl);
                    baseConvTrl(trial,channel,s,:)...
                        = 10^(baseConvTrl(trial,channel,freq,:) ./ 10);

                    baseConvInd(trial,channel,freq,:)...
                        = 10*log10(activityInd ./ baselineInd);
                    baseConvInd(trial,channel,s,:)...
                        = 10^(baseConvInd(trial,channel,freq,:) ./ 10);
                otherwise
                    disp('No baseline type selected.')     
            end
        end
    end
end

% Create the fieldtrip structure.
baseConvData                = activityFreqData;
baseConvData.baseConvTrl    = baseConvTrl;
baseConvData.baseConvInd    = baseConvInd;

% Average over repeats.
cfg                         = [];
cfg.avgoverrpt              = 'yes';
cfg.nanmean                 = 'yes';
baseConvData                = ft_selectdata(cfg,baseConvData);

% Plot to check
% cfg = [];
% 
% cfg.layout              = antLayoutFile;
% cfg.ylim                = [3 35];
% cfg.xlim                = [-3 7];
% cfg.zlim                = 'maxabs';
% % cfg.zlim                = [-5 5];
% cfg.colormap            = jet;
% cfg.colorbar            = 'yes';
% cfg.masknans            = 'yes';
% cfg.baseline            = 'no';
% cfg.channel             = allExceptMastoids;
% 
% figure(1)
% cfg.parameter           = 'baselineConversionTrial';
% ft_multiplotTFR(cfg,baselineConversionData)
% 
% figure(2)
% cfg.parameter           = 'baselineConversionIndividual';
% ft_multiplotTFR(cfg,baselineConversionData)
% 
% 































































