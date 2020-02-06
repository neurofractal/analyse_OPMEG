%% 
% Script to test ft_opm_create and some basic processing on sample AEF data
%%

% %% NIC
% % (1.2) This shouldn't need muching changing is Box Sync is working. 
% pcName          = getenv('COMPUTERNAME');
% if strcmp(pcName,'PRINCE')
%     matlabDir       = 'D:\MATLAB';
% else
%     error('Error: Unknown PC')
% end
% 
% fieldtripDir    = strcat(matlabDir,'\Toolboxes\fieldtrip-20200121');
% workingDir      = strcat(matlabDir,'\Analysis');
% functionsDir    = strcat(workingDir,'\Functions');
% templatesDir    = strcat(workingDir,'\Templates');
% outputsDir      = strcat(workingDir,'\Outputs');
% % dataDir         = strcat(workingDir,strcat('\Data\',ncfg.fileinfo.folder));

%% Paths (RS)
fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
data_dir        = 'D:\data\AEF2\sub-001\ses-001\meg';
save_dir        = 'D:\data\AEF2\sub-001\ses-001\meg';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));

% cd to save dir
cd(save_dir)

% % Or common average ref
% ncfg.sample                 = 1000;
% ncfg.filters.bpfilter       = 'yes';
% ncfg.filters.bpfreq         = [0.5 70];
% ncfg.filters.dftfilter      = 'yes';
% ncfg.filters.demean         = 'no';

%% (2) Start preprocessing.
% Read in the raw data.
cfg             = [];
cfg.precision   = 'single';
cfg.data        = 'sub-001_ses-001_task-beep_run-001_meg.bin';
rawData         = ft_opm_create(cfg);

% Read in the raw data using BIDS
cfg             = [];
cfg.folder      = data_dir;
cfg.bids.task   = 'beep';
cfg.bids.sub    = '001';
cfg.bids.ses    = '001';
cfg.bids.run    = '001';
rawData2         = ft_opm_create(cfg);
clear rawData2

% Plot using ft_databrowser
ft_databrowser([],rawData);

%% Select the sensors of interest
% Select only sensors measuring radial fields
disp('Selecting ensors measuring radial fields only');
cfg             = [];
cfg.channel     = vertcat('-G2-N1-RAD','-G2-MT-RAD','-G2-OJ-RAD',...
    rawData.label(contains(rawData.hdr.fieldori,'RAD')));
rawData_rad     = ft_selectdata(cfg, rawData);

% Denoise
cfg                     = [];
cfg.refchannel          = {'G2-N0-RAD','G2-N4-RAD','G2-N3-RAD','G2-MV-RAD'};
%cfg.channel             = [ncfg.chaninfo.channel];
cfg.truncate            = 'no';
cfg.zscore              = 'no';
cfg.updatesens          = 'yes';
denoiseData             = ft_denoise_pca(cfg,rawData_rad);

cfg = [];
cfg.resamplefs = 1000;
[rawData_rad] = ft_resampledata(cfg, rawData_rad);

% Band-pass filter between 0.5-250Hz to help visualisation
cfg = [];
cfg.demean      = 'yes';
cfg.continuous  = 'yes';
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [1 30];
rawData_rad         = ft_preprocessing(cfg,rawData_rad);

% Plot using ft_databrowser
ft_databrowser([],rawData_rad);

cfg = [];
cfg.channel         = vertcat(rawData_rad.label,{'-G2-N0-RAD',...
    '-G2-N4-RAD','-G2-N3-RAD','-G2-MV-RAD','-G2-OF-RAD'}');
rawData_rad_no_ref = ft_selectdata(cfg, rawData_rad);

ft_databrowser([],rawData_rad_no_ref);

%%
cfg = [];
cfg.dataset                 = 'sub-001_ses-001_task-beep_run-001_meg.bin';
cfg.trialdef.trigchan       = 'trigger';
cfg.trialdef.downsample     = 1000;
%cfg.continuous              = 'yes';
cfg.trialdef.prestim        = 0.1;        % pre-stimulus interval
cfg.trialdef.poststim       = 0.4;        % post-stimulus interval
cfg.trialfun                = 'OPM_TrialFun_RS';
banana                    = ft_definetrial(cfg);

%%
% Redefines the filtered data
cfg                         = [];
data                        = ft_redefinetrial(banana,rawData_rad_no_ref);

ft_databrowser([],data);
avg = ft_timelockanalysis([],data);

% figure;
% for i = 1:12
%     
%     cfg = [];
%     cfg.baseline = [-0.1 0];
%     cfg.xlim = [-0.1 0.4];
%     cfg.baselinetype  = 'relative';
%     %cfg.ylim = [-1e-13 3e-13];
%     cfg.channel = avg.label{i};
%     title(avg.label{i});
%     subplot(3,4,i);ft_singleplotER(cfg,avg);
% end
% 
% 
% cfg = [];
%     cfg.baseline = [-0.1 0];
%     cfg.xlim = [-0.1 0.4];
%     cfg.baselinetype  = 'relative';
%     %cfg.ylim = [-1e-13 3e-13];
%     %cfg.channel = avg.label{i};
%     title('ALL');
% figure;ft_singleplotER(cfg,avg);

cfg = [];
cfg.parameter = 'avg';
ft_databrowser(cfg,avg);













% Extract the trigger channel for later. This doesn't work yet. 
cfg                     = [];
cfg.hdr.chantype        = 'trigger'; % Wrong.
triggerData             = ft_selectdata(cfg,rawData);
% And grab the reference channels.
cfg                     = [];
cfg.channel             = ncfg.chaninfo.reference;
referenceData           = ft_selectdata(cfg,rawData);

% And remove it, leaving just the OPMS. 
cfg                     = [];
cfg.channel             = [ncfg.chaninfo.channel; ncfg.chaninfo.reference];
rawOpmData              = ft_selectdata(cfg,rawData);

% (2.2) Denoise using reference.
cfg                     = [];
cfg.refchannel          = ncfg.chaninfo.reference;
cfg.channel             = [ncfg.chaninfo.channel];
cfg.truncate            = 'no';
cfg.zscore              = 'no';
% cfg.pertrial            = 'yes';
cfg.updatesens          = 'yes';
denoiseData             = ft_denoise_pca(cfg,rawOpmData);

% % Skip that?
% cfg                     = [];
% cfg.channel             = ncfg.chaninfo.channel;
% denoiseData             = ft_selectdata(cfg,rawOpmData);

% Make gradiometers from their pairs. Fairly random for now...
gradData                = [];
gradData.hdr            = denoiseData.hdr;
gradData.fsample        = denoiseData.fsample;
gradData.sampleinfo     = denoiseData.sampleinfo;
gradData.cfg            = denoiseData.cfg;
gradData.time           = denoiseData.time;
gradData.trial{1}       = zeros(length(ncfg.grad.main(:,1)),length(denoiseData.trial{1}(1,:)));
for pairIdx = 1:length(ncfg.grad.main(:,1))
    mainIdx                 = find(ismember(ncfg.slotinfo.channel,ncfg.grad.main{pairIdx,1}));
    emptyIdx                = ismember(ncfg.grad.references(pairIdx,:),'00');
    references              = ncfg.grad.references(pairIdx,~emptyIdx);
    refIdx                  = zeros(size(references));
    for refIdxIdx = 1:length(references)
        refIdx(refIdxIdx)   = find(ismember(ncfg.slotinfo.channel,references(refIdxIdx)));
    end
    
    gradData.label{pairIdx,1}       = ncfg.chaninfo.channel{mainIdx};
    gradData.trial{1}(pairIdx,:)    = denoiseData.trial{1}(mainIdx,:) - mean(denoiseData.trial{1}(refIdx,:),1);
end
gradData.hdr.label          = gradData.label;
gradData.hdr.chantype       = cell(length(ncfg.grad.main(:,1)),1);
gradData.hdr.chantype(:)    = {'meggrad'};
gradData.hdr.chanunit       = cell(length(ncfg.grad.main(:,1)),1);
gradData.hdr.chanunit(:)    = {'fT'};

clear *Idx gradName references rawData rawOpmData

% Downsample from 6000Hz
cfg                     = [];
cfg.resamplefs          = ncfg.sample;
cfg.detrend             = 'no';
dsGradData              = ft_resampledata(cfg,gradData);
dsTriggerData           = ft_resampledata(cfg,triggerData);
dsRefData               = ft_resampledata(cfg,referenceData);
dsMagData               = ft_resampledata(cfg,denoiseData);

% (2.2) Apply some filters.
cfg                     = ncfg.filters;
filteredGradData        = ft_preprocessing(cfg,dsGradData);
filteredRefData         = ft_preprocessing(cfg,dsRefData);
filteredMagData         = ft_preprocessing(cfg,dsMagData);




























