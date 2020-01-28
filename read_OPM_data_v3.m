% (1.2) This shouldn't need muching changing is Box Sync is working. 
pcName          = getenv('COMPUTERNAME');
if strcmp(pcName,'PRINCE')
    matlabDir       = 'D:\MATLAB';
else
    error('Error: Unknown PC')
end

fieldtripDir    = strcat(matlabDir,'\Toolboxes\fieldtrip-20200121');
workingDir      = strcat(matlabDir,'\Analysis');
functionsDir    = strcat(workingDir,'\Functions');
templatesDir    = strcat(workingDir,'\Templates');
outputsDir      = strcat(workingDir,'\Outputs');
% dataDir         = strcat(workingDir,strcat('\Data\',ncfg.fileinfo.folder));

addpath(fieldtripDir)
ft_defaults;
cd(workingDir)
addpath(functionsDir,strcat(functionsDir,'\GUI'),templatesDir,workingDir)

clear *Dir pcName



% Notes.
% > Need to add optitrack data as coilpos/coilori?

% Or common average ref
ncfg.grad.references        = repmat(ncfg.slotinfo.channel',[length(ncfg.slotinfo.channel(:,1)) 1]);
ncfg.grad.main              = ncfg.slotinfo.channel;

ncfg.sample                 = 1000;
ncfg.filters.bpfilter       = 'yes';
ncfg.filters.bpfreq         = [0.5 70];
ncfg.filters.dftfilter      = 'yes';
ncfg.filters.demean         = 'no';

%% (2) Start preprocessing.
% (2.1) Read in the raw data.
cfg             = [];
cfg.data        = 'meg.bin';
cfg.coordystem  = 'coordsystem.json';
cfg.positions   = 'positions.tsv';
cfg.channels    = 'channels.tsv';
cfg.meg         = 'meg.json';

rawData         = ft_opm_create(cfg);

% Extract the trigger channel for later. This doesn't work yet. 
cfg                     = [];
cfg.chantype            = {'TRIG'};
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




























