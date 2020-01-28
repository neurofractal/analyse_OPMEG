%% (1) Get all of the directories and references ready
% (1.1) Change depending on what data is being cleaned.
ncfg.fileinfo.sub           = '001';
ncfg.fileinfo.folder        = 'OPM-AEF_Practice';
ncfg.fileinfo.ses           = '001';
ncfg.fileinfo.task          = 'beep';
ncfg.fileinfo.run           = '001';

ncfg.chaninfo.trigger       = 'NI-TRIG';
ncfg.chaninfo.channel       = {'G2-MW-RAD';'G2-OH-RAD';'G2-OF-RAD';'G2-DU-RAD';...
                                'G2-MU-RAD';'G2-MX-RAD';'G2-MY-RAD';'G2-OI-RAD';...
                                'G2-N2-RAD';'G2-MZ-RAD';'G2-OK-RAD';'G2-OG-RAD'};
ncfg.chaninfo.reference     = {'G2-N0-RAD';'G2-N4-RAD';'G2-N3-RAD';'G2-MV-RAD'};
ncfg.chaninfo.off           = {'G2-N1-RAD';'G2-MT-RAD';'G2-OJ-RAD'};

ncfg.slotinfo.channel       = {'44';'51';'53';'45';...
                                '32';'33';'60';'35';...
                                '52';'42';'43';'34'};
ncfg.slotinfo.off           = {'29';'30';'31'};
ncfg.gradpairs              = {'32','33'; '34','35';...
                                '42','43'; '44','45';...
                                '52','53'; '51','60'};

% (1.2) This shouldn't need muching changing is Box Sync is working. 
pcName          = getenv('COMPUTERNAME');
if strcmp(pcName,'PRINCE')
    matlabDir       = 'D:\MATLAB';
else
    error('Error: Unknown PC')
end

fieldtripDir    = strcat(matlabDir,'\Toolboxes\fieldtrip-20190419');
workingDir      = strcat(matlabDir,'\Analysis');
functionsDir    = strcat(workingDir,'\Functions');
templatesDir    = strcat(workingDir,'\Templates');
outputsDir      = strcat(workingDir,'\Outputs');
dataDir         = strcat(workingDir,strcat('\Data\',ncfg.fileinfo.folder));

addpath(fieldtripDir)
ft_defaults;
cd(workingDir)
addpath(functionsDir,strcat(functionsDir,'\GUI'),templatesDir,workingDir)

% (1.3) Find the right files for this participant.
ncfg.fileinfo.filename        = strcat('sub-',ncfg.fileinfo.sub,'_ses-',...
                                    ncfg.fileinfo.ses,'_task-',ncfg.fileinfo.task,...
                                    '_run-',ncfg.fileinfo.run,'_meg.dat');
ncfg.fileinfo.dataset         = strcat(dataDir,'\',ncfg.fileinfo.filename);

% Load layout file;
load(strcat(templatesDir,'\OPM_layouts\AEF_layout.mat'))

% (1.4) Load some templates that will be used throughout.
clear *Dir pcName

%% (2) Start preprocessing.
% (2.1) Read in the raw data.
cfg                     = [];
cfg.continuous          = 'yes';
cfg.datafile            = ncfg.fileinfo.dataset;
cfg.headerfile          = strcat(ncfg.fileinfo.dataset(1:end-3),'mat');
cfg.chantype            = 'megmag';
cfg.channel             = [ncfg.chaninfo.channel; ncfg.chaninfo.reference;...
                            ncfg.chaninfo.trigger];
cfg.dataformat          = 'spmeeg_mat';
rawData                 = ft_preprocessing(cfg);

% Extract the trigger channel for later.
cfg                     = [];
cfg.channel             = ncfg.chaninfo.trigger;
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

% Make gradiometers from their pairs. Fairly random for now...
gradData                = [];
gradData.hdr            = denoiseData.hdr;
gradData.fsample        = denoiseData.fsample;
gradData.sampleinfo     = denoiseData.sampleinfo;
gradData.cfg            = denoiseData.cfg;
gradData.time           = denoiseData.time;
gradData.trial{1}       = zeros(length(ncfg.gradpairs(:,1)),length(denoiseData.trial{1}(1,:)));
for pairIdx = 1:length(ncfg.gradpairs(:,1))
    refIdx                  = find(ismember(ncfg.slotinfo.channel,ncfg.gradpairs(pairIdx,1)));
    outIdx                  = find(ismember(ncfg.slotinfo.channel,ncfg.gradpairs{pairIdx,2}));
    gradName                = strcat(ncfg.chaninfo.channel{outIdx}(end-5:end-4),'-',ncfg.chaninfo.channel{refIdx}(end-5:end-4),'_GRAD');

    gradData.label{pairIdx,1}       = gradName;
    gradData.trial{1}(pairIdx,:)    = denoiseData.trial{1}(outIdx,:) - denoiseData.trial{1}(refIdx,:);
end
gradData.hdr.label          = gradData.label;
gradData.hdr.chantype       = cell(length(ncfg.gradpairs(:,1)),1);
gradData.hdr.chantype(:)    = {'meggrad'};
gradData.hdr.chanunit       = cell(length(ncfg.gradpairs(:,1)),1);
gradData.hdr.chanunit(:)    = {'fT'};

% Downsample from 6000Hz
cfg                     = [];
cfg.resamplefs          = 1000;
cfg.detrend             = 'no';
dsGradData              = ft_resampledata(cfg,gradData);
dsTriggerData           = ft_resampledata(cfg,triggerData);
dsRefData               = ft_resampledata(cfg,referenceData);
dsMagData               = ft_resampledata(cfg,denoiseData);

% (2.2) Apply some filters.
cfg                     = [];
cfg.channel             = 'all';
cfg.bpfilter            = 'yes';
cfg.bpfreq              = [1 30];
% cfg.bpfiltord           = 3;
cfg.demean              = 'no';
filteredGradData        = ft_preprocessing(cfg,dsGradData);
filteredRefData         = ft_preprocessing(cfg,dsRefData);
filteredMagData         = ft_preprocessing(cfg,dsMagData);

% Have a quick look at the filtered data with the trigger channel.
% cfg                     = [];
% combinedData            = ft_appenddata(cfg,filteredGradData,dsTriggerData,filteredRefData);
% 
% cfg                     = [];
% cfg.channel             = [filteredGradData.label; triggerData.label; filteredRefData.label];
% cfg.viewmode            = 'vertical';
% cfg.continuous          = 'no';
% cfg.plotlabels          = 'some';
% cfg.ylim                = [-7000 7000];
% cfg.chanscale           = [ones(length(filteredGradData.label),1);0.01;ones(4,1)];
% cfg.blocksize           = 5;
% ft_databrowser(cfg,combinedData);

% (2.3) Now set the trial definition to select the whole trial.
cfg                     = [];
cfg.dataset             = ncfg.fileinfo.dataset;
cfg.trialdef.prestim    = .1;
cfg.trialdef.poststim   = .3;
cfg.trialdef.trigchan   = ncfg.chaninfo.trigger;
cfg.trialfun            = 'OPM_TrialFun';
cfg.trialdef.downsample = 1000;
trialCfg                = ft_definetrial(cfg);
epochedGradData             = ft_redefinetrial(trialCfg,filteredGradData);
epochedMagData             = ft_redefinetrial(trialCfg,filteredMagData);

% % (3.2) Take a look at the data and flag any trials for removal.
% cfg                     = [];
% cfg.channel             = 'all';
% cfg.viewmode            = 'vertical';
% cfg.continuous          = 'no';
% cfg.plotlabels          = 'some';
% cfg.ylim                = [-7000 7000];
% ft_databrowser(cfg,epochedGradData);
% ft_databrowser(cfg,epochedMagData);


%% Topoplots
% Use the pairs to find the midpoint between sensors for grads.
layout_grad                 = [];
layout_grad.width           = [zeros(length(ncfg.gradpairs(:,1)),1);layout_opm.width(end-1:end)];
layout_grad.height          = [zeros(length(ncfg.gradpairs(:,1)),1);layout_opm.height(end-1:end)];
layout_grad.mask            = layout_opm.mask;
layout_grad.outline         = layout_opm.outline;
layout_grad.cfg             = layout_opm.cfg;
layout_grad.pos             = [zeros(length(ncfg.gradpairs(:,1)),2);layout_opm.pos(end-1:end,1:2)];
layout_grad.label           = [cell(length(ncfg.gradpairs(:,1)),1);{'COMNT';'SCALE'}];
for pairIdx = 1:length(ncfg.gradpairs(:,1))
    refIdx                          = find(ismember(ncfg.slotinfo.channel,ncfg.gradpairs(pairIdx,1)));
    outIdx                          = find(ismember(ncfg.slotinfo.channel,ncfg.gradpairs{pairIdx,2}));
    
    layout_grad.label{pairIdx,1}    = strcat(ncfg.chaninfo.channel{outIdx}(end-5:end-4),'-',ncfg.chaninfo.channel{refIdx}(end-5:end-4),'_GRAD');
    
    posRefIdx                       = find(ismember(layout_opm.label,ncfg.gradpairs(pairIdx,1)));
    posOutIdx                       = find(ismember(layout_opm.label,ncfg.gradpairs(pairIdx,2)));
    
    layout_grad.pos(pairIdx,:)      = mean([layout_opm.pos(posRefIdx,:);layout_opm.pos(posOutIdx,:)]);
    layout_grad.width(pairIdx,1)    = layout_opm.width(1);
    layout_grad.height(pairIdx,1)   = layout_opm.height(1);
end

% And apply the right labels to the mag layout
layout_mag                  = [];
layout_mag.width            = [zeros(length(ncfg.chaninfo.channel(:,1)),1);layout_opm.width(end-1:end)];
layout_mag.height           = [zeros(length(ncfg.chaninfo.channel(:,1)),1);layout_opm.height(end-1:end)];
layout_mag.mask             = layout_opm.mask;
layout_mag.outline          = layout_opm.outline;
layout_mag.cfg              = layout_opm.cfg;
layout_mag.pos              = [zeros(length(ncfg.chaninfo.channel(:,1)),2);layout_opm.pos(end-1:end,1:2)];
layout_mag.label            = [cell(length(ncfg.chaninfo.channel(:,1)),1);{'COMNT';'SCALE'}];
for pairIdx = 1:length(ncfg.chaninfo.channel(:,1))
    layout_mag.label{pairIdx,1}     = ncfg.chaninfo.channel{pairIdx,1};
    
    slotIdx                         = find(ismember(layout_opm.label,ncfg.slotinfo.channel(pairIdx,1)));
    layout_mag.pos(pairIdx,:)       = layout_opm.pos(slotIdx,:);
    layout_mag.width(pairIdx,1)     = layout_opm.width(1);
    layout_mag.height(pairIdx,1)    = layout_opm.height(1);
end

%% (5) Run an ICA
% (5.2) Run ICA on this downsampled data.
cfg                     = [];
cfg.method              = 'runica';
cfg.numcomponent        = 12;
componentsMag           = ft_componentanalysis(cfg,epochedMagData);
cfg.numcomponent        = 6;
componentsGrad          = ft_componentanalysis(cfg,epochedGradData);

% (5.3) Use the unmixing matrix from the downsampled data to estimate the
% identified components timecourse in the original data.
% (5.4) Use the view IC GUI. Requires Matlab 2017b or above.
cfg                     = [];
cfg.layout              =  layout_mag;
cfg.viewmode            = 'component';
ft_databrowser(cfg,componentsMag)
cfg.layout              =  layout_grad;
ft_databrowser(cfg,componentsGrad)

% Record input of components to be removed into the console.
disp('Enter components to remove in the form [1 4 7].');
componentsToRemoveMag  	= input('User input:');
componentsToRemoveGrad  = input('User input:');

% Remove the undesired components components, if any.
cfg                     = [];
cfg.demean              = 'no'; % Note, this has been enabled by default!
cfg.component           = componentsToRemoveMag; 
cleanMag               = ft_rejectcomponent(cfg,componentsMag,epochedMagData);
cfg.component           = componentsToRemoveGrad; 
cleanGrad               = ft_rejectcomponent(cfg,componentsGrad,epochedGradData);


%% Timelock avg
cfg                     = [];
timelockGradData        = ft_timelockanalysis(cfg,cleanGrad);
timelockMagData         = ft_timelockanalysis(cfg,cleanMag);

cfg                     = [];
cfg.avgoverchan         = 'yes';
timelockGradAvg         = ft_selectdata(cfg,timelockGradData);
timelockMagAvg          = ft_selectdata(cfg,timelockMagData);

% Plot
cfg                     = [];
cfg.viewmode            = 'butterfly';
ft_databrowser(cfg,timelockGradData);
ft_databrowser(cfg,timelockMagData);

ft_databrowser([],timelockGradAvg);
ft_databrowser([],timelockMagAvg);

% Plot
cfg                     = [];
cfg.layout              = layout_grad;
cfg.gridscale          = 100;
cfg.colorbar           = 'yes';
cfg.parameter          = 'avg';
cfg.marker             = 'labels';
cfg.xlim               = [.08 .12];

cfg.zlim               = 'maxmin';
% cfg.channel            = 
cfg.baseline           = 'no';
cfg.markersize         = 2;
cfg.interplimits       = 'electrodes';

cfg.comment            = 'auto';
cfg.interactive        = 'yes';
figure
ft_topoplotER(cfg,timelockGradData);
figure
cfg.layout              = layout_mag;
ft_topoplotER(cfg,timelockMagData);
























