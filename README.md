# analyse_OPMEG
## Scripts to analyse data acquired from the OP-MEG Lab, UCL

![opm](./old/opm_image.jpg)

Copyright (C) 2020-21 Wellcome Trust Centre for Neuroimaging

### Authors:  
- Robert Seymour (rob.seymour@ucl.ac.uk)
- Nic Alexander  (n.alexander@ucl.ac.uk);
- Tim West
- George O'Neill
          
### Load in data

Data from the UCL OP-MEG lab is stored in .bin format, with BIDS-compliant descriptor files, for example: `channels.tsv`, `meg.json`, `coordinatesystem.json`, `positions.tsv`.

**Method 1:**
```matlab
%% Load the OPM data
% Read in the BIDS-organised raw data
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir; % This is overall BIDS directory
cfg.precision   = 'single';
cfg.bids.task   = 'aef';
cfg.bids.sub    = '001';
cfg.bids.ses    = '001';
cfg.bids.run    = '001';
rawData         = ft_opm_create(cfg);
```

**Method 2:**
```matlab
%% Load the OPM data
% Read in the BIDS-organised raw data
disp('Loading data...');
cfg             = [];
cfg.data        = path_to_bin_file;
%                 (requires .json and _channels.tsv in same folder
%                 with same name as the .bin file)
rawData         = ft_opm_create(cfg);
```

### Calculate PSD

```matlab
cfg                 = [];
cfg.channel         = ft_channelselection('MEG',rawData,'qzfm_g2');
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0.1 150];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow1 freq]         = ft_opm_psd(cfg,rawData);
```

### Compare PSDs (before and after pre-processing)

```matlab
% Plot PSD
cfg                     = [];
cfg.channel             = ft_channelselection('MEG',rawData,'qzfm_g2');
cfg.trial_length        = 10;
cfg.method              = 'tim';
cfg.foi                 = [0.1 100];
cfg.plot                = 'yes';
cfg.plot_chans          = 'yes';
cfg.plot_ci             = 'no';
cfg.plot_legend         = 'no';
cfg.transparency        = 0.3;
[pow freq]              = ft_opm_psd_compare(cfg,data1,data2);
```


### Implemented Interference Suppression Algorithms

- HFC (Coming Soon)
- [Zapline](./preprocessing/ft_zapline_window.m)
- [Synthetic Gradiometry](./preprocessing/ft_opm_synth_gradiometer_window.m)
- [Regress motion capture data from OPM data](./preprocessing/regress_motive_OPMdata)


### Example Analysis Pipelines

Outlined in Seymour et al., (2021). Interference suppression techniques for OPM-based MEG: Opportunities and challenges. Under Review.

#### 1. Auditory evoked field paradigm during participant movement
- Data: https://doi.org/10.5281/zenodo.5539414
- Script: [pipeline_tutorial_1.m](https://github.com/FIL-OPMEG/tutorials_interference/blob/main/pipeline_tutorial_2.m)

#### 2. Motor-beta power changes during a finger-tapping paradigm
- Data: https://doi.org/10.5281/zenodo.5539414
- Script: [pipeline_tutorial_2.m](https://github.com/FIL-OPMEG/tutorials_interference/blob/main/pipeline_tutorial_2.m)


*Other messier scripts can be found here: https://github.com/FIL-OPMEG/opm_benchmarking_2020*

### Other

- [Extract Sensor Pos and Ori example](./test_scripts/html/extractSensorPositions_Example.html)


