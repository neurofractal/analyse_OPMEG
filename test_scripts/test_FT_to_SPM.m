%% 
% Script for testing Fieldtrip to SPM Pipeline

fieldtripDir    = 'D:\scripts\fieldtrip-20210606';
script_dir      = 'D:\Github\analyse_OPMEG';
path_to_SPM     = 'D:\scripts\spm12';

% Add Fieldtrip to path
disp('Adding the required scripts to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add other scripts to path
addpath(genpath(script_dir));


%% BIDS data directory
data_dir        = 'D:\data\tutorial_OPM_data';

%% Specify Save Directory
save_dir        = fullfile(data_dir,'results','tutorial2');
mkdir(save_dir);
cd(save_dir);


%% Read in the raw BIDS-organised data
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'motor';
cfg.bids.sub    = '001';
cfg.bids.ses    = '001';
cfg.bids.run    = '001';
rawData         = ft_opm_create(cfg);


%% Resample to 1000Hz
cfg                 = [];
cfg.resamplefs      = 1000;
[rawData]           = ft_resampledata(cfg, rawData);

%% Get Fieldtrip path and remove from your path
[~, ft_path] = ft_version;
% Little hack to turn off warnings
id = 'MATLAB:rmpath:DirNotFound';
warning('off',id);
% Remove Fieldtrip path and external/spm12 and external/spm8
rmpath(ft_path);
rmpath(fullfile(ft_path,'external','spm12'));
rmpath(fullfile(ft_path,'external','spm8'));


%% Add SPM
addpath(path_to_SPM);
spm('defaults', 'eeg');

% Convert data from FT to SPM
[D] = spm_eeg_ft2spm(rawData,'data');

S = [];
S.data = rawData;
S.sMRI = 'D:\data\tutorial_OPM_data\sub-001\ses-002\anat\001.nii';
[D] = spm_opm_create_ft(S);


%% Perform Coreg Using spm_opm_create.m code

sMRI = 'D:\data\tutorial_OPM_data\sub-001\ses-002\anat\001.nii';

%initially used inverse normalised meshes
D = spm_eeg_inv_mesh_ui(D,1,sMRI,3);
save(D);

miMat = zeros(3,3);
fiMat = zeros(3,3);
fid=[];
try % to read the coordsystem.json
    coord = spm_load(S.coordsystem);
    fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
    fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
    fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
    miMat(1,:) = coord.AnatomicalLandmarkCoordinates.coil1;
    miMat(2,:) = coord.AnatomicalLandmarkCoordinates.coil2;
    miMat(3,:) = coord.AnatomicalLandmarkCoordinates.coil3;
    fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
    fid.fid.pnt =fiMat;
    fid.pos= []; % headshape field that is left blank (GRB)
    M = fid;
    M.fid.pnt=miMat;
    M.pnt = D.inv{1}.mesh.fid.pnt;
catch
    try % to find the BIDS coordsystem.json
        coord = spm_load(coordFile);
        fiMat(1,:) = coord.HeadCoilCoordinates.coil1;
        fiMat(2,:) = coord.HeadCoilCoordinates.coil2;
        fiMat(3,:) = coord.HeadCoilCoordinates.coil3;
        miMat(1,:) = coord.AnatomicalLandmarkCoordinates.coil1;
        miMat(2,:) = coord.AnatomicalLandmarkCoordinates.coil2;
        miMat(3,:) = coord.AnatomicalLandmarkCoordinates.coil3;
        fid.fid.label = fieldnames(coord.HeadCoilCoordinates);
        fid.fid.pnt =fiMat;
        fid.pos= []; % headshape field that is left blank (GRB)
        M = fid;
        M.fid.pnt=miMat;
        M.pnt = D.inv{1}.mesh.fid.pnt;
    catch % DEFAULT: transform between fiducials and anatomy is identity
        fid.fid.label = {'nas', 'lpa', 'rpa'}';
        fid.fid.pnt = D.inv{1}.mesh.fid.fid.pnt(1:3,:); % George O'Neill
        fid.pos= [];
        M = fid;
        M.pnt = D.inv{1}.mesh.fid.pnt; % George O'Neill
    end
end



D = fiducials(D, fid);
save(D);
f=fiducials(D);
if ~isempty(f.pnt)
    D = spm_eeg_inv_datareg_ui(D,1,f,M,1);
else
    f.pnt = zeros(0,3);
    D = spm_eeg_inv_datareg_ui(D,1,f,M,0);
end



D.inv{1}.forward.voltype = 'Single Shell';
D = spm_eeg_inv_forward(D);
spm_eeg_inv_checkforward(D,1,1);


D.inv{1}.forward.voltype = 'Single Shell';
D = spm_eeg_inv_forward(D);
nverts = length(D.inv{1}.forward.mesh.vert);
spm_eeg_inv_checkforward(D,1,1);
save(D);




