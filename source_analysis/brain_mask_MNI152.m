function mri_interp = brain_mask_MNI152(mri_path)

% Read MRI
mri = ft_read_mri(mri_path);

% Convert to mm
if ~strcmp(mri.unit,'mm')
    mri = ft_convert_units(mri,'mm');
end

% in order to do this we segment the brain from the anatomical template
if ~exist('mri_mask')

    X = load(which('mri_mask.mat'));
    mri_mask = X.mri_mask;

end


% Interpolate
cfg                         = [];
cfg.interpmethod            = 'nearest';
cfg.parameter               = 'anatomy';
mri_interp                = ft_sourceinterpolate(cfg, mri,mri_mask);

% Mask
log_array = ~logical(mri_mask.anatomy);
mri_interp.anatomy(log_array) = NaN;

[filepath, name, ext] = fileparts(mri_path);
cd(filepath);

new_filename = [name '_MNI152' ext];

mri_interp = rmfield(mri_interp,'cfg');

% Output
cfg = [];
cfg.parameter = 'anatomy';
cfg.filename = new_filename;
ft_sourcewrite(cfg,mri_interp);

end