function D = spm_eeg_inv_cortex_custom(D, val, ctx, space)
% Function for supplying custom cortex mesh for source reconstruction
% FORMAT D = spm_eeg_inv_cortex_custom(D, val, mesh, space)
% Inputs
%   D        - input data struct (required)
%   val      - inv structure space number
%   ctx      - path to gifti file
%   space    - space which gifti file is currently in ('mni'|'inv')
% Outputs
%   D        - same data struct including the new cortex
%__________________________________________________________________________

% Copyright (C) 2020 Wellcome Center for Human Neuroimaging

% George O'Neill
% $Id$



switch space
    case 'mni'
        
        def = D.inv{val}.mesh.def;
        Tmesh = spm_swarp(ctx,def);
        
        D.inv{val}.mesh.tess_mni = export(gifti(ctx),'spm');
        
        if ~D.inv{val}.mesh.template
            [path, name] = spm_fileparts(D.inv{val}.mesh.sMRI);
            outname = fullfile(path, [name ctx]); 
            save(gifti(Tmesh), outname);
            D.inv{val}.mesh.tess_ctx = outname;
        else
            D.inv{val}.mesh.tess_ctx = Tmesh;
        end
        
    case 'inv'
        
        defs.comp{1}.inv.comp{1}.def = {D.inv{val}.mesh.def};
        defs.comp{1}.inv.space = {D.inv{val}.mesh.sMRI};
        defs.out{1}.surf.surface = {ctx};
        defs.out{1}.surf.savedir.savesrc = 1;
        out = spm_deformations(defs);
        D.inv{val}.mesh.tess_ctx = ctx;
        D.inv{val}.mesh.tess_mni  = export(gifti(out.surf{1}), 'spm');
        
end