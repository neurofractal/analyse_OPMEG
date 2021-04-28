function convertHCP_MMP_to_MNI_volumetric
% small script that will project HCM-MMP1.0 fsaverage annot files from here
% (https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446)
% onto MNI 152 2009a NLIN GM space.
% (c) 2016 Andreas Horn ?�BIDMC Boston, MA, USA

% to run, you will need SPM8/12 and the freesurfer/matlab directory on your
% ML path. You need to put the following files into the working directory:
% - lh.pial and rh.pial from freesurfer/subjects/fsaverage/surf
% - mni_icbm152_gm_tal_nlin_asym_09a.nii from http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009
% - freesurfer converted annotation files (rh.HCP-MMP1.annot and
%   lh.HCP-MMP1.annot) from https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446
% - if you want to produce an HD version (0.5 0.5 0.5) you also need Lead
%   Neuroimaging Suite (leadsuite.io)

maxDist=4.5; % maximum distance (mm) for projection ? caution, if too large, subcortical structures will also be (falsely) assigned.
visualize=1; % set to 0 if you don't need to visualize
reslice=1; % set to 1 if you want to produce an 0.5 mm version, too.

% load in LH
[~, lhlab,lhctable]=read_annotation('lh.HCP-MMP1.annot');
[lhvtx,lhfaces]=read_surf('lh.pial');
[~,lhix]=ismember(lhlab,lhctable.table(:,5));

% load in RH
[~, rhlab,rhctable]=read_annotation('rh.HCP-MMP1.annot');
[rhvtx,rhfaces]=read_surf('rh.pial');
[~,rhix]=ismember(rhlab,rhctable.table(:,5));

% visualize Surfaces
if visualize
    figure('Name','Data Conversion','NumberTitle','off','color','w');
    subplot(1,2,1)
    title('Surface data');
    patch('Faces',lhfaces+1,'Vertices',lhvtx,'FaceColor','interp','EdgeColor','none','Facevertexcdata',lhix)
    patch('Faces',rhfaces+1,'Vertices',rhvtx,'FaceColor','interp','EdgeColor','none','Facevertexcdata',rhix)
    axis equal
    axis([-80,80,-110,80])
end

% load ICBM152_nlin 2009a GM template (downloaded from here)
% (http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009):
for resl=0:reslice
    if resl==1
        % resample image to 0.5 0.5 0.5 mm space
        ea_reslice_nii('mni_icbm152_gm_tal_nlin_asym_09a.nii','mni_icbm152_gm_tal_nlin_asym_09a_05.nii',[0.5 0.5 0.5]);
        nii=spm_vol('mni_icbm152_gm_tal_nlin_asym_09a_05.nii');
        suffx='_hd';
    else
        nii=spm_vol('mni_icbm152_gm_tal_nlin_asym_09a.nii');
        suffx='';
    end
    nii.img=spm_read_vols(nii);
    
    % find positive GM voxels (>0.5 probability):
    [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>0.5));
    XYZ_vx=[xx,yy,zz,ones(length(xx),1)]';
    XYZ_mm=nii.mat*XYZ_vx;
    
    [assign,D]=knnsearch([lhvtx;rhvtx],XYZ_mm(1:3,:)');
    
    
    % remove voxels too far away from vertices (defined by maxDist);
    xx(D>maxDist)=[];
    yy(D>maxDist)=[];
    zz(D>maxDist)=[];
    assign(D>maxDist)=[];
    
    ldat=[lhix;rhix];
    
    % set all voxels to zero
    nii.img(:)=0;
    % assign voxel intensities based on label table
    nii.img(sub2ind(size(nii.img),xx,yy,zz))=ldat(assign);
    
    nii.img(nii.img>0)=nii.img(nii.img>0)-1; % first entry in table.table is an indexing entry
    
    nii.fname=['HCP-MMP1_on_MNI152_ICBM2009a_nlin',suffx,'.nii'];
    nii.dt=[16,0];
    spm_write_vol(nii,nii.img);
    gzip(['HCP-MMP1_on_MNI152_ICBM2009a_nlin',suffx,'.nii']);
    delete(['HCP-MMP1_on_MNI152_ICBM2009a_nlin',suffx,'.nii']);
    
    % visualize volumetric data
    if visualize && ~resl
        subplot(1,2,2)
        title('Volumetric data');
        for reg=1:181
            fv=isosurface(nii.img==reg,0.5);
            try % swap xy entries for non-empty labels (as per xy swap in isosurface).
                fv.vertices=[fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)];
            end
            afv(reg).vertices=fv.vertices;
            afv(reg).faces=fv.faces;
            afv(reg).facevertexcdata=repmat(reg,length(fv.vertices),1);
        end
        afv=ea_concatfv(afv);
        
        % convert vertices from voxel to mm space
        afv.vertices=[afv.vertices,ones(length(afv.vertices),1)]';
        afv.vertices=nii.mat*afv.vertices;
        afv.vertices=afv.vertices(1:3,:)';
        patch(afv,'EdgeColor','none','FaceColor','interp')
        axis equal
        axis([-80,80,-110,80])
    end
    
    
end



% Finally, export .txt label file
f=fopen('HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt','w');

for reg=2:181
    fprintf(f,'%d %s\n',reg-1,lhctable.struct_names{reg});
end
fclose(f);

function afv=ea_concatfv(fv)
afv.faces=[];
afv.vertices=[];

for f=1:length(fv)
    afv.faces=[afv.faces;fv(f).faces+length(afv.vertices)];
    afv.vertices=[afv.vertices;fv(f).vertices];
end


if isfield(fv,'facevertexcdata')
    afv.facevertexcdata=[];
    for f=1:length(fv)
        afv.facevertexcdata=[afv.facevertexcdata;fv(f).facevertexcdata];
    end
end
