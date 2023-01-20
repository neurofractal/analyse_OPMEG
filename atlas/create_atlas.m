% 
cd('D:\Github\analyse_OPMEG\atlas\local_global')

% Add Freesurfer_tools and SPM to your path
addpath('D:\Github\OPMsurfer\freesurfer_tools');
addpath('D:\scripts\spm12')

maxDist=5.5; % maximum distance (mm) for projection ? caution, if too large, subcortical structures will also be (falsely) assigned.
visualize=1; % set to 0 if you don't need to visualize
reslice=0; % set to 1 if you want to produce an 0.5 mm version, too.

% load in LH
[~, lhlab,lhctable]=read_annotation('lh.Schaefer2018_100Parcels_17Networks_order.annot');
[lhvtx,lhfaces]=read_surf('lh.pial');

% Delete first superfluous row
lhctable.table(1,:) = [];
lhctable.struct_names(1,:) = [];

[~,lhix]=ismember(lhlab,lhctable.table(:,5));

% load in RH
[~, rhlab,rhctable]=read_annotation('rh.Schaefer2018_100Parcels_17Networks_order.annot');
[rhvtx,rhfaces]=read_surf('rh.pial');

% Delete first superfluous row
rhctable.table(1,:) = [];
rhctable.struct_names(1,:) = [];
[~,rhix]=ismember(rhlab,rhctable.table(:,5));

% Compensate for number of parcels in LH
rhix = rhix+max(lhix);


% visualize Surfaces
figure('Name','Data Conversion','NumberTitle','off','color','w');
subplot(1,2,1)
title('Surface data');
patch('Faces',lhfaces+1,'Vertices',lhvtx,'FaceColor','interp','EdgeColor','none','Facevertexcdata',lhix)
patch('Faces',rhfaces+1,'Vertices',rhvtx,'FaceColor','interp','EdgeColor','none','Facevertexcdata',rhix)
axis equal
axis([-80,80,-110,80])

% load SPM Canonical Brain
nii=spm_vol('D:\scripts\fieldtrip-20210606\template\anatomy\single_subj_T1.nii');
suffx='';
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

%nii.img(nii.img>0)=nii.img(nii.img>0)-1; % first entry in table.table is an indexing entry

nii.fname=['Schaefer2018_100Parcels_17Networks',suffx,'.nii'];
%nii.dt=[16,0];
spm_write_vol(nii,nii.img);
gzip(['Schaefer2018_100Parcels_17Networks',suffx,'.nii']);

% visualize volumetric data
subplot(1,2,2)
title('Volumetric data');
for reg=1:max(rhix)
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
axis([-80,80,-110,80]); 
alpha 0.5;

%% Reead the atlas
atlas = ft_read_atlas('Schaefer2018_100Parcels_17Networks.nii');
count = 1;
atlas.coordsys = 'mni';

%% Add tissue labels
% LH
for h = 1:length(lhctable.struct_names)
    atlas.tissuelabel{count,1} = ['LH_' lhctable.struct_names{h}(15:end)];
    atlas.parcellationlabel{count,1} = ['LH_' lhctable.struct_names{h}(15:end)];
    count = count+1;
end

% RH
for h = 1:length(rhctable.struct_names)
    atlas.tissuelabel{count,1} = ['RH_' rhctable.struct_names{h}(15:end)];
    atlas.parcellationlabel{count,1} =  ['RH_' rhctable.struct_names{h}(15:end)];
    count = count+1;
end

%% Save
atlas.tissue = atlas.parcellation;

save atlas_Schaefer_2018_SPM atlas;

%% Try loading in the atlas again
atlas = ft_read_atlas('atlas_Schaefer_2018_SPM.mat');





