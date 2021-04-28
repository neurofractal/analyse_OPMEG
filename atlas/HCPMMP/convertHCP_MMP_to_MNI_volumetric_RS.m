function convertHCP_MMP_to_MNI_volumetric
% small script that will project HCM-MMP1.0 fsaverage annot files from here
% (https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446)
% onto MNI 152 2009a NLIN GM space.
% (c) 2016 Andreas Horn ? BIDMC Boston, MA, USA

% to run, you will need SPM8/12 and the freesurfer/matlab directory on your
% ML path. You need to put the following files into the working directory:
% - lh.pial and rh.pial from freesurfer/subjects/fsaverage/surf
% - mni_icbm152_gm_tal_nlin_asym_09a.nii from http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009
% - freesurfer converted annotation files (rh.HCP-MMP1.annot and
%   lh.HCP-MMP1.annot) from https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446
% - if you want to produce an HD version (0.5 0.5 0.5) you also need Lead
%   Neuroimaging Suite (leadsuite.io)
%
% EDIT 4/12/2020: Robert Seymour for implementing combined version,
% interpolated onto the SPM canonical brain

% Add Freesurfer_tools and SPM to your path
addpath('/Users/rseymoue/Documents/GitHub/OPMsurfer/freesurfer_tools');
addpath('/Users/rseymoue/Documents/scripts/spm12')


maxDist=5.5; % maximum distance (mm) for projection ? caution, if too large, subcortical structures will also be (falsely) assigned.
visualize=1; % set to 0 if you don't need to visualize
reslice=0; % set to 1 if you want to produce an 0.5 mm version, too.

% Load in table for combined ROIs
combined_table = readtable('MMP_groups.csv');

% load in LH
[~, lhlab,lhctable]=read_annotation('lh.HCP-MMP1.annot');
[lhvtx,lhfaces]=read_surf('lh.pial');

% Delete first superfluous row
lhctable.table(1,:) = [];
lhctable.struct_names(1,:) = [];

[~,lhix]=ismember(lhlab,lhctable.table(:,5));

% Remove random stuff from the names
for i = 1:180
    lhctable.struct_names{i,1} = lhctable.struct_names{i,1}(3:end-4);
end

% Replace numbers (1-180) with 1-22 based on region
lhix_out = lhix;

for i = 1:180
    %[~,indx]        = ismember(lhix,i);
    % Remove hippocampal region
    if strcmp('H',lhctable.struct_names{i})
        disp('Removing hippocampal ROI');
        lhix_out(lhix_out==i)  = 0;
    elseif strcmp('V2',lhctable.struct_names{i})
        disp('Removing V2 ROI')
        lhix_out(lhix_out==i)  = 0;
    elseif strcmp('PHA1',lhctable.struct_names{i})
        disp('Removing Parahippocampal ROI')
        lhix_out(lhix_out==i)  = 0;
    elseif strcmp('PHA2',lhctable.struct_names{i})
        lhix_out(lhix_out==i)  = 0;
    elseif strcmp('PHA3',lhctable.struct_names{i})
        lhix_out(lhix_out==i)  = 0;
    else
        r = find(strcmp(combined_table.ROI,lhctable.struct_names{i}));
        if isempty(r)
            warning([atlas_table.struct_names{i} ' IS WRONG']);
            
        else
            value = combined_table.value(r);
            lhix_out(lhix==i)  = value;
        end
    end
end

lhix = lhix_out;


% load in RH
[~, rhlab,rhctable]=read_annotation('rh.HCP-MMP1.annot');
[rhvtx,rhfaces]=read_surf('rh.pial');

% Delete first superfluous row
rhctable.table(1,:) = [];
rhctable.struct_names(1,:) = [];
[~,rhix]=ismember(rhlab,rhctable.table(:,5));

% Remove random stuff from the names
for i = 1:180
    rhctable.struct_names{i,1} = rhctable.struct_names{i,1}(3:end-4);
end

% For RH
for i = 1:length(rhix)
    if rhix(i) > 0
        rhix(i) = rhix(i)+180;
    end
end

% Replace numbers (1-180) with 23-44 based on region
rhix_out = rhix;

for i = 1:180
    %[~,indx]        = ismember(lhix,i);
    if strcmp('H',lhctable.struct_names{i})
        disp('Removing hippocampal ROI RH')
        rhix_out(rhix_out==i+180)  = 0;
    elseif strcmp('V2',lhctable.struct_names{i})
        disp('Removing V2 ROI RH')
        rhix_out(rhix_out==i+180)  = 0;
    elseif strcmp('PHA1',lhctable.struct_names{i})
        disp('Removing PHA1-3 ROIs RH')
        rhix_out(rhix_out==i+180)  = 0;
    elseif strcmp('PHA2',lhctable.struct_names{i})
        rhix_out(rhix_out==i+180)  = 0;
    elseif strcmp('PHA3',lhctable.struct_names{i})
        rhix_out(rhix_out==i+180)  = 0;
    else
        r = find(strcmp(combined_table.ROI,rhctable.struct_names{i}));
        if isempty(r)
            warning([atlas_table.struct_names{i} ' IS WRONG']);
        else
            value = combined_table.value(r);
            rhix_out(rhix==i+180)  = value+22;
        end
    end
end

rhix = rhix_out;




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

nii=spm_vol('/Users/rseymoue/Documents/scripts/fieldtrip-20191213/template/anatomy/single_subj_T1.nii');
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

nii.fname=['HCP-MMP1_combined_on_spm_brain',suffx,'.nii'];
%nii.dt=[16,0];
spm_write_vol(nii,nii.img);
gzip(['HCP-MMP1_combined_on_spm_brain',suffx,'.nii']);
%delete(['HCP-MMP1_on_MNI152_ICBM2009a_nlin',suffx,'.nii']);

% visualize volumetric data
if visualize
    subplot(1,2,2)
    title('Volumetric data');
    figure;
    for reg=1:44
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

% Export text file with names of regions
hemi            = {'LH','RH'}
names_for_export = unique(combined_table.Combined_ROI,'stable');

names_for_export2 = [];
count              = 1

for h = 1:length(hemi)
    for i = 1:22
        str = regexprep(names_for_export{i}, ' ', '_');

    names_for_export2{count} = [hemi{h} '_' str]
    count = count+1;
    end
end
        
% Finally, export .txt label file
f=fopen('HCP-MMP1_combined_on_spm_brain.txt','w');


for reg=1:44
    fprintf(f,'%s\n',names_for_export2{reg});
end
fclose(f);


