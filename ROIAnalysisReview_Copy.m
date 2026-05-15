% Script summary
%
% This script extracts electric field (E-field) values from a SimNIBS FEM
% simulation and computes summary statistics within specific anatomical
% regions defined by a FreeSurfer aparc_DKTatlas+aseg segmentation.
%
% Processing steps:
% 1. Load the SimNIBS head mesh containing tetrahedral elements and E-field data.
% 2. Load the FreeSurfer segmentation volume (aparc_DKTatlas+aseg).
% 3. Convert tetrahedron centroid coordinates from mesh RAS space to
%    voxel coordinates of the segmentation volume.
% 4. Assign an anatomical label to each tetrahedron based on the voxel
%    containing its centroid.
% 5. Extract electric field magnitude values for predefined ROIs
%    (hippocampus, amygdala, and frontal/cingulate cortical regions).
% 6. Compute mean E-field magnitude within each ROI.
% 7. Linearly scale the E-field values to approximate stimulation at
%    800 mA (simulation originally run at 10 mA).
%
% Output:
% ROI structure containing mean E-field values for each anatomical region.
% author : z5171263

%% 


%Dependencies and paths 
addpath('C:\Users\z5171263\SimNIBS-4.5\simnibs_env\Lib\site-packages\simnibs\matlab_tools');
addpath('C:\fieldtrip-20240113\external\freesurfer');

%  Load FreeSurfer segmentation volume
aseg_file = 'C:\Users\z5171263\Downloads\MRI_files\Donel\aparc_DKTatlast+aseg.nii';
vol  = niftiread(aseg_file);
info = niftiinfo(aseg_file);

%  Load SimNIBS FEM head mesh (.msh)
head_mesh = mesh_load_gmsh4( ...
    fullfile('C:\Users\z5171263\OneDrive - UNSW\Desktop\ECT modelling\compressed NIFTI\simnibs_simulation\donel_TDCS_1_scalar.msh'));

%  Construct coordinate transforms
%  NOTE: Transpose validated (0% out-of-bounds vs ~80% without)
vox2ras = info.Transform.T';
ras2vox = inv(vox2ras);


%  Compute tetrahedral element centres
tets  = head_mesh.tetrahedra;    
nodes = head_mesh.nodes;         

tet_centers = ( ...
    nodes(tets(:,1),:) + ...
    nodes(tets(:,2),:) + ...
    nodes(tets(:,3),:) + ...
    nodes(tets(:,4),:) ) / 4;

%  Convert tetrahedron centres to voxel coordinates
tet_h = [tet_centers ones(size(tet_centers,1),1)]';
vox_coords = ras2vox * tet_h;
vox_coords = round(vox_coords(1:3,:)');


%  Assign FreeSurfer atlas label to each tetrahedron
linear_idx = sub2ind(size(vol), ...
    vox_coords(:,1), ...
    vox_coords(:,2), ...
    vox_coords(:,3));
tet_labels = vol(linear_idx);

%  Confirm number of tetrahedron labels matches E-field elements
length(tet_labels) == length(head_mesh.element_data{2,1}.tetdata);

%  Extract electric field magnitude from mesh
magnE = head_mesh.element_data{2,1}.tetdata;

%  ROI label list used for analysis
unique_labels = unique(tet_labels);

roi_labels = [ ...
    17,18,53,54, ...
    1002,1003,1012,1014,1018,1019,1020,1023, ...
    1026,1027,1028, ...
    2003,2012,2018,2019,2020,2026,2027 ...
];

%  ROI Definitions (FreeSurfer aseg + DKT atlas)

ROI = struct();

ROI(1).id = 17; ROI(1).name = 'LH_Hippocampus';
ROI(2).id = 18; ROI(2).name = 'LH_Amygdala';
ROI(3).id = 53; ROI(3).name = 'RH_Hippocampus';
ROI(4).id = 54; ROI(4).name = 'RH_Amygdala';

%% Left Hemisphere Cortical ROIs 

ROI(5).id = 1002; ROI(5).name = 'LH_CaudalAnteriorCingulate';
ROI(6).id = 1003; ROI(6).name = 'LH_CaudalMiddleFrontal';
ROI(7).id = 1012; ROI(7).name = 'LH_LateralOrbitofrontal';
ROI(8).id = 1014; ROI(8).name = 'LH_MedialOrbitofrontal';
ROI(9).id = 1018; ROI(9).name = 'LH_ParsOpercularis';
ROI(10).id = 1019; ROI(10).name = 'LH_ParsOrbitalis';
ROI(11).id = 1020; ROI(11).name = 'LH_ParsTriangularis';
ROI(12).id = 1023; ROI(12).name = 'LH_PosteriorCingulate';
ROI(13).id = 1026; ROI(13).name = 'LH_RostralAnteriorCingulate';
ROI(14).id = 1027; ROI(14).name = 'LH_RostralMiddleFrontal';
ROI(15).id = 1028; ROI(15).name = 'LH_SuperiorFrontal';

% Right Hemisphere Cortical ROIs

ROI(16).id = 2003; ROI(16).name = 'RH_CaudalMiddleFrontal';
ROI(17).id = 2012; ROI(17).name = 'RH_LateralOrbitofrontal';
ROI(18).id = 2018; ROI(18).name = 'RH_ParsOpercularis';
ROI(19).id = 2019; ROI(19).name = 'RH_ParsOrbitalis';
ROI(20).id = 2020; ROI(20).name = 'RH_ParsTriangularis';
ROI(21).id = 2026; ROI(21).name = 'RH_RostralAnteriorCingulate';
ROI(22).id = 2027; ROI(22).name = 'RH_RostralMiddleFrontal';
ROI(23).id = 2028; ROI(23).name = 'RH_SuperiorFrontal';

%  Extract electric field statistics for each ROI
for i = 1:length(ROI)

    mask = tet_labels == ROI(i).id;
    Evals = magnE(mask);

    ROI(i).mean = mean(Evals) * 80; % scaled to ~800 mA
end


%% 


roi_ids = [18, 1035]; % hippocampus, amygdala, frontal

colors = lines(length(roi_ids));

for i = 1:length(roi_ids)
    mask = tet_labels == roi_ids(i);
    
    scatter3(tet_centers(mask,1), ...
             tet_centers(mask,2), ...
             tet_centers(mask,3), ...
             10, colors(i,:), 'filled');
end

legend('All tets','ROI1','ROI2','ROI3')


gm_mask = (tet_labels >= 1000 & tet_labels < 3000);

gm_idx = find(gm_mask);
gm_idx = gm_idx(randperm(length(gm_idx), round(0.1 * length(gm_idx)))); 

figure
hold on
axis equal
view(3)

title('ROIs on Grey Matter')

% Grey matter background
scatter3(tet_centers(gm_idx,1), ...
         tet_centers(gm_idx,2), ...
         tet_centers(gm_idx,3), ...
         2, [0.7 0.7 0.7], 'filled');

roi_id = 18; 
mask = tet_labels == roi_id;

scatter3(tet_centers(mask,1), ...
         tet_centers(mask,2), ...
         tet_centers(mask,3), ...
         10, 'r', 'filled');

legend('Grey Matter','ROI')