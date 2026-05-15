addpath("C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools");
 
%% Paths - single subject test
assemblynet_file = 'C:\Users\z5171263\Downloads\testing_volbrain\AN001_volbrain\5001 SAG FL3D_VIBE_job1980923_archive\native_structures_job1980923.nii.gz';
mesh_file        = 'C:\Users\z5171263\Downloads\testing_volbrain\edited_AN001_FrontoParietal_simulation_results\edited_AN001_TDCS_1_scalar.msh';
subj_id          = 'AN001';
 
%% Load AssemblyNet structures
vol_ab  = niftiread(assemblynet_file);
info_ab = niftiinfo(assemblynet_file);
vox2ras_ab = info_ab.Transform.T';
ras2vox_ab = inv(vox2ras_ab);
 
%% Load simulation mesh
head_mesh = mesh_load_gmsh4(mesh_file);
 
%% Compute tet centres
tets  = head_mesh.tetrahedra;
nodes = head_mesh.nodes;
tet_centers = (nodes(tets(:,1),:) + nodes(tets(:,2),:) + nodes(tets(:,3),:) + nodes(tets(:,4),:)) / 4;
tet_h = [tet_centers ones(size(tet_centers,1),1)]';
 
%% Map tet centres to AssemblyNet voxel space
vox_coords_ab = ras2vox_ab * tet_h;
vox_coords_ab = round(vox_coords_ab(1:3,:)');
 
inside_mask_ab = vox_coords_ab(:,1) >= 1 & vox_coords_ab(:,1) <= size(vol_ab,1) & ...
                 vox_coords_ab(:,2) >= 1 & vox_coords_ab(:,2) <= size(vol_ab,2) & ...
                 vox_coords_ab(:,3) >= 1 & vox_coords_ab(:,3) <= size(vol_ab,3);
 
n_out_ab = sum(~inside_mask_ab);
if n_out_ab > 0
    fprintf('Out-of-bounds tetrahedra (AssemblyNet): %d (%.4f%%)\n', n_out_ab, 100*n_out_ab/length(inside_mask_ab));
end
 
vox_coords_ab = vox_coords_ab(inside_mask_ab,:);
magnE_ab = head_mesh.element_data{2,1}.tetdata;
magnE_ab = magnE_ab(inside_mask_ab);
linear_idx_ab = sub2ind(size(vol_ab), vox_coords_ab(:,1), vox_coords_ab(:,2), vox_coords_ab(:,3));
tet_labels_ab = vol_ab(linear_idx_ab);
 
%% ROI definitions
ROI_ab = struct();
ROI_ab(1).ids  = [48];                 ROI_ab(1).name  = 'LH_Hippocampus';
ROI_ab(2).ids  = [47];                 ROI_ab(2).name  = 'RH_Hippocampus';
ROI_ab(3).ids  = [32];                 ROI_ab(3).name  = 'LH_Amygdala';
ROI_ab(4).ids  = [31];                 ROI_ab(4).name  = 'RH_Amygdala';
ROI_ab(5).ids  = [167];                ROI_ab(5).name  = 'LH_CaudalAnteriorCingulate_Ctx';
ROI_ab(6).ids  = [166];                ROI_ab(6).name  = 'RH_CaudalAnteriorCingulate_Ctx';
ROI_ab(7).ids  = [143];                ROI_ab(7).name  = 'LH_CaudalMiddleFrontal_Ctx';
ROI_ab(8).ids  = [142];                ROI_ab(8).name  = 'RH_CaudalMiddleFrontal_Ctx';
ROI_ab(9).ids  = [137];                ROI_ab(9).name  = 'LH_LateralOrbitofrontal_Ctx';
ROI_ab(10).ids = [136];                ROI_ab(10).name = 'RH_LateralOrbitofrontal_Ctx';
ROI_ab(11).ids = [147];                ROI_ab(11).name = 'LH_MedialOrbitofrontal_Ctx';
ROI_ab(12).ids = [146];                ROI_ab(12).name = 'RH_MedialOrbitofrontal_Ctx';
ROI_ab(13).ids = [163];                ROI_ab(13).name = 'LH_ParsOpercularis';
ROI_ab(14).ids = [162];                ROI_ab(14).name = 'RH_ParsOpercularis';
ROI_ab(15).ids = [165];                ROI_ab(15).name = 'LH_ParsOrbitalis';
ROI_ab(16).ids = [164];                ROI_ab(16).name = 'RH_ParsOrbitalis';
ROI_ab(17).ids = [205];                ROI_ab(17).name = 'LH_ParsTriangularis';
ROI_ab(18).ids = [204];                ROI_ab(18).name = 'RH_ParsTriangularis';
ROI_ab(19).ids = [167];                ROI_ab(19).name = 'LH_PosteriorCingulate_Ctx';
ROI_ab(20).ids = [166];                ROI_ab(20).name = 'RH_PosteriorCingulate_Ctx';
ROI_ab(21).ids = [101];                ROI_ab(21).name = 'LH_RostralAnteriorCingulate_Ctx';
ROI_ab(22).ids = [100];                ROI_ab(22).name = 'RH_RostralAnteriorCingulate_Ctx';
ROI_ab(23).ids = [191];                ROI_ab(23).name = 'LH_SuperiorFrontal_Ctx';
ROI_ab(24).ids = [190];                ROI_ab(24).name = 'RH_SuperiorFrontal_Ctx';
ROI_ab(25).ids = [187];                ROI_ab(25).name = 'LH_SubcallosalCingulate_Gyrus';
ROI_ab(26).ids = [186];                ROI_ab(26).name = 'RH_SubcallosalCingulate_Gyrus';
ROI_ab(27).ids = [183];                ROI_ab(27).name = 'LH_Precentral_Gyrus';
ROI_ab(28).ids = [182];                ROI_ab(28).name = 'RH_Precentral_Gyrus';
ROI_ab(29).ids = [105, 137, 147, 179]; ROI_ab(29).name = 'LH_Orbital_Gyrus';
ROI_ab(30).ids = [104, 136, 146, 178]; ROI_ab(30).name = 'RH_Orbital_Gyrus';
ROI_ab(31).ids = [143];                ROI_ab(31).name = 'LH_MiddleFrontal_Gyrus';
ROI_ab(32).ids = [142];                ROI_ab(32).name = 'RH_MiddleFrontal_Gyrus';
 
%% Extract mean magE per ROI
for i = 1:length(ROI_ab)
    mask = ismember(tet_labels_ab, ROI_ab(i).ids);
    Evals = magnE_ab(mask);
    ROI_ab(i).mean = mean(Evals) * 800;
end
 
%% Build and display summary table
roi_means = [];
roi_names = {};
for i = 1:length(ROI_ab)
    roi_means(end+1) = ROI_ab(i).mean;
    roi_names{end+1} = ROI_ab(i).name;
end
 
summary_table = array2table(roi_means, 'VariableNames', roi_names, 'RowNames', {subj_id});
disp(summary_table);
 
writetable(summary_table, fullfile(fileparts(mesh_file), 'assemblynet_FP_sim_ROI_summary.xls'), 'WriteRowNames', true);