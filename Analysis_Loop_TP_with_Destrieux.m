addpath("C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools");

%% Paths
nii_folder   = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\ramsay_lakeside'; %DKT altas folder
sim_folder   = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\ramsay_lakeside_TP_simulations';
Destrieux_nii = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\ramsay_lakeside_Destrieux';

nii_files = dir(fullfile(nii_folder, '*_aparc.DKTatlas+aseg.nii'));

summary_table = table();

%% ROI IDs (DKT + Destrieux)
roi_ids = [ ...
    17,18,53,54, ...
    1002,1003,1012,1014,1018,1019,1020,1023,1026,1027,1028, ...
    2002,2003,2012,2014,2018,2019,2020,2023,2026,2027,2028, ...
    11132,12132,11129,12129,11124,12124,11115,12115 ...
];

% Loop over subjects
for s = 1:length(nii_files)

    [~, subj_name, ~] = fileparts(nii_files(s).name);
    subj_id = extractBefore(subj_name, '_aparc');

    fprintf('Processing subject: %s\n', subj_id);

    % Load DKT
    aseg_file = fullfile(nii_folder, nii_files(s).name);
    vol  = niftiread(aseg_file);
    info = niftiinfo(aseg_file);

    % Load Destrieux
    destrieux_file = fullfile(Destrieux_nii, [subj_id '_aparc.a2009s+aseg.nii']);
    use_destrieux = isfile(destrieux_file);

    if use_destrieux
        vol_des  = niftiread(destrieux_file);
        info_des = niftiinfo(destrieux_file);
        vox2ras_des = info_des.Transform.T';
        ras2vox_des = inv(vox2ras_des);
    else
        warning('Destrieux file not found for %s', subj_id);
    end

    % Load simulation mesh
    sim_subj_folder = fullfile(sim_folder, [subj_id '_TemporoParietal_simulation_results']);
    mesh_file = fullfile(sim_subj_folder, [subj_id '_TDCS_1_scalar.msh']);

    if ~isfile(mesh_file)
        warning('Mesh file not found for %s, skipping...', subj_id);
        continue
    end

    head_mesh = mesh_load_gmsh4(mesh_file);

    % Transforms
    vox2ras = info.Transform.T';
    ras2vox = inv(vox2ras);

    tets  = head_mesh.tetrahedra;
    nodes = head_mesh.nodes;

    tet_centers = (nodes(tets(:,1),:) + nodes(tets(:,2),:) + nodes(tets(:,3),:) + nodes(tets(:,4),:)) / 4;
    tet_h = [tet_centers ones(size(tet_centers,1),1)]';

    % Map to DKT
    vox_coords = ras2vox * tet_h;
    vox_coords = round(vox_coords(1:3,:)');

    inside_mask = vox_coords(:,1) >= 1 & vox_coords(:,1) <= size(vol,1) & ...
                  vox_coords(:,2) >= 1 & vox_coords(:,2) <= size(vol,2) & ...
                  vox_coords(:,3) >= 1 & vox_coords(:,3) <= size(vol,3);
    n_out = sum(~inside_mask);
    if n_out > 0
        fprintf('Out-of-bounds tetrahedra (DKT): %d (%.4f%%)\n', n_out, 100*n_out/length(inside_mask));
    end

    vox_coords = vox_coords(inside_mask,:);
    magnE = head_mesh.element_data{2,1}.tetdata;
    magnE = magnE(inside_mask);
    linear_idx = sub2ind(size(vol), vox_coords(:,1), vox_coords(:,2), vox_coords(:,3));
    tet_labels = vol(linear_idx);

    % Map to Destrieux
    if use_destrieux
        vox_coords_des = ras2vox_des * tet_h;
        vox_coords_des = round(vox_coords_des(1:3,:)');

        inside_mask_des = vox_coords_des(:,1) >= 1 & vox_coords_des(:,1) <= size(vol_des,1) & ...
                          vox_coords_des(:,2) >= 1 & vox_coords_des(:,2) <= size(vol_des,2) & ...
                          vox_coords_des(:,3) >= 1 & vox_coords_des(:,3) <= size(vol_des,3);
        n_out_des = sum(~inside_mask_des);
        if n_out_des > 0
            fprintf('Out-of-bounds tetrahedra (Destrieux): %d (%.4f%%)\n', n_out_des, 100*n_out_des/length(inside_mask_des));
        end

        vox_coords_des = vox_coords_des(inside_mask_des,:);
        magnE_des = head_mesh.element_data{2,1}.tetdata;
        magnE_des = magnE_des(inside_mask_des);
        linear_idx_des = sub2ind(size(vol_des), vox_coords_des(:,1), vox_coords_des(:,2), vox_coords_des(:,3));
        tet_labels_des = vol_des(linear_idx_des);
    end

    % ROI Definitions
    ROI = struct();

    % Subcortical
    ROI(1).id = 17; ROI(1).name = 'LH_Hippocampus';
    ROI(2).id = 18; ROI(2).name = 'LH_Amygdala';
    ROI(3).id = 53; ROI(3).name = 'RH_Hippocampus';
    ROI(4).id = 54; ROI(4).name = 'RH_Amygdala';
    
    % LH Cortex
    ROI(5).id = 1002; ROI(5).name = 'LH_CaudalAnteriorCingulate_Ctx';
    ROI(6).id = 1003; ROI(6).name = 'LH_CaudalMiddleFrontal_Ctx';
    ROI(7).id = 1005; ROI(7).name = 'LH_Cuneus_Ctx';
    ROI(8).id = 1006; ROI(8).name = 'LH_Entorhinal_Ctx';
    ROI(9).id = 1007; ROI(9).name = 'LH_Fusiform_Ctx';
    ROI(10).id = 1008; ROI(10).name = 'LH_InferiorParietal_Ctx';
    ROI(11).id = 1009; ROI(11).name = 'LH_InferiorTemporal_Ctx';
    ROI(12).id = 1010; ROI(12).name = 'LH_IsthmusCingulate_Ctx';
    ROI(13).id = 1011; ROI(13).name = 'LH_LateralOccipital_Ctx';
    ROI(14).id = 1012; ROI(14).name = 'LH_LateralOrbitofrontal_Ctx';
    ROI(15).id = 1013; ROI(15).name = 'LH_Lingual_Ctx';
    ROI(16).id = 1014; ROI(16).name = 'LH_MedialOrbitofrontal_Ctx';
    ROI(17).id = 1015; ROI(17).name = 'LH_MiddleTemporal_Ctx';
    ROI(18).id = 1016; ROI(18).name = 'LH_Parahippocampal_Ctx';
    ROI(19).id = 1017; ROI(19).name = 'LH_Paracentral_Ctx';
    ROI(20).id = 1018; ROI(20).name = 'LH_ParsOpercularis';
    ROI(21).id = 1019; ROI(21).name = 'LH_ParsOrbitalis';
    ROI(22).id = 1020; ROI(22).name = 'LH_ParsTriangularis';
    ROI(23).id = 1021; ROI(23).name = 'LH_Pericalcarine_Ctx';
    ROI(24).id = 1022; ROI(24).name = 'LH_Postcentral_Ctx';
    ROI(25).id = 1023; ROI(25).name = 'LH_PosteriorCingulate_Ctx';
    ROI(26).id = 1024; ROI(26).name = 'LH_Precentral_Ctx';
    ROI(27).id = 1025; ROI(27).name = 'LH_Precuneus_Ctx';
    ROI(28).id = 1026; ROI(28).name = 'LH_RostralAnteriorCingulate_Ctx';
    ROI(29).id = 1027; ROI(29).name = 'LH_RostralMiddleFrontal_Ctx';
    ROI(30).id = 1028; ROI(30).name = 'LH_SuperiorFrontal_Ctx';
    ROI(31).id = 1029; ROI(31).name = 'LH_SuperiorParietal_Ctx';
    ROI(32).id = 1030; ROI(32).name = 'LH_SuperiorTemporal_Ctx';
    ROI(33).id = 1031; ROI(33).name = 'LH_Supramarginal_Ctx';
    ROI(34).id = 1034; ROI(34).name = 'LH_TransverseTemporal_Ctx';
    ROI(35).id = 1035; ROI(35).name = 'LH_Insula_Ctx';
    
    % RH Cortex
    ROI(36).id = 2002; ROI(36).name = 'RH_CaudalAnteriorCingulate_Ctx';
    ROI(37).id = 2003; ROI(37).name = 'RH_CaudalMiddleFrontal_Ctx';
    ROI(38).id = 2005; ROI(38).name = 'RH_Cuneus_Ctx';
    ROI(39).id = 2006; ROI(39).name = 'RH_Entorhinal_Ctx';
    ROI(40).id = 2007; ROI(40).name = 'RH_Fusiform_Ctx';
    ROI(41).id = 2008; ROI(41).name = 'RH_InferiorParietal_Ctx';
    ROI(42).id = 2009; ROI(42).name = 'RH_InferiorTemporal_Ctx';
    ROI(43).id = 2010; ROI(43).name = 'RH_IsthmusCingulate_Ctx';
    ROI(44).id = 2011; ROI(44).name = 'RH_LateralOccipital_Ctx';
    ROI(45).id = 2012; ROI(45).name = 'RH_LateralOrbitofrontal_Ctx';
    ROI(46).id = 2013; ROI(46).name = 'RH_Lingual_Ctx';
    ROI(47).id = 2014; ROI(47).name = 'RH_MedialOrbitofrontal_Ctx';
    ROI(48).id = 2015; ROI(48).name = 'RH_MiddleTemporal_Ctx';
    ROI(49).id = 2016; ROI(49).name = 'RH_Parahippocampal_Ctx';
    ROI(50).id = 2017; ROI(50).name = 'RH_Paracentral_Ctx';
    ROI(51).id = 2018; ROI(51).name = 'RH_ParsOpercularis';
    ROI(52).id = 2019; ROI(52).name = 'RH_ParsOrbitalis';
    ROI(53).id = 2020; ROI(53).name = 'RH_ParsTriangularis';
    ROI(54).id = 2021; ROI(54).name = 'RH_Pericalcarine_Ctx';
    ROI(55).id = 2022; ROI(55).name = 'RH_Postcentral_Ctx';
    ROI(56).id = 2023; ROI(56).name = 'RH_PosteriorCingulate_Ctx';
    ROI(57).id = 2024; ROI(57).name = 'RH_Precentral_Ctx';
    ROI(58).id = 2025; ROI(58).name = 'RH_Precuneus_Ctx';
    ROI(59).id = 2026; ROI(59).name = 'RH_RostralAnteriorCingulate_Ctx';
    ROI(60).id = 2027; ROI(60).name = 'RH_RostralMiddleFrontal_Ctx';
    ROI(61).id = 2028; ROI(61).name = 'RH_SuperiorFrontal_Ctx';
    ROI(62).id = 2029; ROI(62).name = 'RH_SuperiorParietal_Ctx';
    ROI(63).id = 2030; ROI(63).name = 'RH_SuperiorTemporal_Ctx';
    ROI(64).id = 2031; ROI(64).name = 'RH_Supramarginal_Ctx';
    ROI(65).id = 2034; ROI(65).name = 'RH_TransverseTemporal_Ctx';
    ROI(66).id = 2035; ROI(66).name = 'RH_Insula_Ctx';



    ROI_des = struct();
    ROI_des(1).id = 11132; ROI_des(1).name = 'LH_SubcallosalCingulate_Gyrus';
    ROI_des(2).id = 12132; ROI_des(2).name = 'RH_SubcallosalCingulate_Gyrus';
    ROI_des(3).id = 11129; ROI_des(3).name = 'LH_Precentral_Gyrus';
    ROI_des(4).id = 12129; ROI_des(4).name = 'RH_Precentral_Gyrus';
    ROI_des(5).id = 11124; ROI_des(5).name = 'LH_Orbital_Gyrus';
    ROI_des(6).id = 12124; ROI_des(6).name = 'RH_Orbital_Gyrus';
    ROI_des(7).id = 11115; ROI_des(7).name = 'LH_MiddleFrontal_Gyrus';
    ROI_des(8).id = 12115; ROI_des(8).name = 'RH_MiddleFrontal_Gyrus';

    for i = 1:length(ROI)
        mask = tet_labels == ROI(i).id;
        Evals = magnE(mask);
        ROI(i).mean = mean(Evals) * 800;
    end

    if use_destrieux
        for i = 1:length(ROI_des)
            mask = tet_labels_des == ROI_des(i).id;
            Evals = magnE_des(mask);
            ROI_des(i).mean = mean(Evals) * 800;
        end
    end

    roi_means = [];
    roi_names = {};

    for i = 1:length(ROI)
        roi_means(end+1) = ROI(i).mean;
        roi_names{end+1} = ROI(i).name;
    end

    for i = 1:length(ROI_des)
        if use_destrieux
            roi_means(end+1) = ROI_des(i).mean;
        else
            roi_means(end+1) = NaN;
        end
        roi_names{end+1} = ROI_des(i).name;
    end

    summary_table = [summary_table; array2table(roi_means, 'VariableNames', roi_names, 'RowNames', {subj_id})];
end

disp(summary_table);
%% 

writetable(summary_table, fullfile(sim_folder, 'Augusta_TP_sim_ROI_summary.xls'), 'WriteRowNames', true);