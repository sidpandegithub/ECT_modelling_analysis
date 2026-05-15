addpath("C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools");

%% Paths
nii_folder   = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\albert_road';
sim_folder   = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\albert_road_TP_simulation';

% Get list of all .nii files
nii_files = dir(fullfile(nii_folder, '*_aparc.DKTatlas+aseg.nii'));

%% Initialize summary table
roi_ids = [ ...
    17,18,53,54, ...
    1002,1003,1012,1014,1018,1019,1020,1023,1026,1027,1028, ...
    2002,2003,2012,2014,2018,2019,2020,2023,2026,2027,2028 ...
];

summary_table = table(); % will store mean E-field per ROI per subject

% Loop over subjects
for s = 1:length(nii_files)
    
    % Subject ID parsing 
    [~, subj_name, ~] = fileparts(nii_files(s).name); % ex., AL_001_aparc.DKTatlas+aseg
    subj_id = extractBefore(subj_name, '_aparc');     % ex., AL_001
    
    fprintf('Processing subject: %s\n', subj_id);
    
    %Load FreeSurfer parcellation
    aseg_file = fullfile(nii_folder, nii_files(s).name);
    vol  = niftiread(aseg_file);
    info = niftiinfo(aseg_file);

    % Locate corresponding simulation folder and scalar mesh
    sim_subj_folder = fullfile(sim_folder, [subj_id '_TemporoParietal_simulation_results']);
    mesh_file = fullfile(sim_subj_folder, [subj_id '_TDCS_1_scalar.msh']);

    if ~isfile(mesh_file)
        warning('Mesh file not found for %s, skipping...', subj_id);
        continue
    end
    head_mesh = mesh_load_gmsh4(mesh_file);

    % Construct coordinate transforms
    %Transpose validated (0% out-of-bounds vs ~80% without)
    vox2ras = info.Transform.T';
    ras2vox = inv(vox2ras);

    % Compute tetrahedral element centres
    tets  = head_mesh.tetrahedra;    
    nodes = head_mesh.nodes;         
    tet_centers = ( ...
        nodes(tets(:,1),:) + ...
        nodes(tets(:,2),:) + ...
        nodes(tets(:,3),:) + ...
        nodes(tets(:,4),:) ) / 4;

    % Convert tetrahedron centres to voxel coordinates
    tet_h = [tet_centers ones(size(tet_centers,1),1)]';
    vox_coords = ras2vox * tet_h;
    vox_coords = round(vox_coords(1:3,:)');

       % --- REMOVE out-of-bounds tetrahedra instead of clipping ---
    inside_mask = ...
        vox_coords(:,1) >= 1 & vox_coords(:,1) <= size(vol,1) & ...
        vox_coords(:,2) >= 1 & vox_coords(:,2) <= size(vol,2) & ...
        vox_coords(:,3) >= 1 & vox_coords(:,3) <= size(vol,3);
    
    % (optional debug print)
    n_out = sum(~inside_mask);
    if n_out > 0
        fprintf('Out-of-bounds tetrahedra: %d (%.4f%%)\n', ...
            n_out, 100*n_out/length(inside_mask));
    end

    vox_coords = vox_coords(inside_mask,:);

    magnE = head_mesh.element_data{2,1}.tetdata;
    magnE = magnE(inside_mask);
    % Assign FreeSurfer atlas label to each tetrahedron
    linear_idx = sub2ind(size(vol), ...
        vox_coords(:,1), ...
        vox_coords(:,2), ...
        vox_coords(:,3));
    tet_labels = vol(linear_idx);

    % Confirm number of tetrahedron labels matches E-field elements
    length(tet_labels) == length(head_mesh.element_data{2,1}.tetdata);

    % Extract electric field magnitude from mesh
    magnE = head_mesh.element_data{2,1}.tetdata;

    % ROI label list used for analysis
    unique_labels = unique(tet_labels);

    roi_labels = [ ...
        17,18,53,54, ...
        1002,1003,1012,1014,1018,1019,1020,1023, ...
        1026,1027,1028, ...
        2002,2003,2012,2014,2018,2019,2020,2023,2026,2027 ...
    ];

    % ROI Definitions (FreeSurfer aseg + DKT atlas)
    ROI = struct();
    
    % Subcortical
    ROI(1).id = 17; ROI(1).name = 'LH_Hippocampus';
    ROI(2).id = 18; ROI(2).name = 'LH_Amygdala';
    ROI(3).id = 53; ROI(3).name = 'RH_Hippocampus';
    ROI(4).id = 54; ROI(4).name = 'RH_Amygdala';
    
    % Left Hemisphere Cortical ROIs
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
    ROI(16).id = 2002; ROI(16).name = 'RH_CaudalAnteriorCingulate';
    ROI(17).id = 2003; ROI(17).name = 'RH_CaudalMiddleFrontal';
    ROI(18).id = 2012; ROI(18).name = 'RH_LateralOrbitofrontal';
    ROI(19).id = 2014; ROI(19).name = 'RH_MedialOrbitofrontal';
    ROI(20).id = 2018; ROI(20).name = 'RH_ParsOpercularis';
    ROI(21).id = 2019; ROI(21).name = 'RH_ParsOrbitalis';
    ROI(22).id = 2020; ROI(22).name = 'RH_ParsTriangularis';
    ROI(23).id = 2023; ROI(23).name = 'RH_PosteriorCingulate';
    ROI(24).id = 2026; ROI(24).name = 'RH_RostralAnteriorCingulate';
    ROI(25).id = 2027; ROI(25).name = 'RH_RostralMiddleFrontal';
    ROI(26).id = 2028; ROI(26).name = 'RH_SuperiorFrontal';

    % Extract electric field statistics for each ROI
    for i = 1:length(ROI)
        mask = tet_labels == ROI(i).id;
        Evals = magnE(mask);
        ROI(i).mean = mean(Evals) * 800; % scaled to ~800 mA
    end

    % Append ROI means to summary table
    roi_means = zeros(1, length(roi_ids));
    roi_names = cell(1, length(roi_ids));  % store names for columns
    
    for k = 1:length(roi_ids)
        idx = find([ROI.id] == roi_ids(k), 1);
        if ~isempty(idx)
            roi_means(k) = ROI(idx).mean;
            roi_names{k} = ROI(idx).name; 
        else
            roi_means(k) = NaN; 
            roi_names{k} = ['ROI_' num2str(roi_ids(k))]; 
        end
    end
    
    summary_table = [summary_table; array2table(roi_means, ...
        'VariableNames', roi_names, 'RowNames', {subj_id})];
end

%% Display the summary table
disp(summary_table);

% save as excel file
writetable(summary_table, fullfile(sim_folder, 'albert_raod_TP_sim_ROI_summary.xls'), 'WriteRowNames', true);