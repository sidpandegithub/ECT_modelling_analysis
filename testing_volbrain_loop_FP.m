addpath("C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools");
 
%% Paths
parcellation_root = 'C:\Users\z5171263\Downloads\VolBrain_Parcellations';
sim_root          = 'C:\Users\z5171263\Downloads\Northside_FP_sims';
 
summary_table = table();
 
%% Find all subject parcellation folders
subj_folders = dir(fullfile(parcellation_root, '*_job*_archive'));
subj_folders = subj_folders([subj_folders.isdir]);
 
for s = 1:length(subj_folders)
 
    folder_name = subj_folders(s).name;  % e.g. AN001_job1980923_archive
 
    % Extract subject ID (everything before _job)
    subj_id = extractBefore(folder_name, '_job');
    fprintf('Processing subject: %s\n', subj_id);
 
    % Find native_structures nii.gz inside archive folder (job number varies)
    archive_folder = fullfile(parcellation_root, folder_name);
    nii_search = dir(fullfile(archive_folder, 'native_structures_*.nii.gz'));
 
    if isempty(nii_search)
        warning('No native_structures file found for %s, skipping...', subj_id);
        continue
    end
 
    assemblynet_file = fullfile(archive_folder, nii_search(1).name);
 
    % Find simulation mesh
    mesh_search = dir(fullfile(sim_root, ['edited_' subj_id '_FrontoParietal_simulation_results'], ...
                      ['edited_' subj_id '_TDCS_1_scalar.msh']));
 
    if isempty(mesh_search)
        warning('Mesh file not found for %s, skipping...', subj_id);
        continue
    end
 
    mesh_file = fullfile(mesh_search(1).folder, mesh_search(1).name);
 
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
    % ROI_ab = struct();
    % ROI_ab(1).ids  = [48];                 ROI_ab(1).name  = 'LH_Hippocampus';
    % ROI_ab(2).ids  = [47];                 ROI_ab(2).name  = 'RH_Hippocampus';
    % ROI_ab(3).ids  = [32];                 ROI_ab(3).name  = 'LH_Amygdala';
    % ROI_ab(4).ids  = [31];                 ROI_ab(4).name  = 'RH_Amygdala';
    % ROI_ab(5).ids  = [167];                ROI_ab(5).name  = 'LH_CaudalAnteriorCingulate_Ctx';
    % ROI_ab(6).ids  = [166];                ROI_ab(6).name  = 'RH_CaudalAnteriorCingulate_Ctx';
    % ROI_ab(7).ids  = [143];                ROI_ab(7).name  = 'LH_CaudalMiddleFrontal_Ctx';
    % ROI_ab(8).ids  = [142];                ROI_ab(8).name  = 'RH_CaudalMiddleFrontal_Ctx';
    % ROI_ab(9).ids  = [137];                ROI_ab(9).name  = 'LH_LateralOrbitofrontal_Ctx';
    % ROI_ab(10).ids = [136];                ROI_ab(10).name = 'RH_LateralOrbitofrontal_Ctx';
    % ROI_ab(11).ids = [147];                ROI_ab(11).name = 'LH_MedialOrbitofrontal_Ctx';
    % ROI_ab(12).ids = [146];                ROI_ab(12).name = 'RH_MedialOrbitofrontal_Ctx';
    % ROI_ab(13).ids = [163];                ROI_ab(13).name = 'LH_ParsOpercularis';
    % ROI_ab(14).ids = [162];                ROI_ab(14).name = 'RH_ParsOpercularis';
    % ROI_ab(15).ids = [165];                ROI_ab(15).name = 'LH_ParsOrbitalis';
    % ROI_ab(16).ids = [164];                ROI_ab(16).name = 'RH_ParsOrbitalis';
    % ROI_ab(17).ids = [205];                ROI_ab(17).name = 'LH_ParsTriangularis';
    % ROI_ab(18).ids = [204];                ROI_ab(18).name = 'RH_ParsTriangularis';
    % ROI_ab(19).ids = [167];                ROI_ab(19).name = 'LH_PosteriorCingulate_Ctx';
    % ROI_ab(20).ids = [166];                ROI_ab(20).name = 'RH_PosteriorCingulate_Ctx';
    % ROI_ab(21).ids = [101];                ROI_ab(21).name = 'LH_RostralAnteriorCingulate_Ctx';
    % ROI_ab(22).ids = [100];                ROI_ab(22).name = 'RH_RostralAnteriorCingulate_Ctx';
    % ROI_ab(23).ids = [191];                ROI_ab(23).name = 'LH_SuperiorFrontal_Ctx';
    % ROI_ab(24).ids = [190];                ROI_ab(24).name = 'RH_SuperiorFrontal_Ctx';
    % ROI_ab(25).ids = [187];                ROI_ab(25).name = 'LH_SubcallosalCingulate_Gyrus';
    % ROI_ab(26).ids = [186];                ROI_ab(26).name = 'RH_SubcallosalCingulate_Gyrus';
    % ROI_ab(27).ids = [183];                ROI_ab(27).name = 'LH_Precentral_Gyrus';
    % ROI_ab(28).ids = [182];                ROI_ab(28).name = 'RH_Precentral_Gyrus';
    % ROI_ab(29).ids = [105, 137, 147, 179]; ROI_ab(29).name = 'LH_Orbital_Gyrus';
    % ROI_ab(30).ids = [104, 136, 146, 178]; ROI_ab(30).name = 'RH_Orbital_Gyrus';
    % ROI_ab(31).ids = [143];                ROI_ab(31).name = 'LH_MiddleFrontal_Gyrus';
    % ROI_ab(32).ids = [142];                ROI_ab(32).name = 'RH_MiddleFrontal_Gyrus';
 
    %% ROI definitions
    
    ROI_ab = struct();
    
    ROI_ab(1).ids  = [48]; ROI_ab(1).name  = 'LH_Hippocampus';
    ROI_ab(2).ids  = [47]; ROI_ab(2).name  = 'RH_Hippocampus';
    ROI_ab(3).ids  = [32]; ROI_ab(3).name  = 'LH_Amygdala';
    ROI_ab(4).ids  = [31]; ROI_ab(4).name  = 'RH_Amygdala';
    
    ROI_ab(5).ids  = [37]; ROI_ab(5).name  = 'LH_Caudate';
    ROI_ab(6).ids  = [36]; ROI_ab(6).name  = 'RH_Caudate';
    
    ROI_ab(7).ids  = [39]; ROI_ab(7).name  = 'LH_CerebellumExterior';
    ROI_ab(8).ids  = [38]; ROI_ab(8).name  = 'RH_CerebellumExterior';
    
    ROI_ab(9).ids  = [41]; ROI_ab(9).name  = 'LH_CerebellumWhiteMatter';
    ROI_ab(10).ids = [40]; ROI_ab(10).name = 'RH_CerebellumWhiteMatter';
    
    ROI_ab(11).ids = [45]; ROI_ab(11).name = 'LH_CerebralWhiteMatter';
    ROI_ab(12).ids = [44]; ROI_ab(12).name = 'RH_CerebralWhiteMatter';
    
    ROI_ab(13).ids = [50]; ROI_ab(13).name = 'LH_InfLateralVentricle';
    ROI_ab(14).ids = [49]; ROI_ab(14).name = 'RH_InfLateralVentricle';
    
    ROI_ab(15).ids = [52]; ROI_ab(15).name = 'LH_LateralVentricle';
    ROI_ab(16).ids = [51]; ROI_ab(16).name = 'RH_LateralVentricle';
    
    ROI_ab(17).ids = [56]; ROI_ab(17).name = 'LH_Pallidum';
    ROI_ab(18).ids = [55]; ROI_ab(18).name = 'RH_Pallidum';
    
    ROI_ab(19).ids = [58]; ROI_ab(19).name = 'LH_Putamen';
    ROI_ab(20).ids = [57]; ROI_ab(20).name = 'RH_Putamen';
    
    ROI_ab(21).ids = [60]; ROI_ab(21).name = 'LH_Thalamus';
    ROI_ab(22).ids = [59]; ROI_ab(22).name = 'RH_Thalamus';
    
    ROI_ab(23).ids = [62]; ROI_ab(23).name = 'LH_VentralDC';
    ROI_ab(24).ids = [61]; ROI_ab(24).name = 'RH_VentralDC';
    
    ROI_ab(25).ids = [75]; ROI_ab(25).name = 'LH_BasalForebrain';
    ROI_ab(26).ids = [76]; ROI_ab(26).name = 'RH_BasalForebrain';
    
    ROI_ab(27).ids = [101]; ROI_ab(27).name = 'LH_AnteriorCingulateGyrus';
    ROI_ab(28).ids = [100]; ROI_ab(28).name = 'RH_AnteriorCingulateGyrus';
    
    ROI_ab(29).ids = [103]; ROI_ab(29).name = 'LH_AnteriorInsula';
    ROI_ab(30).ids = [102]; ROI_ab(30).name = 'RH_AnteriorInsula';
    
    ROI_ab(31).ids = [105]; ROI_ab(31).name = 'LH_AnteriorOrbitalGyrus';
    ROI_ab(32).ids = [104]; ROI_ab(32).name = 'RH_AnteriorOrbitalGyrus';
    
    ROI_ab(33).ids = [107]; ROI_ab(33).name = 'LH_AngularGyrus';
    ROI_ab(34).ids = [106]; ROI_ab(34).name = 'RH_AngularGyrus';
    
    ROI_ab(35).ids = [109]; ROI_ab(35).name = 'LH_CalcarineCortex';
    ROI_ab(36).ids = [108]; ROI_ab(36).name = 'RH_CalcarineCortex';
    
    
    ROI_ab(37).ids = [113]; ROI_ab(37).name = 'LH_CentralOperculum';
    ROI_ab(38).ids = [112]; ROI_ab(38).name = 'RH_CentralOperculum';
    
    ROI_ab(39).ids = [115]; ROI_ab(39).name = 'LH_Cuneus';
    ROI_ab(40).ids = [114]; ROI_ab(40).name = 'RH_Cuneus';
    
    ROI_ab(41).ids = [117]; ROI_ab(41).name = 'LH_EntorhinalArea';
    ROI_ab(42).ids = [116]; ROI_ab(42).name = 'RH_EntorhinalArea';
    
    ROI_ab(43).ids = [119]; ROI_ab(43).name = 'LH_FrontalOperculum';
    ROI_ab(44).ids = [118]; ROI_ab(44).name = 'RH_FrontalOperculum';
    
    ROI_ab(45).ids = [121]; ROI_ab(45).name = 'LH_FrontalPole';
    ROI_ab(46).ids = [120]; ROI_ab(46).name = 'RH_FrontalPole';
    
    ROI_ab(47).ids = [123]; ROI_ab(47).name = 'LH_FusiformGyrus';
    ROI_ab(48).ids = [122]; ROI_ab(48).name = 'RH_FusiformGyrus';
    
    ROI_ab(49).ids = [125]; ROI_ab(49).name = 'LH_GyrusRectus';
    ROI_ab(50).ids = [124]; ROI_ab(50).name = 'RH_GyrusRectus';
    
    ROI_ab(51).ids = [129]; ROI_ab(51).name = 'LH_InfOccipitalGyrus';
    ROI_ab(52).ids = [128]; ROI_ab(52).name = 'RH_InfOccipitalGyrus';
    
    ROI_ab(53).ids = [133]; ROI_ab(53).name = 'LH_InfTemporalGyrus';
    ROI_ab(54).ids = [132]; ROI_ab(54).name = 'RH_InfTemporalGyrus';
    
    ROI_ab(55).ids = [135]; ROI_ab(55).name = 'LH_LingualGyrus';
    ROI_ab(56).ids = [134]; ROI_ab(56).name = 'RH_LingualGyrus';
    
    ROI_ab(57).ids = [137]; ROI_ab(57).name = 'LH_LateralOrbitalGyrus';
    ROI_ab(58).ids = [136]; ROI_ab(58).name = 'RH_LateralOrbitalGyrus';
    
    ROI_ab(59).ids = [139]; ROI_ab(59).name = 'LH_MiddleCingulateGyrus';
    ROI_ab(60).ids = [138]; ROI_ab(60).name = 'RH_MiddleCingulateGyrus';
    
    ROI_ab(61).ids = [141]; ROI_ab(61).name = 'LH_MedialFrontalCortex';
    ROI_ab(62).ids = [140]; ROI_ab(62).name = 'RH_MedialFrontalCortex';
    
    
    ROI_ab(63).ids = [143]; ROI_ab(63).name = 'LH_MiddleFrontalGyrus';
    ROI_ab(64).ids = [142]; ROI_ab(64).name = 'RH_MiddleFrontalGyrus';
    
    ROI_ab(65).ids = [145]; ROI_ab(65).name = 'LH_MiddleOccipitalGyrus';
    ROI_ab(66).ids = [144]; ROI_ab(66).name = 'RH_MiddleOccipitalGyrus';
    
    ROI_ab(67).ids = [147]; ROI_ab(67).name = 'LH_MedialOrbitalGyrus';
    ROI_ab(68).ids = [146]; ROI_ab(68).name = 'RH_MedialOrbitalGyrus';
    
    ROI_ab(69).ids = [149]; ROI_ab(69).name = 'LH_PostcentralGyrusMedial';
    ROI_ab(70).ids = [148]; ROI_ab(70).name = 'RH_PostcentralGyrusMedial';
    
    ROI_ab(71).ids = [151]; ROI_ab(71).name = 'LH_PrecentralGyrusMedial';
    ROI_ab(72).ids = [150]; ROI_ab(72).name = 'RH_PrecentralGyrusMedial';
    
    ROI_ab(73).ids = [153]; ROI_ab(73).name = 'LH_SupFrontalGyrusMedial';
    ROI_ab(74).ids = [152]; ROI_ab(74).name = 'RH_SupFrontalGyrusMedial';
    
    ROI_ab(75).ids = [155]; ROI_ab(75).name = 'LH_MiddleTemporalGyrus';
    ROI_ab(76).ids = [154]; ROI_ab(76).name = 'RH_MiddleTemporalGyrus';
    
    ROI_ab(77).ids = [157]; ROI_ab(77).name = 'LH_OccipitalPole';
    ROI_ab(78).ids = [156]; ROI_ab(78).name = 'RH_OccipitalPole';
    
    ROI_ab(79).ids = [161]; ROI_ab(79).name = 'LH_OccipitalFusiformGyrus';
    ROI_ab(80).ids = [160]; ROI_ab(80).name = 'RH_OccipitalFusiformGyrus';
    
    ROI_ab(81).ids = [163]; ROI_ab(81).name = 'LH_OpercularInfFrontalGyrus';
    ROI_ab(82).ids = [162]; ROI_ab(82).name = 'RH_OpercularInfFrontalGyrus';
    
    ROI_ab(83).ids = [165]; ROI_ab(83).name = 'LH_OrbitalInfFrontalGyrus';
    ROI_ab(84).ids = [164]; ROI_ab(84).name = 'RH_OrbitalInfFrontalGyrus';
    
    ROI_ab(85).ids = [167]; ROI_ab(85).name = 'LH_PosteriorCingulateGyrus';
    ROI_ab(86).ids = [166]; ROI_ab(86).name = 'RH_PosteriorCingulateGyrus';
    
    ROI_ab(87).ids = [169]; ROI_ab(87).name = 'LH_Precuneus';
    ROI_ab(88).ids = [168]; ROI_ab(88).name = 'RH_Precuneus';
    
    ROI_ab(89).ids = [171]; ROI_ab(89).name = 'LH_ParahippocampalGyrus';
    ROI_ab(90).ids = [170]; ROI_ab(90).name = 'RH_ParahippocampalGyrus';
    
    ROI_ab(91).ids = [173]; ROI_ab(91).name = 'LH_PosteriorInsula';
    ROI_ab(92).ids = [172]; ROI_ab(92).name = 'RH_PosteriorInsula';
    
    ROI_ab(93).ids = [175]; ROI_ab(93).name = 'LH_ParietalOperculum';
    ROI_ab(94).ids = [174]; ROI_ab(94).name = 'RH_ParietalOperculum';
    
    ROI_ab(95).ids = [177]; ROI_ab(95).name = 'LH_PostcentralGyrus';
    ROI_ab(96).ids = [176]; ROI_ab(96).name = 'RH_PostcentralGyrus';
    
    ROI_ab(97).ids = [179]; ROI_ab(97).name = 'LH_PosteriorOrbitalGyrus';
    ROI_ab(98).ids = [178]; ROI_ab(98).name = 'RH_PosteriorOrbitalGyrus';
    
    ROI_ab(99).ids = [181]; ROI_ab(99).name = 'LH_PlanumPolare';
    ROI_ab(100).ids = [180]; ROI_ab(100).name = 'RH_PlanumPolare';
    
    ROI_ab(101).ids = [183]; ROI_ab(101).name = 'LH_PrecentralGyrus';
    ROI_ab(102).ids = [182]; ROI_ab(102).name = 'RH_PrecentralGyrus';
    
    ROI_ab(103).ids = [185]; ROI_ab(103).name = 'LH_PlanumTemporale';
    ROI_ab(104).ids = [184]; ROI_ab(104).name = 'RH_PlanumTemporale';
    
    ROI_ab(105).ids = [187]; ROI_ab(105).name = 'LH_SubcallosalArea';
    ROI_ab(106).ids = [186]; ROI_ab(106).name = 'RH_SubcallosalArea';
    
    ROI_ab(107).ids = [191]; ROI_ab(107).name = 'LH_SupFrontalGyrus';
    ROI_ab(108).ids = [190]; ROI_ab(108).name = 'RH_SupFrontalGyrus';
    
    ROI_ab(109).ids = [193]; ROI_ab(109).name = 'LH_SupplementaryMotorCortex';
    ROI_ab(110).ids = [192]; ROI_ab(110).name = 'RH_SupplementaryMotorCortex';
    
    ROI_ab(111).ids = [195]; ROI_ab(111).name = 'LH_SupramarginalGyrus';
    ROI_ab(112).ids = [194]; ROI_ab(112).name = 'RH_SupramarginalGyrus';
    
    ROI_ab(113).ids = [197]; ROI_ab(113).name = 'LH_SupOccipitalGyrus';
    ROI_ab(114).ids = [196]; ROI_ab(114).name = 'RH_SupOccipitalGyrus';
    
    ROI_ab(115).ids = [199]; ROI_ab(115).name = 'LH_SupParietalLobule';
    ROI_ab(116).ids = [198]; ROI_ab(116).name = 'RH_SupParietalLobule';
    
    ROI_ab(117).ids = [201]; ROI_ab(117).name = 'LH_SupTemporalGyrus';
    ROI_ab(118).ids = [200]; ROI_ab(118).name = 'RH_SupTemporalGyrus';
    
    ROI_ab(119).ids = [203]; ROI_ab(119).name = 'LH_TemporalPole';
    ROI_ab(120).ids = [202]; ROI_ab(120).name = 'RH_TemporalPole';
    
    ROI_ab(121).ids = [205]; ROI_ab(121).name = 'LH_TriangularInfFrontalGyrus';
    ROI_ab(122).ids = [204]; ROI_ab(122).name = 'RH_TriangularInfFrontalGyrus';
    
    ROI_ab(123).ids = [207]; ROI_ab(123).name = 'LH_TransverseTemporalGyrus';
    ROI_ab(124).ids = [206]; ROI_ab(124).name = 'RH_TransverseTemporalGyrus';

    
    %% Extract mean magE per ROI
    for i = 1:length(ROI_ab)
        mask = ismember(tet_labels_ab, ROI_ab(i).ids);
        Evals = magnE_ab(mask);
        ROI_ab(i).mean = mean(Evals) * 800;
    end
 
    %% Build summary row
    roi_means = [];
    roi_names = {};
    for i = 1:length(ROI_ab)
        roi_means(end+1) = ROI_ab(i).mean;
        roi_names{end+1} = ROI_ab(i).name;
    end
 
    summary_table = [summary_table; array2table(roi_means, 'VariableNames', roi_names, 'RowNames', {subj_id})];
 
end
 
disp(summary_table);
writetable(summary_table, fullfile(sim_root, 'assemblynet_FP_sim_ROI_summary.xls'), 'WriteRowNames', true);