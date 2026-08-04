addpath("C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools");

mesh_file    = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\volume\m2m_AL002\AL002.msh';
dkt_nii_file = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\ramsay_lakeside\AL_002_aparc.DKTatlas+aseg.nii';
vb_file      = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\analysis\volume\AL_002_job2010908_archive\native_structures_job2010908.nii.gz';


head_mesh = mesh_load_gmsh4(mesh_file);
nodes    = head_mesh.nodes;
tris     = head_mesh.triangles;
tri_regs = head_mesh.triangle_regions;
 
disp('Unique triangle region tags in this mesh:'); disp(unique(tri_regs));
GM_SURFACE_TAG = 1002;   % <-- confirm against list printed above
gm_tris = tris(tri_regs == GM_SURFACE_TAG, :);
fprintf('GM surface: %d triangles\n', size(gm_tris,1));
 
tri_centers = ( nodes(gm_tris(:,1),:) + nodes(gm_tris(:,2),:) + nodes(gm_tris(:,3),:) ) / 3;
tri_h = [tri_centers ones(size(tri_centers,1),1)]';
 

%  DKT LABELS (nearest-labeled-voxel search, restricted to the 66 valid

vol_dkt  = niftiread(dkt_nii_file);
info_dkt = niftiinfo(dkt_nii_file);
ras2vox_dkt = inv(info_dkt.Transform.T');
vox_dkt = (ras2vox_dkt * tri_h)';
vox_dkt = vox_dkt(:,1:3);
 
dkt_valid_ids = [17,18,53,54,1002,1003,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1034,1035,2002,2003,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2034,2035];
 
dkt_mask = ismember(vol_dkt, dkt_valid_ids);
[ix, iy, iz] = ind2sub(size(vol_dkt), find(dkt_mask));
dkt_labeled_coords = [ix, iy, iz];
dkt_labeled_ids    = vol_dkt(sub2ind(size(vol_dkt), ix, iy, iz));
 
fprintf('DKT labeled voxels available for matching: %d\n', numel(dkt_labeled_ids));
 
inside_dkt = vox_dkt(:,1)>=1 & vox_dkt(:,1)<=size(vol_dkt,1) & ...
             vox_dkt(:,2)>=1 & vox_dkt(:,2)<=size(vol_dkt,2) & ...
             vox_dkt(:,3)>=1 & vox_dkt(:,3)<=size(vol_dkt,3);
 
dkt_labels = nan(size(gm_tris,1),1);
MAX_SEARCH_DIST_VOX_DKT = 3;
[nn_idx_dkt, nn_dist_dkt] = knnsearch(dkt_labeled_coords, vox_dkt(inside_dkt,:));
good_match_dkt = nn_dist_dkt <= MAX_SEARCH_DIST_VOX_DKT;
idx_inside_dkt = find(inside_dkt);
dkt_labels(idx_inside_dkt(good_match_dkt)) = dkt_labeled_ids(nn_idx_dkt(good_match_dkt));
 
%% 
%  VOLBRAIN LABELS (nearest-labeled-voxel search - fills small label gaps)
vol_vb  = niftiread(vb_file);
info_vb = niftiinfo(vb_file);
ras2vox_vb = inv(info_vb.Transform.T');
vox_vb = (ras2vox_vb * tri_h)';
vox_vb = vox_vb(:,1:3);
 
CORTEX_MIN = 100; CORTEX_MAX = 207;
[ix, iy, iz] = ind2sub(size(vol_vb), find(vol_vb >= CORTEX_MIN & vol_vb <= CORTEX_MAX));
labeled_voxel_coords = [ix, iy, iz];
labeled_voxel_ids    = vol_vb(sub2ind(size(vol_vb), ix, iy, iz));
 
inside_vb = vox_vb(:,1)>=1 & vox_vb(:,1)<=size(vol_vb,1) & ...
            vox_vb(:,2)>=1 & vox_vb(:,2)<=size(vol_vb,2) & ...
            vox_vb(:,3)>=1 & vox_vb(:,3)<=size(vol_vb,3);
 
vb_labels = nan(size(gm_tris,1),1);
MAX_SEARCH_DIST_VOX = 3;
[nn_idx, nn_dist] = knnsearch(labeled_voxel_coords, vox_vb(inside_vb,:));
good_match = nn_dist <= MAX_SEARCH_DIST_VOX;
idx_inside = find(inside_vb);
vb_labels(idx_inside(good_match)) = labeled_voxel_ids(nn_idx(good_match));
 
%% CROSSWALK: DKT id -> set of "agreeing" volBrain ids

crosswalk = containers.Map('KeyType','double','ValueType','any');
% Direct 1:1
crosswalk(1005)=[115]; crosswalk(2005)=[114];   % Cuneus
crosswalk(1006)=[117]; crosswalk(2006)=[116];   % Entorhinal
crosswalk(1007)=[123]; crosswalk(2007)=[122];   % Fusiform
crosswalk(1009)=[133]; crosswalk(2009)=[132];   % InferiorTemporal
crosswalk(1012)=[137]; crosswalk(2012)=[136];   % LateralOrbitofrontal
crosswalk(1013)=[135]; crosswalk(2013)=[134];   % Lingual
crosswalk(1014)=[147]; crosswalk(2014)=[146];   % MedialOrbitofrontal
crosswalk(1015)=[155]; crosswalk(2015)=[154];   % MiddleTemporal
crosswalk(1016)=[171]; crosswalk(2016)=[170];   % Parahippocampal
crosswalk(1018)=[163]; crosswalk(2018)=[162];   % ParsOpercularis
crosswalk(1019)=[165]; crosswalk(2019)=[164];   % ParsOrbitalis
crosswalk(1020)=[205]; crosswalk(2020)=[204];   % ParsTriangularis
crosswalk(1021)=[109]; crosswalk(2021)=[108];   % Pericalcarine
crosswalk(1022)=[177]; crosswalk(2022)=[176];   % Postcentral
crosswalk(1023)=[167]; crosswalk(2023)=[166];   % PosteriorCingulate
crosswalk(1024)=[183]; crosswalk(2024)=[182];   % Precentral
crosswalk(1025)=[169]; crosswalk(2025)=[168];   % Precuneus
crosswalk(1029)=[199]; crosswalk(2029)=[198];   % SuperiorParietal
crosswalk(1030)=[201]; crosswalk(2030)=[200];   % SuperiorTemporal
crosswalk(1031)=[195]; crosswalk(2031)=[194];   % Supramarginal
crosswalk(1034)=[207]; crosswalk(2034)=[206];   % TransverseTemporal
crosswalk(17)=[48];    crosswalk(53)=[47];      % Hippocampus
crosswalk(18)=[32];    crosswalk(54)=[31];      % Amygdala
% Composite / merged (DKT id -> multiple acceptable volBrain ids)
crosswalk(1002)=[101]; crosswalk(1026)=[101];   % Caudal+RostralACC (LH) -> AnteriorCingulateGyrus
crosswalk(2002)=[100]; crosswalk(2026)=[100];   % (RH)
crosswalk(1003)=[143]; crosswalk(1027)=[143];   % Caudal+RostralMFG (LH) -> MiddleFrontalGyrus
crosswalk(2003)=[142]; crosswalk(2027)=[142];   % (RH)
crosswalk(1028)=[191,153];  crosswalk(2028)=[190,152];   % SuperiorFrontal -> lateral+medial
crosswalk(1035)=[103,173];  crosswalk(2035)=[102,172];   % Insula -> anterior+posterior
crosswalk(1017)=[151,149];  crosswalk(2017)=[150,148];   % Paracentral -> medial pre+postcentral
crosswalk(1011)=[129,145,197]; crosswalk(2011)=[128,144,196]; % LateralOccipital -> inf+mid+sup occipital
crosswalk(1008)=[107];      crosswalk(2008)=[106];       % InferiorParietal ~ Angular (approx)
% No volBrain counterpart -> isthmus cingulate excluded (left out of map = "no correspondence")
 
%% BUILD DIFFERENCE (AGREEMENT) LABELS PER TRIANGLE
%  1 = agree, 0 = disagree, NaN = no correspondence / missing data
agreement = nan(size(gm_tris,1),1);
for t = 1:size(gm_tris,1)
    d = dkt_labels(t); v = vb_labels(t);
    if isnan(d) || isnan(v), continue; end
    if ~isKey(crosswalk, d), continue; end   % e.g. isthmus cingulate - no defined correspondence
    expected_vb_ids = crosswalk(d);
    agreement(t) = double(ismember(v, expected_vb_ids));
end
 
n_valid = sum(~isnan(agreement));
n_agree = sum(agreement==1);
fprintf('Triangles with a defined correspondence: %d (%.1f%% of GM surface)\n', ...
    n_valid, 100*n_valid/size(gm_tris,1));
fprintf('Agreement: %d / %d = %.1f%%\n', n_agree, n_valid, 100*n_agree/n_valid);
 
%%  COLORS
% DKT: distinct color per region (consistent, fixed order over 1000-2035)
dkt_ids_all = [17,18,53,54,1002,1003,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1034,1035,2002,2003,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2034,2035];
dkt_cmap = hsv(numel(dkt_ids_all));
dkt_id_to_color = containers.Map(dkt_ids_all, num2cell(dkt_cmap,2));
dkt_face_color = repmat([0.6 0.6 0.6], size(gm_tris,1), 1);
for t = 1:size(gm_tris,1)
    if ~isnan(dkt_labels(t)) && isKey(dkt_id_to_color, dkt_labels(t))
        dkt_face_color(t,:) = dkt_id_to_color(dkt_labels(t));
    end
end
 
% volBrain: fixed colormap over cortical range (as before)
vb_ids_all = CORTEX_MIN:CORTEX_MAX;
vb_cmap = hsv(numel(vb_ids_all));
vb_id_to_color = containers.Map(vb_ids_all, num2cell(vb_cmap,2));
vb_face_color = repmat([0.6 0.6 0.6], size(gm_tris,1), 1);
valid_vb = ~isnan(vb_labels);
vb_face_color(valid_vb,:) = cell2mat(arrayfun(@(id) vb_id_to_color(id), vb_labels(valid_vb), 'UniformOutput', false));
 
% Difference: green = agree, red = disagree, grey = no data/correspondence
diff_face_color = repmat([0.75 0.75 0.75], size(gm_tris,1), 1);  % grey default
diff_face_color(agreement==1,:) = repmat([0.10 0.75 0.10], sum(agreement==1), 1);  % green
diff_face_color(agreement==0,:) = repmat([0.85 0.10 0.10], sum(agreement==0), 1);  % red
 
%% ID -> NAME LOOKUPS (for legend text)
dkt_id_to_name = containers.Map('KeyType','double','ValueType','char');
dkt_id_to_name(17)='LH_Hippocampus'; dkt_id_to_name(53)='RH_Hippocampus';
dkt_id_to_name(18)='LH_Amygdala';    dkt_id_to_name(54)='RH_Amygdala';
dkt_names_ctx = {'CaudalAnteriorCingulate','CaudalMiddleFrontal','Cuneus','Entorhinal','Fusiform', ...
    'InferiorParietal','InferiorTemporal','IsthmusCingulate','LateralOccipital','LateralOrbitofrontal', ...
    'Lingual','MedialOrbitofrontal','MiddleTemporal','Parahippocampal','Paracentral','ParsOpercularis', ...
    'ParsOrbitalis','ParsTriangularis','Pericalcarine','Postcentral','PosteriorCingulate','Precentral', ...
    'Precuneus','RostralAnteriorCingulate','RostralMiddleFrontal','SuperiorFrontal','SuperiorParietal', ...
    'SuperiorTemporal','Supramarginal','TransverseTemporal','Insula'};
dkt_ids_ctx_lh = [1002,1003,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1034,1035];
dkt_ids_ctx_rh = [2002,2003,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2034,2035];
for i = 1:length(dkt_names_ctx)
    dkt_id_to_name(dkt_ids_ctx_lh(i)) = ['LH_' dkt_names_ctx{i}];
    dkt_id_to_name(dkt_ids_ctx_rh(i)) = ['RH_' dkt_names_ctx{i}];
end
 
vb_id_to_name = containers.Map('KeyType','double','ValueType','char');
vb_pairs = {101,'LH_AnteriorCingulateGyrus';100,'RH_AnteriorCingulateGyrus';103,'LH_AnteriorInsula';102,'RH_AnteriorInsula'; ...
    105,'LH_AnteriorOrbitalGyrus';104,'RH_AnteriorOrbitalGyrus';107,'LH_AngularGyrus';106,'RH_AngularGyrus'; ...
    109,'LH_CalcarineCortex';108,'RH_CalcarineCortex';113,'LH_CentralOperculum';112,'RH_CentralOperculum'; ...
    115,'LH_Cuneus';114,'RH_Cuneus';117,'LH_EntorhinalArea';116,'RH_EntorhinalArea'; ...
    119,'LH_FrontalOperculum';118,'RH_FrontalOperculum';121,'LH_FrontalPole';120,'RH_FrontalPole'; ...
    123,'LH_FusiformGyrus';122,'RH_FusiformGyrus';125,'LH_GyrusRectus';124,'RH_GyrusRectus'; ...
    129,'LH_InfOccipitalGyrus';128,'RH_InfOccipitalGyrus';133,'LH_InfTemporalGyrus';132,'RH_InfTemporalGyrus'; ...
    135,'LH_LingualGyrus';134,'RH_LingualGyrus';137,'LH_LateralOrbitalGyrus';136,'RH_LateralOrbitalGyrus'; ...
    139,'LH_MiddleCingulateGyrus';138,'RH_MiddleCingulateGyrus';141,'LH_MedialFrontalCortex';140,'RH_MedialFrontalCortex'; ...
    143,'LH_MiddleFrontalGyrus';142,'RH_MiddleFrontalGyrus';145,'LH_MiddleOccipitalGyrus';144,'RH_MiddleOccipitalGyrus'; ...
    147,'LH_MedialOrbitalGyrus';146,'RH_MedialOrbitalGyrus';149,'LH_PostcentralGyrusMedial';148,'RH_PostcentralGyrusMedial'; ...
    151,'LH_PrecentralGyrusMedial';150,'RH_PrecentralGyrusMedial';153,'LH_SupFrontalGyrusMedial';152,'RH_SupFrontalGyrusMedial'; ...
    155,'LH_MiddleTemporalGyrus';154,'RH_MiddleTemporalGyrus';157,'LH_OccipitalPole';156,'RH_OccipitalPole'; ...
    161,'LH_OccipitalFusiformGyrus';160,'RH_OccipitalFusiformGyrus';163,'LH_OpercularInfFrontalGyrus';162,'RH_OpercularInfFrontalGyrus'; ...
    165,'LH_OrbitalInfFrontalGyrus';164,'RH_OrbitalInfFrontalGyrus';167,'LH_PosteriorCingulateGyrus';166,'RH_PosteriorCingulateGyrus'; ...
    169,'LH_Precuneus';168,'RH_Precuneus';171,'LH_ParahippocampalGyrus';170,'RH_ParahippocampalGyrus'; ...
    173,'LH_PosteriorInsula';172,'RH_PosteriorInsula';175,'LH_ParietalOperculum';174,'RH_ParietalOperculum'; ...
    177,'LH_PostcentralGyrus';176,'RH_PostcentralGyrus';179,'LH_PosteriorOrbitalGyrus';178,'RH_PosteriorOrbitalGyrus'; ...
    181,'LH_PlanumPolare';180,'RH_PlanumPolare';183,'LH_PrecentralGyrus';182,'RH_PrecentralGyrus'; ...
    185,'LH_PlanumTemporale';184,'RH_PlanumTemporale';187,'LH_SubcallosalArea';186,'RH_SubcallosalArea'; ...
    191,'LH_SupFrontalGyrus';190,'RH_SupFrontalGyrus';193,'LH_SupplementaryMotorCortex';192,'RH_SupplementaryMotorCortex'; ...
    195,'LH_SupramarginalGyrus';194,'RH_SupramarginalGyrus';197,'LH_SupOccipitalGyrus';196,'RH_SupOccipitalGyrus'; ...
    199,'LH_SupParietalLobule';198,'RH_SupParietalLobule';201,'LH_SupTemporalGyrus';200,'RH_SupTemporalGyrus'; ...
    203,'LH_TemporalPole';202,'RH_TemporalPole';205,'LH_TriangularInfFrontalGyrus';204,'RH_TriangularInfFrontalGyrus'; ...
    207,'LH_TransverseTemporalGyrus';206,'RH_TransverseTemporalGyrus'};
for i = 1:size(vb_pairs,1)
    vb_id_to_name(vb_pairs{i,1}) = vb_pairs{i,2};
end
 
%%  PLOTTING (DKT and volBrain get a legend panel; Agreement does not)

views = struct('name', {'Left_lateral','Right_lateral','Frontal'}, ...
                'az',   {180, 0, 90}, 'el', {0, 0, 0}, ...
                'hemi', {'LH','RH','BOTH'});
 
datasets = struct('prefix', {'DKT','volBrain','Agreement'}, ...
                   'color',  {dkt_face_color, vb_face_color, diff_face_color}, ...
                   'id_map', {dkt_id_to_name, vb_id_to_name, containers.Map()}, ...
                   'ids_present', {unique(dkt_labels(~isnan(dkt_labels))), unique(vb_labels(~isnan(vb_labels))), []}, ...
                   'colormap_fn', {dkt_id_to_color, vb_id_to_color, containers.Map()});
 
all_figs = gobjects(length(datasets), length(views));
 
for d = 1:length(datasets)
    add_legend = ~strcmp(datasets(d).prefix, 'Agreement');
    for v = 1:length(views)
        all_figs(d,v) = figure('Name', [datasets(d).prefix '_' views(v).name], 'Color','w', ...
                                'Position', [100 100 1100 700]);
        if add_legend
            ax_brain = axes('Position', [0.02 0.05 0.68 0.9]);
        else
            ax_brain = axes('Position', [0.05 0.05 0.9 0.9]);
        end
        patch('Faces', gm_tris, 'Vertices', nodes, 'FaceVertexCData', datasets(d).color, ...
              'FaceColor', 'flat', 'EdgeColor', 'none');
        axis equal off
        camlight('headlight'); lighting gouraud
        view(views(v).az, views(v).el);
        title([datasets(d).prefix ' - ' strrep(views(v).name,'_',' ')], 'Interpreter','none');
 
        if add_legend
            % Filter legend entries to the hemisphere(s) visible in this view
            ids = datasets(d).ids_present;
            names = cell(size(ids));
            for k = 1:length(ids)
                if isKey(datasets(d).id_map, ids(k))
                    names{k} = datasets(d).id_map(ids(k));
                else
                    names{k} = sprintf('ID_%d', ids(k));
                end
            end
            switch views(v).hemi
                case 'LH'
                    keep = startsWith(names, 'LH_');
                case 'RH'
                    keep = startsWith(names, 'RH_');
                otherwise
                    keep = true(size(names));
            end
            legend_ids = ids(keep);
            legend_names = names(keep);
 
            ax_legend = axes('Position', [0.72 0.03 0.26 0.94]);
            axis(ax_legend, [0 1 0 1]); axis(ax_legend, 'off'); hold(ax_legend, 'on');
            n = length(legend_ids);
            row_h = 1 / max(n,1);
            for k = 1:n
                y_top = 1 - (k-1)*row_h;
                col = datasets(d).colormap_fn(legend_ids(k));
                rectangle(ax_legend, 'Position', [0, y_top-row_h*0.8, 0.08, row_h*0.7], ...
                          'FaceColor', col, 'EdgeColor', 'none');
                text(ax_legend, 0.11, y_top - row_h*0.45, strrep(legend_names{k}, '_', ' '), ...
                     'FontSize', 6.5, 'Interpreter', 'none', 'VerticalAlignment','middle');
            end
            hold(ax_legend, 'off');
        end
    end
end
 