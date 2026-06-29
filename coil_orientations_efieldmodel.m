addpath('C:\Users\z5171263\SimNIBS-4.5\simnibs_env\Lib\site-packages\simnibs\matlab_tools');
addpath('C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools');

coord_MNI = [-39.1, 51.6, 30.2];
subdir = 'C:\Users\z5171263\Downloads\E_field_modelling\m2m\m2m_vivien';

coord_subjectspace = mni2subject_coords(coord_MNI, subdir);
head_mesh = mesh_load_gmsh4(fullfile(subdir, 'vivien.msh'));

gm_label = 1002;
gm_triangles = head_mesh.triangles(head_mesh.triangle_regions == gm_label, :);
gm_nodes_idx = unique(gm_triangles(:));
gm_nodes = head_mesh.nodes(gm_nodes_idx, :);

% Find closest GM node
d = vecnorm(gm_nodes - coord_subjectspace, 2, 2);
[~, idx_local] = min(d);
target_gm = gm_nodes(idx_local, :);
target_node_index = gm_nodes_idx(idx_local);

% Surface normal
tri_idx = any(gm_triangles == target_node_index, 2);
local_tris = gm_triangles(tri_idx, :);
normals = zeros(size(local_tris,1), 3);
for i = 1:size(local_tris,1)
    tri = local_tris(i,:);
    v1 = head_mesh.nodes(tri(1),:);
    v2 = head_mesh.nodes(tri(2),:);
    v3 = head_mesh.nodes(tri(3),:);
    n = cross(v2-v1, v3-v1);
    normals(i,:) = n / norm(n);
end
normal = mean(normals, 1);
normal = normal / norm(normal);

% Ensure normal points outward
head_centre = mean(head_mesh.nodes, 1);
if dot(normal, target_gm - head_centre) < 0
    normal = -normal;
end

% Project midline into cortical tangent plane
midline = [0 1 0];
midline_tangent = midline - dot(midline, normal) * normal;
midline_tangent = midline_tangent / norm(midline_tangent);

% Second tangent axis — negative cross product (key fix for correct rotation direction)
cross_term = -cross(normal, midline_tangent);
cross_term = cross_term / norm(cross_term);

%% Compute ref points for each angle
angles = [-45 -22.5 0 22.5 45 67.5 90 112.5 135];
offset = 20; % mm

fprintf('\n%-15s  %-35s  %-35s\n', 'Condition', 'pos (coil centre)', 'pos_ref (direction reference)');
fprintf('%s\n', repmat('-', 1, 90));

for i = 1:length(angles)
    a = angles(i);
    coil_vec = cosd(a)*midline_tangent + sind(a)*cross_term;
    coil_vec = coil_vec / norm(coil_vec);
    ref_point = target_gm + offset * coil_vec;
    fprintf('CT_A_%-8g   pos=[%6.2f %6.2f %6.2f]   ref=[%6.2f %6.2f %6.2f]\n', ...
        a, ...
        target_gm(1), target_gm(2), target_gm(3), ...
        ref_point(1), ref_point(2), ref_point(3));
end