addpath('C:\Users\z5171263\SimNIBS-4.5\simnibs_env\Lib\site-packages\simnibs\matlab_tools');
addpath('C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools');

% coord_MNI = [-39.1, 51.6, 30.2];%vivien
% coord_MNI = [-44, 42.4, 26.5]; %donel
% coord_MNI = [-39.8, 47.1, 27.0]; %steve
% coord_MNI = [-41.8, 46.6, 31.2]; %cesar
% coord_MNI = [-41.8, 42.5, 27.6]; %adam
% coord_MNI = [-41.6, 50.8, 41.4];% nahian
% coord_MNI = [-43.9, 44.0, 30.2];%zhaoxia
% coord_MNI = [-40.0, 43.9, 28.3];%yann
% coord_MNI = [-41.7, 46.4, 33.9]%anthony
coord_MNI = [-42.6, 48.5, 29.3]%albert


subdir  = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\E_field_modelling\m2m\m2m_albert';
% subdir = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\E_field_modelling\m2m\m2m_donel';

coord_subjectspace = mni2subject_coords(coord_MNI, subdir);
head_mesh = mesh_load_gmsh4(fullfile(subdir,'albert.msh'));
% head_mesh = mesh_load_gmsh4(fullfile(subdir, 'donel.msh'));

gm_label = 1002;

gm_triangles = head_mesh.triangles(head_mesh.triangle_regions == gm_label, :);
gm_nodes_idx = unique(gm_triangles(:));
gm_nodes = head_mesh.nodes(gm_nodes_idx, :);

d = vecnorm(gm_nodes - coord_subjectspace, 2, 2);
[~, idx_local] = min(d);

target_gm = gm_nodes(idx_local, :);
target_node_index = gm_nodes_idx(idx_local);

tri_idx = any(gm_triangles == target_node_index, 2);
local_tris = gm_triangles(tri_idx, :);

normals = zeros(size(local_tris,1),3);

for i = 1:size(local_tris,1)
    tri = local_tris(i,:);
    v1 = head_mesh.nodes(tri(1),:);
    v2 = head_mesh.nodes(tri(2),:);
    v3 = head_mesh.nodes(tri(3),:);
    n = cross(v2-v1, v3-v1);
    normals(i,:) = n / norm(n);
end

normal = mean(normals,1);
normal = normal / norm(normal);

midline = [0 1 0];
midline_tangent = midline - dot(midline,normal)*normal;

if norm(midline_tangent) < 1e-6
    midline = [1 0 0];
    midline_tangent = midline - dot(midline,normal)*normal;
end

midline_tangent = midline_tangent / norm(midline_tangent);

angles = [-45 -22.5 0 22.5 45 67.5 90 112.5 135];
%% 

figure
hold on
axis equal
axis off
view(3)

trisurf(gm_triangles, ...
    head_mesh.nodes(:,1), ...
    head_mesh.nodes(:,2), ...
    head_mesh.nodes(:,3), ...
    'FaceColor',[0.6 0.6 0.6], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.9);

plot3(target_gm(1),target_gm(2),target_gm(3), ...
    'ko','MarkerSize',8,'LineWidth',2)

scale = 10;

quiver3(target_gm(1),target_gm(2),target_gm(3), ...
    normal(1)*scale,normal(2)*scale,normal(3)*scale, ...
    'b','LineWidth',3,'HandleVisibility','off');

cmap = [
    hex2rgb('#d73027');
    hex2rgb('#f46d43');
    hex2rgb('#fdae61');
    hex2rgb('#dfc27d');
    hex2rgb('#c7eae5');
    hex2rgb('#80cdc1');
    hex2rgb('#35978f');
    hex2rgb('#01665e');
    hex2rgb('#003c30')];

h = gobjects(length(angles),1);

for i = 1:length(angles)
    a = angles(i);
    coil_vec = cosd(a)*midline_tangent + sind(a)*cross(normal,midline_tangent);
    coil_vec = coil_vec / norm(coil_vec);
    h(i) = quiver3(target_gm(1),target_gm(2),target_gm(3), ...
        coil_vec(1)*scale, ...
        coil_vec(2)*scale, ...
        coil_vec(3)*scale, ...
        'Color',cmap(i,:), ...
        'LineWidth',2);
end

theta = linspace(0,2*pi,100);
radius = scale * 0.8;

% circle_pts = target_gm + ...
%     radius*(cos(theta')*midline_tangent + sin(theta')*cross(normal,midline_tangent));
% 
% plot3(circle_pts(:,1), circle_pts(:,2), circle_pts(:,3), 'k--','LineWidth',1);

camlight
lighting gouraud

labels = arrayfun(@(x) sprintf('%g°',x), angles, 'UniformOutput', false);
legend(h, labels, 'Location','northeast')



%% Find coil angle most perpendicular to local gyrus direction

perpendicularity = zeros(length(angles),1);

for i = 1:length(angles)

    a = angles(i);

    coil_vec = cosd(a)*midline_tangent + ...
               sind(a)*cross(normal,midline_tangent);

    coil_vec = coil_vec / norm(coil_vec);

    % Angle between coil and gyrus (90 deg = perpendicular)
    perpendicularity(i) = acosd(abs(dot(coil_vec,gyrus_direction)));

end

% Find closest to 90 degrees
[~,best_idx] = max(perpendicularity);
best_angle = angles(best_idx);

fprintf('Most perpendicular coil angle: %.1f degrees\n', best_angle);
fprintf('Angle relative to gyrus: %.1f degrees\n', perpendicularity(best_idx));
%% Estimate gyrus direction using a wider local patch

radius = 8; % mm - tune this (5-10mm is typical for gyral width/length scale)

% Find all GM nodes within radius of target (Euclidean, not just 1-ring)
d_all = vecnorm(gm_nodes - target_gm, 2, 2);
patch_mask = d_all <= radius;
local_coords = gm_nodes(patch_mask, :);

fprintf('Gyrus direction estimated from %d nodes within %.1f mm\n', ...
    size(local_coords,1), radius);

% PCA on the wider patch
centered = local_coords - mean(local_coords,1);
[~,S,V] = svd(centered,'econ');

% Sanity check: how elongated is the patch? (ratio of 1st to 2nd singular value)
elongation_ratio = S(1,1) / S(2,2);
fprintf('Elongation ratio (1st/2nd component): %.2f\n', elongation_ratio);

gyrus_direction = V(:,1)';
gyrus_direction = gyrus_direction / norm(gyrus_direction);

% Project onto cortical tangent plane (same as before)
gyrus_direction = gyrus_direction - dot(gyrus_direction,normal)*normal;
gyrus_direction = gyrus_direction / norm(gyrus_direction);

%% Plot gyrus direction (wide patch) and best perpendicular coil angle

figure
hold on
axis equal
axis off
view(3)

trisurf(gm_triangles, ...
    head_mesh.nodes(:,1), ...
    head_mesh.nodes(:,2), ...
    head_mesh.nodes(:,3), ...
    'FaceColor',[0.6 0.6 0.6], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.8);

% Highlight the local patch used for the PCA fit
plot3(local_coords(:,1), local_coords(:,2), local_coords(:,3), ...
    'y.', 'MarkerSize', 14);

plot3(target_gm(1),target_gm(2),target_gm(3), ...
    'ko','MarkerSize',8,'LineWidth',2);

scale = 15;

% Plot gyrus orientation (now from the wide patch)
quiver3(target_gm(1),target_gm(2),target_gm(3), ...
    gyrus_direction(1)*scale, ...
    gyrus_direction(2)*scale, ...
    gyrus_direction(3)*scale, ...
    'g','LineWidth',3);

% Calculate best coil vector
coil_vec_best = cosd(best_angle)*midline_tangent + ...
                sind(best_angle)*cross(normal,midline_tangent);
coil_vec_best = coil_vec_best / norm(coil_vec_best);

% Plot best coil orientation
quiver3(target_gm(1),target_gm(2),target_gm(3), ...
    coil_vec_best(1)*scale, ...
    coil_vec_best(2)*scale, ...
    coil_vec_best(3)*scale, ...
    'r','LineWidth',3);

legend({'Cortex','Local patch (radius fit)','Target','Gyrus direction', ...
        sprintf('Best coil angle %.1f°',best_angle)}, ...
        'Location','northeast')

title(sprintf('Radius = %.1f mm, elongation ratio = %.2f', ...
    radius, elongation_ratio))

camlight
lighting gouraud