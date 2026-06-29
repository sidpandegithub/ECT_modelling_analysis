addpath('C:\Users\z5171263\SimNIBS-4.5\simnibs_env\Lib\site-packages\simnibs\matlab_tools');
addpath('C:\Users\z5171263\SimNIBS-4.6\simnibs_env\Lib\site-packages\simnibs\matlab_tools');

coord_MNI = [-39.1, 51.6, 30.2]; % vivien
% coord_MNI = [-44, 42.4, 26.5]; % donel

% subdir = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\ECT modelling\vivien MRI\m2m_vivien';
% subdir = 'C:\Users\z5171263\OneDrive - UNSW\Desktop\ECT modelling\compressed NIFTI\m2m_donel';
% subdir  = 'C:\Users\z5171263\Downloads\E_field_modelling\m2m\m2m_vivien';
subdir = 'C:\Users\z5171263\Downloads\m2m_UA019';

coord_subjectspace  = mni2subject_coords(coord_MNI,subdir);

% head_mesh = mesh_load_gmsh4(fullfile(subdir,'vivien.msh'));
head_mesh = mesh_load_gmsh4(fullfile(subdir,'donel.msh'));

gm_label = 1002;

gm_triangles = head_mesh.triangles(head_mesh.triangle_regions == gm_label, :);
gm_nodes_idx = unique(gm_triangles(:));
gm_nodes = head_mesh.nodes(gm_nodes_idx, :);

% Find closest GM node to target
d = vecnorm(gm_nodes - coord_subjectspace, 2, 2);
[~, idx_local] = min(d);

target_gm = gm_nodes(idx_local, :);
target_node_index = gm_nodes_idx(idx_local);

% Find triangles touching the target node
tri_idx = any(gm_triangles == target_node_index, 2);
local_tris = gm_triangles(tri_idx, :);

% Compute surface normal
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

%% Define midline

midline = [0 1 0];

% Project midline into cortical plane
midline_tangent = midline - dot(midline,normal)*normal;
midline_tangent = midline_tangent / norm(midline_tangent);

%% Coil angles used in experiment
angles = [-45 -22.5 0 22.5 45 67.5 90 112.5 135];

% Define midline
midline = [0 1 0];

% Project midline into cortical plane
midline_tangent = midline - dot(midline,normal)*normal;
midline_tangent = midline_tangent / norm(midline_tangent);

%% Plot

figure
hold on
axis equal
view(3)

xlabel('X')
ylabel('Y')
zlabel('Z')

title('GM Surface + Coil Orientations')

trisurf(gm_triangles, ...
        head_mesh.nodes(:,1), ...
        head_mesh.nodes(:,2), ...
        head_mesh.nodes(:,3), ...
        'FaceColor',[0.8 0.8 0.8], ...
        'EdgeColor','none', ...
        'FaceAlpha',0.9);

% Plot target
plot3(target_gm(1),target_gm(2),target_gm(3), ...
      'ko','MarkerSize',8,'LineWidth',2)

scale = 10;

% Surface normal
h_normal = quiver3(target_gm(1),target_gm(2),target_gm(3), ...
        normal(1)*scale,normal(2)*scale,normal(3)*scale, ...
        'b','LineWidth',3, 'HandleVisibility','off');

% Generate colours
cmap = lines(length(angles));

% Store handles for legend
h = gobjects(length(angles),1);

% Plot coil orientations
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

camlight
lighting gouraud

% Build legend labels
labels = arrayfun(@(x) sprintf('%g°',x), angles, 'UniformOutput', false);

legend(h, labels, 'Location','bestoutside')

%% zoomed view 
figure
hold on
axis equal
axis off
view(3)

xlabel('X')
ylabel('Y')
zlabel('Z')

title('GM Surface + Coil Orientations')

% Plot GM surface
trisurf(gm_triangles, ...
        head_mesh.nodes(:,1), ...
        head_mesh.nodes(:,2), ...
        head_mesh.nodes(:,3), ...
        'FaceColor',[0.8 0.8 0.8], ...
        'EdgeColor','none', ...
        'FaceAlpha',0.9);

% Plot target
plot3(target_gm(1),target_gm(2),target_gm(3), ...
      'ko','MarkerSize',8,'LineWidth',2)

scale = 10;

% Surface normal
h_normal = quiver3(target_gm(1),target_gm(2),target_gm(3), ...
        normal(1)*scale,normal(2)*scale,normal(3)*scale, ...
        'b','LineWidth',3, 'HandleVisibility','off');

% Custom colours for each angle
cmap = [
    hex2rgb('#d73027');  % -45
    hex2rgb('#f46d43');  % -22.5
    hex2rgb('#fdae61');  % 0
    hex2rgb('#dfc27d');  % 22.5
    hex2rgb('#c7eae5');  % 45
    hex2rgb('#80cdc1');  % 67.5
    hex2rgb('#35978f');  % 90
    hex2rgb('#01665e');  % 112.5
    hex2rgb('#003c30')]; % 135

% Store handles for legend
h = gobjects(length(angles),1);

% Plot coil orientations
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

camlight
lighting gouraud

% Build legend labels
labels = arrayfun(@(x) sprintf('%g°',x), angles, 'UniformOutput', false);
legend(h, labels, 'Location','bestoutside')

% Zoom around target
zoom_range = 20; % adjust this to zoom in/out
xlim([target_gm(1)-zoom_range, target_gm(1)+zoom_range])
ylim([target_gm(2)-zoom_range, target_gm(2)+zoom_range])
zlim([target_gm(3)-zoom_range, target_gm(3)+zoom_range])

%% 

%% Plot

figure
hold on
axis equal
view(3)

% Remove axes
axis off

xlabel('X')
ylabel('Y')
zlabel('Z')

title('GM Surface + Coil Orientations')

trisurf(gm_triangles, ...
        head_mesh.nodes(:,1), ...
        head_mesh.nodes(:,2), ...
        head_mesh.nodes(:,3), ...
        'FaceColor',[0.8 0.8 0.8], ...
        'EdgeColor','none', ...
        'FaceAlpha',0.9);

% Plot target
plot3(target_gm(1),target_gm(2),target_gm(3), ...
      'ko','MarkerSize',8,'LineWidth',2)

scale = 10;

% Surface normal
h_normal = quiver3(target_gm(1),target_gm(2),target_gm(3), ...
        normal(1)*scale,normal(2)*scale,normal(3)*scale, ...
        'b','LineWidth',3, 'HandleVisibility','off');

% Generate colours
cmap = [
    hex2rgb('#d73027');  % -45
    hex2rgb('#f46d43');  % -22.5
    hex2rgb('#fdae61');  % 0
    hex2rgb('#dfc27d');  % 22.5
    hex2rgb('#c7eae5');  % 45
    hex2rgb('#80cdc1');  % 67.5
    hex2rgb('#35978f');  % 90
    hex2rgb('#01665e');  % 112.5
    hex2rgb('#003c30')]; % 135

% Store handles for legend
h = gobjects(length(angles),1);

% Plot coil orientations
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

camlight
lighting gouraud

% Build legend labels
labels = arrayfun(@(x) sprintf('%g°',x), angles, 'UniformOutput', false);

% Place legend inside figure and closer to plot
legend(h, labels, 'Location','northeast')  % you can also try 'northeastoutside', 'northwest', etc.
