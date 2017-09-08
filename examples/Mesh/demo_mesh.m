clc; clear;
%% Mesh
mesh = TriangleMesh();
% setup
mesh.set_points_tri([0 0 1 0 0.5 0.7]');
mesh.set_facets_tri([0 1 1 2 2 0]');

mesh = mesh.build_tri(); % not ready.

mesh = mesh.refine_tri('q34.0a0.0000125');
mesh.getInfo_tri();

C = mesh.connectivity();
[pts,~, elem, ~, ~] = mesh.getData_tri();

K= 8;
map = MetisPartition('PartGraphRecursive', C, K);

figure; hold on;
for val = 0:K-1
v = elem(:, map==val);
p = pts(:, unique(v));
scatter(p(1,:), p(2,:));
end
hold off;
