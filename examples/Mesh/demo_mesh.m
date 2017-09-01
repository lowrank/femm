clc; clear;
%% Mesh
mesh = TriangleMesh();
% setup
mesh.set_points_tri([0 0 1 0 0.5 0.7]');
mesh.set_facets_tri([0 1 1 2 2 0]');

mesh = mesh.build_tri(); % not ready.

mesh = mesh.refine_tri('q34.0a0.0000125');
mesh.getInfo_tri();