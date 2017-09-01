clc; clear;
%% create Mesh
mesh = TriangleMesh();
% setup
N = 256;
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));

nodes = [cos(theta); sin(theta)];
edges = [0:(N-2); 1:(N-1)]; 
edges(2,end) = 0;

%%
mesh.set_points_tri(reshape(nodes, size(nodes, 2) * 2, 1));
mesh.set_facets_tri(reshape(edges, size(edges, 2) * 2, 1));

mesh = mesh.build_tri(); % not ready.

mesh = mesh.refine_tri('q34.0a0.00125');
% mesh.getInfo_tri();

%% Function Space

f = @(x)(sin(2*pi*x(1,:).*sin(2*pi*x(2,:))));
V = FunctionSpace(mesh, 2);
u = Function(V, f); % allocation.
v = Function(V);

trisurf(V.elems(1:3, :)', V.nodes(1,:), V.nodes(2,:), v, 'EdgeColor', 'none');shading interp;axis equal;view(2);


