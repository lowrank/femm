clc; clear;
global logger_level
logger_level = 0;
%% create Mesh
mesh = TriangleMesh();

%% setup
N = 256;
theta = linspace(0, 2*pi, N);
theta = theta(1:(N-1));

nodes = [cos(theta); sin(theta)];
edges = [0:(N-2); 1:(N-1)]; 
edges(2,end) = 0;

%% re-format input.
mesh.set_points_tri(reshape(nodes, size(nodes, 2) * 2, 1));
mesh.set_facets_tri(reshape(edges, size(edges, 2) * 2, 1));
mesh = mesh.build_tri(); % not ready.
mesh = mesh.refine_tri('q34.0a0.00125');  % ready now.


%% 2nd order FE Space, mesh must be refined before applied.
V = FunctionSpace(mesh, 2, 4); % partition is applied.

%% 4th order quadrature
fb = FormBuilder(V, 4);

