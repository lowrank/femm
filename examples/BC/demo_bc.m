clc;clear;

opt = struct('deg', 3, 'qdeg', 6, 'min_area', 4e-3, 'edge', [0 1 1 0; 0 0 1 1]);
V = femm(opt);

S = V.build('s', 1);
Q = V.build('e', 1, 'all');

f =@(x)(6 * x(1,:)); 
load_f = f(V.quad2d);
L = V.build('l', load_f);

bcs = BC('robin');
bcs.set_constraint('x-1');
bcs.set_constraint('x');
bcs.set_constraint('y-1');
bcs.set_constraint('y');
[e1, e2, e3, e4] = bcs.get_boundary(V.space.edges, V.space.nodes, 4);

q1 = V.builder.qnodes1D(V.space.nodes, V.edge.qnodes, e1);
q2 = V.builder.qnodes1D(V.space.nodes, V.edge.qnodes, e2);
q3 = V.builder.qnodes1D(V.space.nodes, V.edge.qnodes, e3);
q4 = V.builder.qnodes1D(V.space.nodes, V.edge.qnodes, e4);

g1 = @(x)(4 * (x(1,:) == 1));
g2 = @(x)(0 * (x(1,:) == 0));
g3 = @(x)(x(1,:).^3 .* (x(2,:) == 1));
g4 = @(x)(x(1,:).^3 .* (x(2,:) == 0));

load_g1 = g1(q1);
load_g2 = g2(q2);
load_g3 = g3(q3);
load_g4 = g4(q4);

R1 = V.build('g', load_g1, e1);
R2 = V.build('g', load_g2, e2);
R3 = V.build('g', load_g3, e3);
R4 = V.build('g', load_g4, e4);

A = -S - Q;
b = L  - R1 - R2 - R3 - R4;

x = A\b;
trisurf(V.space.elems(1:3,:)', V.space.nodes(1,:), V.space.nodes(2,:), x);

fprintf('L2 error is %2.6e\n', sqrt((x' - V.space.nodes(1,:).^3)* V.build('m', 1) * (x' - V.space.nodes(1,:).^3)'));