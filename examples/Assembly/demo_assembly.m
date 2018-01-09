clc;clear;

opt = struct('deg', 3, 'qdeg', 6, 'min_area', 1e-3, 'edge', [0 1 1 0; 0 0 1 1]);
V = femm(opt);

S = V.build('s', 1);
Q = V.build('e', 1, 'all');

f =@(x)(6 * x(1,:)); 
load_f = f(V.quad2d);
L = V.build('l', load_f);

g =@(x)(x(1,:).^3 .* (x(2,:) == 0) + x(1,:).^3 .* (x(2,:) == 1) + 4 * (x(1,:) == 1));
load_g = g(V.quad1d);
R = V.build('g', load_g, 'all');

A = -S - Q;
b = L  - R;

x = A\b;
trisurf(V.space.elems(1:3,:)', V.space.nodes(1,:), V.space.nodes(2,:), x);

[DX, DY] = V.builder.reference_grad(V.rffs.nodes);
[gradX, gradY] = V.gradient(x, DX, DY);


fprintf('L2 error is %2.6e\n', sqrt((x' - V.space.nodes(1,:).^3)* V.build('m', 1) * (x' - V.space.nodes(1,:).^3)'));