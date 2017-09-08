
%% 2d quadrature.
deg = 15;
qm = QuadMode(2, deg);
[q, w] = qm.export();
 error = 0;
%%  test performance.
for i = 0:deg
    for j = 0:(deg-i)
        val = factorial(i) * factorial(j)/factorial(i+j+2);
        out = (q(1,:).^i .* q(2,:).^j)*w/2.0 ;
        error = max(error, abs((out -val)/val));
    end
end

logger('*cyan', sprintf('     maximum relative error is %4.6e\n', error));

deg = 13;
qm = QuadMode(1, deg);
[q,w] = qm.export();
error = 0.;
%%  test performance.
for i = 0:deg / 2
    val = 2.0/(2 * i+1);
    out = (q.^(2*i))*w;
    error = max(error, abs((out -val)/val));
end

logger('*cyan', sprintf('     maximum relative error is %4.6e\n', error));