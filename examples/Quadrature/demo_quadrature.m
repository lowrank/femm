clc; clear;
%% 2d quadrature.
qm = QuadMode(2);
deg = 25;
[qx, qy, qw] = qm.get_vr_quad(deg); % degree which is accurate upto.

 error = 0;
%%  test performance.
for i = 0:deg
    for j = 0:(deg-i)
        val = factorial(i) * factorial(j)/factorial(i+j+2);
        out = (qx.^i .* qy.^j)'*qw ;
        error = max(error, abs((out -val)/val));
    end
end

logger('*cyan', sprintf('     maximum relative error is %4.6e\n', error));


qm = QuadMode(1);
deg = 13;
[qx, qw] = qm.get_vr_quad(deg);
error = 0.;
%%  test performance.
for i = 0:deg / 2
    val = 2.0/(2 * i+1);
    out = (qx.^(2*i))'*qw;
    error = max(error, abs((out -val)/val));
end

logger('*cyan', sprintf('     maximum relative error is %4.6e\n', error));