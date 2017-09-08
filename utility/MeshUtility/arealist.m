function ret = arealist(p, t)
n = size(t, 2);
ret = zeros(n, 1);
for i = 1:n
    u = t(1, i);
    v = t(2, i);
    w = t(3, i);
    det = (p(1, v) - p(1, u)) * (p(2, w) - p(2, u)) - (p(2, v) - p(2, u))*(p(1, w) - p(1, u));
    ret(i) = 0.5 * det;
end
end