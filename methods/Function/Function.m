function ret = Function(V, handle)
assert((nargin == 2 || nargin == 1), 'Error: WrongNumberInput');
if (nargin == 1)
    ret = zeros(size(V.nodes, 2), 1);
else
    ret = handle(V.nodes);
end
end

