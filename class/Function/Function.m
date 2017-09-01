classdef Function < handle
    properties(Access = private)
    end
    
    properties (Access = public)
        deg
        data
    end
    
    methods 
        function this = Function(V, f)
            if (isa(V, 'FunctionSpace') == 1) 
                this.deg = V.deg;

                if (nargin == 2)
                    this.data = f(V.nodes'); % caution: nodes is 2 x n matrix.
                else
                    this.data = zeros(size(V.nodes, 2) , 1);
                end
            else
                this.deg = V.deg;
                this.data = zeros(V.n, 1);
            end
        end
        
        function ret = plus(obj1, obj2)
            assert(isa(obj1, 'Function') == 1);
            assert(isa(obj2, 'Function') == 1);
            n1 = size(obj1.data, 1);
            n2 = size(obj2.data, 1);
            
            assert(n1 == n2);
            ret = Function(struct('deg', obj1.deg, 'n', n1));
            ret.data = obj1.data + obj2.data;
            ret.deg = max(obj1.deg, obj2.deg);
        end
        
    end
    
end

