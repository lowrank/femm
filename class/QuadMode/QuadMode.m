classdef QuadMode < handle

    properties (Access = private)
        address
        dimension
    end

    methods
        function this = QuadMode(dim, degree)
            this.address = ModeWrapper('new', dim, degree);
            this.dimension = dim;
            logger('*Green', sprintf('Initializing QuadMode at %d.\n', this.address));
        end

        function delete(this)
            logger('*Red', sprintf('  Destroying QuadMode at %d.\n', this.address));
            ModeWrapper('delete', this.address);    
        end

        function [output] = getAddress(this)
            output = this.address;
        end
        
        function [qpts, qwts] = export(this)
  	        [qpts, qwts] = ModeWrapper('export', this.address);
        end
    end
end 