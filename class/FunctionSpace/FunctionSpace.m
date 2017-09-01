classdef FunctionSpace < handle

    properties
        address
        
        nodes
        elems
        edges
        edgeIds;
        neighbors;
    end

    methods
        function this = FunctionSpace(tm, deg)
            this.address = FunctionSpaceWrapper('new');   
            % get all information from triangle mesh.
            [p, s, t, e, n] = tm.getData_tri();
            % build up function space.
            logger('*Green', sprintf('Initializing FunctionSpace at %d.\n', this.address));
            build(this, p, s, t, e, n, deg);
            [this.nodes, this.elems, this.edges, this.edgeIds, this.neighbors] = getData(this);
            
        end

        function delete(this)
            logger('*Red', sprintf('    Destroy FunctionSpace at %d.\n', this.address));
            FunctionSpaceWrapper('delete', this.address);
        end
        
        function build(this, p, s, t, e, n, deg)
            % taking care of starting zero.
            s = s-1; t = t-1;
            e = e-1; n = n-1;
            FunctionSpaceWrapper('build',this.address, p,s,t,e,n, deg);
            logger('*cyan', sprintf('    Building FunctionSpace at %d.\n', this.address));
        end

        function [fp, ft, fe, feId, fn] = getData(this)
            [fp, ft, fe, feId, fn] = FunctionSpaceWrapper('getData',this.address);
            logger('*Magenta', sprintf('     Extract FunctionSpace at %d.\n', this.address));
        end
    end
end
