classdef TriangleMesh < handle
    
    properties 
        address
    end
    
    methods
        function this = TriangleMesh()
            this.address = TriangleWrapper('new');
            logger('*[0, 1, 0]',sprintf( 'Initializing TriangleMesh at %d.\n', this.address));
        end
        
        function delete(this)
            logger('*[1 ,0, 0]', sprintf( '    Destory TriangleMesh at %d.\n', this.address));
            TriangleWrapper('delete', this.address);
        end
        
        function getInfo_tri(this)
            TriangleWrapper('getInfo_tri', this.address);
        end
        
        function [output] = getAddress(this)
            output = this.address;
        end
        
        function [p, s, t, e, n] = getData_tri(this)
            [p ,s ,t, e, n] = TriangleWrapper('getData_tri', this.address);
        end

        function set_points_tri(this, points)
            TriangleWrapper('set_points_tri', this.address, points);
        end

        function set_facets_tri(this, facets)
            TriangleWrapper('set_facets_tri', this.address, facets);
        end

        function that = build_tri(this)
            that = TriangleMesh();
            TriangleWrapper('build_tri', this.address, that.address);
            logger('*[1, 1, 0.2]', sprintf('    Building TriangleMesh at %d.\n', that.address));
            
        end
        
        function that = refine_tri(this, switches)
            that = TriangleMesh();
            TriangleWrapper('refine_tri', this.address, that.address, switches);
            logger('*[0, 1., 1.]',sprintf('    Refining TriangleMesh at %d.\n', that.address));
            
        end
        
        function conn = connectivity(this)
            [I, J, V] = TriangleWrapper('connectivity', this.address);
            conn = sparse(I,J,V);
        end
        

    end
    
    methods (Static)
        function about_tri()
            TriangleWrapper('verionInfo_tri');
        end
    end
    
end
