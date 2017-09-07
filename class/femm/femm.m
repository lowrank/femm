classdef femm < handle    
    properties
        facet
        edge
        space
        quadrature
        builder
        data
    end
    
    methods
        function obj = femm(opt)
            assert(isfield(opt, 'kway'));
            assert(isfield(opt, 'edge'));
            assert(size(opt.edge, 1) == 2);
            assert(isfield(opt, 'min_area'));
            if isfield(opt, 'qdeg')
                d = max(opt.deg * 2, opt.qdeg);
            else
                d = opt.deg * 2;
            end
            % now builds a reference mesh for all preliminary information.
            % this information should be built by other means in theory.
            % Now just uses brute-force.
            rftm = TriangleMesh();
            rftm.set_points_tri([0 0 1 0 0 1]');   
            rftm.set_facets_tri([0 1 1 2 2 0]');
            rftm = rftm.build_tri();
            rftm = rftm.refine_tri('q34a1');
            
            obj.builder = FormBuilder();
            % reference function space, also by brute-force. Partition is
            % not allowed here.
            rffs = FunctionSpace(rftm, opt.deg);
            %% Geometry
            % fill the information of facet.
            [obj.facet.qnodes, obj.facet.weights] = ...
                QuadMode(2, d).export();
            [obj.facet.ref, obj.facet.refx, ...
                obj.facet.refy] = obj.builder.reference2D(rffs.nodes,...
                obj.facet.qnodes);
            
            % fill the information of edge.
            [obj.edge.qnodes, obj.edge.weigths] = ...
                QuadMode(1,d).export();
            
            [obj.edge.ref, ~] = obj.builder.reference1D(opt.deg, obj.edge.qnodes);
            
            % iniitialize mesh.
            mesh = TriangleMesh();
            hull = reshape(opt.edge, 2 * size(opt.edge, 2), 1);
            idx  = circshift(reshape(repmat(0:size(opt.edge, 2)-1, 2, 1),...
                2 * size(opt.edge, 2), 1), [-1,1]);
            mesh.set_points_tri(hull);
            mesh.set_facets_tri(idx);
            mesh = mesh.build_tri(); % not ready to go
            mesh = mesh.refine_tri(sprintf('q34.0a%f', opt.min_area));  % ready.

            %% Function Space
            obj.space = FunctionSpace(mesh, opt.deg, opt.kway);
            obj.quadrature = obj.builder.qnodes2D(obj.space.nodes, obj.facet.qnodes, obj.space.elems);            
            
            %% Data : det(Area), Jacobian
            obj.data.area = arealist(obj.space.nodes, obj.space.elems);
            assert(all(obj.data.area > 0)); % orientation. 
        end        
    end
    
end

