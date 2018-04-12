classdef femm < handle    
    properties
        facet
        edge
        space
        quad2d
        quad1d
        builder
        data
        rffs
    end
    
    methods
        function obj = femm(opt)
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
            obj.rffs = rffs;
            %% Geometry
            % fill the information of facet.
            [obj.facet.qnodes, obj.facet.weights] = ...
                QuadMode(2, d).export();
            [obj.facet.ref, obj.facet.refx, ...
                obj.facet.refy] = obj.builder.reference2D(rffs.nodes,...
                obj.facet.qnodes);
            
            % fill the information of edge.
            [obj.edge.qnodes, obj.edge.weights] = ...
                QuadMode(1,d).export();
            
            [obj.edge.ref, ~] = obj.builder.reference1D(opt.deg, obj.edge.qnodes);
            
            % iniitialize mesh.
            if ~isfield(opt, 'mesh')
                mesh = TriangleMesh();
                hull = reshape(opt.edge, 2 * size(opt.edge, 2), 1);
                idx  = circshift(reshape(repmat(0:size(opt.edge, 2)-1, 2, 1),...
                    2 * size(opt.edge, 2), 1), [-1,1]);
                mesh.set_points_tri(hull);
                mesh.set_facets_tri(idx);
                mesh = mesh.build_tri(); % not ready to go
                mesh = mesh.refine_tri(sprintf('q34.0a%f', opt.min_area));  % ready.
            else
                mesh = opt.mesh;
            end

            %% Function Space
            if (isfield(opt, 'kway') && opt.kway > 1)
                obj.space = FunctionSpace(mesh, opt.deg, opt.kway);
            else
                obj.space = FunctionSpace(mesh, opt.deg);
            end
            obj.quad2d = obj.builder.qnodes2D(obj.space.nodes, obj.facet.qnodes, obj.space.elems);   
            obj.quad1d = obj.builder.qnodes1D(obj.space.nodes, obj.edge.qnodes, obj.space.edges);
            
            %% Data : det(Area), Jacobian
            obj.data.area = arealist(obj.space.nodes, obj.space.elems);
            assert(all(obj.data.area > 0)); % orientation. 
            
        end
        
        %% set boundaries
        
        
        
        %% build matrix and vector
        function ret = build(obj, flag, F, bc)
            if flag == 'm'
                [I, J, V] = obj.builder.assema(obj.space.nodes, obj.space.elems,...
                    obj.facet.ref, obj.facet.weights, F, obj.data.area);
                ret = sparse(I,J,V);
            elseif flag == 's'
                [I, J, V] = obj.builder.assems(obj.space.nodes, obj.space.elems,...
                    obj.facet.refx, obj.facet.refy, obj.facet.weights, F, obj.data.area);
                ret = sparse(I,J, V);
            elseif flag == 'e'
                if strcmp(bc, 'all')
                    [I, J, V] = obj.builder.assemlbc(obj.space.nodes, obj.space.edges, ...
                        obj.edge.ref, obj.edge.weights, F);
                    n = size(obj.space.nodes, 2);
                    ret = sparse(I, J, V, n,n);
                else
                    [I, J, V] = obj.builder.assemlbc(obj.space.nodes, bc, ...
                    obj.edge.ref, obj.edge.weights, F);
                    n = size(obj.space.nodes, 2);
                    ret = sparse(I, J, V, n,n);
                end
            elseif flag == 'l'
                ret = obj.builder.asseml(obj.space.nodes, obj.facet.qnodes, obj.space.elems, obj.facet.ref, obj.facet.weights, F, obj.data.area);
            elseif flag == 'g'
                if strcmp(bc, 'all')
                    ret = obj.builder.assemrbc(obj.space.nodes, obj.edge.qnodes, obj.space.edges, ...
                        obj.edge.ref, obj.edge.weights, F);
                else
                    ret = obj.builder.assemrbc(obj.space.nodes, obj.edge.qnodes, bc, ...
                        obj.edge.ref, obj.edge.weights, F);
                end
            else
                error('FEMM:BUILD:Unrecongnized flag');
            end
                
        end
        
        function  [gx, gy] = gradient(obj, u, DX, DY)
            [gx, gy] = obj.builder.assemble_grad(u, obj.space.nodes,...
                obj.space.elems, DX, DY);
        end
        
        %% gradient use.
        function w = adj(obj, u, v, Fs, Fa)
            w = obj.builder.assemble_node(obj.space.nodes, obj.space.elems, ...
                obj.facet.ref, obj.facet.refx, obj.facet.refy,...
                obj.facet.weights, Fs, Fa, u, v);
        end
    end
    
end

