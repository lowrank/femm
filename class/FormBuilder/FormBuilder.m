classdef FormBuilder < handle
    
    properties (Access = public)
        address
        r
        rx
        ry
        w
        re
        we
    end
    
    methods 
        function this = FormBuilder(fs, qdeg)
            this.address = FormWrapper('new');
            qm = QuadMode(2);
            [x_, y_, this.w] = qm.get_vr_quad(qdeg); 
            q = [x_, y_]';
          
            % brute-force for reference element. do not change.
            reftm = TriangleMesh();
            reftm.set_points_tri([0 0 1 0 0 1]');
            reftm.set_facets_tri([0 1 1 2 2 0]');
            reftm = reftm.build_tri();
            reftm = reftm.refine_tri('q34a1');
            
            [refp, refs, reft, refe, refn] = reftm.getData_tri();
            
            reffs = FunctionSpace(reftm, fs.deg);
            reffs.build(refp, refs, reft, refe, refn, fs.deg);
            
            p = reffs.nodes;
            [this.r, this.rx, this.ry] = ref2d(this, p, q);
            
            qm = QuadMode(1);
            
            pe = linspace(-1, 1, fs.deg + 1);
            pe = [pe(1), pe(end), pe(2:end-1)];
            
            [x_, this.we] =  qm.get_vr_quad(qdeg);
            qe = x_';
            
            [this.re] = ref1d(this, pe, qe);
            
            reffs.delete();
            reftm.delete();
                
        end
        
        function delete(this)
            FormWrapper('delete', this.address);
        end
        
        function [r, rx, ry] = ref2d(this, p, q)
            [r, rx, ry] = FormWrapper('ref2d', this.address, p, q);
        end
        
        function [r] = ref1d(this,p, q)
            r = FormWrapper('ref1d', this.address, p,q);
        end
    end
    
end