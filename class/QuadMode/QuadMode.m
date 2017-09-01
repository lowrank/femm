classdef QuadMode < handle

    properties (Access = private)
        address
        dimension
    end

    methods
        function this = QuadMode(dim)
            this.address = ModeWrapper('new', dim);
            this.dimension = dim;
            logger('*Green', sprintf('Initializing QuadMode at %d.\n', this.address));
        end

        function delete(this)
            if (this.dimension == 2)
                logger('*Red', sprintf('  Destroying QuadMode at %d.\n', this.address));
                ModeWrapper('delete_tri', this.address);
            else
                logger('*Red', sprintf('  Destroying QuadMode at %d.\n', this.address));
                ModeWrapper('delete_tet', this.address);
            end
        end

        function [output] = getAddress(this)
            output = this.address;
        end

        function [varargout] = get_vr_quad(this, degree)
            if (this.dimension == 2)
                varargout = cell(3, 1);
                [varargout{1}, varargout{2}, varargout{3}] = ModeWrapper('load_vr_tri', this.address, degree);
                [varargout{1}, varargout{2}, varargout{3}] = this.affine_tri(varargout{1}, varargout{2}, varargout{3});
            else
                varargout = cell(4, 1);
                [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = ModeWrapper('load_vr_tet', this.address, degree);
                [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = this.affine_tet(varargout{1}, varargout{2}, varargout{3}, varargout{4});
            end
        end
    end

    methods (Static)
        function help()
            ModeWrapper('help');
            cprintf('*Cyan', '- help : get help string. \n'); 
            cprintf('*Cyan', '- get_vr_quad : get unscaled quadrature nodes and weights. \n');
            cprintf('*Cyan', '- getAddress : get address of pointer. \n');
            cprintf('*Cyan', '- delete : call destructor automatically. \n');
            cprintf('*Cyan', '- affine : map quadrature points to reference triangle. \n');
        end

        % todo : change following function adapt to 3d.



        function [x, y, w] = affine_tri(x, y, w)
            A = [1./2. -1./sqrt(3)/2.; 0. 1.0/sqrt(3.0)];
            b = [1./3.; 1./3.];
            ref = 0.5;
            
            xy = A * [x y]' + repmat(b, 1, size(x, 1));
            w = ref * w/(sum(w));
            x = xy(1, :)';
            y = xy(2, :)';
        end

        function [x,y,z,w] = affine_tet(x,y,z,w)
            
            A = [1.  -1./sqrt(3) -1./sqrt(6); 0. 2.0/sqrt(3) -1./sqrt(6); 0. 0. sqrt(6)/2.0];
            b = [-1./2.; -1./2.; -1./2.];
            
            ref = 1.0/6.0;
            
            xyz = A*[x y z]' + repmat(b, 1, size(x, 1));
            w = ref * w/(sum(w));
            x = (xyz(1,:)' + 1)./2.;
            y = (xyz(2,:)' + 1)./2.;
            z = (xyz(3,:)' + 1)./2.;
        end
    end

end 