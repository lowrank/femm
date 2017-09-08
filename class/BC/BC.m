classdef BC < handle

properties (Access = private)
  address % ID of the session instance.
end

methods
  function this = BC(varargin)
  	if nargin == 1
    	this.address = BCWrapper('new', varargin{1});
    elseif nargin == 0
    	this.address = BCWrapper('placeholder');
    end
  end

  function setDirichlet(this, edges)
      BCWrapper('set_dirichlet', this.address, edges);
  end
  
  function delete(this)
      BCWrapper('delete', this.address);
  end
    
  
  function report(this)
      BCWrapper('report', this.address);
  end
  
  function [dof, ndof] = dofs(this, N)
      [dof, ndof]  = BCWrapper('dofs', this.address, N);
  end
  
  function [ndof] = getDirichlet(this)
      ndof = BCWrapper('get_dirichlet', this.address);
  end
  function [] = set_boundary(this, expr)
  	BCWrapper('set_boundary', this.address, expr);
  end
  
  function [varargout] = get_boundary(this, edges, nodes, NargOut)
     varargout = cell(NargOut, 1);
 
    [varargout{:}] = BCWrapper('get_boundary', this.address, edges, nodes);
  end
end
end
