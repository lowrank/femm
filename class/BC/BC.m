classdef BC < handle

properties (Access = private)
  address % ID of the session instance.
end

methods
  function this = BC(id)
      this.address = BCWrapper('new', id);
  end

  function delete(this)
      BCWrapper('delete', this.address);
  end
  
  function set_constraint(this, str)
      BCWrapper('push_expr', this.address, str);
  end
    
  function [varargout] = get_boundary(this, edges, nodes, Nout) 
      varargout = cell(Nout, 1);
      [varargout{:}] = BCWrapper('filter_segments', this.address, edges, nodes);
  end
end
end
