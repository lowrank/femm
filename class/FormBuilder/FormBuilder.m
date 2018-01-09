classdef FormBuilder < handle

properties (Access = public)
  address 
end

methods
  function this = FormBuilder()
    this.address = FormWrapper('new');
  end

  function delete(this)
    FormWrapper('delete', this.address);
  end
  
  function [F, DX, DY] = reference2D(this, nodes, qnodes)
  	[F, DX, DY] = FormWrapper('reference2D', this.address, nodes, qnodes);
  end
  
  function [DX, DY] = reference_grad(this, nodes)
      [DX, DY] = FormWrapper('reference_grad', this.address, nodes);
  end
  
  function [F, DX] = reference1D(this, deg, qnodes)
  	[F, DX] = FormWrapper('reference1D', this.address, deg, qnodes);
  end
  
  function [I, J, V] = assema(this, pnodes, pelems, ref_fnk, weights, extern, Arealist)
  	[I, J, V] = FormWrapper('assema', this.address, pnodes, pelems, ref_fnk, weights, extern, Arealist);
  end
  
  function [I, J, V] = assems(this, pnodes, pelems, ref_gradx, ref_grady, weights, extern, Arealist)
  	[I, J, V] = FormWrapper('assems', this.address, pnodes, pelems, ref_gradx, ref_grady, weights, extern, Arealist);
  end
  
  function [L] = asseml(this, pnodes, qnodes, pelems, ref, weights, extern, Arealist) 
    [L] = FormWrapper('asseml', this.address, pnodes, qnodes, pelems, ref, weights, extern, Arealist);
  end
  
  function [L] = assemrbc(this, pnodes, qnodes, pedges, ref, weights, extern) 
    [L] = FormWrapper('assemrbc', this.address, pnodes, qnodes, pedges, ref, weights, extern);
  end
  
  function [I, J, V] = assemlbc(this, pnodes, pedges, ref, weights, extern) 
    [I, J, V] = FormWrapper('assemlbc', this.address, pnodes, pedges, ref, weights, extern);
  end
  
  function [C] = qnodes2D(this, pnodes, qnodes, pelems) 
  	[C] = FormWrapper('qnodes2D', this.address, pnodes, qnodes, pelems);
  end
  
  function [C] = qnodes1D(this, pnodes, qnodes, pedges) 
  	[C] = FormWrapper('qnodes1D', this.address, pnodes, qnodes, pedges);
  end
  
  function [I, J ,V] = assemble_grad_x_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern)
     [I, J, V] =  FormWrapper('assemex_gradfunc_x',  this.address, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern);
  end
  
  function [I, J ,V] = assemble_grad_y_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern)
     [I, J, V] =  FormWrapper('assemex_gradfunc_y',  this.address, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern);
  end
  
  function [I, J ,V, W] = assemble_grad_xy_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern_x, extern_y)
     [I, J, V, W] =  FormWrapper('assemex_gradfunc_xy',  this.address, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern_x, extern_y);
  end
  
  function [I, J ,V] = assemble_load_matrix(this, pnodes, pelems, ref, weights, extern)
     [I, J, V] =  FormWrapper('assemex_lm',  this.address, pnodes, pelems, ref, weights, extern);
  end
  
  function [w] = assemble_elem(this, pnodes, pelems, ref, refx, refy, weights, extern_s, extern_a, u, v)
    w = FormWrapper('assemloe', this.address, pnodes, pelems, ref, refx, refy, weights, extern_s, extern_a, u, v);
  end
  
  function [w] = assemble_node(this, pnodes, pelems, ref, refx, refy, weights, extern_s, extern_a, u, v)
    w = FormWrapper('assemlon', this.address, pnodes, pelems, ref, refx, refy, weights, extern_s, extern_a, u, v);
  end
  
  function [gx, gy] = assemble_grad(this, u, pnodes, pelems, refx, refy)
      [gx, gy] = FormWrapper('grad', this.address, u, pnodes, pelems, refx, refy);
  end
  
  % Other methods...
end
end
