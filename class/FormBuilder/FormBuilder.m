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
        
        function [r, rx, ry] = ref(this, p, q)
            [r, rx, ry] = FormWrapper('ref', this.address, p, q);
        end
    end
    
end