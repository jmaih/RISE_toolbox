classdef rfvar < abstvar
    
    properties(Constant,Hidden)
        
        optimize = false
        
    end
    
    methods(Access=private)
        
        varargout=read_identification_restrictions(varargin)
        
    end
    
    methods
        
        function self=rfvar(varargin)
            
            self=self@abstvar(varargin{:});
            
            if nargin>0
                
                self=abstvar.recreate_parameters(self,1);
                
            end
            
        end
                
        varargout=identification(varargin)
                
        varargout=solve(varargin)
                
        varargout=structural_shocks(varargin)
        
        varargout=bootstrap(varargin)
        
    end
    
end
