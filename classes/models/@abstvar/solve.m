function [sol,M]=solve(self,param)

if nargin<2,param=[]; end

if isempty(param),param=self.estim_param; end

if isstruct(param)
    
    if isfield(param,{'B','S','Q'})
        % these are the core fields. There will be others in a SVAR
        
        sol=param;
        
        M=[];
        
        return
        
    else
        
        error('wrong field names in parameters')
        
    end
    
end

[sol,M]=vartools.solve(self.mapping,param);

end

