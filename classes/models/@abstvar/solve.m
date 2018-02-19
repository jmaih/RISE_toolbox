function [sol,retcode]=solve(self,param,yt)

if nargin<3
    
    yt=[]; % data needed for endogenous probabilities
    
    if nargin<2
        
        param=[];
        
    end
    
end

retcode=0;
        
if isempty(param),param=self.estim_.estim_param; end

if isstruct(param)
    
    if isfield(param,{'B','S','Q'})
        % these are the core fields. There will be others in a SVAR
        
        sol=param;
        
        return
        
    else
        
        error('wrong field names in parameters')
        
    end
    
end

if isempty(yt)
    
    first_obs=vartools.xpand_panel(self.estim_.Y(:,1,:));
    
    yt=first_obs;
    
end

M=vartools.estim2states(param,...
    self.estim_.links.theMap,...
    self.mapping.nparams,...
    self.mapping.nregimes);
        
[sol,retcode]=vartools.solve(self.mapping,M,yt);

end

