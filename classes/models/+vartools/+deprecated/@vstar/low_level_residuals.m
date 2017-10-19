function [resid,G,retcode,obj]=low_level_residuals(obj)

if isempty(obj)
    
    resid=struct();
    
    return
    
end

[G,retcode,obj]=low_level_transitions(obj);

resid=[];

if retcode
    
    return
    
end

variables_locations_in_data=obj.variables_locations_in_data;

endo_data=obj.data(variables_locations_in_data.endo_id,:);

exo_data=obj.data(variables_locations_in_data.det_id,:);

[y,X,~,T]=vartools.set_y_and_x(endo_data,exo_data,obj.nlags,obj.constant);
        
p=obj.endogenous.number;

resid=zeros(p,T); 

s=obj.solution;

for t=1:T
    
    Bt=s.B(:,:,1);
    
    for ireg=2:size(s.B,3)
        
        Gt=G(ireg-1,t);
        
        Bt=Bt+Gt*s.B(:,:,ireg);
        
    end
    
    ut=y(:,t)-Bt*X(:,t);
    
    resid(:,t)=ut;
    
end


end