function v=cov_endo_priors(obj)

% this file demonstrates how to setup a simple endogenous prior problem

myList={'GHAT','PAIHAT','RHAT'};

if nargin==0
    
    sd=struct();
    
    sd.GHAT=0.006584682257472;
    
    sd.PAIHAT=0.002593276205803;
    
    sd.RHAT=0.006167923408467;
    
    n=numel(myList);
    
    v=cell(n,1);
    
    for ii=1:n
        
        vname=myList{ii};
        
        v{ii}={sd.(vname),0.0015,'gamma'};
        
    end
    
    
    v = struct('priors',{v},...
        'kf_filtering_level',0);
    
else
    
    pos=locate_variables(myList,obj.endogenous.name);
    
    % setting the number of autocovariances to 0 such that only the
    % variances are computed.
    ac=theoretical_autocovariances(obj,'autocov_ar',0);
    
    variances=diag(ac(pos,pos,1));
    
    mystdevs=sqrt(variances);
    
    v = mystdevs;
    
end


end