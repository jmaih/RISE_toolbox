function v=simple_endo_priors(obj,filtration) %#ok<INUSD>

% this file demonstrates how to setup a simple endogenous prior problem

nconstr=5;

is_initial=nargin==0;

if is_initial
    
    v = cell(nconstr,1);
    
else
    
    myirfs=irf(obj);
    
    C_A=double(myirfs.EPS_A.C);
    
    v = zeros(nconstr,1);
    
end

pointer=0;

for ii=1:nconstr
    
    pointer=pointer+1;
    
    if is_initial

        v{pointer}={0.005,2*0.005,'gamma'};
        
    else
        
        v(pointer)=C_A(pointer+1);
        
    end
    
end

if is_initial
    
    v = struct('priors',{v},...
        'kf_filtering_level',0);
    
end

end