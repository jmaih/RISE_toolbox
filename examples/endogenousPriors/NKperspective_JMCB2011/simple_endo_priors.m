function v=simple_endo_priors(obj,filtration) %#ok<INUSD>

% this file demonstrates how to setup a simple endogenous prior problem

nconstr=5;

if nargin==1
    
    v = cell(nconstr,1);

	for ii=1:nconstr
	    	    
	    v{ii}={0.005,2*0.005,'gamma'};
	            
	end

    v = struct('priors',{v},...
        'kf_filtering_level',0);
    
else
    
    myirfs=irf(obj);
    
    C_A=double(myirfs.EPS_A.C);
    
    v = zeros(nconstr,1);

	pointer=0;
	
	for ii=1:nconstr
	    
	    pointer=pointer+1;
	    
	    v(pointer)=C_A(pointer+1);
	            
	end
    
end


end